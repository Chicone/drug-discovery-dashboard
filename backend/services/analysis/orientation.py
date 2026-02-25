from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import mdtraj as md
import numpy as np


def _pick_traj_file(job_dir: Path) -> Path:
    for name in [
        "md_nojump.xtc",
        "md_fit.xtc",
        "md_centered.xtc",
        "md_pbc.xtc",
        "md_unwrapped.xtc",
        "md.xtc",
    ]:
        p = job_dir / "out" / name
        if p.exists() and p.stat().st_size > 0:
            return p
    raise FileNotFoundError("No trajectory found in job_dir/out/.")


def _detect_ligand_resname(system_gro: Path) -> str:
    # Keep it simple: you said your Martini ligands are ORT.
    # If you want auto-detect later, we can improve.
    return "ORT"


def _principal_axis(coords: np.ndarray) -> np.ndarray:
    """
    coords: (n_atoms, 3) for one frame
    Returns unit vector (3,) corresponding to largest-variance axis.
    """
    x = coords - coords.mean(axis=0, keepdims=True)
    cov = (x.T @ x) / max(1, x.shape[0] - 1)
    vals, vecs = np.linalg.eigh(cov)
    v = vecs[:, np.argmax(vals)]
    v = v / (np.linalg.norm(v) + 1e-12)
    return v


def _angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    c = float(np.clip(np.dot(u, v), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def compute_ligand_orientation_single(
    job_dir: Path,
    ligand_resname: Optional[str] = None,
    ref_resids: Tuple[int, int] = (168, 253),
    contact_cutoff_nm: float = 0.6,
) -> Dict:
    """
    Metrics:
      - ligand principal-axis angle vs Z (0..90 deg, symmetrized)
      - ligand COM z (and relative to protein COM z)
      - contacts to residues 168 and 253 (distance min per frame)
    """
    system = job_dir / "out" / "system.gro"
    if not system.exists():
        return {"error": f"Missing file: {system}"}

    try:
        traj_path = _pick_traj_file(job_dir)
    except FileNotFoundError as e:
        return {"error": str(e)}

    if ligand_resname is None:
        ligand_resname = _detect_ligand_resname(system)

    t = md.load(str(traj_path), top=str(system))
    top = t.topology

    lig_idx = top.select(f"resname {ligand_resname}")
    if lig_idx.size == 0:
        return {"error": f"No atoms found for ligand resname={ligand_resname}"}

    prot_idx = top.select("protein")
    if prot_idx.size == 0:
        return {"error": "No protein atoms found (top.select('protein') empty)"}

    # Residue selections for contacts (by residue number in the topology)
    r1, r2 = ref_resids
    r1_idx = top.select(f"protein and resSeq {r1}")
    r2_idx = top.select(f"protein and resSeq {r2}")
    # r1_idx = top.select(f"protein and resid {r1}")
    # r2_idx = top.select(f"protein and resid {r2}")

    if r1_idx.size == 0:
        return {"error": f"No atoms for protein resid {r1} (check numbering)"}
    if r2_idx.size == 0:
        return {"error": f"No atoms for protein resid {r2} (check numbering)"}

    z = np.array([0.0, 0.0, 1.0], dtype=float)

    n = t.n_frames
    times_ps = t.time.astype(float)  # mdtraj time is in ps
    times_ns = (times_ps / 1000.0).tolist()

    angles = []
    lig_com_z = []
    lig_minus_prot_com_z = []

    d168 = []
    d253 = []

    # Precompute per-frame COMs cheaply
    lig_xyz = t.xyz[:, lig_idx, :]  # (n, n_lig, 3)
    prot_xyz = t.xyz[:, prot_idx, :]

    lig_com = lig_xyz.mean(axis=1)   # (n, 3)
    prot_com = prot_xyz.mean(axis=1) # (n, 3)

    # Distances: compute min distance between ligand atoms and residue atoms
    r1_xyz = t.xyz[:, r1_idx, :]
    r2_xyz = t.xyz[:, r2_idx, :]

    for i in range(n):
        v = _principal_axis(lig_xyz[i])
        a = _angle_deg(v, z)
        # Symmetrize axis sign (v and -v are same axis):
        a = min(a, 180.0 - a)
        # Often you want 0..90:
        if a > 90.0:
            a = 180.0 - a

        angles.append(a)

        lig_com_z.append(float(lig_com[i, 2]))
        lig_minus_prot_com_z.append(float(lig_com[i, 2] - prot_com[i, 2]))

        # min pairwise distance for contacts (nm)
        d1 = np.linalg.norm(
            lig_xyz[i][:, None, :] - r1_xyz[i][None, :, :],
            axis=2,
        ).min()
        d2 = np.linalg.norm(
            lig_xyz[i][:, None, :] - r2_xyz[i][None, :, :] ,
            axis=2,
        ).min()
        d168.append(float(d1))
        d253.append(float(d2))

    angles_np = np.array(angles)
    d168_np = np.array(d168)
    d253_np = np.array(d253)

    contacts_168 = float(np.mean(d168_np < contact_cutoff_nm))
    contacts_253 = float(np.mean(d253_np < contact_cutoff_nm))

    return {
        "times_ns": times_ns,
        "angle_deg": angles,
        "lig_com_z_nm": lig_com_z,
        "lig_minus_prot_com_z_nm": lig_minus_prot_com_z,
        "d168_nm": d168,
        "d253_nm": d253,
        "summary": {
            "angle_mean_deg": float(angles_np.mean()),
            "angle_std_deg": float(angles_np.std()),
            "contact_frac_168": contacts_168,
            "contact_frac_253": contacts_253,
            "contact_cutoff_nm": contact_cutoff_nm,
        },
        "meta": {
            "traj_file": str(traj_path.name),
            "ligand_resname": ligand_resname,
            "ref_resids": [r1, r2],
        },
    }

import numpy as np

def compute_ligand_orientation_multi(job_dirs):
    results = {}
    angle_series = []

    for job_dir in job_dirs:
        job_id = job_dir.name
        single = compute_ligand_orientation_single(job_dir)
        results[job_id] = single


        if "angle_deg" in single:
            angle_series.append(np.array(single["angle_deg"]))

    # ---------------------------------------------
    # Build aggregate if possible
    # ---------------------------------------------
    aggregate = None

    if len(angle_series) > 0:
        min_len = min(len(a) for a in angle_series)
        trimmed = np.array([a[:min_len] for a in angle_series])

        aggregate = {
            "mean": trimmed.mean(axis=0).tolist(),
            "min": trimmed.min(axis=0).tolist(),
            "max": trimmed.max(axis=0).tolist(),
        }

    print("AGGREGATE:", aggregate is not None)
    if aggregate:
        print("Aggregate keys:", aggregate.keys())
        print("Length mean:", len(aggregate["mean"]))

    return {
        "replicas": [
            {
                "job_id": job_id,
                "angle_deg": data["angle_deg"],
                "times_ns": data["times_ns"],
            }
            for job_id, data in results.items()
            if "angle_deg" in data
        ],
        "aggregate": aggregate,
    }

