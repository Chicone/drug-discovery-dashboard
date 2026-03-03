# backend/services/analysis/activation.py

import pathlib
from typing import Dict, Iterable, List, Tuple

import mdtraj as md
import numpy as np


def _select_residues_atoms(
    traj: md.Trajectory,
    resids: Iterable[int],
) -> np.ndarray:
    """
    Return atom indices for a list of residue indices.

    Note:
        This assumes mdtraj "resSeq" / "resid" matches your numbering.
        You may need to adapt the selection string to your topology
        (eg using resSeq instead of resid).
    """
    top = traj.topology
    # Example for resid selection. Change to resSeq if needed.
    sel_terms = [f"resid {r}" for r in resids]
    atom_indices = top.select(" or ".join(sel_terms))

    if atom_indices.size == 0:
        raise ValueError(
            f"No atoms found for residues {list(resids)}. "
            "Check numbering (resid vs resSeq) and selection string."
        )

    return atom_indices


def compute_segment_com(
    traj: md.Trajectory,
    resids: Iterable[int],
) -> np.ndarray:
    """
    Compute COM per frame for all atoms in the given residue set.

    Returns:
        np.ndarray of shape (n_frames, 3)
    """
    atom_indices = _select_residues_atoms(traj, resids)

    masses = np.array(
        [atom.element.mass for atom in traj.topology.atoms]
    )[atom_indices]

    coords = traj.xyz[:, atom_indices, :]  # (n_frames, n_atoms, 3)

    mass_sum = masses.sum()
    com = (coords * masses[None, :, None]).sum(axis=1) / mass_sum

    return com


def compute_helix_pair_distance(
    traj: md.Trajectory,
    helix_a_resids: Iterable[int],
    helix_b_resids: Iterable[int],
) -> np.ndarray:
    """
    Distance between COMs of two helix IC segments per frame.

    This gives you a classic TM3-TM6 IC distance trace.
    """
    com_a = compute_segment_com(traj, helix_a_resids)
    com_b = compute_segment_com(traj, helix_b_resids)

    dist = np.linalg.norm(com_a - com_b, axis=1)
    return dist  # shape (n_frames,)


def compute_segment_displacement_along_axis(
    traj: md.Trajectory,
    resids: Iterable[int],
    axis: Tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> np.ndarray:
    """
    Project COM displacement of a helix segment along a given axis.

    Typical use:
        - axis ≈ membrane normal (0, 0, 1) in OPM-oriented systems
        - negative values = "downwards" if you define axis consistently.

    Returns:
        np.ndarray of shape (n_frames,) with displacement in nm.
    """
    axis_vec = np.asarray(axis, dtype=float)
    axis_vec /= np.linalg.norm(axis_vec)

    com = compute_segment_com(traj, resids)  # (n_frames, 3)
    com0 = com[0]                            # reference frame

    disp = com - com0                        # (n_frames, 3)
    proj = disp @ axis_vec                   # (n_frames,)

    return proj


def aggregate_replicas_1d(
    replicas: List[Dict[str, List[float]]],
    key: str = "values",
) -> Dict[str, List[float]]:
    """
    Aggregate 1D traces across replicas, using min length across them.

    Expects each replica to be of the form:
        { "name": "...", "values": [float, ...] }
    or "key" set to something else if needed.
    """
    if not replicas:
        raise ValueError("No replicas provided for aggregation.")

    min_len = min(len(r[key]) for r in replicas)

    trimmed = [
        np.asarray(r[key][:min_len], dtype=float) for r in replicas
    ]

    stacked = np.vstack(trimmed)  # (n_reps, min_len)

    mean = stacked.mean(axis=0)
    min_v = stacked.min(axis=0)
    max_v = stacked.max(axis=0)

    return {
        "mean": mean.tolist(),
        "min": min_v.tolist(),
        "max": max_v.tolist(),
    }

import mdtraj as md
import numpy as np
from typing import Iterable, List, Dict, Tuple


# -------------------------
# Default correct resSeq ranges for A2A AR IC tips
# (You can override via parameters)
# -------------------------
A2A_TM3_DEFAULT = list(range(126, 136))   # around R3.50 (IC tip of TM3)
A2A_TM6_DEFAULT = list(range(224, 236))   # around L6.30 (IC tip of TM6)


# -------------------------
# Helpers
# -------------------------

def _select_by_resSeq(top: md.Topology, resSeqs: Iterable[int]) -> np.ndarray:
    resSeqs = set(resSeqs)
    idx = [atom.index for atom in top.atoms if atom.residue.resSeq in resSeqs]
    if not idx:
        raise ValueError(f"No atoms found for resSeq values {resSeqs}")
    return np.array(idx, dtype=int)


def _mass_com(coords, masses):
    total = masses.sum()
    return (coords * masses[:, None]).sum(axis=0) / total


def _axis_from_box(traj: md.Trajectory):
    box = traj.unitcell_vectors.mean(axis=0)
    lengths = np.linalg.norm(box, axis=1)
    axis = box[np.argmin(lengths)]
    axis = axis / np.linalg.norm(axis)
    return axis


# -------------------------
# MAIN FUNCTION (same signature)
# -------------------------

def compute_activation_metrics_multi(
    traj_paths,
    tm3_ic_resids=None,
    tm6_ic_resids=None,
    axis=(0, 0, 1),
):

    results = []

    for job_dir in traj_paths:
        job_id = job_dir.name
        out = job_dir / "out"

        # topology
        system = out / "npt.gro"
        if not system.exists():
            raise RuntimeError(f"[{job_id}] Missing npt.gro")

        # detect trajectory
        traj_file = None

        parts = sorted(out.glob("md.part*.xtc"))
        if parts:
            traj_file = parts[-1]

        if traj_file is None:
            for name in [
                "md_nojump.xtc", "md_centered.xtc", "md_fit.xtc",
                "md_pbc.xtc", "md_unwrapped.xtc", "md.xtc",
            ]:
                p = out / name
                if p.exists() and p.stat().st_size > 0:
                    traj_file = p
                    break

        if traj_file is None:
            raise RuntimeError(f"[{job_id}] No trajectory found in {out}")

        print(f"[DEBUG] Loading {traj_file}")

        traj = md.load(str(traj_file), top=str(system))
        times_ps = traj.time.copy()

        # Align on protein
        prot_idx = traj.top.select("protein")
        traj.superpose(traj, 0, atom_indices=prot_idx)

        # Axis
        if axis == (0, 0, 1):
            axis_vec = _axis_from_box(traj)
        else:
            axis_vec = np.array(axis, dtype=float)
            axis_vec /= np.linalg.norm(axis_vec)

        print(f"[DEBUG] Axis used for {job_id}: {axis_vec}")

        # Determine TM3/TM6 resSeq lists
        tm3_res = tm3_ic_resids if tm3_ic_resids is not None else A2A_TM3_DEFAULT
        tm6_res = tm6_ic_resids if tm6_ic_resids is not None else A2A_TM6_DEFAULT

        print(f"[DEBUG] TM3 residues for {job_id}: {tm3_res}")
        print(f"[DEBUG] TM6 residues for {job_id}: {tm6_res}")

        # Select atoms (resSeq!)
        tm3_idx = _select_by_resSeq(traj.top, tm3_res)
        tm6_idx = _select_by_resSeq(traj.top, tm6_res)

        masses = np.array([atom.element.mass for atom in traj.top.atoms])

        # COM per frame
        tm3_com = np.array([
            _mass_com(traj.xyz[i, tm3_idx, :], masses[tm3_idx])
            for i in range(traj.n_frames)
        ])
        tm6_com = np.array([
            _mass_com(traj.xyz[i, tm6_idx, :], masses[tm6_idx])
            for i in range(traj.n_frames)
        ])

        # IC distance
        ic_dist = np.linalg.norm(tm6_com - tm3_com, axis=1)

        # Displacement
        disp = np.dot(tm6_com - tm6_com[0], axis_vec)

        results.append({
            "job_id": job_id,
            "tm3_tm6_distance": ic_dist.tolist(),
            "tm6_displacement": disp.tolist(),
            "times_ps": times_ps.tolist(),
            "timestep_ps": float(times_ps[1] - times_ps[0]) if traj.n_frames > 1 else 0.0,
        })

    # Aggregate
    min_len = min(len(r["tm3_tm6_distance"]) for r in results)
    dist = np.array([r["tm3_tm6_distance"][:min_len] for r in results])
    disp = np.array([r["tm6_displacement"][:min_len] for r in results])

    aggregate = {
        "tm3_tm6_distance": {
            "mean": dist.mean(axis=0).tolist(),
            "min":  dist.min(axis=0).tolist(),
            "max":  dist.max(axis=0).tolist(),
        },
        "tm6_displacement": {
            "mean": disp.mean(axis=0).tolist(),
            "min":  disp.min(axis=0).tolist(),
            "max":  disp.max(axis=0).tolist(),
        },
    }

    return {"replicas": results, "aggregate": aggregate}


def compute_activation_metrics_multi_old(
    traj_paths,
    tm3_ic_resids=None,
    tm6_ic_resids=None,
    axis=(0, 0, 1),
):

    # Default GPCR IC residues (last 10 of TM3, first 10 of TM6)
    if tm3_ic_resids is None:
        tm3_ic_resids = list(range(99, 109))    # example only; adjust correctly
    if tm6_ic_resids is None:
        tm6_ic_resids = list(range(217, 227))   # example only; adjust correctly

    axis = np.array(axis) / np.linalg.norm(axis)

    replicas = []

    for job_dir in traj_paths:
        job_id = job_dir.name

        # -------------------------
        # Load trajectory robustly
        # -------------------------
        traj, times_ps = load_job_traj(job_dir)

        # Align traj here
        traj.superpose(traj, 0)


        # -------------------------
        # Compute TM3-TM6 IC distance
        # -------------------------
        tm3_idx = traj.top.select("resid " + " ".join(map(str, tm3_ic_resids)))
        tm6_idx = traj.top.select("resid " + " ".join(map(str, tm6_ic_resids)))

        tm3_com = traj.xyz[:, tm3_idx, :].mean(axis=1)
        tm6_com = traj.xyz[:, tm6_idx, :].mean(axis=1)

        ic_dist = np.linalg.norm(tm6_com - tm3_com, axis=1)

        # -------------------------
        # TM6 IC displacement along axis
        # -------------------------
        proj = np.dot(tm6_com - tm6_com[0], axis)

        replicas.append({
            "job_id": job_id,
            "tm3_tm6_distance": ic_dist.tolist(),
            "tm6_displacement": proj.tolist(),
            "times_ps": times_ps.tolist(),
            "timestep_ps": float(times_ps[1] - times_ps[0]) if len(times_ps) > 1 else 0.0
        })

    # -------------------------
    # Aggregate: align length
    # -------------------------
    min_len = min(len(r["tm3_tm6_distance"]) for r in replicas)
    tm3tm6_trimmed = [r["tm3_tm6_distance"][:min_len] for r in replicas]

    # -------------------------
    # Aggregate: align length
    # -------------------------
    min_len = min(len(r["tm3_tm6_distance"]) for r in replicas)

    dist_trimmed = np.array([r["tm3_tm6_distance"][:min_len] for r in replicas])
    disp_trimmed = np.array([r["tm6_displacement"][:min_len] for r in replicas])

    aggregate = {
        "tm3_tm6_distance": {
            "mean": dist_trimmed.mean(axis=0).tolist(),
            "min": dist_trimmed.min(axis=0).tolist(),
            "max": dist_trimmed.max(axis=0).tolist(),
        },
        "tm6_displacement": {
            "mean": disp_trimmed.mean(axis=0).tolist(),
            "min": disp_trimmed.min(axis=0).tolist(),
            "max": disp_trimmed.max(axis=0).tolist(),
        }
    }

    return {
        "replicas": replicas,
        "aggregate": aggregate,
    }

import numpy as np
import mdtraj as md

def load_job_traj(job_dir):
    """
    Identical trajectory detection logic used in COM / RMSD analysis.
    Returns:
        traj (md.Trajectory),
        times_ps (ndarray)
    """

    job_id = job_dir.name

    # -------------------------
    # Topology (npt.gro or system.gro)
    # -------------------------
    system = job_dir / "out" / "npt.gro"
    if not system.exists():
        raise RuntimeError(f"[{job_id}] Missing npt.gro")

    # -------------------------
    #  Trajectory detection (same as compute_ligand_com_multi)
    # -------------------------
    traj_file = None

    # 1) continuation files
    part_files = sorted((job_dir / "out").glob("md.part*.xtc"))
    if part_files:
        traj_file = part_files[-1]

    # 2) fallback standard names
    if traj_file is None:
        for name in [
            "md_nojump.xtc",
            "md_centered.xtc",
            "md_fit.xtc",
            "md_pbc.xtc",
            "md_unwrapped.xtc",
            "md.xtc",
        ]:
            p = job_dir / "out" / name
            if p.exists() and p.stat().st_size > 0:
                traj_file = p
                break

    if traj_file is None:
        raise RuntimeError(f"[{job_id}] No trajectory found in {job_dir}/out")

    print(f"[DEBUG] Loading traj for job {job_id}: {traj_file.name}")

    traj = md.load(str(traj_file), top=str(system))
    return traj, traj.time