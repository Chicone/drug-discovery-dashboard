import numpy as np
POCKET_RESIDS = [
    84,85,86,88,89,
    168,169,177,181,
    249,253,270,274,
    277,278
]

def compute_com(traj, atom_indices):
    """
    Compute center of mass per frame for selected atoms.
    Martini beads all have mass defined in topology,
    so we use real masses from mdtraj.
    """

    masses = np.array([atom.element.mass for atom in traj.topology.atoms])[atom_indices]
    coords = traj.xyz[:, atom_indices, :]  # (n_frames, n_atoms, 3)

    mass_sum = masses.sum()
    com = (coords * masses[None, :, None]).sum(axis=1) / mass_sum

    return com  # shape (n_frames, 3)

def select_ligand_atoms(traj, ligand_resname="ORT"):
    top = traj.topology
    atom_indices = top.select(f"resname {ligand_resname}")

    if len(atom_indices) == 0:
        raise ValueError(f"Ligand {ligand_resname} not found.")

    return atom_indices

def select_martini_pocket_sc_beads(traj, pocket_resids):
    top = traj.topology

    atom_indices = []

    for atom in top.atoms:
        if atom.residue.resSeq in pocket_resids:
            if atom.name.startswith("SC"):
                atom_indices.append(atom.index)

    if len(atom_indices) == 0:
        raise ValueError("No SC beads found for given pocket residues.")

    return np.array(atom_indices)

def compute_ligand_pocket_com_distance(
    traj,
    pocket_resids,
    ligand_resname="ORT"
):
    """
    Returns distance (nm) per frame between:
    - Ligand COM
    - Orthosteric pocket COM (SC beads only)
    """

    print("DEBUG compute_ligand_pocket_com_distance entered")

    ligand_atoms = select_ligand_atoms(traj, ligand_resname)
    pocket_atoms = select_martini_pocket_sc_beads(traj, pocket_resids)

    ligand_com = compute_com(traj, ligand_atoms)
    pocket_com = compute_com(traj, pocket_atoms)

    distances = np.linalg.norm(ligand_com - pocket_com, axis=1)

    return distances # nm

from pathlib import Path
from typing import List, Dict, Union
import mdtraj as md

from pathlib import Path
import numpy as np
import mdtraj as md

def compute_ligand_com_multi(job_dirs):

    replicas = []

    for job_dir in job_dirs:
        job_id = job_dir.name

        system = job_dir / "out" / "npt.gro"
        if not system.exists():
            raise RuntimeError(f"Missing npt.gro for job {job_id}")

        # ---- Trajectory detection ----
        traj_file = None
        out_dir = job_dir / "out"

        # 1 prefer processed trajectories first
        for name in [
            "md_fit.xtc",
            "md_centered.xtc",
            "md_nojump.xtc",
            "md_whole.xtc",
            "md.xtc",
        ]:
            p = out_dir / name
            if p.exists() and p.stat().st_size > 0:
                traj_file = p
                break

            # 2 fallback to raw continuation files only if needed
            if traj_file is None:
                part_files = sorted(out_dir.glob("md.part*.xtc"))
                if part_files:
                    traj_file = part_files[-1]

        if traj_file is None:
            raise RuntimeError(f"No trajectory found for job {job_id}")

        # ---- Load trajectory ----
        traj = md.load(str(traj_file), top=str(system))

        # ---- Compute COM distance ----
        distances = compute_ligand_pocket_com_distance(
            traj,
            POCKET_RESIDS,
            ligand_resname="ORT"
        )

        times_ps = traj.time

        replicas.append({
            "job_id": job_id,
            "values": distances.tolist(),
            "times_ps": times_ps.tolist(),
            "timestep_ps": float(times_ps[1] - times_ps[0]) if len(times_ps) > 1 else 0.0
        })

    # -------- Aggregate --------
    min_len = min(len(r["values"]) for r in replicas)

    trimmed = [r["values"][:min_len] for r in replicas]

    mean = np.mean(trimmed, axis=0)
    min_v = np.min(trimmed, axis=0)
    max_v = np.max(trimmed, axis=0)

    return {
        "replicas": replicas,
        "aggregate": {
            "mean": mean.tolist(),
            "min": min_v.tolist(),
            "max": max_v.tolist()
        }
    }






