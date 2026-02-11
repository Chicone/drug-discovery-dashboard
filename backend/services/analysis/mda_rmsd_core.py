# backend/services/analysis/mda_rmsd_core.py

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

def compute_ligand_rmsd(
    gro_path: str,
    xtc_path: str,
    ligand_resname: str = "ORT"
):
    """
    Compute ligand RMSD through a trajectory, aligned on the protein backbone (CA or BB).
    Works for Martini or AA MD.

    Returns:
        dict: { "frames": [...], "rmsd": [...], "mean": float }
    """
    # Load
    u = mda.Universe(gro_path, xtc_path)

    # --- selections ---
    # For Martini: BB is backbone bead
    # For AA: CA
    protein_sel = None
    if len(u.select_atoms("protein and name CA")) > 0:
        protein_sel = "protein and name CA"
    else:
        protein_sel = "protein and name BB"

    protein = u.select_atoms(protein_sel)
    ligand = u.select_atoms(f"resname {ligand_resname}")

    if len(ligand) == 0:
        raise ValueError(f"Ligand with resname '{ligand_resname}' not found.")

    # Save reference positions AFTER loading gro
    ref_protein = protein.positions.copy()
    ref_ligand = ligand.positions.copy()

    rmsd_list = []
    frame_ids = []

    for ts in u.trajectory:
        frame_ids.append(ts.frame)

        # Rigid superposition of the entire frame onto reference protein
        align.alignto(u, ref_protein, select=protein_sel, weights=None)

        # RMSD of ligand
        delta = ligand.positions - ref_ligand
        val = np.sqrt((delta * delta).sum(axis=1).mean())
        rmsd_list.append(float(val))

    return {
        "frames": frame_ids,
        "rmsd": rmsd_list,
        "mean": float(np.mean(rmsd_list)),
    }

