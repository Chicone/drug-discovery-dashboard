# backend/routers/analysis.py
from fastapi import APIRouter, Query
from typing import List, Optional
from backend.services.analysis.ligand_rmsd import compute_ligand_rmsd_multi
from backend.services.md_helpers import MD_RUNS_DIR
import mdtraj as md

from backend.services.analysis.distance import (
    compute_ligand_com_multi,
)

router = APIRouter()

@router.get("/ligand_rmsd")
def ligand_rmsd(
    job_ids: list[str] = Query(
        None,
        description="Repeat job_ids=... or provide a single comma-separated value",
    ),
):

    # If someone sends ?job_ids=a,b as a single item, split it.
    if not job_ids:
        return {"error": "job_ids is required"}

    ids: List[str] = []
    for item in job_ids:
        ids.extend([x.strip() for x in item.split(",") if x.strip()])

    job_dirs = [MD_RUNS_DIR / job_id for job_id in ids]

    # -----------------------------------------------------
    # DEBUG: show exactly which directories are being used
    # -----------------------------------------------------
    print("DEBUG ligand_rmsd(): job_ids parsed =", ids)
    print("DEBUG ligand_rmsd(): job_dirs =", job_dirs)

    return compute_ligand_rmsd_multi(job_dirs)


@router.get("/ligand_com_distance")
def ligand_com_distance(
    job_ids: list[str] = Query(
        None,
        description="Repeat job_ids=... or provide a single comma-separated value",
    ),
):
    if not job_ids:
        return {"error": "job_ids is required"}

    ids: List[str] = []
    for item in job_ids:
        ids.extend([x.strip() for x in item.split(",") if x.strip()])

    job_dirs = [MD_RUNS_DIR / job_id for job_id in ids]

    return compute_ligand_com_multi(job_dirs)