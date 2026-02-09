from fastapi import APIRouter, UploadFile, File, Form, Query
from fastapi.responses import JSONResponse, FileResponse

from backend.services import docking_jobs

router = APIRouter(prefix="/api/docking", tags=["Docking"])


@router.post("/vina")
async def dock_vina(
    receptor: UploadFile = File(...),
    ligand: UploadFile = File(...),
    center_x: float | None = Form(None),
    center_y: float | None = Form(None),
    center_z: float | None = Form(None),
    size_x: float = Form(22.0),
    size_y: float = Form(22.0),
    size_z: float = Form(22.0),
    exhaustiveness: int = Form(8),
    num_modes: int = Form(9),
    seed: int | None = Form(None),
):
    return await docking_jobs.run_vina_docking(
        receptor=receptor,
        ligand=ligand,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        seed=seed,
    )


@router.get("/runs")
def list_runs(limit: int = Query(0, ge=0)):
    return docking_jobs.list_runs(limit)


@router.get("/runs/{run_id}")
def get_run(run_id: str):
    return docking_jobs.get_run(run_id)


@router.get("/runs/{run_id}/analyze/{pose_idx}")
def analyze_pose(run_id: str, pose_idx: int, cutoff: float = 4.0):
    return docking_jobs.analyze_pose(run_id, pose_idx, cutoff)


@router.get("/runs/{run_id}/complex/{mode}")
def get_complex(run_id: str, mode: int):
    return docking_jobs.get_complex(run_id, mode)


@router.post("/runs/{run_id}/save/{mode}")
def save_pose(run_id: str, mode: int):
    return docking_jobs.save_pose(run_id, mode)


@router.delete("/runs/{run_id}")
def delete_run(run_id: str):
    return docking_jobs.delete_run(run_id)
