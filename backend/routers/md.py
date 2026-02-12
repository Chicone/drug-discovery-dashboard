from fastapi import APIRouter, UploadFile, File, Form
from typing import Optional

from backend.services.md_jobs import create_md_job_service
# Add these imports at the top if missing:
from fastapi.responses import JSONResponse, FileResponse
from fastapi import Query
from pathlib import Path

from backend.services.md_jobs import (
    list_jobs,
    get_job_status,
    read_log_chunk,
    delete_job,
    stop_job,
    _md_job_dir,
)
from backend.services.md_jobs import run_cgmd_setup_job
from backend.services.md_helpers import MD_RUNS_DIR

# ------------------------------------------------------
# Models
# ------------------------------------------------------
from pydantic import BaseModel

router = APIRouter()

class MartiniSetupRequest(BaseModel):
    workdir: str
    aa_pdb: str
    name: str
    martini_ff: str
    nt: int = 1
    do_em: bool = False
    do_nvt: bool = False
    do_npt: bool = False
    do_md: bool = False

@router.post("/run-cgmd-setup")
def run_cgmd_setup(request: MartiniSetupRequest):
    return run_cgmd_setup_job(request)

@router.post("/jobs")
async def create_md_job(
    protein_pdb: UploadFile = File(...),
    preset: str = Form("cg_popc_50ns"),
    scenario: str = Form("protein_only"),
    environment: str = Form("membrane"),  # "membrane" | "solvated_box"
    workflow: str = Form("build_only"),  # "build_only" | "equilibrate" | "full"
    parent_job_id: Optional[str] = Form(None),
    md_ns: Optional[float] = Form(None),
    nt: int = Form(1),
    orthosteric_ligand: Optional[UploadFile] = File(None),
    allosteric_pose: Optional[UploadFile] = File(None),
):
    return create_md_job_service(
        protein_pdb,
        preset,
        scenario,
        environment,
        workflow,
        parent_job_id,
        md_ns,
        nt,
        orthosteric_ligand,
        allosteric_pose,
    )




@router.get("/jobs")
def list_md_jobs(limit: int = Query(20, ge=0)):
    return list_jobs(limit)


from backend.services.md_jobs import get_job
@router.get("/jobs/{job_id}")
def get_md_job(job_id: str):
    data = get_job(job_id)
    if data is None:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    return data


from fastapi.responses import Response
from backend.services.md_jobs import read_log_window
@router.get("/jobs/{job_id}/log")
def get_md_job_log(
    job_id: str,
    offset: int = Query(0, ge=0),
    chunk_size: int = Query(200_000, ge=1, le=1_000_000),
):
    data, new_offset, window_size = read_log_window(job_id, offset, chunk_size)

    if data is None:
        return Response(
            "",
            media_type="text/plain",
            headers={
                "X-Log-Offset": "0",
                "X-Log-Size": "0",
            },
        )

    return Response(
        content=data,
        media_type="text/plain",
        headers={
            "X-Log-Offset": str(new_offset),
            "X-Log-Size": str(window_size),
        },
    )


from backend.services.md_jobs import list_job_files

@router.get("/jobs/{job_id}/files")
def list_md_job_files(job_id: str):
    files = list_job_files(job_id)
    if files is None:
        return JSONResponse(status_code=404, content={"error": "Job not found"})
    return {"files": files}


from backend.services.md_jobs import create_job_zip

@router.get("/jobs/{job_id}/download")
def download_md_job(job_id: str):
    zip_path = create_job_zip(job_id)
    if zip_path is None:
        return JSONResponse(status_code=404, content={"error": "Job not found"})

    return FileResponse(
        path=str(zip_path),
        filename=f"md_{job_id}.zip",
        media_type="application/zip",
    )

@router.delete("/jobs/{job_id}")
def delete_md_job(job_id: str):
    ok = delete_job(job_id)

    if not ok:
        raise HTTPException(status_code=404, detail="Job not found")

    return {"status": "deleted", "job_id": job_id}


from fastapi import HTTPException
from backend.services.md_jobs import stop_job

@router.post("/jobs/{job_id}/stop")
def stop_md_job(job_id: str):
    ok = stop_job(job_id)

    if not ok:
        return JSONResponse(
            status_code=400,
            content={"error": "No running process found"}
        )

    return {"status": "stopped", "job_id": job_id}

@router.get("/md/log_tail")
def get_log_tail(job_id: str, offset: int = 0):
    job_dir = MD_RUNS_DIR / job_id
    log_path = job_dir / "log.txt"

    if not log_path.exists():
        return {"chunk": "", "offset": offset}

    size = log_path.stat().st_size

    if offset > size:
        # File was truncated or restarted
        offset = 0

    with open(log_path, "rb") as f:
        f.seek(offset)
        chunk = f.read()

    return {
        "chunk": chunk.decode(errors="replace"),
        "offset": offset + len(chunk)
    }







