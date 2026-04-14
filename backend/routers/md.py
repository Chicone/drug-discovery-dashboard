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
from backend.services.md_jobs import create_backmap_job

import subprocess


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

class BackmapRequest(BaseModel):
    time_ps: float = 0.0
    run_em: bool = True
    run_nvt: bool = True
    run_npt: bool = True

@router.post("/run-cgmd-setup")
def run_cgmd_setup(request: MartiniSetupRequest):
    return run_cgmd_setup_job(request)

@router.post("/jobs")
async def create_md_job(
    protein_pdb: Optional[UploadFile] = File(None),
    preset: str = Form("cg_popc_50ns"),
    scenario: str = Form("protein_only"),
    workflow: str = Form("build_only"),
    parent_job_id: Optional[str] = Form(None),
    md_ns: Optional[float] = Form(None),
    nt: int = Form(1),
    orthosteric_ligand: Optional[UploadFile] = File(None),
    orthosteric_smiles: str = Form(None),
    ligand_case: str | None = Form(None),
    allosteric_pose: Optional[UploadFile] = File(None),
    comment: Optional[str] = Form(None),
    ligand_rtf: Optional[UploadFile] = File(None),
    ligand_g_rtf: Optional[UploadFile] = File(None),
    ligand_prm: Optional[UploadFile] = File(None),
    ligand_itp: Optional[UploadFile] = File(None),
    ligand_coord_template: Optional[UploadFile] = File(None),
    ligand_toppar: Optional[UploadFile] = File(None),
):
    return create_md_job_service(
        protein_pdb,
        preset,
        scenario,
        workflow,
        parent_job_id,
        md_ns,
        nt,
        orthosteric_ligand,
        orthosteric_smiles,
        allosteric_pose,
        ligand_case,
        comment,
        ligand_rtf,
        ligand_g_rtf,
        ligand_prm,
        ligand_itp,
        ligand_coord_template,
        ligand_toppar,
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

@router.get("/md/jobs/{job_id}/log")
def get_log(job_id: str, offset: int = 0):

    job_dir = MD_RUNS_DIR / job_id
    log_path = job_dir / "log.txt"

    if not log_path.exists():
        return Response("", headers={"X-Log-Offset": "0"})

    size = log_path.stat().st_size

    # 🔥 Always keep last 2MB window
    MAX_WINDOW = 2_000_000

    window_start = max(0, size - MAX_WINDOW)
    window_size = size - window_start

    # Tail mode
    if offset > window_size:
        offset = window_size

    with open(log_path, "rb") as f:
        f.seek(window_start + offset)
        chunk = f.read()

    return Response(
        chunk,
        media_type="text/plain",
        headers={"X-Log-Offset": str(offset + len(chunk))}
    )

@router.post("/{job_id}/open-chimerax")
def open_chimerax(job_id: str):
    import re
    import subprocess
    from pathlib import Path

    base = _md_job_dir(job_id)
    out_dir = base / "out"

    if not out_dir.exists():
        raise HTTPException(status_code=404, detail="Job out/ directory not found")

    # ---------------------------------------------------------
    # 1. Decide which structure to open
    #    Prefer AA backmapped structure if present
    # ---------------------------------------------------------
    topology = None
    structure_kind = None

    aa_backmapped_pdb = out_dir / "backmapped.pdb"
    aa_backmapped_gro = out_dir / "backmapped.gro"
    aa_frame_gro = out_dir / "frame.gro"
    cg_system_pdb = out_dir / "system.pdb"

    if aa_backmapped_pdb.exists():
        topology = aa_backmapped_pdb
        structure_kind = "aa_backmapped_pdb"

    elif aa_backmapped_gro.exists():
        # Convert GRO -> PDB for ChimeraX
        aa_backmapped_pdb = out_dir / "backmapped.pdb"

        subprocess.run(
            [
                "gmx", "editconf",
                "-f", str(aa_backmapped_gro),
                "-o", str(aa_backmapped_pdb),
            ],
            check=True,
        )

        topology = aa_backmapped_pdb
        structure_kind = "aa_backmapped_gro"

    elif aa_frame_gro.exists():
        aa_frame_pdb = out_dir / "frame.pdb"

        subprocess.run(
            [
                "gmx", "editconf",
                "-f", str(aa_frame_gro),
                "-o", str(aa_frame_pdb),
            ],
            check=True,
        )

        topology = aa_frame_pdb
        structure_kind = "aa_frame_gro"

    elif cg_system_pdb.exists():
        topology = cg_system_pdb
        structure_kind = "cg_system_pdb"

    else:
        raise HTTPException(
            status_code=404,
            detail="No structure file found (expected backmapped.gro/pdb, frame.gro, or system.pdb)",
        )

    # ---------------------------------------------------------
    # 2. Look for trajectory if present
    #    This is optional for AA jobs
    # ---------------------------------------------------------
    part_files = []
    part_files += list(out_dir.glob("md.part*.xtc"))
    part_files += list(out_dir.glob("npt.part*.xtc"))
    part_files += list(out_dir.glob("nvt.part*.xtc"))
    part_files += list(out_dir.glob("em.part*.xtc"))

    chosen = None

    if part_files:
        def part_index(path):
            m = re.search(r"\.part(\d+)\.", path.name)
            return int(m.group(1)) if m else -1

        chosen = max(part_files, key=part_index)

    else:
        candidates = [
            "md.xtc",
            "npt.xtc",
            "nvt.xtc",
            "em.xtc",
        ]

        for name in candidates:
            p = out_dir / name
            if p.exists() and p.stat().st_size > 0:
                chosen = p
                break

    # ---------------------------------------------------------
    # 3. Launch ChimeraX
    # ---------------------------------------------------------
    script = Path("~/PyCharm/ddd/backend/data/chimera/view.cxc").expanduser()

    cmd = [
        "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX",
        str(topology),
    ]

    # Only add script + trajectory if trajectory exists
    if chosen is not None:
        cmd += [str(script), str(chosen)]

    subprocess.Popen(cmd)

    return {
        "status": "launched",
        "structure": topology.name,
        "structure_kind": structure_kind,
        "trajectory": chosen.name if chosen is not None else None,
    }

@router.post("/{job_id}/backmap")
def backmap_md_job(job_id: str, request: BackmapRequest):
    try:
        return create_backmap_job(
            parent_job_id=job_id,
            time_ps=request.time_ps,
            run_em=request.run_em,
            run_nvt=request.run_nvt,
            run_npt=request.run_npt,
        )
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))






