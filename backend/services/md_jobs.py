# backend/services/md_jobs.py

from __future__ import annotations

import shutil
import re
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

from zoneinfo import ZoneInfo
from fastapi.responses import JSONResponse

import sys
import threading
from backend.services.md_helpers import MD_RUNS_DIR

from backend.pipelines.runtime import RUNNING_JOBS


from backend.services.md_helpers import (
    md_job_dir,
    safe_mkdir,
    write_json,
    read_json,
    md_status_from_dir,
)

# ============================================================
# PATHS
# ============================================================
PROJECT_ROOT = next(
    p for p in Path(__file__).resolve().parents
    if (p / "backend").exists()
)
PIPELINE_BY_WORKFLOW = {
    "build_only": "backend.pipelines.pipeline",
    "build_and_equilibrate": "backend.pipelines.pipeline",
    "run_md": "backend.pipelines.pipeline",
    "build_and_run_full": "backend.pipelines.pipeline",
}

# ============================================================
# MAIN JOB CREATION SERVICE
# ============================================================

def create_md_job_service(
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
    ligand_case=None,
    comment=None,
):
    """
    This function contains EXACTLY the logic currently in md.py:create_md_job,
    without modifying behavior.
    """
    # print("DEBUG scenario =", scenario)
    #
    # print("DEBUG protein_pdb =", protein_pdb)
    # print("DEBUG orthosteric_ligand =", orthosteric_ligand)
    # print("DEBUG scenario =", scenario)

    if scenario != "ligand_water" and protein_pdb is None:
        return JSONResponse(
            status_code=400,
            content={"error": "Protein PDB required for this scenario"},
        )

    if ligand_case is None:
        ligand_case = "default"

    try:
        allowed_scenarios = {
            # Protein only
            "protein_only_membrane",
            "protein_only_water",

            # Protein + orthosteric
            "protein_plus_orthosteric_membrane",
            "protein_plus_orthosteric_water",

            # Protein + orthosteric + allosteric
            "protein_plus_orthosteric_plus_allosteric_membrane",
            "protein_plus_orthosteric_plus_allosteric_water",

            # Ligand only
            "ligand_water",

            # Dimer
            "dimer_membrane",
        }

        if scenario not in allowed_scenarios:
            return JSONResponse(
                status_code=400,
                content={"error": f"Invalid scenario: {scenario}"},
            )

        allowed_workflow = {
            "build_only",
            "build_and_equilibrate",
            "run_md",
            "build_and_run_full",
        }
        if workflow not in allowed_workflow:
            return JSONResponse(status_code=400,
                                content={"error": f"Invalid workflow: {workflow}"})

        if workflow == "run_md" and not parent_job_id:
            return JSONResponse(
                status_code=400,
                content={"error": "run_md requires parent_job_id"},
            )

        if md_ns is not None:
            try:
                md_ns = float(md_ns)
            except ValueError:
                return JSONResponse(
                    status_code=400,
                    content={"error": "md_ns must be a number (float, in ns)"},
                )

            if md_ns <= 0:
                return JSONResponse(
                    status_code=400,
                    content={"error": "md_ns must be > 0"},
                )

        if workflow == "run_md":
            parent_dir = _md_job_dir(parent_job_id)
            if not parent_dir.exists():
                return JSONResponse(
                    status_code=404,
                    content={"error": "parent_job_id not found"},
                )
        if workflow == "build_and_run_full":
            # Treat it like build_and_equilibrate…
            # …but tell the container it must continue into production
            preset = "m3_popc_full"  # define a new preset in your preset list
            parent_job_id = None  # fresh system

        # Validate scenario requirements (Option A: orthosteric can be extracted)
        if scenario == "protein_plus_orthosteric" and orthosteric_ligand is None:
            # allowed: we'll try to extract after saving protein.pdb
            pass

        if scenario == "protein_plus_orthosteric_plus_allosteric":
            if allosteric_pose is None:
                return JSONResponse(
                    status_code=400,
                    content={"error": "Scenario requires allosteric_pose"},
                )
            # orthosteric_ligand may be None -> extracted from protein.pdb

        MD_RUNS_DIR.mkdir(parents=True, exist_ok=True)

        job_id = str(uuid.uuid4())
        job_dir = _md_job_dir(job_id)
        job_dir.mkdir(exist_ok=False)

        input_dir = job_dir / "input"
        out_dir = job_dir / "out"
        _safe_mkdir(input_dir)
        _safe_mkdir(out_dir)

        # -------------------------------------------------
        # Reuse system from ANY parent job (build / equil / md)
        # -------------------------------------------------
        start_gro = None
        start_cpt = None
        if workflow == "run_md":
            parent_dir = _md_job_dir(parent_job_id)
            parent_out = parent_dir / "out"

            # system.top is mandatory
            top_src = parent_out / "system.top"
            if not top_src.exists():
                return JSONResponse(
                    status_code=400,
                    content={"error": "Missing system.top in parent job"},
                )
            shutil.copy2(top_src, out_dir / "system.top")

            # Copy all topology includes
            for itp in parent_out.glob("*.itp"):
                shutil.copy2(itp, out_dir / itp.name)

            # Ensure Protein.itp is present (old system compatibility)
            protein0 = out_dir / "Protein_0.itp"
            protein = out_dir / "Protein.itp"

            if protein0.exists() and not protein.exists():
                shutil.copy2(protein0, protein)

            # First find BEST available starting structure (newest logic)
            candidate_gros = [
                "npt.gro", "md.gro", "nvt.gro", "em.gro", "system.gro",
            ]

            for fname in candidate_gros:
                f = parent_out / fname
                if f.exists():
                    start_gro = fname
                    shutil.copy2(f, out_dir / fname)
                    break

            if start_gro is None:
                return JSONResponse(
                    status_code=400,
                    content={"error": "No usable .gro found in parent job"},
                )

            # NEW: Ensure every MD job has a system.gro for analysis (RMSD, etc.)
            # ----------------------------------------------------------------------
            sys_src = parent_out / "system.gro"
            if sys_src.exists():
                shutil.copy2(sys_src, input_dir / "system.gro")
            else:
                # fallback: clone the chosen start_gro as system.gro
                fallback_src = parent_out / start_gro
                shutil.copy2(fallback_src, input_dir / "system.gro")

            # Second: find checkpoint (.cpt) — optional
            # Prefer MD checkpoint first
            candidate_cpts = [
                "md.cpt",
                "md_prev.cpt",  # optional: if you want to support fallback
                "npt.cpt",
                "nvt.cpt",
                "em.cpt",
            ]

            for fname in candidate_cpts:
                f = parent_out / fname
                if f.exists():
                    start_cpt = fname
                    shutil.copy2(f, out_dir / fname)
                    break

            if start_cpt is None:
                return JSONResponse(
                    status_code=400,
                    content={"error": "No usable .cpt found in parent job"},
                )

        if scenario != "ligand_water" and protein_pdb is None:
            return JSONResponse(
                status_code=400,
                content={"error": "Protein PDB required for this scenario"},
            )

        # --- Save protein ---
        pdb_path = None
        if protein_pdb is not None:
            pdb_path = input_dir / "protein.pdb"
            protein_pdb.file.seek(0)
            with open(pdb_path, "wb") as f:
                shutil.copyfileobj(protein_pdb.file, f)

        orth_path = None
        orth_extracted = False
        extract_summary = None

        needs_orth = scenario.startswith("protein_plus_orthosteric")

        # Option A: extract orthosteric from uploaded protein.pdb if not provided
        if needs_orth and orthosteric_ligand is None:
            original_with_ligands = input_dir / "protein_with_ligands.pdb"
            shutil.copyfile(pdb_path, original_with_ligands)

            extracted_orth = input_dir / "orthosteric_extracted.pdb"

            extract_summary = _extract_orthosteric_from_pdb(
                pdb_in=original_with_ligands,
                protein_out=pdb_path,  # overwrite protein.pdb as protein-only
                orth_out=extracted_orth,
            )

            if extract_summary.get("orth_hetatm_lines", 0) <= 0:
                return JSONResponse(
                    status_code=400,
                    content={
                        "error": (
                            "No orthosteric ligand found in uploaded protein PDB. "
                            "Upload orthosteric_ligand or provide a receptor PDB "
                            "containing the orthosteric as HETATM."
                        )
                    },
                )

            orth_path = extracted_orth
            orth_extracted = True

        # If orthosteric ligand explicitly uploaded, save it and prefer it
        if orthosteric_ligand is not None:
            orth_name = Path(orthosteric_ligand.filename).name
            orth_path = input_dir / orth_name
            orthosteric_ligand.file.seek(0)
            with open(orth_path, "wb") as f:
                shutil.copyfileobj(orthosteric_ligand.file, f)

        allo_path = None
        if allosteric_pose is not None:
            allo_name = Path(allosteric_pose.filename).name
            allo_path = input_dir / allo_name
            allosteric_pose.file.seek(0)
            with open(allo_path, "wb") as f:
                shutil.copyfileobj(allosteric_pose.file, f)

        created_at = datetime.now(
            ZoneInfo("Europe/Paris")
        ).isoformat(timespec="seconds")

        run_record = {
            "job_id": job_id,
            "created_at": created_at,
            "preset": preset,
            "scenario": scenario,
            "workflow": workflow,
            "parent_job_id": parent_job_id,
            "protein_filename": (
                protein_pdb.filename if protein_pdb is not None else None
            ),
            "orthosteric_filename": (
                orthosteric_ligand.filename
                if orthosteric_ligand is not None
                else ("orthosteric_extracted.pdb" if orth_extracted else None)
            ),
            "allosteric_pose_filename": (
                allosteric_pose.filename
                if allosteric_pose is not None else None
            ),
            "orthosteric_extracted": orth_extracted,
            "orthosteric_extract_summary": extract_summary,
            "ligand_case": ligand_case,
        }

        _write_json(job_dir / "run.json", run_record)

        # Initialize status + log so frontend can read immediately
        _write_json(job_dir / "status.json", {"status": "queued"})
        (job_dir / "log.txt").write_text("", encoding="utf-8")

        # Save params for container (future-proof for your md-runner)
        params = {
            "preset": preset,
            "scenario": scenario,
            "workflow": workflow,
            "parent_job_id": parent_job_id,
            "md_ns": md_ns,
            "start_gro": start_gro,
            "start_cpt": start_cpt,
            "orthosteric_smiles": orthosteric_smiles,
            "ligand_case": ligand_case,
            "gmx": {
                "nt": int(nt),
                "ntmpi": 1
            },
            "files": {
                "protein_pdb": "protein.pdb",
                "orthosteric_ligand": orth_path.name if orth_path else None,
                "allosteric_pose": allo_path.name if allo_path else None,
            },
            "comment": comment
        }
        _write_json(input_dir / "params.json", params)

        import traceback

        try:
            _launch_md_local(job_dir)
            _write_json(job_dir / "status.json", {"status": "running"})
        except Exception as e:
            tb = traceback.format_exc()
            # write to job log so you can see it from the UI
            (job_dir / "log.txt").write_text(
                "FAILED TO LAUNCH LOCAL MD JOB\n\n"
                f"Exception: {repr(e)}\n\n"
                f"Traceback:\n{tb}\n",
                encoding="utf-8",
            )
            _write_json(job_dir / "status.json", {"status": "error", "error": str(e)})
            return JSONResponse(
                status_code=500,
                content={"error": f"Failed to launch local MD job: {e}"},
            )

        # try:
        #     _launch_md_docker(job_dir)
        #     _write_json(job_dir / "status.json", {"status": "running"})
        # except Exception:
        #     # _launch_md_docker already wrote error + status
        #     return JSONResponse(status_code=500, content={"error": "Failed to launch Docker"})

        return JSONResponse({"job_id": job_id})

    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})
        return JSONResponse(status_code=500, content={"error": str(e)})


# ============================================================
# DIRECTORIES & CONSTANTS
# ============================================================

# Same directory your main.py uses
BASE_DIR = Path(__file__).resolve().parents[2]
MD_RUNS_DIR = BASE_DIR / "backend" / "data" / "md_runs"
MD_RUNS_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# GENERIC HELPERS
# ============================================================

def _safe_mkdir(p: Path):
    """Create directory if not exists."""
    p.mkdir(parents=True, exist_ok=True)


def _write_json(path: Path, obj):
    """Write JSON in UTF-8."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def _read_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _md_job_dir(job_id: str) -> Path:
    return MD_RUNS_DIR / job_id


# ============================================================
# LOG UTILITIES
# ============================================================

def read_log_chunk(job_id: str, offset: int) -> tuple[str, int]:
    """
    Return (chunk, new_offset).
    """
    job_dir = _md_job_dir(job_id)
    log_file = job_dir / "log.txt"

    if not log_file.exists():
        return "", offset

    text = log_file.read_text(encoding="utf-8")
    if offset >= len(text):
        return "", len(text)

    return text[offset:], len(text)


# ============================================================
# JOB STATUS UTILITIES
# ============================================================

def get_job_status(job_id: str) -> dict:
    job_dir = _md_job_dir(job_id)
    status_file = job_dir / "status.json"

    if not status_file.exists():
        return {"status": "unknown"}

    return _read_json(status_file)


def list_jobs(limit: Optional[int] = None) -> list[dict]:
    """
    Returns a list of { job_id, status, scenario, preset, ... } from run.json.
    Sorted by creation date (descending).
    """
    jobs = []

    for job_dir in MD_RUNS_DIR.iterdir():
        if not job_dir.is_dir():
            continue

        run_json = job_dir / "run.json"
        status_json = job_dir / "status.json"

        if not run_json.exists():
            continue

        try:
            run = _read_json(run_json)
        except Exception:
            continue

        status = "unknown"
        if status_json.exists():
            try:
                st = _read_json(status_json)
                status = st.get("status", "unknown")
            except Exception:
                pass

        jobs.append({
            "job_id": run["job_id"],
            "created_at": run.get("created_at"),
            "scenario": run.get("scenario"),
            "preset": run.get("preset"),
            "status": status,
        })

    # Sort by creation date descending
    jobs.sort(key=lambda j: j.get("created_at", ""), reverse=True)

    if limit and limit > 0:
        jobs = jobs[:limit]

    return jobs


# ============================================================
# ORTHOSTERIC EXTRACTION (placeholder)
# ============================================================

def _extract_orthosteric_from_pdb(
    pdb_in: Path,
    protein_out: Path,
    orth_out: Path,
) -> dict:
    """
    Split a PDB into:
      - protein_out: ATOM records (+ TER/END)
      - orth_out:    HETATM records for "ligand-like" residues
    Excludes common waters/ions from the ligand file.

    Returns a small summary dict.
    """
    exclude_resnames = {
        "HOH", "WAT", "DOD",
        "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU",
        "CO", "NI", "CD", "HG", "BR", "I",
        "SO4", "PO4",
        # add more if needed
    }

    atom_lines = []
    het_lines = []
    conect_lines = []

    for line in pdb_in.read_text(encoding="utf-8", errors="replace").splitlines(True):
        rec = line[:6].strip()

        if rec == "ATOM":
            atom_lines.append(line)
            continue

        if rec == "HETATM":
            # PDB fixed columns: resname at 18-20 (0-based 17:20)
            resname = line[17:20].strip().upper()
            if resname and resname not in exclude_resnames:
                het_lines.append(line)
            continue

        if rec == "CONECT":
            conect_lines.append(line)
            continue

        # Keep TER/END in protein file if present
        if rec in {"TER", "END", "ENDMDL"}:
            atom_lines.append(line)
            continue

    protein_out.write_text("".join(atom_lines) + "\n", encoding="utf-8")

    # Include CONECT if you want (harmless, sometimes helpful)
    orth_text = "".join(het_lines + conect_lines).strip() + "\n"
    orth_out.write_text(orth_text, encoding="utf-8")

    return {
        "protein_atoms_lines": len(atom_lines),
        "orth_hetatm_lines": len(het_lines),
        "conect_lines": len(conect_lines),
    }


# ============================================================
# MD LAUNCHING (placeholder)
# ============================================================

def _launch_md_local(job_dir: Path) -> None:

    print(">>> RUNTIME: _launch_md_local() CALLED")

    input_dir = job_dir / "input"
    params_path = input_dir / "params.json"

    if not params_path.exists():
        raise RuntimeError("params.json missing")

    params = json.loads(params_path.read_text(encoding="utf-8"))

    workflow = params.get("workflow", "full")

    module = PIPELINE_BY_WORKFLOW.get(workflow)
    if module is None:
        raise RuntimeError(f"Unknown workflow: {workflow}")

    job_id = job_dir.name


    input_dir = job_dir / "input"
    cmd = [
        sys.executable,
        "-m", module,
        "--workdir", str(job_dir / "out"),
        "--aa_pdb", str(input_dir / params["files"]["protein_pdb"]),
    ]

    nt = params.get("gmx", {}).get("nt", 1)
    cmd += ["--nt", str(nt)]

    orth = params["files"].get("orthosteric_ligand")
    if orth:
        cmd += ["--orthosteric_pdb", str(input_dir / orth)]

    orth_smiles = params.get("orthosteric_smiles")
    if orth_smiles:
        cmd += ["--orthosteric_smiles", orth_smiles]

    env = os.environ.copy()
    env["PYTHONPATH"] = str(PROJECT_ROOT)
    env["MARTINI_FF"] = "/Users/luiscamara/miniforge3/envs/drugdash/martini_v300"

    if "MARTINI_FF" in env:
        cmd += ["--martini_ff", env["MARTINI_FF"]]

    if workflow == "build_and_run_full":
        cmd.extend(["--do-em", "--do-nvt", "--do-npt"])

    preset = params.get("preset", "")

    if preset.endswith("_build"):
        pass
    elif preset.endswith("_em"):
        cmd.append("--do-em")
    elif preset.endswith("_eq"):
        cmd.extend(["--do-em", "--do-nvt", "--do-npt"])
    # elif workflow == "run_md":
    elif preset.endswith("_prod") or workflow == "run_md":
        cmd.append("--do-md")

    _write_json(job_dir / "status.json", {"status": "running"})

    def worker():
        try:
            code = _run_and_stream(
                job_dir,
                cmd,
                cwd=PROJECT_ROOT,
                env=env,
            )

            _write_json(job_dir / "status.json", {
                "status": "done" if code == 0 else "error",
                "exit_code": code,
            })

            _append_log(job_dir, f"\n[MD] Finished with exit code {code}\n")

        except Exception as e:
            _append_log(job_dir, f"\n[MD] EXCEPTION: {e}\n")
            _write_json(job_dir / "status.json", {
                "status": "error",
                "exception": str(e),
            })

    threading.Thread(target=worker, daemon=True).start()

def _append_log(job_dir: Path, line: str) -> None:
    with open(job_dir / "log.txt", "a", encoding="utf-8") as f:
        f.write(line)

from collections import deque
import os
import signal
import subprocess
from pathlib import Path

def _run_and_stream(
    job_dir: Path,
    cmd: list[str],
    cwd: Path,
    env=None,
) -> int:

    job_id = job_dir.name

    _append_log(job_dir, "\n$ " + " ".join(cmd) + "\n")

    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
        preexec_fn=os.setsid,  # start new process group
    )

    RUNNING_JOBS[job_id] = proc
    _write_json(job_dir / "pid.json", {"pid": proc.pid})

    # Rolling buffer of recent lines
    last_lines = deque(maxlen=300)

    # Diagnostics trackers
    last_energy = None
    last_temp = None
    last_pressure = None
    last_step_line = None

    fatal_detected = False
    fatal_trigger_line = None

    assert proc.stdout is not None

    for line in iter(proc.stdout.readline, ""):

        last_lines.append(line)

        # Track useful diagnostic info
        if "Epot" in line:
            last_energy = line
        if "Temperature" in line:
            last_temp = line
        if "Pressure" in line:
            last_pressure = line
        if "Step" in line:
            last_step_line = line

        lower = line.lower()

        # Only detect real fatal error
        if "fatal error" in lower:
            fatal_detected = True
            fatal_trigger_line = line

        # Write line to log
        _append_log(job_dir, line)

    # Process finished naturally
    code = proc.wait()

    if code != 0:
        _append_log(job_dir, "\n\n### MD TERMINATED WITH ERROR ###\n")
        _append_log(job_dir, f"Return code: {code}\n")

        if fatal_detected and fatal_trigger_line:
            _append_log(job_dir, "\nFatal trigger line:\n")
            _append_log(job_dir, fatal_trigger_line)

        if last_step_line:
            _append_log(job_dir, "\nLast step line:\n")
            _append_log(job_dir, last_step_line)

        if last_energy:
            _append_log(job_dir, "\nLast energy line:\n")
            _append_log(job_dir, last_energy)

        if last_temp:
            _append_log(job_dir, "\nLast temperature line:\n")
            _append_log(job_dir, last_temp)

        if last_pressure:
            _append_log(job_dir, "\nLast pressure line:\n")
            _append_log(job_dir, last_pressure)

        _append_log(job_dir, "\n--- Last 300 log lines before termination ---\n")
        _append_log(job_dir, "".join(last_lines))

    RUNNING_JOBS.pop(job_id, None)
    return code



# ============================================================
# JOB CONTROL
# ============================================================

import json
def stop_job(job_id: str):
    job_dir = _md_job_dir(job_id)
    pid_file = job_dir / "pid.json"

    if not pid_file.exists():
        return False

    try:
        data = json.loads(pid_file.read_text())
        pid = data["pid"]
    except Exception:
        return False

    try:
        # Kill the whole process group
        os.killpg(os.getpgid(pid), signal.SIGTERM)
    except ProcessLookupError:
        pass

    _write_json(job_dir / "status.json", {"status": "stopped"})
    return True



def delete_job(job_id: str) -> bool:
    job_dir = _md_job_dir(job_id)
    if not job_dir.exists():
        return False

    shutil.rmtree(job_dir)
    return True

def get_job(job_id: str) -> dict | None:
    """
    Return full run.json + computed status.
    Equivalent to old get_md_job endpoint.
    """
    job_dir = _md_job_dir(job_id)
    run_json = job_dir / "run.json"

    if not run_json.exists():
        return None

    try:
        data = _read_json(run_json)
    except Exception:
        return None

    # Determine status: read status.json if present
    status_json = job_dir / "status.json"
    status = "unknown"
    if status_json.exists():
        try:
            st = _read_json(status_json)
            status = st.get("status", "unknown")
        except Exception:
            pass

    data["status"] = status

    # --- Load comment from params.json ---
    params_json = job_dir / "input" / "params.json"
    if params_json.exists():
        try:
            params = _read_json(params_json)
            data["comment"] = params.get("comment", "")
        except Exception:
            data["comment"] = ""
    else:
        data["comment"] = ""

    return data

def read_log_window(job_id: str, offset: int, chunk_size: int):
    job_dir = _md_job_dir(job_id)
    log_path = job_dir / "log.txt"

    if not log_path.exists():
        return None, 0, 0

    file_size = log_path.stat().st_size

    offset = min(offset, file_size)

    with open(log_path, "rb") as f:
        f.seek(offset)
        data = f.read(chunk_size)

    new_offset = offset + len(data)

    return data, new_offset, file_size



def list_job_files(job_id: str):
    """
    Returns a list of {"name": path, "size": bytes} inside out/.
    Returns None if job does not exist.
    """
    job_dir = _md_job_dir(job_id)
    out_dir = job_dir / "out"

    if not out_dir.exists():
        return None

    files = []
    for p in out_dir.rglob("*"):
        if p.is_file():
            rel = p.relative_to(out_dir).as_posix()
            files.append({
                "name": rel,
                "size": p.stat().st_size,
            })

    files.sort(key=lambda x: x["name"])
    return files

import zipfile
def create_job_zip(job_id: str):
    """
    Creates (or overwrites) out.zip inside the job directory.
    Returns Path(zipfile) or None if job not found.
    """
    job_dir = _md_job_dir(job_id)
    out_dir = job_dir / "out"

    if not out_dir.exists():
        return None

    zip_path = job_dir / "out.zip"

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in out_dir.rglob("*"):
            if p.is_file():
                z.write(p, arcname=p.relative_to(out_dir).as_posix())

    return zip_path

def _launch_md_docker(job_dir: Path):
    """
    Placeholder for future Docker-based MD launching.
    Currently unimplemented.
    """
    raise NotImplementedError(
        "Docker-based MD launching is not implemented yet. "
        "Use _launch_md_local instead."
    )



def delete_job(job_id: str) -> bool:
    job_dir = _md_job_dir(job_id)
    if not job_dir.exists():
        return False

    shutil.rmtree(job_dir)
    return True

from backend.pipelines.pipeline import main as run_pipeline

def run_cgmd_setup_job(request):
    print("REQUEST BODY:", request.model_dump())

    args = [
        "--workdir", request.workdir,
        "--aa_pdb", request.aa_pdb,
        "--name", request.name,
        "--martini_ff", request.martini_ff,
        "--nt", str(request.nt),
    ]

    # Optional orthosteric ligand PDB
    if getattr(request, "orthosteric_pdb", None) is not None:
        args += ["--orthosteric_pdb", request.orthosteric_pdb]

    # Pipeline stages
    if request.do_em:
        args.append("--do-em")
    if request.do_nvt:
        args.append("--do-nvt")
    if request.do_npt:
        args.append("--do-npt")
    if request.do_md:
        args.append("--do-md")

    run_pipeline(args)
    return {"status": "ok"}








