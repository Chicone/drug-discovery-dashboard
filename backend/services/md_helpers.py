from pathlib import Path
import json
from typing import Optional


# Where MD runs are stored
MD_RUNS_DIR = Path(__file__).resolve().parents[1] / "data" / "md_runs"


def md_job_dir(job_id: str) -> Path:
    return MD_RUNS_DIR / job_id


def safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_json(p: Path, obj: dict) -> None:
    p.write_text(json.dumps(obj, indent=2), encoding="utf-8")


def read_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def md_status_from_dir(job_dir: Path) -> str:
    status_file = job_dir / "status.json"
    if status_file.exists():
        try:
            return read_json(status_file).get("status", "unknown")
        except Exception:
            return "unknown"

    # Infer status from marker files
    if (job_dir / "out" / "done.marker").exists():
        return "done"
    if (job_dir / "out" / "error.marker").exists():
        return "error"

    return "running"
