import subprocess
from pathlib import Path

from .paths import (
    cg_xtc,
    cg_tpr,
    extracted_frame_path,
    aa_prep_dir,
)

import subprocess
from pathlib import Path

from .paths import cg_xtc, cg_tpr, extracted_frame_path


def extract_frame(
    job_id: str,
    time_ps: float,
    traj_path: Path | None = None,
    tpr_path: Path | None = None,
) -> Path:
    """
    Extract one frame from a CG trajectory at time_ps.

    If traj_path or tpr_path are given, use them directly.
    Otherwise fall back to the default path helpers.
    """

    traj = traj_path if traj_path is not None else cg_xtc(job_id)
    tpr = tpr_path if tpr_path is not None else cg_tpr(job_id)
    out = extracted_frame_path(job_id)

    if not traj.exists():
        raise RuntimeError(f"Trajectory not found: {traj}")

    if not tpr.exists():
        raise RuntimeError(f"TPR not found: {tpr}")

    out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "gmx", "trjconv",
        "-f", str(traj),
        "-s", str(tpr),
        "-o", str(out),
        "-dump", str(time_ps),
    ]

    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    stdout, _ = proc.communicate("0\n")

    if proc.returncode != 0:
        raise RuntimeError(
            f"gmx trjconv failed with code {proc.returncode}\n{stdout}"
        )

    if not out.exists():
        raise RuntimeError(f"Frame extraction failed, output not created: {out}")

    return out
