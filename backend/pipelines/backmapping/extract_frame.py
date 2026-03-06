import subprocess
from pathlib import Path

from .paths import (
    cg_xtc,
    cg_tpr,
    extracted_frame_path,
    aa_prep_dir,
)


def extract_frame(job_id: str, time_ps: float) -> Path:
    """
    Extract a frame from a Martini trajectory using gmx trjconv.
    """

    traj = cg_xtc(job_id)
    tpr = cg_tpr(job_id)

    if not traj.exists():
        raise RuntimeError(f"Trajectory not found: {traj}")

    if not tpr.exists():
        raise RuntimeError(f"TPR not found: {tpr}")

    out_dir = aa_prep_dir(job_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    frame_out = extracted_frame_path(job_id)

    cmd = [
        "gmx",
        "trjconv",
        "-f", str(traj),
        "-s", str(tpr),
        "-o", str(frame_out),
        "-dump", str(time_ps),
    ]

    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        text=True,
    )

    # select "System" group automatically
    proc.communicate("0\n")

    if proc.returncode != 0:
        raise RuntimeError("gmx trjconv failed")

    return frame_out
