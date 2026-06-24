from pathlib import Path

# Order matters: newest / best first
RESTART_GRO_CANDIDATES = [
    "npt.gro",
    "md.gro",       # finished MD
    "nvt.gro",
    "em.gro",
    "system.gro",
]

def find_restart_structure(job_dir: Path) -> Path:
    out_dir = job_dir / "out"
    if not out_dir.exists():
        raise FileNotFoundError("Parent job has no out/ directory")

    for name in RESTART_GRO_CANDIDATES:
        p = out_dir / name
        if p.exists():
            return p

    raise FileNotFoundError(
        "No restartable GRO found. Expected one of: "
        + ", ".join(RESTART_GRO_CANDIDATES)
    )
