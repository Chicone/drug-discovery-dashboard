from pathlib import Path

# project root
BASE_DIR = Path(__file__).resolve().parents[3]

DATA_DIR = BASE_DIR / "backend" / "data"
MD_RUNS_DIR = DATA_DIR / "md_runs"


def cg_run_dir(job_id: str) -> Path:
    """
    Directory of an existing Martini MD run.
    """
    return MD_RUNS_DIR / job_id


def cg_out_dir(job_id: str) -> Path:
    """
    'out' directory containing trajectory files.
    """
    return cg_run_dir(job_id) / "out"


def aa_prep_dir(job_id: str) -> Path:
    """
    Directory where backmapping outputs will be written.
    """
    return MD_RUNS_DIR / f"{job_id}_aa_prep"


def extracted_frame_path(job_id: str) -> Path:
    return aa_prep_dir(job_id) / "frame.gro"


def backmapped_structure_path(job_id: str) -> Path:
    return aa_prep_dir(job_id) / "backmapped.gro"


def cg_xtc(job_id: str) -> Path:
    print(cg_out_dir(job_id) / "md.part0001.xtc")
    return cg_out_dir(job_id) / "md.part0001.xtc"


def cg_tpr(job_id: str) -> Path:
    return cg_out_dir(job_id) / "md.tpr"