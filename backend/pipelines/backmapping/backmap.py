from pathlib import Path
import subprocess

from .paths import (
    extracted_frame_path,
    backmapped_structure_path,
    aa_prep_dir,
)
from .generate_mapping import generate_backward_mapping

# Locate Backward inside the project
BACKWARD = (
    Path(__file__).resolve()
    .parents[3] / "tools" / "backward" / "backward.py"
)

from .paths import cg_out_dir
def run_backmapping(job_id: str) -> Path:
    """
    Convert a Martini frame to an atomistic structure
    using Backward.
    """

    itp = cg_out_dir(job_id) / "Orthosteric.itp"

    if not itp.exists():
        raise RuntimeError(f"Ligand topology not found: {itp}")

    # Generate mapping automatically
    residue = detect_residue_name(itp)

    mapping_file = (
        Path(__file__).resolve().parents[3]
        / "tools"
        / "backward"
        / "Mapping"
        / f"{residue.lower()}.charmm36.map"
    )

    if not mapping_file.exists():
        generate_backward_mapping(itp, mapping_file)

    frame = extracted_frame_path(job_id)

    if not frame.exists():
        raise RuntimeError(f"Frame not found: {frame}")

    out_dir = aa_prep_dir(job_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    output = backmapped_structure_path(job_id)

    cmd = [
        "python",
        str(BACKWARD),
        "-f", str(frame),
        "-o", str(output),
        "-to", "charmm36",
    ]

    print("Running Backward:")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

    return output

def detect_residue_name(itp_file):
    for line in itp_file.read_text().splitlines():

        parts = line.split()

        if len(parts) >= 5 and parts[0].isdigit():
            return parts[3]   # residue column

    raise RuntimeError("Could not detect residue name")