from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from .backmapping.extract_frame import extract_frame
from .backmapping.backmap import run_backmapping
from .backmapping.paths import extracted_frame_path


def backmap_frame_to_aa(
    job_id: str,
    time_ps: float | None = None,
    input_gro: Path | None = None,
    run_em: bool = True,
    run_nvt: bool = True,
    run_npt: bool = True,
    aa_job_out_dir: Path | None = None,
    traj_path: Path | None = None,
    tpr_path: Path | None = None,
) -> dict:
    """
    Minimal AA reconstruction wrapper.

    Current scope:
    - if input_gro is provided, use it directly as frame.gro
    - else extract a frame from the CG trajectory
    - run Backward to produce backmapped.gro

    Not implemented yet:
    - AA topology preparation
    - EM / NVT / NPT
    """

    if input_gro is None and time_ps is None:
        raise ValueError("Either time_ps or input_gro must be provided")

    # ---------------------------------------------------------
    # 1. Prepare frame.gro
    # ---------------------------------------------------------
    if input_gro is not None:
        dst = extracted_frame_path(job_id)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_gro, dst)
        frame_path = dst
        print(f"[aa] Using static CG structure as frame: {frame_path}")
    else:
        if traj_path is None:
            raise ValueError(
                "traj_path must be provided when extracting a trajectory frame"
            )

        extract_frame(job_id, float(time_ps), traj_path=traj_path, tpr_path=tpr_path,)
        frame_path = extracted_frame_path(job_id)
        print(f"[aa] Extracted trajectory frame from {traj_path}: {frame_path}")

    if not frame_path.exists():
        raise FileNotFoundError(f"Frame not found: {frame_path}")

    # ---------------------------------------------------------
    # 2. Run Backward
    # ---------------------------------------------------------
    backmapped_path = run_backmapping(job_id)

    if not backmapped_path.exists():
        raise FileNotFoundError(
            f"Backmapping completed but output not found: {backmapped_path}"
        )

    # ---------------------------------------------------------
    # 3. Optionally copy final outputs into this AA job's out/
    # ---------------------------------------------------------
    copied_backmapped = None
    copied_frame = None

    if aa_job_out_dir is not None:
        aa_job_out_dir.mkdir(parents=True, exist_ok=True)

        copied_frame = aa_job_out_dir / "frame.gro"
        copied_backmapped = aa_job_out_dir / "backmapped.gro"

        # Avoid SameFileError if source and destination are already the same
        if frame_path.resolve() != copied_frame.resolve():
            shutil.copy2(frame_path, copied_frame)
            print(f"[aa] Copied frame.gro -> {copied_frame}")
        else:
            copied_frame = frame_path
            print(f"[aa] frame.gro already in output dir: {copied_frame}")

        if backmapped_path.resolve() != copied_backmapped.resolve():
            shutil.copy2(backmapped_path, copied_backmapped)
            print(f"[aa] Copied backmapped.gro -> {copied_backmapped}")
        else:
            copied_backmapped = backmapped_path
            print(f"[aa] backmapped.gro already in output dir: {copied_backmapped}")

    # ---------------------------------------------------------
    # 4. Placeholder for future AA relaxation
    # ---------------------------------------------------------
    if run_em or run_nvt or run_npt:
        print(
            "[aa] NOTE: run_em/run_nvt/run_npt were requested, "
            "but AA relaxation is not implemented yet in aa_pipeline.py"
        )

    return {
        "job_id": job_id,
        "frame_gro": str(copied_frame or frame_path),
        "backmapped_gro": str(copied_backmapped or backmapped_path),
        "run_em_requested": run_em,
        "run_nvt_requested": run_nvt,
        "run_npt_requested": run_npt,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.workdir).resolve()
    job_dir = out_dir.parent
    input_dir = job_dir / "input"

    params_path = input_dir / "params.json"
    if not params_path.exists():
        raise FileNotFoundError(f"Missing params.json: {params_path}")

    params = json.loads(params_path.read_text())

    parent_job_id = params["parent_job_id"]
    time_ps = params.get("time_ps", 0.0)
    run_em = params.get("run_em", True)
    run_nvt = params.get("run_nvt", True)
    run_npt = params.get("run_npt", True)

    source_xtc = params.get("source_xtc")
    source_static_gro = params.get("source_static_gro")

    print(f"[aa] parent_job_id = {parent_job_id}")
    print(f"[aa] source_xtc = {source_xtc}")
    print(f"[aa] source_static_gro = {source_static_gro}")

    if source_xtc:
        traj_path = out_dir / source_xtc
        if not traj_path.exists():
            raise FileNotFoundError(f"Trajectory source not found: {traj_path}")

        xtc_name = traj_path.name
        if xtc_name.startswith("nvt"):
            tpr_name = "nvt.tpr"
        elif xtc_name.startswith("npt"):
            tpr_name = "npt.tpr"
        elif xtc_name.startswith("md"):
            tpr_name = "md.tpr"
        else:
            raise RuntimeError(f"Cannot infer matching TPR for trajectory: {xtc_name}")

        tpr_path = out_dir / tpr_name
        if not tpr_path.exists():
            raise FileNotFoundError(f"Matching TPR not found: {tpr_path}")

        result = backmap_frame_to_aa(
            job_id=parent_job_id,
            time_ps=float(time_ps),
            input_gro=None,
            run_em=run_em,
            run_nvt=run_nvt,
            run_npt=run_npt,
            aa_job_out_dir=out_dir,
            traj_path=traj_path,
            tpr_path=tpr_path,
        )
    else:
        if not source_static_gro:
            raise RuntimeError(
                "No trajectory source and no static structure source were provided"
            )

        static_gro = out_dir / source_static_gro
        if not static_gro.exists():
            raise FileNotFoundError(f"Static source gro not found: {static_gro}")

        result = backmap_frame_to_aa(
            job_id=parent_job_id,
            time_ps=None,
            input_gro=static_gro,
            run_em=run_em,
            run_nvt=run_nvt,
            run_npt=run_npt,
            aa_job_out_dir=out_dir,
            traj_path=None,
        )

    print("[aa] Reconstruction finished")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()