from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from .backmapping.extract_frame import extract_frame
from .backmapping.backmap import run_backmapping
from .backmapping.paths import extracted_frame_path
from .backmapping.prepare_aa_system import (
    prepare_aa_system,
    extract_ligand_from_complex_pdb,
    prepare_ligand_topology,
    integrate_ligand_into_aa_complex,
    prepare_aa_equilibration_system,
)

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

    enable_ligand = False

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
            raise RuntimeError(
                f"Cannot infer matching TPR for trajectory: {xtc_name}"
            )

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

    protein_pdb = Path(params["protein_pdb"])
    complex_pdb = Path(params["complex_pdb"])

    extracted_protein_pdb = out_dir / "aa_build" / "protein_from_complex.pdb"
    extracted_ligand_pdb = out_dir / "aa_build" / "ligand_from_complex.pdb"

    extracted_protein_pdb.parent.mkdir(parents=True, exist_ok=True)

    extract_info = extract_ligand_from_complex_pdb(
        complex_pdb=complex_pdb,
        protein_out=extracted_protein_pdb,
        ligand_out=extracted_ligand_pdb,
        ligand_resname=params.get("ligand_resname"),
    )

    result["extracted_protein_pdb"] = str(extracted_protein_pdb)
    result["extracted_ligand_pdb"] = str(extracted_ligand_pdb)
    result["ligand_extract_info"] = extract_info

    aa_info = prepare_aa_system(
        aa_job_out_dir=out_dir,
        protein_pdb=extracted_protein_pdb,
        forcefield=params.get("forcefield", "charmm36-jul2022"),
        water_model=params.get("water_model", "tip3p"),
    )

    result["aa_build_dir"] = str(aa_info["aa_build_dir"])
    result["protein_processed_gro"] = str(aa_info["protein_processed_gro"])
    result["topol_top"] = str(aa_info["topol_top"])
    result["posre_itp"] = str(aa_info["posre_itp"])

    ligand_smiles = Path(params["ligand_smiles"])
    ligand_itp = Path(params["ligand_itp"]) if params.get("ligand_itp") else None
    ligand_rtf = Path(params["ligand_rtf"]) if params.get("ligand_rtf") else None
    ligand_g_rtf = (
        Path(params["ligand_g_rtf"]) if params.get("ligand_g_rtf") else None
    )
    ligand_prm = Path(params["ligand_prm"]) if params.get("ligand_prm") else None
    ligand_coord_template = (
        Path(params["ligand_coord_template"])
        if params.get("ligand_coord_template") else None
    )
    ligand_toppar = Path(params["ligand_toppar"]) if params.get("ligand_toppar") else None

    print(f"[aa] ligand_itp = {ligand_itp}")
    print(f"[aa] ligand_coord_template = {ligand_coord_template}")
    print(
        f"[aa] ligand_coord_template exists = "
        f"{ligand_coord_template.exists() if ligand_coord_template else None}"
    )

    ligand_info = prepare_ligand_topology(
        aa_job_out_dir=out_dir,
        ligand_pdb=extracted_ligand_pdb,
        ligand_smiles=ligand_smiles,
        ligand_name=params.get("ligand_resname", "UNL"),
        ligand_itp=ligand_itp,
        ligand_rtf=ligand_rtf,
        ligand_g_rtf=ligand_g_rtf,
        ligand_prm=ligand_prm,
        ligand_coord_template=ligand_coord_template,
        ligand_toppar=ligand_toppar,
    )


    result["ligand_input_pdb"] = str(ligand_info["ligand_input_pdb"])
    result["ligand_input_smi"] = str(ligand_info["ligand_input_smi"])
    result["ligand_from_smiles_sdf"] = str(ligand_info["ligand_from_smiles_sdf"])
    result["ligand_parametrization_ready"] = (
        ligand_info["ligand_parametrization_ready"]
    )

    if ligand_info["ligand_itp"] is not None:
        result["ligand_itp"] = str(ligand_info["ligand_itp"])

    if ligand_info["ligand_rtf"] is not None:
        result["ligand_rtf"] = str(ligand_info["ligand_rtf"])

    if ligand_info["ligand_g_rtf"] is not None:
        result["ligand_g_rtf"] = str(ligand_info["ligand_g_rtf"])

    if ligand_info["ligand_prm"] is not None:
        result["ligand_prm"] = str(ligand_info["ligand_prm"])

    if ligand_info.get("ligand_toppar") is not None:
        result["ligand_toppar"] = str(ligand_info["ligand_toppar"])

    if ligand_info["ligand_coord_template"] is not None:
        result["ligand_coord_template"] = str(
            ligand_info["ligand_coord_template"]
        )

    if enable_ligand and ligand_info["ligand_itp"] is not None:
        if ligand_info["ligand_coord_template"] is None:
            raise RuntimeError(
                "Ligand topology exists, but no matching "
                "ligand_coord_template was provided. Please provide a "
                "46-atom PDB/MOL2/GRO that matches ligand.itp exactly."
            )

        integration_info = integrate_ligand_into_aa_complex(
            aa_job_out_dir=out_dir,
            protein_gro=Path(aa_info["protein_processed_gro"]),
            topol_top=Path(aa_info["topol_top"]),
            ligand_coord_template=Path(ligand_info["ligand_coord_template"]),
            ligand_itp=Path(ligand_info["ligand_itp"]),
            ligand_toppar=(
                Path(ligand_info["ligand_toppar"])
                if ligand_info.get("ligand_toppar") is not None else None
            ),
            ligand_resname=params.get("ligand_resname", "UNL"),
        )

        result["complex_gro"] = str(integration_info["complex_gro"])
        result["ligand_gro"] = str(integration_info["ligand_gro"])
        result["topol_top"] = str(integration_info["topol_top"])

        complex_gro_for_eq = Path(integration_info["complex_gro"])
        topol_top_for_eq = Path(integration_info["topol_top"])

    else:
        print("[aa] Skipping ligand integration (protein-only mode)")
        result["ligand_gro"] = None
        complex_gro_for_eq = Path(aa_info["protein_processed_gro"])
        topol_top_for_eq = Path(aa_info["topol_top"])
        result["complex_gro"] = str(complex_gro_for_eq)
        result["topol_top"] = str(topol_top_for_eq)

    eq_info = prepare_aa_equilibration_system(
        aa_job_out_dir=out_dir,
        complex_gro=complex_gro_for_eq,
        topol_top=topol_top_for_eq,
        run_em=run_em,
        run_nvt=run_nvt,
        run_npt=run_npt,
    )

    result["boxed_gro"] = str(eq_info["boxed_gro"])
    result["solvated_gro"] = str(eq_info["solvated_gro"])
    result["ionized_gro"] = str(eq_info["ionized_gro"])
    result["em_gro"] = str(eq_info["em_gro"])
    result["nvt_gro"] = str(eq_info["nvt_gro"])
    result["npt_gro"] = str(eq_info["npt_gro"])

    print("[aa] Reconstruction finished")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()