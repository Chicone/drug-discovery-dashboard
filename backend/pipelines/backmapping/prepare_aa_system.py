from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from rdkit import Chem

THIS_DIR = Path(__file__).resolve().parent
AA_MDP_DIR = THIS_DIR / "aa_mdp_templates"

IONS_MDP = AA_MDP_DIR / "ions.mdp"
EM_MDP = AA_MDP_DIR / "em.mdp"
NVT_MDP = AA_MDP_DIR / "nvt.mdp"
NPT_MDP = AA_MDP_DIR / "npt.mdp"

def _run(
    cmd: list[str],
    cwd: Path | None = None,
    stdin_text: str | None = None,
) -> None:
    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        input=stdin_text,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}"
        )


def prepare_aa_system(
    aa_job_out_dir: Path,
    protein_pdb: Path,
    forcefield: str = "charmm36-jul2022",
    water_model: str = "tip3p",
) -> dict[str, Path]:
    aa_build_dir = aa_job_out_dir / "aa_build"
    aa_build_dir.mkdir(parents=True, exist_ok=True)

    protein_input = aa_build_dir / "protein_input.pdb"
    shutil.copy2(protein_pdb, protein_input)

    protein_processed_gro = aa_build_dir / "protein_processed.gro"
    topol_top = aa_build_dir / "topol.top"
    posre_itp = aa_build_dir / "posre.itp"

    cmd = [
        "gmx",
        "pdb2gmx",
        "-f", str(protein_input),
        "-o", str(protein_processed_gro),
        "-p", str(topol_top),
        "-i", str(posre_itp),
        "-ff", forcefield,
        "-water", water_model,
        "-ignh",
        "-ter",
    ]

    _run(cmd, cwd=aa_build_dir, stdin_text="1\n0\n")

    return {
        "aa_build_dir": aa_build_dir,
        "protein_input": protein_input,
        "protein_processed_gro": protein_processed_gro,
        "topol_top": topol_top,
        "posre_itp": posre_itp,
    }


def extract_ligand_from_complex_pdb(
    complex_pdb: Path,
    protein_out: Path,
    ligand_out: Path,
    ligand_resname: str | None = None,
) -> dict:
    exclude_resnames = {
        "HOH", "WAT", "DOD",
        "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU",
        "CO", "NI", "CD", "HG", "BR", "I",
        "SO4", "PO4",
    }

    atom_lines = []
    het_lines = []

    for line in complex_pdb.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines(True):
        rec = line[:6].strip()

        if rec == "ATOM":
            atom_lines.append(line)
            continue

        if rec == "HETATM":
            resname = line[17:20].strip().upper()

            if not resname or resname in exclude_resnames:
                continue

            if ligand_resname is None or resname == ligand_resname.upper():
                het_lines.append(line)

            continue

        if rec in {"TER", "END", "ENDMDL"}:
            atom_lines.append(line)
            continue

    protein_out.write_text("".join(atom_lines) + "\n", encoding="utf-8")
    ligand_out.write_text("".join(het_lines).strip() + "\n", encoding="utf-8")

    if len(het_lines) == 0:
        raise RuntimeError(
            f"No ligand atoms extracted from {complex_pdb}. "
            f"Check ligand_resname."
        )

    return {
        "protein_atoms_lines": len(atom_lines),
        "ligand_hetatm_lines": len(het_lines),
    }


def prepare_ligand_topology(
    aa_job_out_dir: Path,
    ligand_pdb: Path,
    ligand_smiles: Path,
    ligand_name: str = "UNL",
    ligand_itp: Path | None = None,
    ligand_rtf: Path | None = None,
    ligand_g_rtf: Path | None = None,
    ligand_prm: Path | None = None,
    ligand_coord_template: Path | None = None,
    ligand_toppar: Path | None = None,
) -> dict[str, Path | bool | None]:
    """
    Prepare ligand inputs for the AA reconstruction stage.

    Current scope:
    - copy the extracted docked ligand coordinates
    - copy the ligand SMILES
    - generate a simple SDF from the SMILES for bookkeeping
    - prefer a GROMACS ligand .itp if available
    - otherwise stage CHARMM-GUI files (RTF / GRTF / PRM)

    Notes:
    - If ligand_itp is present, that is the preferred route for later
      integration into topol.top.
    """

    aa_build_dir = aa_job_out_dir / "aa_build"
    aa_build_dir.mkdir(parents=True, exist_ok=True)

    ligand_input_pdb = aa_build_dir / "ligand_input.pdb"
    ligand_input_smi = aa_build_dir / "ligand_input.smi"
    ligand_from_smiles_sdf = aa_build_dir / "ligand_from_smiles.sdf"

    shutil.copy2(ligand_pdb, ligand_input_pdb)
    shutil.copy2(ligand_smiles, ligand_input_smi)

    smiles_line = ligand_smiles.read_text(encoding="utf-8").strip()
    if not smiles_line:
        raise RuntimeError(f"Empty ligand SMILES file: {ligand_smiles}")

    smiles = smiles_line.split()[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"Could not parse ligand SMILES: {smiles}")

    mol.SetProp("_Name", ligand_name)
    writer = Chem.SDWriter(str(ligand_from_smiles_sdf))
    writer.write(mol)
    writer.close()

    copied_ligand_itp = None
    copied_ligand_rtf = None
    copied_ligand_g_rtf = None
    copied_ligand_prm = None
    copied_ligand_toppar = None
    parametrization_ready = False

    if ligand_itp is not None:
        copied_ligand_itp = aa_build_dir / "ligand.itp"
        shutil.copy2(ligand_itp, copied_ligand_itp)
        parametrization_ready = True

    if ligand_rtf is not None:
        copied_ligand_rtf = aa_build_dir / "ligand.rtf"
        shutil.copy2(ligand_rtf, copied_ligand_rtf)

    if ligand_g_rtf is not None:
        copied_ligand_g_rtf = aa_build_dir / "ligand_g.rtf"
        shutil.copy2(ligand_g_rtf, copied_ligand_g_rtf)

    if ligand_prm is not None:
        copied_ligand_prm = aa_build_dir / "ligand.prm"
        shutil.copy2(ligand_prm, copied_ligand_prm)

    if ligand_toppar is not None:
        copied_ligand_toppar = aa_build_dir / "toppar.str"
        shutil.copy2(ligand_toppar, copied_ligand_toppar)

    copied_ligand_coord_template = None
    if ligand_coord_template is not None:
        suffix = ligand_coord_template.suffix.lower()
        copied_ligand_coord_template = aa_build_dir / f"ligand_template{suffix}"
        shutil.copy2(ligand_coord_template, copied_ligand_coord_template)

    return {
        "ligand_input_pdb": ligand_input_pdb,
        "ligand_input_smi": ligand_input_smi,
        "ligand_from_smiles_sdf": ligand_from_smiles_sdf,
        "ligand_itp": copied_ligand_itp,
        "ligand_rtf": copied_ligand_rtf,
        "ligand_g_rtf": copied_ligand_g_rtf,
        "ligand_prm": copied_ligand_prm,
        "ligand_parametrization_ready": parametrization_ready,
        "ligand_coord_template": copied_ligand_coord_template,
    }

def _insert_ligand_include_into_topol(
    topol_top: Path,
    ligand_itp_name: str = "ligand.itp",
    ligand_toppar_name: str | None = None,
) -> None:
    """
    Ensure topol.top uses a local posre.itp include and includes
    ligand parameter files before ligand.itp, all before [ system ].
    """
    lines = topol_top.read_text(encoding="utf-8").splitlines()

    normalized = []
    found_posre = False
    found_ligand = False
    found_toppar = False

    for line in lines:
        stripped = line.strip()

        if "posre.itp" in stripped and stripped.startswith("#include"):
            normalized.append('#include "posre.itp"')
            found_posre = True
            continue

        if stripped == f'#include "{ligand_itp_name}"':
            found_ligand = True

        if ligand_toppar_name is not None and stripped == f'#include "{ligand_toppar_name}"':
            found_toppar = True

        normalized.append(line)

    lines = normalized

    # Ensure local posre.itp include exists
    if not found_posre:
        inserted = False
        for i, line in enumerate(lines):
            if line.strip() == "#ifdef POSRES":
                lines.insert(i + 1, '#include "posre.itp"')
                inserted = True
                break

        if not inserted:
            for i, line in enumerate(lines):
                if line.strip().lower() == "[ system ]":
                    lines.insert(i, '#include "posre.itp"')
                    inserted = True
                    break

    # Build include block to insert before [ system ]
    include_block = []

    if ligand_toppar_name is not None and not found_toppar:
        include_block.append(f'#include "{ligand_toppar_name}"')

    if not found_ligand:
        include_block.append(f'#include "{ligand_itp_name}"')

    if include_block:
        inserted = False
        for i, line in enumerate(lines):
            if line.strip().lower() == "[ system ]":
                for j, inc in enumerate(include_block):
                    lines.insert(i + j, inc)
                inserted = True
                break

        if not inserted:
            raise RuntimeError(
                f"Could not find [ system ] block in {topol_top} "
                f"to insert ligand include(s)"
            )

    topol_top.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _append_ligand_to_molecules_block(
    topol_top: Path,
    ligand_resname: str = "UNL",
) -> None:
    """
    Ensure the ligand appears once in the [ molecules ] block.
    """
    lines = topol_top.read_text(encoding="utf-8").splitlines()

    molecules_idx = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            molecules_idx = i
            break

    if molecules_idx is None:
        raise RuntimeError(f"[ molecules ] block not found in {topol_top}")

    # Check if already present
    for line in lines[molecules_idx + 1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.split()[0] == ligand_resname:
            topol_top.write_text("\n".join(lines) + "\n", encoding="utf-8")
            return

    # Append at end of molecules block
    lines.append(f"{ligand_resname:<16s} 1")
    topol_top.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _pdb_to_gro(
    input_pdb: Path,
    output_gro: Path,
) -> None:
    cmd = [
        "gmx",
        "editconf",
        "-f", str(input_pdb),
        "-o", str(output_gro),
    ]
    _run(cmd)


def _merge_gro_files(
    protein_gro: Path,
    ligand_gro: Path,
    output_gro: Path,
) -> None:
    protein_lines = protein_gro.read_text(encoding="utf-8").splitlines()
    ligand_lines = ligand_gro.read_text(encoding="utf-8").splitlines()

    if len(protein_lines) < 3 or len(ligand_lines) < 3:
        raise RuntimeError("Invalid .gro file(s) for merge")

    protein_title = protein_lines[0]
    protein_n = int(protein_lines[1].strip())
    ligand_n = int(ligand_lines[1].strip())

    protein_atoms = protein_lines[2:-1]
    ligand_atoms = ligand_lines[2:-1]
    box_line = protein_lines[-1]

    merged = []
    merged.append(f"{protein_title} + ligand")
    merged.append(str(protein_n + ligand_n))
    merged.extend(protein_atoms)
    merged.extend(ligand_atoms)
    merged.append(box_line)

    output_gro.write_text("\n".join(merged) + "\n", encoding="utf-8")

def integrate_ligand_into_aa_complex(
    aa_job_out_dir: Path,
    protein_gro: Path,
    topol_top: Path,
    ligand_coord_template: Path,
    ligand_itp: Path,
    ligand_toppar: Path | None = None,
    ligand_resname: str = "UNL",
) -> dict[str, Path]:
    """
    Integrate a ligand into the AA protein topology/coordinates.

    Requires a ligand coordinate template that matches ligand.itp exactly.
    """

    aa_build_dir = aa_job_out_dir / "aa_build"
    aa_build_dir.mkdir(parents=True, exist_ok=True)

    staged_ligand_itp = aa_build_dir / "ligand.itp"
    if ligand_itp.resolve() != staged_ligand_itp.resolve():
        shutil.copy2(ligand_itp, staged_ligand_itp)

    staged_ligand_toppar = None
    if ligand_toppar is not None:
        staged_ligand_toppar = aa_build_dir / "toppar.str"
        if ligand_toppar.resolve() != staged_ligand_toppar.resolve():
            shutil.copy2(ligand_toppar, staged_ligand_toppar)

    ligand_gro = aa_build_dir / "ligand.gro"
    complex_gro = aa_build_dir / "complex.gro"

    _insert_ligand_include_into_topol(
        topol_top,
        ligand_itp_name="ligand.itp",
        ligand_toppar_name=(
            staged_ligand_toppar.name if staged_ligand_toppar is not None else None
        ),
    )
    _append_ligand_to_molecules_block(topol_top, ligand_resname=ligand_resname)

    _ligand_structure_to_gro(ligand_coord_template, ligand_gro)

    def _count_gro_atoms(gro_path: Path) -> int:
        lines = gro_path.read_text().splitlines()
        if len(lines) < 2:
            raise RuntimeError(f"Invalid GRO file: {gro_path}")
        return int(lines[1].strip())

    def _count_itp_atoms(itp_path: Path) -> int:
        in_atoms = False
        count = 0

        for raw_line in itp_path.read_text().splitlines():
            line = raw_line.strip()

            if not line or line.startswith(";"):
                continue

            if line.startswith("["):
                in_atoms = (line == "[ atoms ]")
                continue

            if in_atoms:
                count += 1

        if count == 0:
            raise RuntimeError(f"No [ atoms ] entries found in {itp_path}")

        return count

    gro_atoms = _count_gro_atoms(ligand_gro)
    itp_atoms = _count_itp_atoms(staged_ligand_itp)

    if gro_atoms != itp_atoms:
        raise RuntimeError(
            f"Ligand coordinate/topology mismatch: ligand.gro has {gro_atoms} atoms "
            f"but ligand.itp expects {itp_atoms} atoms. "
            f"Template used: {ligand_coord_template}"
        )

    _merge_gro_files(protein_gro, ligand_gro, complex_gro)

    return {
        "ligand_itp": staged_ligand_itp,
        "ligand_toppar": staged_ligand_toppar,
        "ligand_gro": ligand_gro,
        "complex_gro": complex_gro,
        "topol_top": topol_top,
    }

def prepare_aa_equilibration_system(
    aa_job_out_dir: Path,
    complex_gro: Path,
    topol_top: Path,
    run_em: bool = True,
    run_nvt: bool = True,
    run_npt: bool = True,
    box_distance_nm: float = 1.0,
    box_type: str = "cubic",
) -> dict[str, Path]:
    _check_aa_mdp_templates()

    if run_npt and not run_nvt:
        raise ValueError("run_npt=True requires run_nvt=True")

    eq_dir = aa_job_out_dir / "aa_equil"
    eq_dir.mkdir(parents=True, exist_ok=True)

    boxed_gro = eq_dir / "boxed.gro"
    solvated_gro = eq_dir / "solvated.gro"
    ions_tpr = eq_dir / "ions.tpr"
    ionized_gro = eq_dir / "ionized.gro"

    em_tpr = eq_dir / "em.tpr"
    nvt_tpr = eq_dir / "nvt.tpr"
    npt_tpr = eq_dir / "npt.tpr"

    _run([
        "gmx", "editconf",
        "-f", str(complex_gro),
        "-o", str(boxed_gro),
        "-c",
        "-d", str(box_distance_nm),
        "-bt", box_type,
    ])

    _run([
        "gmx", "solvate",
        "-cp", str(boxed_gro),
        "-cs", "spc216.gro",
        "-o", str(solvated_gro),
        "-p", str(topol_top),
    ])

    _run([
        "gmx", "grompp",
        "-f", str(IONS_MDP),
        "-c", str(solvated_gro),
        "-p", str(topol_top),
        "-o", str(ions_tpr),
        "-maxwarn", "1",
    ])

    _run(
        [
            "gmx", "genion",
            "-s", str(ions_tpr),
            "-o", str(ionized_gro),
            "-p", str(topol_top),
            "-pname", "NA",
            "-nname", "CL",
            "-neutral",
            "-conc", "0.15",
        ],
        stdin_text="SOL\n",
    )

    em_gro = ionized_gro
    nvt_gro = em_gro
    npt_gro = nvt_gro

    if run_em:
        _run([
            "gmx", "grompp",
            "-f", str(EM_MDP),
            "-c", str(ionized_gro),
            "-p", str(topol_top),
            "-o", str(em_tpr),
        ])
        _run(["gmx", "mdrun", "-deffnm", str(eq_dir / "em")])
        em_gro = eq_dir / "em.gro"

    if run_nvt:
        _run([
            "gmx", "grompp",
            "-f", str(NVT_MDP),
            "-c", str(em_gro),
            "-r", str(em_gro),
            "-p", str(topol_top),
            "-o", str(nvt_tpr),
        ])
        _run(["gmx", "mdrun", "-deffnm", str(eq_dir / "nvt")])
        nvt_gro = eq_dir / "nvt.gro"

    if run_npt:
        _run([
            "gmx", "grompp",
            "-f", str(NPT_MDP),
            "-c", str(nvt_gro),
            "-r", str(nvt_gro),
            "-t", str(eq_dir / "nvt.cpt"),
            "-p", str(topol_top),
            "-o", str(npt_tpr),
        ])
        _run(["gmx", "mdrun", "-deffnm", str(eq_dir / "npt")])
        npt_gro = eq_dir / "npt.gro"

    return {
        "boxed_gro": boxed_gro,
        "solvated_gro": solvated_gro,
        "ionized_gro": ionized_gro,
        "em_gro": em_gro,
        "nvt_gro": nvt_gro,
        "npt_gro": npt_gro,
        "topol_top": topol_top,
    }

def _check_aa_mdp_templates() -> None:
    for p in [IONS_MDP, EM_MDP, NVT_MDP, NPT_MDP]:
        if not p.exists():
            raise FileNotFoundError(f"Missing MDP template: {p}")


def _ligand_structure_to_gro(
    input_structure: Path,
    output_gro: Path,
) -> None:
    suffix = input_structure.suffix.lower()

    if suffix == ".gro":
        shutil.copy2(input_structure, output_gro)
        return

    if suffix in {".pdb", ".mol2"}:
        _run([
            "gmx", "editconf",
            "-f", str(input_structure),
            "-o", str(output_gro),
        ])
        return

    raise RuntimeError(
        f"Unsupported ligand coordinate template format: {input_structure}"
    )

