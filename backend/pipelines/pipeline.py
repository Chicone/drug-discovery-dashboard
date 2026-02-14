#!/usr/bin/env python3
"""
FULL PIPELINE: OPM-oriented A2A GPCR -> Martini3 CG -> POPC membrane
-> EM -> (optional) NVT -> (optional) NPT -> (optional) production MD.

Based on the reference pipeline steps in:
A2A_OPM_Martini_Pipeline.txt

Requirements (must be on PATH):
  - martinize2
  - insane
  - gmx  (GROMACS)

Environment:
  - MARTINI_FF: path to your Martini v3 directory containing .itp files
    Example:
      export MARTINI_FF=/.../martini_v300

Typical usage:
  python build_cg_membrane_md.py \
    --workdir ~/PyCharm/ddd/pdb_membrane \
    --aa_pdb A2AR_AF_OPM.pdb \
    --name Protein \
    --do-em --do-nvt --do-npt --do-md
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
import json

MDP_DIR = Path(__file__).resolve().parent / "mdp_templates"

# ------------------------------- Config ------------------------------


@dataclass(frozen=True)
class MartiniIncludes:
    base_itp: str = "martini_v3.0.0.itp"
    lipids_itp: str = "martini_v3.0.0_phospholipids_v1.itp"
    solvents_itp: str = "martini_v3.0.0_solvents_v1.itp"
    ions_itp: str = "martini_v3.0.0_ions_v1.itp"


@dataclass(frozen=True)
class PipelineConfig:
    workdir: Path
    aa_pdb: Path
    protein_name: str
    martini_ff: Path

    # martinize2
    cg_pdb: Path
    cg_top: Path
    itp_src: Path
    itp_dst: Path

    # insane
    system_gro: Path
    system_top: Path

    # md
    em_mdp: Path
    em_tpr: Path
    em_deffnm: str

    nvt_mdp: Path
    nvt_tpr: Path
    nvt_deffnm: str

    npt_mdp: Path
    npt_tpr: Path
    npt_deffnm: str

    md_mdp: Path
    md_tpr: Path
    md_deffnm: str

    orth_itp: Optional[Path] = None
    allo_itp: Optional[Path] = None


# ------------------------------ Helpers ------------------------------


def run(cmd: List[str], cwd: Path, env: Optional[dict] = None) -> None:
    """Run a command, stream output, raise on error."""
    print("\n$ " + " ".join(cmd))
    try:
        subprocess.run(
            cmd,
            cwd=str(cwd),
            env=env,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Command failed with exit code {e.returncode}: {cmd}")


def must_exist(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")
    print(f"[write] {path}")


def replace_in_file(path: Path, pattern: str, repl: str) -> int:
    """Regex replace in file; returns number of substitutions."""
    data = path.read_text(encoding="utf-8")
    new_data, n = re.subn(pattern, repl, data, flags=re.MULTILINE)
    if n:
        path.write_text(new_data, encoding="utf-8")
        print(f"[patch] {path} ({n} substitutions)")
    return n


def patch_gro_ions(gro_path: Path) -> None:
    """
    INSANE may output CL- and NA+ in .gro. Martini usually expects CL and NA.
    We'll do safe text replacement (exact tokens).
    """
    data = gro_path.read_text(encoding="utf-8")

    # Replace only in the atom-name/resname fields by simple token matches.
    # This is pragmatic: in .gro, 'CL-' and 'NA+' appear as residue names.
    data2 = data.replace("CL-", "CL ").replace("NA+", "NA ")
    if data2 != data:
        gro_path.write_text(data2, encoding="utf-8")
        print(f"[patch] Ion naming fixed in {gro_path}")
    else:
        print(f"[patch] No ion naming changes needed in {gro_path}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def extract_protein_only(src_pdb: Path, dst_pdb: Path) -> None:
    """
    Keep only ATOM records (protein), drop HETATM.
    """
    lines = []
    for line in src_pdb.read_text().splitlines():
        if line.startswith("ATOM"):
            lines.append(line)
    dst_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")

def extract_ligand(
    complex_pdb: Path,
    dst_pdb: Path,
    resname: Optional[str] = None,
) -> None:
    """
    Extract ligand from a docked complex.
    If resname is None, extract all HETATM.
    """
    lines = []
    for line in complex_pdb.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        if resname is not None and line[17:20].strip() != resname:
            continue
        lines.append(line)

    if not lines:
        raise RuntimeError(f"No ligand atoms found in {complex_pdb}")

    dst_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")

def patch_small_beads(itp_path: Path):
    """
    Shrink the only bead in the ligand that can be made smaller:
    SC5 (Small apolar bead) -> TC5 (Tiny apolar bead).

    All other beads in the ligand (TN5a, TN4, TN3, TN5d, TC3, TC5)
    are ALREADY Tiny beads in Martini 3 and must not be modified.
    """
    if not itp_path.exists():
        raise FileNotFoundError(f"Cannot patch missing ITP: {itp_path}")

    # Only valid shrinkable mapping for your ligand:
    replacements = {
        " SC5 ": " TC4 ",
    }

    text = itp_path.read_text()

    changed = 0
    for old, new in replacements.items():
        n = text.count(old)
        if n:
            text = text.replace(old, new)
            changed += n

    itp_path.write_text(text)
    print(f"[patch] Applied {changed} SC5→TC5 replacements to {itp_path}")


from pathlib import Path

def patch_ligand_20pc(itp_path: Path):
    """
    Reduce ONLY dihedral force constants to 10% of original.
    Do NOT modify angles.
    Do NOT modify equilibrium angles.
    """

    if not itp_path.exists():
        raise FileNotFoundError(f"Cannot patch missing ITP: {itp_path}")

    lines = itp_path.read_text().splitlines()
    new_lines = []

    section = None

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("["):
            section = stripped.lower()

        # --- Dihedrals only ---
        if section == "[dihedrals]" and stripped and not stripped.startswith(";"):
            parts = re.split(r"\s+", stripped)
            if len(parts) >= 7:
                try:
                    k = float(parts[6])
                    parts[6] = f"{k * 0.2:.2f}"  # 70% softer
                    line = "   " + "   ".join(parts)
                except ValueError:
                    pass

        new_lines.append(line)

    itp_path.write_text("\n".join(new_lines) + "\n")
    print("[patch] Reduced dihedral force constants to 10%")

import re
from pathlib import Path


def patch_ligand(
    itp_path: Path,
    *,
    bond_factor: float = 1.0,
    angle_factor: float = 1.0,
    dihedral_factor: float = 0.15,
    bond_k_min: float = 0.0,
    angle_k_min: float = 0.0,
    dihedral_k_min: float = 50.0,
) -> None:
    """
    Independently scale bonded force constants in:

        - [ bonds ]     (k column 5)
        - [ angles ]    (k column 6)
        - [ dihedrals ] (k column 7)

    Only scales force constants >= corresponding *_k_min.
    """

    if not itp_path.exists():
        raise FileNotFoundError(f"Cannot patch missing ITP: {itp_path}")

    lines = itp_path.read_text().splitlines()
    out = []
    section = None
    changed_b = 0
    changed_a = 0
    changed_d = 0

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("["):
            section = stripped.lower()

        if not stripped or stripped.startswith(";"):
            out.append(line)
            continue

        parts = re.split(r"\s+", stripped)

        try:
            # -------------------------
            # BONDS
            # i j func r0 k
            # -------------------------
            if section == "[bonds]" and len(parts) >= 5:
                k_old = float(parts[4])
                if k_old >= bond_k_min and bond_factor != 1.0:
                    parts[4] = f"{k_old * bond_factor:.2f}"
                    line = "   " + "   ".join(parts)
                    changed_b += 1

            # -------------------------
            # ANGLES
            # i j k func theta0 k
            # -------------------------
            elif section == "[angles]" and len(parts) >= 6:
                k_old = float(parts[5])
                if k_old >= angle_k_min and angle_factor != 1.0:
                    parts[5] = f"{k_old * angle_factor:.2f}"
                    line = "   " + "   ".join(parts)
                    changed_a += 1

            # -------------------------
            # DIHEDRALS
            # i j k l func phi0 k mult
            # -------------------------
            elif section == "[dihedrals]" and len(parts) >= 7:
                k_old = float(parts[6])
                if k_old >= dihedral_k_min and dihedral_factor != 1.0:
                    parts[6] = f"{k_old * dihedral_factor:.2f}"
                    line = "   " + "   ".join(parts)
                    changed_d += 1

        except ValueError:
            pass

        out.append(line)

    itp_path.write_text("\n".join(out) + "\n")

    print(
        f"[patch] Bonds: {changed_b}, "
        f"Angles: {changed_a}, "
        f"Dihedrals: {changed_d}"
    )


# ------------------------------ Steps --------------------------------



def step_martinize(cfg: PipelineConfig) -> PipelineConfig:
    print("\n=== 1) MARTINIZE2 (AA -> CG) ===")
    must_exist(cfg.aa_pdb, "AA PDB input")

    cg_name = "Protein_cg.pdb"
    top_name = "Protein.top"
    itp_name = "Protein_0.itp"

    cmd = [
        "martinize2",
        "-f", cfg.aa_pdb.name,
        "-x", cg_name,
        "-o", top_name,
        "-ff", "martini3001",
        "-ss", "C",
        "-name", cfg.protein_name,
        "-p", "backbone",
        "-pf", "1000",
        "-maxwarn", "30",
    ]
    run(cmd, cwd=cfg.workdir)

    cg_path = cfg.workdir / cg_name
    top_path = cfg.workdir / top_name
    itp_src = cfg.workdir / itp_name

    must_exist(cg_path, "CG PDB output")
    must_exist(top_path, "CG TOP output")
    must_exist(itp_src, "martinize2 ITP output (Protein_0.itp)")

    shutil.copyfile(itp_src, cfg.itp_dst)
    print(f"[copy] {itp_src.name} → {cfg.itp_dst.name}")

    # return the updated config
    return cfg.__class__(**{
        **cfg.__dict__,
        "cg_pdb": cg_path,
        "cg_top": top_path,
        "itp_src": itp_src,
    })


def step_insane(cfg: PipelineConfig) -> None:
    print("\n=== 2) INSANE (build POPC bilayer + solvent + ions) ===")
    must_exist(cfg.cg_pdb, "CG PDB input for insane")

    cmd = [
        "insane",
        "-f", str(cfg.cg_pdb),
        "-o", str(cfg.system_gro),
        "-p", str(cfg.system_top),
        "-l", "POPC",
        "-sol", "W",
        "-salt", "0.02",
        "-center",
        # "-dm",  "0",
        "-ring",
        "-d", "1.2",
        "-box", "14.0,14.0,20.0",
        # "-box", "16.0,16.0,22.0",
         # "-box","11.0,11.0,15.0",
    ]
    run(cmd, cwd=cfg.workdir)

    must_exist(cfg.system_gro, "system.gro output")
    must_exist(cfg.system_top, "system.top output")

def step_insane_water(cfg: PipelineConfig, include_protein: bool) -> None:
    """
    Build a water-only system (no membrane).
    Uses INSANE but without -l POPC.
    """

    if include_protein:
        must_exist(cfg.cg_pdb, "CG PDB input for insane")
        input_pdb = str(cfg.cg_pdb)
    else:
        ligand_pdb = cfg.workdir / "orthosteric_cg.pdb"
        must_exist(ligand_pdb, "orthosteric_cg.pdb for ligand-only system")

        centered = cfg.workdir / "orthosteric_cg_centered.pdb"

        run([
            "gmx", "editconf",
            "-f", str(ligand_pdb),
            "-o", str(centered),
            "-c"
        ], cwd=cfg.workdir)

        input_pdb = str(centered)

    cmd = [
        "insane",
        "-f", input_pdb,
        "-o", str(cfg.system_gro),
        "-p", str(cfg.system_top),
        "-sol", "W",
        "-salt", "0.15",
        "-box", "12.0,12.0,12.0",
    ]

    run(cmd, cwd=cfg.workdir)

    must_exist(cfg.system_gro, "system.gro output")
    must_exist(cfg.system_top, "system.top output")


def step_patch_system_top(
    cfg: PipelineConfig,
    inc: MartiniIncludes,
    scenario: str,
) -> None:
    """
    Patch INSANE's system.top:
      - remove legacy martini.itp include
      - prepend Martini 3 includes + Protein.itp + Orthosteric.itp
      - ensure [ molecules ] contains ORT 1 if orth_itp exists
      - fix NA+/CL- -> NA/CL
    """
    must_exist(cfg.system_top, "system.top to patch")
    if scenario != "ligand_water":
        must_exist(cfg.itp_dst, "Protein.itp to include")

    orth_itp = getattr(cfg, "orth_itp", None)
    allo_itp = getattr(cfg, "allo_itp", None)

    print("[patch] orth_itp =", orth_itp)
    print("[patch] allo_itp =", allo_itp)

    lines = cfg.system_top.read_text(encoding="utf-8",
                                    errors="replace").splitlines()

    # ------------------------------------------------------------
    # 1) Remove ONLY the specific legacy includes (do not filter
    #    the rest of the file, or you'll drop ORTHO etc.)
    # ------------------------------------------------------------
    kept = []
    for l in lines:
        s = l.strip()
        if s.startswith("#include") and "martini.itp" in s:
            continue
        if s.startswith("#include") and 'Protein.top' in s:
            continue
        kept.append(l)

    # ------------------------------------------------------------
    # 2) Fix ion names everywhere (topology)
    # ------------------------------------------------------------
    kept = [l.replace("NA+", "NA").replace("CL-", "CL") for l in kept]

    # ------------------------------------------------------------
    # REMOVE Protein molecule entry for ligand-only systems
    # ------------------------------------------------------------
    if scenario == "ligand_water":
        cleaned = []
        inside_mol = False

        for line in kept:
            s = line.strip().lower()

            if s == "[ molecules ]":
                inside_mol = True
                cleaned.append(line)
                continue

            if inside_mol:
                if s.startswith("[") and s.endswith("]"):
                    inside_mol = False
                    cleaned.append(line)
                    continue

                # Skip Protein line
                if s.startswith("protein"):
                    print("[patch] Removed Protein from [ molecules ]")
                    continue

            cleaned.append(line)

        kept = cleaned

    # ------------------------------------------------------------
    # 3) Prepend Martini 3 includes if not already present
    # ------------------------------------------------------------
    def inc_line(fname: str) -> str:
        return f'#include "{(cfg.martini_ff / fname).as_posix()}"'

    include_block = [
        inc_line(inc.base_itp),
        inc_line(inc.lipids_itp),
        inc_line(inc.solvents_itp),
        inc_line(inc.ions_itp),
    ]

    # Only include Protein.itp if this scenario has protein
    if scenario != "ligand_water":
        include_block.append('#include "Protein.itp"')

    if orth_itp is not None:
        must_exist(orth_itp, "Orthosteric ligand itp")
        include_block.append(f'#include "{orth_itp.name}"')

    if allo_itp is not None:
        must_exist(allo_itp, "Allosteric ligand itp")
        include_block.append(f'#include "{allo_itp.name}"')

    already_patched = any("martini_v3.0.0" in l for l in kept)
    if not already_patched:
        kept = include_block + [""] + kept

    # ------------------------------------------------------------
    # 4) Ensure [ molecules ] contains ORTHO / ALLO counts
    #    (section-aware; won't clobber anything)
    # ------------------------------------------------------------
    entries = []
    if orth_itp is not None:
        entries.append(("ORT", 1))
    if allo_itp is not None:
        entries.append(("ALL", 1))

    if entries:
        # Find [ molecules ]
        mol_i = None
        for i, l in enumerate(kept):
            if l.strip().lower() == "[ molecules ]":
                mol_i = i
                break
        if mol_i is None:
            raise RuntimeError("system.top missing [ molecules ] section")

        # Determine end of molecules block (next [ section ] or EOF)
        mol_end = len(kept)
        for j in range(mol_i + 1, len(kept)):
            s = kept[j].strip()
            if s.startswith("[") and s.endswith("]"):
                mol_end = j
                break

        block = kept[mol_i:mol_end]
        block_stripped = [b.strip() for b in block]

        def has_mol(name: str) -> bool:
            # match first token in a molecules line
            for b in block_stripped:
                if not b or b.startswith(";"):
                    continue
                toks = b.split()
                if toks and toks[0] == name:
                    return True
            return False

        # Build *new* entries that need to be inserted
        to_add = []
        for name, count in entries:
            if not has_mol(name):
                to_add.append(f"{name:<8} {count}")

        if to_add:
            if scenario == "ligand_water":
                # ORT must be first molecule (matches system.gro)
                insert_at = mol_i + 1
                print("[patch] inserted at top (ligand_water):", ", ".join(to_add))
            else:
                # Insert after protein
                insert_at = None
                j = mol_i + 1
                while j < mol_end:
                    s = kept[j].strip()
                    if s and not s.startswith(";"):
                        insert_at = j + 1
                        break
                    j += 1
                if insert_at is None:
                    insert_at = mol_end
                print("[patch] inserted after protein:", ", ".join(to_add))

            kept = kept[:insert_at] + to_add + kept[insert_at:]
        else:
            print("[patch] [ molecules ] already contains ligand entries")

    # ------------------------------------------------------------
    # FIX: Ensure protein moleculetype matches martinize2 output
    # martinize2 always creates Protein_0.itp with:
    #   [ moleculetype ]
    #   Protein_0    1
    #
    # INSANE writes "Protein" in system.top, so we replace it.
    # ------------------------------------------------------------
    if scenario != "ligand_water":
        for idx, line in enumerate(kept):
            if line.strip().lower() == "[ molecules ]":
                # find first real molecules entry
                j = idx + 1
                while j < len(kept) and (
                    kept[j].strip() == "" or kept[j].strip().startswith(";")
                ):
                    j += 1

                # Replace the name only (first column)
                parts = kept[j].split()
                if parts:
                    old = parts[0]
                    parts[0] = "Protein_0"
                    kept[j] = f"{parts[0]:<15} {parts[1]}"
                    print(f"[patch] Molecule name fixed: {old} -> Protein_0")
                break


    # ------------------------------------------------------------
    # 5) Write once (avoid multi-pass overwrite bugs)
    # ------------------------------------------------------------
    cfg.system_top.write_text("\n".join(kept) + "\n", encoding="utf-8")

    # ------------------------------------------------------------
    # 6) Patch ions in .gro too (if present)
    # ------------------------------------------------------------
    if hasattr(cfg, "system_gro") and cfg.system_gro.exists():
        gro = cfg.system_gro.read_text(encoding="utf-8", errors="replace")
        gro = gro.replace("NA+", "NA ").replace("CL-", "CL ")
        cfg.system_gro.write_text(gro, encoding="utf-8")

def step_fix_ions(cfg: PipelineConfig) -> None:
    print("\n=== 3) FIX ION NAMING (CL-, NA+) ===")
    must_exist(cfg.system_gro, "system.gro for ion fix")
    patch_gro_ions(cfg.system_gro)

def step_write_mdps(cfg: PipelineConfig, md_ns_override: Optional[float] = None) -> None:
    print("\n=== 4) WRITE MDP FILES ===")
    print("[DEBUG] md_ns_override =", md_ns_override)

    # -------------------------------------------------
    # Select MDP template set based on scenario
    # -------------------------------------------------
    params_file = cfg.workdir.parent / "input" / "params.json"

    scenario = None
    if params_file.exists():
        params = json.loads(params_file.read_text())
        scenario = params.get("scenario")

    if scenario == "ligand_water":
        template_dir = MDP_DIR / "ligand_water"
        print("[mdp] Using ligand_water MDP templates")

    elif scenario.endswith("_water"):
        template_dir = MDP_DIR / "protein_only"
        print("[mdp] Using protein_water MDP templates")

    elif scenario.endswith("_membrane"):
        template_dir = MDP_DIR / "full"
        print("[mdp] Using membrane/full MDP templates")

    else:
        raise ValueError(f"Unknown scenario for MDP selection: {scenario}")

    def cp(name, out_path):
        src = template_dir / name
        if not src.exists():
            raise FileNotFoundError(f"Missing MDP template: {src}")
        write_text(out_path, src.read_text())

    # Copy default templates
    cp("em.mdp",  cfg.em_mdp)
    cp("nvt.mdp", cfg.nvt_mdp)
    cp("npt.mdp", cfg.npt_mdp)
    cp("md.mdp",  cfg.md_mdp)

    # -------------------------------------------------
    # Apply MD duration override (md_ns)
    # -------------------------------------------------
    if md_ns_override is not None:
        md_text = cfg.md_mdp.read_text()

        # convert ns → nsteps
        # dt is read from mdp
        dt = None
        for line in md_text.splitlines():
            if line.strip().startswith("dt"):
                dt = float(line.split("=")[1])
                break

        if dt is None:
            raise RuntimeError("md.mdp missing 'dt' parameter")

        nsteps = int((md_ns_override * 1000) / dt)

        md_text = re.sub(
            r"^nsteps\s*=\s*\d+",
            f"nsteps = {nsteps}",
            md_text,
            flags=re.MULTILINE
        )

        cfg.md_mdp.write_text(md_text)
        print(f"[patch] Applied md_ns override → nsteps = {nsteps}")



def step_grompp_mdrun_em(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 5-6) EM: grompp + mdrun ===")
    must_exist(cfg.em_mdp, "em.mdp")
    must_exist(cfg.system_gro, "system.gro")
    must_exist(cfg.system_top, "system.top")

    run(
        [
            "gmx",
            "grompp",
            "-f", str(cfg.em_mdp),
            "-c", str(cfg.system_gro),
            "-p", str(cfg.system_top),
            "-o", str(cfg.em_tpr),
            "-maxwarn", "1",
        ],
        cwd=cfg.workdir,
    )
    must_exist(cfg.em_tpr, "em.tpr")

    run(
        [
            "gmx", "mdrun",
            "-ntmpi", "1",
            "-nt", str(nt),
            "-v",
            "-deffnm", cfg.em_deffnm,
        ],
        cwd=cfg.workdir,
    )

    must_exist(cfg.workdir / f"{cfg.em_deffnm}.gro", "em.gro output")


def step_grompp_mdrun_nvt(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 8.1) NVT: grompp + mdrun ===")
    c_in = cfg.workdir / f"{cfg.em_deffnm}.gro"
    must_exist(c_in, "em.gro input for NVT")

    run(
        [
            "gmx", "grompp",
            "-f", str(cfg.nvt_mdp),
            "-c", str(c_in),
            "-r", str(c_in),
            "-p", str(cfg.system_top),
            "-o", str(cfg.nvt_tpr),
        ],
        cwd=cfg.workdir,
    )
    must_exist(cfg.nvt_tpr, "nvt.tpr")

    run(
        [
            "gmx",
            "mdrun",
            "-ntmpi", "1",
            "-nt", str(nt),
            "-v",
            "-deffnm", cfg.nvt_deffnm,
        ],
        cwd=cfg.workdir,
    )

    must_exist(cfg.workdir / f"{cfg.nvt_deffnm}.gro", "nvt.gro output")


def step_grompp_mdrun_npt(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 8.2) NPT: grompp + mdrun ===")

    workdir = cfg.workdir

    # ------------------------------------------------
    # Is this a continuation-run or a normal pipeline?
    # ------------------------------------------------
    params_file = workdir.parent / "input" / "params.json"
    if not params_file.exists():
        raise FileNotFoundError(f"Missing params.json: {params_file}")
    params = json.loads(params_file.read_text())

    continuation_mode = params.get("workflow") == "run_md"

    # ------------------------------------------------
    # Determine input structure for grompp
    # ------------------------------------------------
    if continuation_mode:
        start_gro = params.get("start_gro")
        if not start_gro:
            raise RuntimeError("params.json must define 'start_gro' for continuation")

        gro_in = workdir / start_gro

        start_cpt = params.get("start_cpt")
        cpt_in = workdir / start_cpt if start_cpt else None

        if cpt_in is not None and not cpt_in.exists():
            cpt_in = None

    else:
        gro_in = workdir / f"{cfg.nvt_deffnm}.gro"
        cpt_in = workdir / f"{cfg.nvt_deffnm}.cpt"
        if not cpt_in.exists():
            cpt_in = None

    # ------------------------------------------------
    # Run grompp (must always get a .gro, never .cpt)
    # ------------------------------------------------
    cmd = [
        "gmx", "grompp",
        "-f", str(cfg.npt_mdp),
        "-c", str(gro_in),
        "-r", str(gro_in),
        "-p", str(cfg.system_top),
        "-o", str(cfg.npt_tpr),
        "-maxwarn", "1",
    ]

    if cpt_in:
        cmd.extend(["-t", str(cpt_in)])

    run(cmd, cwd=workdir)
    must_exist(cfg.npt_tpr, "npt.tpr")

    # ------------------------------------------------
    # Run mdrun (can use -cpi checkpoint)
    # ------------------------------------------------
    mdrun_cmd = [
        "gmx", "mdrun",
        "-noappend",
        "-ntmpi", "1",
        "-nt", str(nt),
        "-v",
        "-deffnm", cfg.npt_deffnm,
        # "-c", f"{cfg.npt_deffnm}.gro", # <<< FORCE OUTPUT NAME
    ]

    # if cpt_in:
    #     mdrun_cmd.insert(3, "-cpi")
    #     mdrun_cmd.insert(4, str(cpt_in))

    run(mdrun_cmd, cwd=workdir)

    # Normalize output names: npt.part0001.* → npt.*
    suffixes = [".gro", ".log", ".edr", ".cpt", ".xtc"]

    for suf in suffixes:
        part = workdir / f"{cfg.npt_deffnm}.part0001{suf}"
        final = workdir / f"{cfg.npt_deffnm}{suf}"
        if part.exists() and not final.exists():
            part.rename(final)

    must_exist(workdir / f"{cfg.npt_deffnm}.gro", "npt.gro output")


def step_grompp_mdrun_md(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 8.3) PRODUCTION MD: grompp + mdrun ===")

    workdir = cfg.workdir

    # ------------------------------------------------
    # Load job parameters
    # ------------------------------------------------
    params_file = workdir.parent / "input" / "params.json"
    if not params_file.exists():
        raise FileNotFoundError(f"Missing params.json: {params_file}")

    params = json.loads(params_file.read_text())

    # ------------------------------------------------
    # Determine the correct starting structure
    # ------------------------------------------------
    start_gro = params.get("start_gro")
    if not start_gro:
        start_gro = "npt.gro"

    # start_gro = params.get("start_gro", "npt.gro")
    start_cpt = params.get("start_cpt")

    gro_in = workdir / start_gro
    if not gro_in.exists():
        raise FileNotFoundError(f"start_gro '{start_gro}' not found")

    # ------------------------------------------------
    # Correct continuation logic
    # ------------------------------------------------

    md_cpt = workdir / "md.cpt"
    if md_cpt.exists():
        # Case 2 or 3 — MD had started before and crashed or continuing
        cpt_in = md_cpt
    else:
        # No md.cpt exists → first MD after NPT
        # Do NOT use NPT checkpoint
        cpt_in = None

    # ------------------------------------------------
    # Generate MD TPR (uses md.mdp, NOT npt.mdp)
    # ------------------------------------------------
    run([
        "gmx", "grompp",
        "-f", str(cfg.md_mdp),
        "-c", str(gro_in),
        "-p", str(cfg.system_top),
        "-o", str(cfg.md_tpr),
        "-maxwarn", "1",
    ], cwd=workdir)

    # ------------------------------------------------
    # Run MD
    # ------------------------------------------------
    mdrun_cmd = [
        "gmx", "mdrun",
        "-noappend",
        "-ntmpi", "1",
        "-nt", str(nt),
        "-v",
        "-deffnm", cfg.md_deffnm,
    ]

    if cpt_in:
        mdrun_cmd.insert(3, "-cpi")
        mdrun_cmd.insert(4, str(cpt_in))

    run(mdrun_cmd, cwd=workdir)

    # ------------------------------------------------------------
    # Normalize output names (md.partXXXX.* → md.*)
    # Always pick the highest part number (latest continuation)
    # ------------------------------------------------------------
    import re

    # Find all part files for this MD
    part_files = list(workdir.glob(f"{cfg.md_deffnm}.part*.xtc"))

    if part_files:
        # Extract numeric suffix from filenames
        def part_index(path):
            m = re.search(r"\.part(\d+)\.", path.name)
            return int(m.group(1)) if m else -1

        # Pick the highest index → latest continuation
        latest_part = max(part_files, key=part_index)
        latest_idx = part_index(latest_part)
        tag = f".part{latest_idx:04d}"

        suffixes = [".gro", ".log", ".edr", ".cpt", ".xtc"]

        for suf in suffixes:
            src = workdir / f"{cfg.md_deffnm}{tag}{suf}"
            dst = workdir / f"{cfg.md_deffnm}{suf}"
            if src.exists():
                # overwrite the destination, always choose latest
                shutil.copy2(src, dst)
                print(f"[normalize] {src.name} -> {dst.name}")
    else:
        print("[normalize] No md.partXXXX.* files found")

def merge_cg_pdbs(out_pdb: Path, *pdbs: Path) -> None:
    """
    Merge CG PDBs into a strict PDB (INSANE-friendly).
    Protein stays ATOM.
    Ligand detected by resname ORT/ORTH or chain L, becomes HETATM ORT L 1.
    Atom serials are renumbered consecutively.
    """

    def infer_elem(atom_name: str) -> str:
        s = re.sub(r"[^A-Za-z]", "", atom_name).upper()
        if s.startswith(("CL", "BR")):
            return s[:2]
        return (s[:1] or "C")

    def parse_line(line: str):
        # Prefer fixed columns when present
        if len(line) >= 54:
            rec = line[0:6].strip()
            name = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21].strip() or "A"
            resid_txt = line[22:26].strip() or "1"
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occ = float(line[54:60].strip() or "1.00") if len(line) >= 60 else 1.00
            bfac = float(line[60:66].strip() or "0.00") if len(line) >= 66 else 0.00
            elem = line[76:78].strip().upper() if len(line) >= 78 else ""
            resid = int(resid_txt)
            if not elem:
                elem = infer_elem(name)
            return rec, name, resn, chain, resid, x, y, z, occ, bfac, elem

        # Fallback: whitespace parsing (handles odd lines)
        parts = line.split()
        # ATOM serial name resname chain resid x y z [occ] [bfac] [elem]
        rec = parts[0]
        name = parts[2]
        resn = parts[3]
        chain = parts[4]
        resid = int(parts[5])
        x, y, z = map(float, parts[6:9])
        occ = float(parts[9]) if len(parts) > 9 else 1.00
        bfac = float(parts[10]) if len(parts) > 10 else 0.00
        elem = parts[-1].upper() if re.fullmatch(r"[A-Za-z]{1,2}", parts[-1]) else infer_elem(name)
        return rec, name, resn, chain, resid, x, y, z, occ, bfac, elem

    out = []
    serial = 1

    for pdb in pdbs:
        for line in pdb.read_text().splitlines():
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            rec, name, resn, chain, resid, x, y, z, occ, bfac, elem = parse_line(line)

            is_ligand = (resn.upper() in {"ORT", "ORTH"} or chain.upper() == "L")
            if is_ligand:
                rec = "HETATM"
                resn = "ORT"
                chain = "L"
                resid = 1
            else:
                rec = "ATOM"
                resn = (resn[:3] or "RES")
                chain = (chain[:1] or "A")

            out.append(
                f"{rec:<6s}{serial:5d} "
                f"{name:>4s} "
                f"{resn:>3s} "
                f"{chain:1s}{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{occ:6.2f}{bfac:6.2f}"
                f"{'':10s}{elem:>2s}"
            )
            serial += 1

    out.append("END")
    out_pdb.write_text("\n".join(out) + "\n")

import hashlib

def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()

def compute_cache_key(protein_aa: Path,
                      ligand_aa: Optional[Path],
                      smiles: str) -> str:

    h = hashlib.sha256()
    h.update(protein_aa.read_bytes())

    if ligand_aa and ligand_aa.exists():
        h.update(ligand_aa.read_bytes())
    else:
        h.update(b"NO-LIGAND")

    h.update(smiles.encode("utf8"))

    return h.hexdigest()[:16]


def cache_check_before_martinize(cfg, args, cache_root: Path, scenario: str):
    """
    Checks whether a cached martinized complex exists.
    If yes:
        - restores all CG files into cfg.workdir
        - returns (True, updated_cfg, cache_dir)
    If no:
        - returns (False, cfg, cache_dir)
    """

    # Build simple cache key from protein + ligand AA PDBs
    import hashlib
    hasher = hashlib.sha256()

    if scenario != "ligand_water":
        hasher.update(cfg.aa_pdb.read_bytes())
    else:
        hasher.update(b"NO-PROTEIN")

    # Ligand AA PDB (if present)
    if args.orthosteric_pdb:
        lig = Path(args.orthosteric_pdb)
        if lig.exists():
            hasher.update(lig.read_bytes())

    # Scenario must affect cache key
    hasher.update(scenario.encode("utf-8"))

    key = hasher.hexdigest()[:16]
    cache_dir = cache_root / key

    # If cache does not exist → miss
    if not cache_dir.exists():
        print(f"\n[CACHE MISS] No cached martinized complex ({key})")
        return False, cfg, cache_dir

    # Otherwise → HIT
    print(f"\n[CACHE HIT] Using cached martinized complex: {cache_dir}")

    # Restore important files
    restore_list = [
        "Protein.top",
        "Protein_0.itp",
        "Protein.itp",
        "Protein_cg.pdb",
        "orthosteric_aa.pdb",
        "orthosteric_cg.pdb",
        "Orthosteric.itp",
        "cg_complex.pdb",
    ]

    for fname in restore_list:
        src = cache_dir / fname
        if src.exists():
            shutil.copy2(src, cfg.workdir / fname)

    # Update cfg with restored paths
    new_cfg = cfg.__class__(**{
        **cfg.__dict__,
        "cg_pdb": cfg.workdir / "cg_complex.pdb",
        "cg_top": cfg.workdir / "Protein.top",
        "itp_src": cfg.workdir / "Protein_0.itp",
        "itp_dst": cfg.workdir / "Protein.itp",
        "orth_itp": (cfg.workdir / "Orthosteric.itp")
                    if (cache_dir / "Orthosteric.itp").exists()
                    else None,
    })

    return True, new_cfg, cache_dir

def cache_save_after_martinize(cfg, cache_dir: Path):
    """
    Save the martinized protein + ligand CG complex into cache_dir.
    Only files that matter before INSANE are saved.
    """
    restore_list = [
        "Protein.top",
        "Protein_0.itp",
        "Protein.itp",
        "Protein_cg.pdb",
        "orthosteric_aa.pdb",
        "orthosteric_cg.pdb",
        "Orthosteric.itp",
        "cg_complex.pdb",
    ]

    cache_dir.mkdir(parents=True, exist_ok=True)

    for fname in restore_list:
        src = cfg.workdir / fname
        if src.exists():
            shutil.copy2(src, cache_dir / fname)

    print(f"[CACHE SAVED] Martinized complex stored in: {cache_dir}")




# ------------------------------ Main ---------------------------------


def build_config(args: argparse.Namespace) -> PipelineConfig:
    workdir = Path(args.workdir).expanduser().resolve()
    aa_pdb = (workdir / args.aa_pdb).resolve()
    martini_ff = Path(args.martini_ff).expanduser().resolve()

    protein_name = args.name

    cg_pdb = workdir / "A2AR_AF_cg.pdb"
    cg_top = workdir / "A2AR_AF.top"
    itp_src = workdir / "molecule_0.itp"
    itp_dst = workdir / "Protein.itp"

    system_gro = workdir / "system.gro"
    system_top = workdir / "system.top"

    em_mdp = workdir / "em.mdp"
    em_tpr = workdir / "em.tpr"
    em_deffnm = "em"

    nvt_mdp = workdir / "nvt.mdp"
    nvt_tpr = workdir / "nvt.tpr"
    nvt_deffnm = "nvt"

    npt_mdp = workdir / "npt.mdp"
    npt_tpr = workdir / "npt.tpr"
    npt_deffnm = "npt"

    md_mdp = workdir / "md.mdp"
    md_tpr = workdir / "md.tpr"
    md_deffnm = "md"

    return PipelineConfig(
        workdir=workdir,
        aa_pdb=aa_pdb,
        protein_name=protein_name,
        martini_ff=martini_ff,
        cg_pdb=cg_pdb,
        cg_top=cg_top,
        itp_src=itp_src,
        itp_dst=itp_dst,
        system_gro=system_gro,
        system_top=system_top,
        em_mdp=em_mdp,
        em_tpr=em_tpr,
        em_deffnm=em_deffnm,
        nvt_mdp=nvt_mdp,
        nvt_tpr=nvt_tpr,
        nvt_deffnm=nvt_deffnm,
        npt_mdp=npt_mdp,
        npt_tpr=npt_tpr,
        npt_deffnm=npt_deffnm,
        md_mdp=md_mdp,
        md_tpr=md_tpr,
        md_deffnm=md_deffnm,
    )


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", required=True)
    p.add_argument("--aa_pdb", required=True,
                   help="AA oriented PDB inside workdir")
    p.add_argument("--orthosteric_pdb", default=None,
                   help="AA complex PDB (or ligand PDB) used to extract orthosteric ligand", )

    p.add_argument("--name", default="Protein",
                   help="Protein name used by martinize2/topology")
    p.add_argument("--martini_ff", default=os.environ.get("MARTINI_FF", ""),
                   help="Path to Martini v3 directory (or set MARTINI_FF)")
    p.add_argument("--nt", type=int, default=1,
                   help="OpenMP threads for mdrun (-nt)")
    p.add_argument("--do-em", action="store_true")
    p.add_argument("--do-nvt", action="store_true")
    p.add_argument("--do-npt", action="store_true")
    p.add_argument("--do-md", action="store_true")
    p.add_argument("--skip-martinize", action="store_true")
    p.add_argument("--skip-insane", action="store_true")
    p.add_argument("--skip-patch-top", action="store_true")
    p.add_argument("--skip-ion-fix", action="store_true")
    p.add_argument("--skip-mdp-write", action="store_true")
    return p.parse_args(argv)

def main(argv=None) -> None:
    args = parse_args(argv)
    cfg = build_config(args)

    CG_CACHE = Path("backend/data/cg_cache")
    CG_CACHE.mkdir(parents=True, exist_ok=True)

    # ===============================================================
    # Read scenario from params.json (affects caching and build)
    # ===============================================================
    params_file = cfg.workdir.parent / "input" / "params.json"

    if not params_file.exists():
        raise RuntimeError("Missing params.json: scenario must be defined")

    params = json.loads(params_file.read_text())

    scenario = params.get("scenario")
    if scenario is None:
        raise RuntimeError("params.json missing 'scenario'")

    # ============================================
    # MD-ONLY CONTINUATION MODE (SKIP FULL PIPELINE)
    # ============================================
    print("ARE YOU DOING MD? ", args.do_md)

    if args.do_md:
        print("MD-only mode detected -> skipping ALL build steps.")

        # Load MD override from params.json
        params_file = cfg.workdir.parent / "input" / "params.json"
        md_ns = None
        if params_file.exists():
            params = json.loads(params_file.read_text())
            md_ns = params.get("md_ns")

        # 1) Write MD parameters (md.mdp with correct nsteps)
        step_write_mdps(cfg, md_ns_override=md_ns)

        # -------------------------------------------------------------
        # 2) NEW — Always generate system.pdb for visualization (Chimera)
        # -------------------------------------------------------------
        # md continuation jobs DO have system.gro in OUT from parent job,
        # because md_jobs.py copies it into input/system.gro
        system_gro = cfg.workdir / "system.gro"
        system_pdb = cfg.workdir / "system.pdb"

        if system_gro.exists():
            run([
                "gmx", "editconf",
                "-f", str(system_gro),
                "-o", str(system_pdb),
            ], cwd=cfg.workdir)
            print("[MD-only] Generated system.pdb for visualization.")
        else:
            print("[MD-only] WARNING: system.gro not found, could not create system.pdb.")

        # 3) Continue MD
        step_grompp_mdrun_md(cfg, nt=args.nt)

        return


    ensure_dir(cfg.workdir)

    if not cfg.martini_ff or not cfg.martini_ff.exists():
        raise SystemExit(
            "Martini force-field path missing.\n"
            "Pass --martini_ff /path/to/martini_v300 or export MARTINI_FF."
        )

    inc = MartiniIncludes()

    print(f"[cfg] workdir    = {cfg.workdir}")
    print(f"[cfg] aa_pdb     = {cfg.aa_pdb}")
    print(f"[cfg] martini_ff = {cfg.martini_ff}")
    print(f"[cfg] nt         = {args.nt}")

    # ===============================================================
    # 1) Extract protein-only AA PDB  (monomer or dimer preserved)
    # ===============================================================
    if scenario != "ligand_water":
        protein_aa_clean = cfg.workdir / "protein_aa_clean.pdb"
        extract_protein_only(cfg.aa_pdb, protein_aa_clean)

        cfg = cfg.__class__(**{**cfg.__dict__, "aa_pdb": protein_aa_clean})

    # ===============================================================
    # SPECIAL CASE: Dimer scenario (placeholder)
    # ===============================================================
    params_file = cfg.workdir.parent / "input" / "params.json"
    if params_file.exists():
        params = json.loads(params_file.read_text())
        if params.get("scenario") == "dimer_two_ligands":
            raise RuntimeError(
                "Dimer build scenario detected. "
                "Dimer construction is not implemented yet."
            )

    # ===============================================================
    # 2) Process orthosteric ligand *if provided*
    # ===============================================================
    ortho_aa = None
    if args.orthosteric_pdb:
        ortho_aa = cfg.workdir / "orthosteric_aa.pdb"
        extract_ligand(Path(args.orthosteric_pdb), ortho_aa)

    # CACHE CHECK
    hit, cfg, cache_dir = cache_check_before_martinize(
        cfg, args, CG_CACHE, scenario
    )

    SKIP_LIG = hit  # skip ligand build if cache hit
    args.skip_martinize = hit

    # ===============================================================
    # 3) Martinize PROTEIN (skip for ligand-only case)
    # ===============================================================
    if not args.skip_martinize and scenario != "ligand_water":
        step_martinize(cfg)

    # --- Fix config paths after martinize2 ---
    if scenario != "ligand_water":
        cfg = cfg.__class__(**{
            **cfg.__dict__,
            "cg_pdb": cfg.workdir / "Protein_cg.pdb",
            "cg_top": cfg.workdir / "Protein.top",
            "itp_src": cfg.workdir / "Protein_0.itp",
            "itp_dst": cfg.workdir / "Protein.itp",
        })

    # ===============================================================
    # 4) Martinize LIGAND (if provided)
    # ===============================================================
    orth_itp = None
    orth_cg = None

    if args.orthosteric_pdb and not SKIP_LIG:
        # Fresh ligand build
        from .ligand_cg_builder import build_ligand

        SMILES_ORTHO = "n1c(N2CCCCC2)c(C#N)c(c3ccccc3)c(C#N)c(N)1"

        orth_cg = cfg.workdir / "orthosteric_cg.pdb"
        orth_itp = cfg.workdir / "Orthosteric.itp"

        build_ligand(
            SMILES_ORTHO,
            complex_pdb=Path(args.orthosteric_pdb),
            cg_pdb_out=orth_cg,
            itp_out=orth_itp,
        )

        # patch_small_beads(orth_itp)
        # patch_ligand(orth_itp, factor=0.0, k_min=10)

        if scenario != "ligand_water":
            # Combine protein + ligand into cg_complex.pdb
            cg_complex = cfg.workdir / "cg_complex.pdb"
            merge_cg_pdbs(cg_complex, cfg.cg_pdb, orth_cg)

            # Update cfg
            cfg = cfg.__class__(**{
                **cfg.__dict__,
                "cg_pdb": cg_complex,
                "orth_itp": orth_itp,
            })
        else:
            # Ligand-only system: use ligand CG directly
            cfg = cfg.__class__(**{
                **cfg.__dict__,
                "cg_pdb": orth_cg,
                "orth_itp": orth_itp,
            })

        # Save to cache since ligand was built now
        cache_save_after_martinize(cfg, cache_dir)

    elif args.orthosteric_pdb and SKIP_LIG:
        orth_itp = cfg.workdir / "Orthosteric.itp"
        orth_cg = cfg.workdir / "orthosteric_cg.pdb"

        if orth_itp.exists():
            # patch_small_beads(orth_itp)
            patch_ligand(
                orth_itp,
                bond_factor=1,
                angle_factor=1,
                dihedral_factor=1,
                bond_k_min=0,
                angle_k_min=0,
                dihedral_k_min=0,
            )
            print("[CACHE] Patched Orthosteric.itp from cache.")

        if scenario != "ligand_water":
            cfg = cfg.__class__(**{
                **cfg.__dict__,
                "orth_itp": orth_itp,
                "cg_pdb": cfg.workdir / "cg_complex.pdb",
            })
        else:
            cfg = cfg.__class__(**{
                **cfg.__dict__,
                "orth_itp": orth_itp,
                "cg_pdb": orth_cg,
            })


    elif not args.orthosteric_pdb:
        # No ligand in this project
        cfg = cfg.__class__(**{
            **cfg.__dict__,
            "orth_itp": None,
        })

    # ===============================================================
    # 5) INSANE Build environment depending on scenario
    # ===============================================================
    if not args.skip_insane:

        # Membrane systems
        if scenario.endswith("_membrane"):
            print(f"\n[SCENARIO] {scenario} → membrane system")
            step_insane(cfg)

        # Water-only systems
        elif scenario.endswith("_water"):
            print(f"\n[SCENARIO] {scenario} → solvated box")

            if scenario == "ligand_water":
                step_insane_water(cfg, include_protein=False)
            else:
                step_insane_water(cfg, include_protein=True)

        else:
            raise ValueError(f"Unknown scenario: {scenario}")

    if not args.skip_patch_top:
        step_patch_system_top(cfg, inc, scenario)

    if not args.skip_ion_fix:
        step_fix_ions(cfg)

    if not args.skip_mdp_write:
        # Load MD duration override
        params_file = cfg.workdir.parent / "input" / "params.json"
        md_ns = None
        if params_file.exists():
            params = json.loads(params_file.read_text())
            md_ns = params.get("md_ns")

        step_write_mdps(cfg, md_ns_override=md_ns)

    # Create pdb of the full system
    run([
        "gmx", "editconf",
        "-f", str(cfg.workdir / "system.gro"),
        "-o", str(cfg.workdir / "system.pdb"),
    ], cwd=cfg.workdir)

    # ===============================================================
    # NEW: preset-controlled pipeline
    # ===============================================================
    params_file = cfg.workdir.parent / "input" / "params.json"
    if params_file.exists():
        params = json.loads(params_file.read_text())
        preset = params.get("preset")
    else:
        preset = None

    if preset == "m3_popc_full":
        print("\n=== FULL PRESET: EM → NVT → NPT → MD ===")
        step_grompp_mdrun_em(cfg, nt=args.nt)
        step_grompp_mdrun_nvt(cfg, nt=args.nt)
        step_grompp_mdrun_npt(cfg, nt=args.nt)
        step_grompp_mdrun_md(cfg, nt=args.nt)
        print("\nFULL PRESET FINISHED.")
        return

    if args.do_em:
        step_grompp_mdrun_em(cfg, nt=args.nt)

    if args.do_nvt:
        step_grompp_mdrun_nvt(cfg, nt=args.nt)

    if args.do_npt:
        step_grompp_mdrun_npt(cfg, nt=args.nt)

    # if args.do_md:
    #     step_grompp_mdrun_md(cfg, nt=args.nt)

    print("\nDONE.")
    print("Visualize with:")
    print("  vmd system.gro")
    if (cfg.workdir / "em.gro").exists():
        print("  vmd em.gro")
    if (cfg.workdir / "npt.gro").exists():
        print("  vmd npt.gro")


if __name__ == "__main__":
    main()
