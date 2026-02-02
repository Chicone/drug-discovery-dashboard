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


# ------------------------------ MDPs ---------------------------------


EM_MDP = """\
integrator    = steep
emtol         = 1000.0
emstep        = 0.01
nsteps        = 50000

cutoff-scheme = Verlet
nstlist       = 20
rlist         = 1.1
coulombtype   = Reaction-Field
rcoulomb      = 1.1
vdwtype       = Cut-off
rvdw          = 1.1

pbc           = xyz
tcoupl        = no
pcoupl        = no
constraints   = none
"""

# Minimal placeholders. Tune as you like, but these get you a runnable baseline.
# For Martini 3 GPCRs you may want protein position restraints initially.
NVT_MDP = """\
integrator               = md
dt                       = 0.02
nsteps                   = 25000        ; 0.5 ns

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.1
coulombtype              = Reaction-Field
rcoulomb                 = 1.1
vdwtype                  = Cut-off
rvdw                     = 1.1

tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0
ref_t                    = 310

pcoupl                   = no
pbc                      = xyz

constraints              = none
"""

NPT_MDP = """\
integrator               = md
dt                       = 0.02
nsteps                   = 50000        ; 1.0 ns

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.1
coulombtype              = Reaction-Field
rcoulomb                 = 1.1
vdwtype                  = Cut-off
rvdw                     = 1.1

tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0
ref_t                    = 310

pcoupl                   = berendsen
pcoupltype               = semiisotropic
tau_p                    = 12.0
ref_p                    = 1.0 1.0
compressibility          = 3e-4 3e-4

pbc                      = xyz
constraints              = none
"""

MD_MDP = """\
integrator               = md
dt                       = 0.02
nsteps                   = 500000       ; 10 ns (increase to 5e6 for 100 ns)

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.1
coulombtype              = Reaction-Field
rcoulomb                 = 1.1
vdwtype                  = Cut-off
rvdw                     = 1.1

tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0
ref_t                    = 310

pcoupl                   = parrinello-rahman
pcoupltype               = semiisotropic
tau_p                    = 12.0
ref_p                    = 1.0 1.0
compressibility          = 3e-4 3e-4

pbc                      = xyz
constraints              = none


nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000   ; save frame every 5000 steps
"""


# ------------------------------ Steps --------------------------------


def step_martinize(cfg: PipelineConfig) -> None:
    print("\n=== 1) MARTINIZE2 (AA -> CG) ===")
    must_exist(cfg.aa_pdb, "AA PDB input")

    cmd = [
        "martinize2",
        "-f",
        str(cfg.aa_pdb),
        "-x",
        str(cfg.cg_pdb),
        "-o",
        str(cfg.cg_top),
        "-ff",
        "martini3001",
        "-ss",
        "C",
        "-name",
        cfg.protein_name,
        "-p",
        "backbone",
        "-pf",
        "1000",
        "-maxwarn",
        "30",
    ]
    run(cmd, cwd=cfg.workdir)

    must_exist(cfg.cg_pdb, "CG PDB output")
    must_exist(cfg.cg_top, "CG TOP output")
    must_exist(cfg.itp_src, "martinize2 ITP output (molecule_0.itp)")

    # Rename molecule_0.itp -> Protein.itp
    shutil.copyfile(cfg.itp_src, cfg.itp_dst)
    print(f"[copy] {cfg.itp_src.name} -> {cfg.itp_dst.name}")


def step_insane(cfg: PipelineConfig) -> None:
    print("\n=== 2) INSANE (build POPC bilayer + solvent + ions) ===")
    must_exist(cfg.cg_pdb, "CG PDB input for insane")

    cmd = [
        "insane",
        "-f",
        str(cfg.cg_pdb),
        "-o",
        str(cfg.system_gro),
        "-p",
        str(cfg.system_top),
        "-l",
        "POPC",
        "-sol",
        "W",
        "-salt",
        "0.15",
        "-center",
        "-dm",
        "0",
        "-box",
        "12.0,12.0,10.0",
    ]
    run(cmd, cwd=cfg.workdir)

    must_exist(cfg.system_gro, "system.gro output")
    must_exist(cfg.system_top, "system.top output")


def step_patch_system_top(cfg: PipelineConfig,
                          inc: MartiniIncludes) -> None:
    print("\n=== 2b) PATCH system.top (#includes + protein itp) ===")
    must_exist(cfg.system_top, "system.top to patch")
    must_exist(cfg.itp_dst, "Protein.itp to include")

    # Build include block using MARTINI_FF path.
    def inc_line(fname: str) -> str:
        return f'#include "{(cfg.martini_ff / fname).as_posix()}"'

    includes = "\n".join(
        [
            inc_line(inc.base_itp),
            inc_line(inc.lipids_itp),
            inc_line(inc.solvents_itp),
            inc_line(inc.ions_itp),
            '#include "Protein.itp"',
        ]
    )

    data = cfg.system_top.read_text(encoding="utf-8")

    # If includes already present, don't duplicate.
    if 'Protein.itp' in data and inc.base_itp in data:
        print("[patch] system.top already seems to have includes. Skipping.")
    else:
        new_data = includes + "\n\n" + data
        cfg.system_top.write_text(new_data, encoding="utf-8")
        print("[patch] Preprended Martini includes to system.top")

    # INSANE often writes molecule_0 in [ molecules ].
    replace_in_file(cfg.system_top, r"^\s*molecule_0\s+1\s*$",
                    "Protein   1")


def step_fix_ions(cfg: PipelineConfig) -> None:
    print("\n=== 3) FIX ION NAMING (CL-, NA+) ===")
    must_exist(cfg.system_gro, "system.gro for ion fix")
    patch_gro_ions(cfg.system_gro)


def step_write_mdps(cfg: PipelineConfig) -> None:
    print("\n=== 4) WRITE MDP FILES ===")
    write_text(cfg.em_mdp, EM_MDP)
    write_text(cfg.nvt_mdp, NVT_MDP)
    write_text(cfg.npt_mdp, NPT_MDP)
    write_text(cfg.md_mdp, MD_MDP)


def step_grompp_mdrun_em(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 5-6) EM: grompp + mdrun ===")
    must_exist(cfg.em_mdp, "em.mdp")
    must_exist(cfg.system_gro, "system.gro")
    must_exist(cfg.system_top, "system.top")

    run(
        [
            "gmx",
            "grompp",
            "-f",
            str(cfg.em_mdp),
            "-c",
            str(cfg.system_gro),
            "-p",
            str(cfg.system_top),
            "-o",
            str(cfg.em_tpr),
            "-maxwarn",
            "1",
        ],
        cwd=cfg.workdir,
    )
    must_exist(cfg.em_tpr, "em.tpr")

    run(
        [
            "gmx",
            "mdrun",
            "-ntmpi",
            "1",
            "-nt",
            str(nt),
            "-v",
            "-deffnm",
            cfg.em_deffnm,
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
            "gmx",
            "grompp",
            "-f",
            str(cfg.nvt_mdp),
            "-c",
            str(c_in),
            "-p",
            str(cfg.system_top),
            "-o",
            str(cfg.nvt_tpr),
        ],
        cwd=cfg.workdir,
    )
    must_exist(cfg.nvt_tpr, "nvt.tpr")

    run(
        [
            "gmx",
            "mdrun",
            "-ntmpi",
            "1",
            "-nt",
            str(nt),
            "-v",
            "-deffnm",
            cfg.nvt_deffnm,
        ],
        cwd=cfg.workdir,
    )

    must_exist(cfg.workdir / f"{cfg.nvt_deffnm}.gro", "nvt.gro output")


def step_grompp_mdrun_npt(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 8.2) NPT: grompp + mdrun ===")
    c_in = cfg.workdir / f"{cfg.nvt_deffnm}.gro"
    must_exist(c_in, "nvt.gro input for NPT")

    run(
        [
            "gmx",
            "grompp",
            "-f",
            str(cfg.npt_mdp),
            "-c",
            str(c_in),
            "-r",
            str(c_in),
            "-p",
            str(cfg.system_top),
            "-o",
            str(cfg.npt_tpr),
            "-maxwarn",
            "1",
        ],
        cwd=cfg.workdir,
    )
    must_exist(cfg.npt_tpr, "npt.tpr")

    run(
        [
            "gmx",
            "mdrun",
            "-ntmpi",
            "1",
            "-nt",
            str(nt),
            "-v",
            "-deffnm",
            cfg.npt_deffnm,
        ],
        cwd=cfg.workdir,
    )

    must_exist(cfg.workdir / f"{cfg.npt_deffnm}.gro", "npt.gro output")


def step_grompp_mdrun_md(cfg: PipelineConfig, nt: int) -> None:
    print("\n=== 8.3) PRODUCTION MD: grompp + mdrun ===")
    c_in = cfg.workdir / f"{cfg.npt_deffnm}.gro"
    must_exist(c_in, "npt.gro input for MD")

    run(
        [
            "gmx",
            "grompp",
            "-f",
            str(cfg.md_mdp),
            "-c",
            str(c_in),
            "-p",
            str(cfg.system_top),
            "-o",
            str(cfg.md_tpr),
        ],
        cwd=cfg.workdir,
    )
    must_exist(cfg.md_tpr, "md.tpr")

    run(
        [
            "gmx",
            "mdrun",
            "-ntmpi",
            "1",
            "-nt",
            str(nt),
            "-v",
            "-deffnm",
            cfg.md_deffnm,
        ],
        cwd=cfg.workdir,
    )


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

    if not args.skip_martinize:
        step_martinize(cfg)

    if not args.skip_insane:
        step_insane(cfg)

    if not args.skip_patch_top:
        step_patch_system_top(cfg, inc)

    if not args.skip_ion_fix:
        step_fix_ions(cfg)

    if not args.skip_mdp_write:
        step_write_mdps(cfg)

    if args.do_em:
        step_grompp_mdrun_em(cfg, nt=args.nt)

    if args.do_nvt:
        step_grompp_mdrun_nvt(cfg, nt=args.nt)

    if args.do_npt:
        step_grompp_mdrun_npt(cfg, nt=args.nt)

    if args.do_md:
        step_grompp_mdrun_md(cfg, nt=args.nt)

    print("\nDONE.")
    print("Visualize with:")
    print("  vmd system.gro")
    if (cfg.workdir / "em.gro").exists():
        print("  vmd em.gro")
    if (cfg.workdir / "npt.gro").exists():
        print("  vmd npt.gro")


if __name__ == "__main__":
    main()
