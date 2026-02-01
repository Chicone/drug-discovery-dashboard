# backend/md/pipelines/core.py

import json
import shutil
import subprocess
from pathlib import Path

from pathlib import Path
MDP_DIR = Path(__file__).resolve().parent / "mdp"

# -------------------------------------------------------
# Utilities
# -------------------------------------------------------

def load_params(job_dir: Path):
    params_path = job_dir / "input" / "params.json"
    return json.loads(params_path.read_text())

def append_log(job_dir: Path, line: str):
    with open(job_dir / "log.txt", "a", encoding="utf-8") as f:
        f.write(line)

def run_cmd(job_dir: Path, cmd, cwd=None):
    append_log(job_dir, f"\n$ {' '.join(cmd)}\n")

    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert proc.stdout
    for line in proc.stdout:
        append_log(job_dir, line)

    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"Command failed (rc={rc}): {' '.join(cmd)}")
    return rc


# -------------------------------------------------------
# Build steps
# -------------------------------------------------------

def build_membrane_system(job_dir: Path,
                          workdir: Path,
                          params: dict,
                          cfg: dict) -> None:
    steps = cfg["steps"]

    workdir.mkdir(parents=True, exist_ok=True)

    in_pdb = job_dir / "input" / "protein.pdb"
    out_pdb = workdir / "protein.pdb"
    if not in_pdb.exists():
        raise FileNotFoundError(f"Missing input PDB: {in_pdb}")
    shutil.copyfile(in_pdb, out_pdb)

    martini_dir = (
        Path.home()
        / "miniforge3"
        / "envs"
        / "drugdash"
        / "martini_v300"
    )

    def _remove_include_lines(text: str, needle: str) -> str:
        out = []
        for line in text.splitlines(keepends=True):
            s = line.strip()
            if s.startswith("#include") and needle in s:
                continue
            out.append(line)
        return "".join(out)

    def _remove_exact_include(text: str, target: str) -> str:
        out = []
        for line in text.splitlines(keepends=True):
            if line.strip() == target:
                continue
            out.append(line)
        return "".join(out)

    # ------------------------------------------------------------
    # 1) martinize2
    # ------------------------------------------------------------
    if "martinize" in steps:
        run_cmd(
            job_dir,
            [
                "martinize2",
                "-f", str(out_pdb),
                "-o", "Protein.top",
                "-x", "Protein_cg.pdb",
                "-name", "Protein",
                "-ff", "martini3001",
                "-p", "backbone",
                "-pf", "1000",
                "-maxwarn", "30",
            ],
            cwd=workdir,
        )

        # martinize2 emits Protein_0.itp when -name Protein is used
        p0 = workdir / "Protein_0.itp"
        p1 = workdir / "Protein.itp"
        if not p0.exists():
            raise FileNotFoundError(f"Missing martinize output: {p0}")

        # Rename Protein_0.itp -> Protein.itp
        p0.replace(p1)

        # Replace moleculetype Protein_0 -> Protein in Protein.itp
        txt = p1.read_text(encoding="utf-8", errors="replace")
        txt = txt.replace("Protein_0", "Protein")
        txt = _remove_include_lines(txt, "martini.itp")
        p1.write_text(txt, encoding="utf-8")

        # Patch Protein.top:
        # - include Protein.itp (not Protein_0.itp)
        # - replace Protein_0 -> Protein
        # - remove any martini.itp include
        pt = workdir / "Protein.top"
        if not pt.exists():
            raise FileNotFoundError(f"Missing martinize output: {pt}")

        t = pt.read_text(encoding="utf-8", errors="replace")
        t = t.replace("Protein_0.itp", "Protein.itp")
        t = t.replace("Protein_0", "Protein")
        t = _remove_include_lines(t, "martini.itp")
        pt.write_text(t, encoding="utf-8")

    # ------------------------------------------------------------
    # 2) insane
    # ------------------------------------------------------------
    if "insane" in steps:
        lipid_cfg = cfg.get("lipids", {"POPC": 1.0})
        lipids = ",".join(lipid_cfg.keys())

        run_cmd(
            job_dir,
            [
                "insane",
                "-f", "Protein_cg.pdb",
                "-o", "system.gro",
                "-p", "system.top",
                "-l", lipids,
                "-sol", "W",
                "-salt", "0.15",
                "-center",
                "-box", "12,12,10",
            ],
            cwd=workdir,
        )

        # Fix ions in insane outputs (per your working txt)
        gro = workdir / "system.gro"
        top = workdir / "system.top"

        if gro.exists():
            s = gro.read_text(encoding="utf-8", errors="replace")
            s = s.replace("CL-", "CL ")
            s = s.replace("NA+", "NA ")
            gro.write_text(s, encoding="utf-8")

        if top.exists():
            s = top.read_text(encoding="utf-8", errors="replace")
            s = s.replace("CL-", "CL")
            s = s.replace("NA+", "NA")
            top.write_text(s, encoding="utf-8")

    # ------------------------------------------------------------
    # 3) patch system.top
    # ------------------------------------------------------------
    if "patch_top" in steps:
        system_top = workdir / "system.top"
        if not system_top.exists():
            raise FileNotFoundError(
                f"Expected insane output missing: {system_top}"
            )

        orig = system_top.read_text(encoding="utf-8", errors="replace")

        # Remove legacy Martini 2 include and any include of Protein.top
        orig = _remove_include_lines(orig, "martini.itp")
        orig = _remove_exact_include(orig, '#include "Protein.top"')

        # Ensure names are consistent after rename
        orig = orig.replace("Protein_0.itp", "Protein.itp")
        orig = orig.replace("Protein_0", "Protein")

        includes = (
            f'#include "{martini_dir / "martini_v3.0.0.itp"}"\n'
            f'#include "{martini_dir / "martini_v3.0.0_phospholipids_v1.itp"}"\n'
            f'#include "{martini_dir / "martini_v3.0.0_solvents_v1.itp"}"\n'
            f'#include "{martini_dir / "martini_v3.0.0_ions_v1.itp"}"\n'
            '#include "Protein.itp"\n\n'
        )

        # Prepend includes once, but always persist cleaned content
        if "martini_v3.0.0.itp" not in orig:
            system_top.write_text(includes + orig, encoding="utf-8")
        else:
            system_top.write_text(orig, encoding="utf-8")

    # ------------------------------------------------------------
    # 4) equilibration / production
    # ------------------------------------------------------------



def build_solvated_system(job_dir: Path, workdir: Path, params: dict, cfg: dict):
    steps = cfg["steps"]

    workdir.mkdir(parents=True, exist_ok=True)

    # Always copy AA input into workdir
    in_pdb = job_dir / "input" / "protein.pdb"
    out_pdb = workdir / "protein.pdb"
    if not in_pdb.exists():
        raise FileNotFoundError(f"Missing input PDB: {in_pdb}")
    shutil.copyfile(in_pdb, out_pdb)

    # 1) martinize
    if "martinize" in steps:
        run_cmd(job_dir, [
            "martinize2",
            "-f", str(out_pdb),
            "-o", "Protein.top",
            "-x", "Protein_cg.pdb",
            "-name", "Protein",
            "-ff", "martini3001",
            "-p", "backbone",
            "-pf", "1000",
            "-maxwarn", "30",
        ], cwd=workdir)

    # 2) solvated box (NO membrane) via insane
    if "insane" in steps:
        # Here we build water + ions around the protein, but no lipids
        # insane supports no lipids by omitting -l (or setting -l to empty),
        # depending on version. The robust way is: keep -sol W and specify box.
        run_cmd(job_dir, [
            "insane",
            "-f", "Protein_cg.pdb",
            "-o", "system.gro",
            "-p", "system.top",
            "-sol", "W",
            "-salt", "0.15",
            "-center",
            "-box", "12,12,12",
        ], cwd=workdir)

    # 3) patch system.top includes (same as membrane case)
    if "patch_top" in steps:
        top = workdir / "system.top"
        if not top.exists():
            raise FileNotFoundError(f"Expected insane output missing: {top}")

        orig = top.read_text(encoding="utf-8")

        includes = (
            '#include "martini_v3.0.0.itp"\n'
            '#include "martini_v3.0.0_solvents_v1.itp"\n'
            '#include "martini_v3.0.0_ions_v1.itp"\n'
            '#include "Protein.top"\n\n'
        )

        if "martini_v3.0.0.itp" not in orig:
            top.write_text(includes + orig, encoding="utf-8")


# -------------------------------------------------------
# Equilibration (EM/NVT/NPT)
# -------------------------------------------------------

def run_equilibration(job_dir: Path, workdir: Path, params: dict, cfg: dict):
    steps = cfg["steps"]

    system_gro = workdir / "system.gro"
    system_top = workdir / "system.top"

    if any(s in steps for s in ("em", "nvt", "npt")):
        if not system_gro.exists():
            raise FileNotFoundError(
                f"Missing {system_gro}. Preset must include 'insane'."
            )
        if not system_top.exists():
            raise FileNotFoundError(
                f"Missing {system_top}. Preset must include 'patch_top'."
            )

    # Helper: choose threads if present
    gmx_cfg = params.get("gmx", {})
    ntmpi = str(gmx_cfg.get("ntmpi", 1))
    nt = str(gmx_cfg.get("nt", 4))

    # ---------------- EM ----------------
    if "em" in steps:
        run_cmd(job_dir, [
            "gmx", "grompp",
            "-f", str(MDP_DIR / "em.mdp"),
            "-c", "system.gro",
            "-p", "system.top",
            "-o", "em.tpr",
            "-maxwarn", "1",
        ], cwd=workdir)

        run_cmd(job_dir, [
            "gmx", "mdrun",
            "-ntmpi", ntmpi,
            "-nt", nt,
            "-v",
            "-deffnm", "em",
        ], cwd=workdir)

    # ---------------- NVT ----------------
    if "nvt" in steps:
        # input is em.gro if EM ran, otherwise system.gro
        nvt_in = "em.gro" if (workdir / "em.gro").exists() else "system.gro"

        run_cmd(job_dir, [
            "gmx", "grompp",
            "-f", str(MDP_DIR / "nvt.mdp"),
            "-c", nvt_in,
            "-r", nvt_in,
            "-p", "system.top",
            "-o", "nvt.tpr",
            "-maxwarn", "1",
        ], cwd=workdir)

        run_cmd(job_dir, [
            "caffeinate",
            "-dimsu",
            "gmx", "mdrun",
            "-ntmpi", ntmpi,
            "-nt", nt,
            "-v",
            "-deffnm", "nvt",
        ], cwd=workdir)

    # ---------------- NPT ----------------
    if "npt" in steps:
        # input is nvt.gro if NVT ran, otherwise em.gro/system.gro
        if (workdir / "nvt.gro").exists():
            npt_in = "nvt.gro"
        elif (workdir / "em.gro").exists():
            npt_in = "em.gro"
        else:
            npt_in = "system.gro"

        run_cmd(job_dir, [
            "gmx", "grompp",
            "-f", str(MDP_DIR / "npt.mdp"),
            "-c", npt_in,
            "-r", npt_in,
            "-p", "system.top",
            "-o", "npt.tpr",
            "-maxwarn", "1",
        ], cwd=workdir)

        run_cmd(job_dir, [
            "gmx", "mdrun",
            "-ntmpi", ntmpi,
            "-nt", nt,
            "-v",
            "-deffnm", "npt",
        ], cwd=workdir)

# -------------------------------------------------------
# Production
# -------------------------------------------------------

def run_production(job_dir: Path, workdir: Path, params: dict, cfg: dict):
    steps = cfg["steps"]
    if "md" not in steps:
        return

    system_top = workdir / "system.top"
    if not system_top.exists():
        raise FileNotFoundError(
            f"Missing {system_top}. Preset must include build + patch_top."
        )

    # Choose starting structure for production
    start_gro = params.get("start_gro")

    if not start_gro:
        raise ValueError(
            "start_gro not provided in params. "
            "MD runs must explicitly specify a starting structure."
        )

    start_gro_path = workdir / start_gro
    if not start_gro_path.exists():
        raise FileNotFoundError(
            f"Requested start_gro '{start_gro}' not found in {workdir}"
        )

    # Threads (optional)
    gmx_cfg = params.get("gmx", {})
    ntmpi = str(gmx_cfg.get("ntmpi", 1))
    nt = str(gmx_cfg.get("nt", 4))

   # Start from preset default
    md_ns = cfg.get("md_ns", 50)

    # Override from params if provided
    md_ns_param = params.get("md_ns")
    if md_ns_param is not None:
        md_ns = float(md_ns_param)

    # Safety
    if md_ns <= 0:
        raise ValueError("md_ns must be > 0")

    # Time step (ps)
    dt_ps = float(cfg.get("md_dt_ps", 0.01))

    total_ps = float(md_ns) * 1000.0
    nsteps = int(round(total_ps / dt_ps))

    md_mdp = workdir / "md.mdp"

    mdp_lines = [
        "; MARTINI 3 — PRODUCTION RUN (GPCR)",
        "integrator          = md",
        f"dt                  = {dt_ps}",
        f"nsteps              = {nsteps}",
        "",
        "; Velocity generation (replica diversification)",
        "gen_vel             = yes",
        "gen_temp            = 310",
        "gen_seed            = -1",
        "",
        "; Output control",
        "nstxout             = 0",
        "nstvout             = 0",
        "nstfout             = 0",
        "nstenergy           = 1000",
        "nstlog              = 1000",
        "nstxout-compressed  = 2000",
        "",
        "; Neighbour searching",
        "cutoff-scheme       = Verlet",
        "nstlist             = 20",
        "rlist               = 1.1",
        "rvdw                = 1.1",
        "rcoulomb            = 1.1",
        "",
        "; Electrostatics (Martini 3)",
        "coulomb-type        = reaction-field",
        "epsilon_rf          = 15",
        "",
        "; Temperature coupling",
        "tcoupl              = v-rescale",
        "tc-grps             = Protein POPC W ION",
        "tau-t               = 1.0 1.0 1.0 1.0",
        "ref-t               = 310 310 310 310",
        "",
        "; Pressure coupling",
        "pcoupl              = c-rescale",
        "pcoupltype          = semiisotropic",
        "tau-p               = 12.0",
        "ref-p               = 1.0 1.0",
        "compressibility     = 3e-4 3e-4",
        "",
        "refcoord_scaling    = all",
        "pbc                 = xyz",
        "constraints         = none",
        "",
    ]

    md_mdp.write_text("\n".join(mdp_lines), encoding="utf-8")

    # Build tpr
    run_cmd(job_dir, [
        "gmx", "grompp",
        "-f", "md.mdp",
        "-c", start_gro,
        "-p", "system.top",
        "-o", "md.tpr",
        "-maxwarn", "1",
    ], cwd=workdir)

    # Run production
    run_cmd(job_dir, [
        "gmx", "mdrun",
        "-ntmpi", ntmpi,
        "-nt", nt,
        "-v",
        "-deffnm", "md",
    ], cwd=workdir)

