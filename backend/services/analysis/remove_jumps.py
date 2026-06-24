#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def run_gmx(cmd, selection_text, cwd):
    print("\n$ " + " ".join(str(x) for x in cmd))
    result = subprocess.run(
        [str(x) for x in cmd],
        input=selection_text,
        text=True,
        capture_output=True,
        cwd=str(cwd),
    )

    print(result.stdout)
    if result.stderr:
        print(result.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"GROMACS command failed with exit code {result.returncode}"
        )


def generate_analysis_trajectories(job_dir: Path, gmx: str = "gmx"):
    out_dir = job_dir / "out"
    if not out_dir.exists():
        raise RuntimeError(f"Missing out/ directory: {out_dir}")

    tpr = out_dir / "md.tpr"
    if not tpr.exists():
        raise RuntimeError(f"Missing md.tpr: {tpr}")

    xtc = out_dir / "md.xtc"
    if not xtc.exists() or xtc.stat().st_size == 0:
        parts = sorted(out_dir.glob("md.part*.xtc"))
        if not parts:
            raise RuntimeError(
                f"No md.xtc or md.part*.xtc found in {out_dir}"
            )
        xtc = parts[-1]

    print(f"[INFO] Using input trajectory: {xtc.name}")

    whole_xtc = out_dir / "md_whole.xtc"
    nojump_xtc = out_dir / "md_nojump.xtc"
    centered_xtc = out_dir / "md_centered.xtc"
    fit_xtc = out_dir / "md_fit.xtc"

    run_gmx(
        [
            gmx, "trjconv",
            "-s", tpr,
            "-f", xtc,
            "-o", whole_xtc,
            "-pbc", "whole",
        ],
        selection_text="0\n",
        cwd=out_dir,
    )

    run_gmx(
        [
            gmx, "trjconv",
            "-s", tpr,
            "-f", whole_xtc,
            "-o", nojump_xtc,
            "-pbc", "nojump",
        ],
        selection_text="0\n",
        cwd=out_dir,
    )

    run_gmx(
        [
            gmx, "trjconv",
            "-s", tpr,
            "-f", nojump_xtc,
            "-o", centered_xtc,
            "-center",
            "-pbc", "mol",
            "-ur", "compact",
        ],
        selection_text="1\n0\n",
        cwd=out_dir,
    )

    run_gmx(
        [
            gmx, "trjconv",
            "-s", tpr,
            "-f", centered_xtc,
            "-o", fit_xtc,
            "-fit", "rot+trans",
        ],
        selection_text="1\n0\n",
        cwd=out_dir,
    )

    print("\n[OK] Generated:")
    print(f"  {whole_xtc.name}")
    print(f"  {nojump_xtc.name}")
    print(f"  {centered_xtc.name}")
    print(f"  {fit_xtc.name}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "job_dir",
        type=Path,
        help="Path to the job directory containing out/",
    )
    parser.add_argument(
        "--gmx",
        default="gmx",
        help="GROMACS executable (default: gmx)",
    )
    args = parser.parse_args()

    generate_analysis_trajectories(args.job_dir, gmx=args.gmx)


if __name__ == "__main__":
    main()