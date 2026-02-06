# backend/pipelines/full.py

import argparse
from pathlib import Path

from .core import (
    load_params,
    build_membrane_system,
    build_solvated_system,
    run_equilibration,
    run_production,
)
from .presets import PRESETS


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", required=True)
    args = parser.parse_args()

    job_dir = Path(args.job_dir)
    workdir = job_dir / "out"

    params = load_params(job_dir)

    preset = params.get("preset")
    if not preset:
        raise ValueError("params.json missing 'preset'")

    cfg = PRESETS.get(preset)
    if cfg is None:
        raise ValueError(
            f"Unknown preset: {preset}. Allowed: {', '.join(PRESETS.keys())}"
        )

    # 1) build system
    if cfg.get("environment", params.get("environment")) == "membrane":
        build_membrane_system(job_dir, workdir, params, cfg)
    else:
        build_solvated_system(job_dir, workdir, params, cfg)

    # 2) equilibration (only runs steps present in cfg["steps"])
    run_equilibration(job_dir, workdir, params, cfg)

    # 3) production (only runs if "md" in steps)
    run_production(job_dir, workdir, params, cfg)


if __name__ == "__main__":
    main()

