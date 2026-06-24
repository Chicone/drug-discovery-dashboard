# backend/md/pipelines/build_system.py

import argparse
from pathlib import Path

from .core import load_params, build_membrane_system, build_solvated_system
from .presets import PRESETS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", required=True)
    args = parser.parse_args()

    job_dir = Path(args.job_dir)
    workdir = job_dir / "out"
    workdir.mkdir(parents=True, exist_ok=True)

    params = load_params(job_dir)

    preset = params.get("preset", "m3_popc_build")
    cfg = PRESETS.get(preset)
    if cfg is None:
        raise ValueError(
            f"Unknown preset: {preset}. Allowed: {', '.join(sorted(PRESETS))}"
        )

    env = params["environment"]

    if env == "membrane":
        build_membrane_system(job_dir, workdir, params, cfg)
    else:
        build_solvated_system(job_dir, workdir, params, cfg)


if __name__ == "__main__":
    main()

