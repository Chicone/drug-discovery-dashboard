# backend/pipelines/build_and_equilibrate.py

import argparse
from pathlib import Path

import backend.pipelines.presets as presets_module

from backend.pipelines.core import (
    load_params,
    build_membrane_system,
    run_equilibration,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", required=True)
    args = parser.parse_args()

    job_dir = Path(args.job_dir)
    workdir = job_dir / "out"

    params = load_params(job_dir)

    preset = params.get("preset")
    if not preset:
        raise ValueError("params.json missing required key: 'preset'")

    PRESETS = presets_module.PRESETS
    cfg = PRESETS.get(preset)
    if cfg is None:
        raise ValueError(
            f"Unknown preset: {preset}. Allowed: {', '.join(sorted(PRESETS))}"
        )

    # 1) Build prerequisites (martinize/insane/patch_top/fix_ions) as defined
    build_membrane_system(job_dir, workdir, params, cfg)

    # 2) Run EM/NVT/NPT as defined
    run_equilibration(job_dir, workdir, params, cfg)


if __name__ == "__main__":
    main()


