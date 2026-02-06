# backend/md/pipelines/production.py

import argparse
from pathlib import Path

from .core import load_params, run_production
from . import presets as presets_module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", required=True)
    args = parser.parse_args()

    job_dir = Path(args.job_dir)
    workdir = job_dir / "out"
    workdir.mkdir(parents=True, exist_ok=True)

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

    # 🚀 Production MD only
    run_production(job_dir, workdir, params, cfg)


if __name__ == "__main__":
    main()
