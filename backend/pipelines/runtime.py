# backend/md/runtime.py

from typing import Dict
import subprocess

# job_id -> subprocess.Popen
RUNNING_JOBS: Dict[str, subprocess.Popen] = {}
