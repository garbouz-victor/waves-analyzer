#!/usr/bin/env python3
"""Minimal PW1 CLI. Reference output NEVER qualifies as a target result."""
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "sloshing_visualization/src"))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pinned-wetting-matplotlib")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

from sloshing.pinned_wetting.mission import main

if __name__ == "__main__":
    sys.exit(main(ROOT))
