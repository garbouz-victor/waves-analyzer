#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sloshing.validation.contact import contact_mesh_study

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("validation_results/contact"))
    parser.add_argument("--time-output", type=Path, default=Path("validation_results/time"))
    args = parser.parse_args()
    contact_mesh_study(args.output, args.time_output)
