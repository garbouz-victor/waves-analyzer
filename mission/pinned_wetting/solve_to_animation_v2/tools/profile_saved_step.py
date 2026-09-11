"""Read-only one-step replay and direct-ordering diagnostic, never a native run."""
import argparse
import cProfile
from dataclasses import replace
import io
import json
from pathlib import Path
import pstats
import time

import h5py
import numpy as np
from scipy.sparse.linalg import splu

from sloshing.pinned_wetting import free_boundary as fb
from sloshing.pinned_wetting.mission import job, read_json, sha256, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native", type=Path)
    parser.add_argument("--step", type=int)
    parser.add_argument("--orderings", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    output = root / "sloshing_visualization/output/pinned_wetting/free_boundary_fixed_conditions"
    import yaml
    cfg = yaml.safe_load((root / "mission/pinned_wetting/CONFIG.yaml").read_text())
    with job(root, "PW2 read-only saved-step cost diagnostic", cfg["resources"]) as (_, heartbeat):
        before = sha256(args.native)
        with h5py.File(args.native, "r") as h:
            n = args.step or int(h.attrs["last_accepted_step"])
            if n < 1:
                raise ValueError("A noninitial accepted state is required")
            controls = fb.Controls(**json.loads(h.attrs["numerical_identity"])["controls"])
            current, previous = h["states"][f"{n:08d}"], h["states"][f"{n-1:08d}"]
            row = json.loads(current.attrs["diagnostics"])
            if "rezoning_before_step" in current:
                transfer = current["rezoning_before_step"]
                print("saved transfer datasets:", list(transfer), flush=True)
                # These are the committed starting fields, not a new transfer.
                coordinates, velocity = transfer["geometry"][:], transfer["velocity"][:]
            else:
                coordinates, velocity = previous["geometry"][:], previous["velocity"][:]
            expected = [current[key][:] for key in ("geometry", "velocity", "pressure_midpoint")]
            mesh = replace(fb.initial_mesh(controls), doflocs=coordinates)
            orientation = h["topology/orientation"][:]
        capture = {}
        original = fb.PicardLinearSolver.solve
        def observed_solve(self, matrix, rhs, velocity_size):
            capture.update(matrix=matrix.copy(), rhs=rhs.copy())
            return original(self, matrix, rhs, velocity_size)
        profile = cProfile.Profile()
        started = time.monotonic()
        try:
            fb.PicardLinearSolver.solve = observed_solve
            answer = profile.runcall(fb.step, mesh, velocity, row["dt_s"], controls, orientation)
        finally:
            fb.PicardLinearSolver.solve = original
        replay_wall = time.monotonic() - started
        stream = io.StringIO()
        pstats.Stats(profile, stream=stream).strip_dirs().sort_stats("cumulative").print_stats(35)
        print(stream.getvalue(), flush=True)
        output.mkdir(parents=True, exist_ok=True)
        stem = f"profile_{args.native.parent.name}_{n}_{sha256(Path(fb.__file__))[:8]}"
        profile.dump_stats(str(output / (stem + ".prof")))
        report = {"scope": "READ_ONLY_REPLAY_NOT_NEW_NATIVE_DATA", "native": str(args.native.resolve()),
                  "step": n, "replay_wall_s": replay_wall, "profile": stream.getvalue(),
                  "replay_max_abs_differences": {key: float(np.max(np.abs(actual-want)))
                      for key, actual, want in zip(("geometry", "velocity", "pressure"),
                                                   (answer[0].p, answer[1], answer[2]), expected)},
                  "diagnostics": answer[3], "orderings": {}}
        A, rhs = capture["matrix"], capture["rhs"]
        baseline = fb.spsolve(A, rhs) if args.orderings else None
        for ordering in (("COLAMD", "MMD_AT_PLUS_A", "MMD_ATA") if args.orderings else ()):
            started = time.monotonic()
            factor = splu(A, permc_spec=ordering)
            solution = factor.solve(rhs)
            report["orderings"][ordering] = {"wall_s": time.monotonic()-started,
                "factor_nnz": factor.L.nnz+factor.U.nnz,
                "max_abs_residual": float(np.max(np.abs(A@solution-rhs))),
                "max_abs_solution_difference": float(np.max(np.abs(solution-baseline)))}
            print(ordering, report["orderings"][ordering], flush=True)
            heartbeat()
        report.update(native_sha256=before, native_unchanged=before == sha256(args.native),
                      replay_source_sha256=sha256(Path(fb.__file__)))
        write_json(output/(stem+".json"), report)
        print(json.dumps({k:v for k,v in report.items() if k != "profile"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
