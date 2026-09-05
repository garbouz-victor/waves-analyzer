#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sloshing.validation.artifacts import pyplot, save_figure, write_tables
from sloshing.validation.stiffness import scalar_rows, scalar_time_rows, stiff_fem_comparison


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("validation_results/time"))
    parser.add_argument("--unit-only", action="store_true")
    parser.add_argument("--physical-only", action="store_true")
    parser.add_argument("--workers", type=int, choices=(1, 2, 3), default=1)
    args = parser.parse_args()
    if args.physical_only:
        from sloshing.validation.physical import physical_time_study
        physical_time_study(args.output, workers=args.workers)
        return
    for name, rows in (("scalar_stiff", scalar_rows()), ("scalar_time", scalar_time_rows()),
                       ("fem_stiff", stiff_fem_comparison())):
        write_tables(args.output, name, rows)
        for row in rows:
            print(name, row, flush=True)
        plt = pyplot()
        fig, ax = plt.subplots()
        if name == "scalar_stiff":
            for method in ("midpoint", "sdirk2", "exact"):
                ax.semilogx([-r["lambda_dt"] for r in rows], [r[method] for r in rows], "o-", label=method)
            ax.set(xlabel="|lambda dt|", ylabel="amplification")
        else:
            for method in ("midpoint", "sdirk2"):
                group = [r for r in rows if r["integrator"] == method]
                if name == "scalar_time":
                    ax.loglog([r["dt"] for r in group], [r["error"] for r in group], "o-", label=method)
                    ax.set(xlabel="dt (s)", ylabel="absolute error at t=1")
                else:
                    ax.plot([r["step"] for r in group], [r["amplitude"] for r in group], "o-", label=method)
                    ax.set(xlabel="step", ylabel="signed FEM mode amplitude")
        ax.legend()
        ax.grid(True)
        save_figure(fig, args.output, name)
        plt.close(fig)
    if not args.unit_only:
        from sloshing.validation.physical import physical_time_study
        physical_time_study(args.output, workers=args.workers)


if __name__ == "__main__":
    main()
