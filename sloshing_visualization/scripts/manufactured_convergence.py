#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sloshing.validation.spatial import convergence_rows, manufactured_run
from sloshing.validation.artifacts import pyplot, save_figure, write_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("validation_results/manufactured"))
    args = parser.parse_args()
    rows = convergence_rows()
    half = manufactured_run(32, rows[-1]["dt"]/2)
    contamination = {key: abs(half[key]-rows[-1][key])/half[key]
                     for key in ("velocity_l2", "velocity_h1", "pressure_l2", "eta_l2")}
    print("dt/2 check", contamination, flush=True)
    write_tables(args.output, "spatial", rows, {"halved_dt": half, "relative_error_change": contamination})
    plt = pyplot()
    fig, ax = plt.subplots()
    for key in ("velocity_l2", "velocity_h1", "pressure_l2", "eta_l2", "divergence_l2"):
        ax.loglog([r["h"] for r in rows], [r[key] for r in rows], "o-", label=key)
    ax.set(xlabel="h (m)", ylabel="FEM quadrature error", title="Analytic MMS, uniform mesh")
    ax.legend()
    ax.grid(True)
    save_figure(fig, args.output, "spatial")


if __name__ == "__main__":
    main()
