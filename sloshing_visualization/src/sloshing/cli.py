"""CLI for a reproducible solve; rendering will be a separate stage."""

import argparse
from dataclasses import asdict
from pathlib import Path
import time

from .config import MESH_PRESETS, SimulationConfig, VISCOSITY_PRESETS
from .diagnostics import PINNED_MESSAGE
from .solver import SloshingSolver
from .storage import save_run


def main(argv=None):
    parser = argparse.ArgumentParser(description="Viscous linear sloshing with strict no-slip and pinned contacts")
    parser.add_argument("--config", type=Path, help="JSON config; explicit CLI arguments override its values")
    for name in ("a", "d", "g", "alpha-deg", "nu", "t-end", "dt", "snapshot-dt", "x-grading", "z-grading"):
        parser.add_argument(f"--{name}", type=float)
    parser.add_argument("--nu-preset", type=float, choices=VISCOSITY_PRESETS)
    parser.add_argument("--mesh", choices=tuple(MESH_PRESETS))
    parser.add_argument("--integrator", choices=("midpoint", "sdirk2"))
    for name in ("nx", "nz", "visualization-nx", "visualization-nz"):
        parser.add_argument(f"--{name}", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-visualization-grid", action="store_true")
    args = parser.parse_args(argv)
    base = SimulationConfig.from_json(args.config) if args.config else SimulationConfig()
    values = asdict(base)
    values.update({key: value for key, value in vars(args).items() if key in values and value is not None})
    if args.nu_preset is not None:
        if args.nu is not None:
            parser.error("Use either --nu or --nu-preset")
        values["nu"] = args.nu_preset
    config = SimulationConfig(**values)
    output = args.output or Path("results") / (config.run_name + ".h5")
    if output.exists() or (output.with_suffix("") / "diagnostics.csv").exists():
        raise FileExistsError(f"Run already exists: {output}; choose --output")
    print(PINNED_MESSAGE, flush=True)
    print(f"nu={config.nu:g} m²/s, alpha={config.alpha_deg:g}°, mesh={config.resolution}, "
          f"dt={config.dt:g} s, end={config.t_end:g} s", flush=True)
    start = time.perf_counter()
    solver = SloshingSolver(config)
    print(f"FEM: {solver.fem.mesh.nelements} triangles, "
          f"{solver.fem.velocity.N} velocity and {solver.fem.pressure.N} pressure dofs; "
          f"assembly/factorization {time.perf_counter()-start:.2f} s", flush=True)

    def progress(state, row):
        print(f"t={state.time:6.3f} E={row['total_energy']:.8e} "
              f"div_L2={row['divergence_l2']:.3e} weak_div={row['weak_divergence_l2']:.3e} "
              f"volume_error={row['volume_error']:.3e} slope={row['max_slope']:.3g} "
              f"energy_budget={row['energy_balance_relative']:.3e}", flush=True)

    save_run(solver, output, visualization_grid=not args.no_visualization_grid, progress=progress)
    print(f"Saved {output} and {output.with_suffix('') / 'diagnostics.csv'} "
          f"in {time.perf_counter()-start:.2f} s", flush=True)
    return 0
