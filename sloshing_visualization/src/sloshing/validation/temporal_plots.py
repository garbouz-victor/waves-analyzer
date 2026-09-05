"""Error plots reveal differences concealed by almost coincident trajectories."""

import json
from pathlib import Path

import h5py
import numpy as np

from .artifacts import pyplot, save_figure, write_tables


def internal_refinement(output):
    """Norm-curve differences at every internal step, distinct from field norms."""
    output = Path(output)
    rows = []
    for nu in (.01, .1, 1.):
        for method in ("midpoint", "sdirk2"):
            for dt in (.005, .0025):
                root = output/"runs"
                with h5py.File(root/f"nu{nu:g}_{method}_dt{dt:g}.h5") as a, h5py.File(root/f"nu{nu:g}_{method}_dt{dt/2:g}.h5") as b:
                    ta, tb = a["steps/time"][:], b["steps/time"][:]
                    row = {"nu": nu, "integrator": method, "dt": dt, "reference_dt": dt/2}
                    for key in ("omega_L2", "omega_wall_left_L2", "omega_wall_right_L2", "omega_surface_layer_L2", "omega_max"):
                        exact_times = np.interp(ta, tb, b[f"steps/{key}"][:])
                        delta = abs(a[f"steps/{key}"][:]-exact_times)
                        row[key+"_norm_curve_difference_max"] = float(np.max(delta))
                        row[key+"_norm_curve_relative_max"] = float(np.max(delta)/max(np.max(b[f"steps/{key}"][:]), 1e-30))
                    for component in ("u", "w"):
                        for i in range(7):
                            key = f"probe_{component}_{i}"
                            delta = abs(a[f"steps/{key}"][:]-np.interp(ta, tb, b[f"steps/{key}"][:]))
                            row[key+"_error_max"] = float(np.max(delta))
                    row["omega_probe_symmetry_max"] = float(max(np.max(abs(a["steps/omega_probe_0"][:]-a["steps/omega_probe_1"][:])),
                                                                   np.max(abs(a["steps/omega_probe_2"][:]-a["steps/omega_probe_3"][:]))))
                    corner_error = a["steps/omega_probe_2"][:]-np.interp(ta, tb, b["steps/omega_probe_2"][:])
                    early = (ta > 0) & (ta <= .05+1e-12)
                    # A smaller dt reveals additional startup times. Compare
                    # maxima on the SAME early grid for both refinement pairs.
                    shared = early & np.isclose(ta/.005, np.round(ta/.005), atol=1e-12, rtol=0)
                    row["corner_early_error_max"] = float(np.max(abs(corner_error[early])))
                    row["corner_early_error_on_common_0_005_grid"] = float(np.max(abs(corner_error[shared])))
                    rows.append(row)
    write_tables(output, "internal_refinement", rows,
                 {"warning": "These are differences OF norms at every dt; physical_refinement stores norms OF field differences at snapshots."})
    return rows


def temporal_error_plots(output, nus=(.01, .1, 1.)):
    output = Path(output)
    with (output/"physical_refinement.json").open() as handle:
        rows = json.load(handle)["rows"]
    plt = pyplot()
    for nu in nus:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for method in ("midpoint", "sdirk2"):
            group = [r for r in rows if r["nu"] == nu and r["integrator"] == method]
            for ax, key in zip(axes, ("eta_relative_max", "wall_left_relative_max")):
                ax.loglog([r["dt"] for r in group], [r[key] for r in group], "o-", label=method)
                ax.set(xlabel="dt (s)", ylabel=key)
                ax.grid(True)
            axes[0].legend()
        fig.suptitle(f"Errors versus dt/2, nu={nu:g}")
        fig.tight_layout()
        save_figure(fig, output, f"refinement_errors_nu{nu:g}")
        plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        for ax, method in zip(axes, ("midpoint", "sdirk2")):
            with h5py.File(output/"runs"/f"nu{nu:g}_{method}_dt0.00125.h5") as fine:
                tf = fine["steps/time"][:]
                wf = fine["steps/omega_probe_2"][:]
            for dt in (.005, .0025):
                with h5py.File(output/"runs"/f"nu{nu:g}_{method}_dt{dt:g}.h5") as h:
                    t = h["steps/time"][:]
                    e = h["steps/omega_probe_2"][:]-np.interp(t, tf, wf)
                    select = t <= .05+1e-12
                    ax.plot(t[select], e[select], "o-", ms=3, label=f"dt={dt:g}")
            ax.set(xlabel="t (s)", ylabel="omega error at (-0.995,-0.005) (1/s)", title=method)
            ax.legend()
            ax.grid(True)
        fig.suptitle(f"Internal-step errors versus dt=0.00125, nu={nu:g}")
        fig.tight_layout()
        save_figure(fig, output, f"internal_step_errors_nu{nu:g}")
        plt.close(fig)
