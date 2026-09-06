"""Scientific qualification figures, not production animation assets."""

import json
from pathlib import Path

import h5py
import numpy as np

from ..artifacts import pyplot, save_figure, write_tables
from ..contact import surface_value
from .dataset import validate_cache
from .fem_evaluation import FEMGeometry, evaluate_map
from .signals import horizontal_quadrature, modal_events


def history_plots(path, output):
    validate_cache(path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = pyplot()
    with h5py.File(path) as h:
        c = json.loads(h.attrs["config_json"])
        t = h["times"][:]
        d = h["diagnostics"]
        x = h["mesh/surface_x"][:]
        fig, ax = plt.subplots(figsize=(8, 4.8))
        query = np.unique(np.r_[np.linspace(-c["a"], c["a"], 1001), x])
        selected = [v for v in (0, .1, .25, .5, 1, 2, 3, 4, 5) if v <= t[-1]+1e-12]
        for when in selected:
            i = np.argmin(abs(t-when))
            ax.plot(query, surface_value(x, h["fields/eta"][i], query), label=f"t={t[i]:g}")
        ax.plot(x[[0, -1]], h["fields/eta"][0][[0, -1]], "ko", ms=4, label="pinned endpoints")
        ax.set(xlabel="x (m)", ylabel="eta (m)", title=f"Actual P2 surfaces: {c['mesh']}, alpha={c['alpha_deg']}°")
        ax.legend(ncol=3, fontsize=8)
        ax.grid(True)
        fig.tight_layout()
        save_figure(fig, output, "surface_profiles")
        plt.close(fig)
        specifications = {
            "slope_history": (["slope_global", "slope_bulk", "slope_intermediate", "slope_contact"], "max |eta_x| (dimensionless)"),
            "energy_history": (["kinetic_energy", "gravitational_potential_energy", "total_energy"], "energy / density / span (m^4/s^2)"),
            "vorticity_history": (["omega_L2", "omega_wall_left_L2", "omega_wall_right_L2", "omega_surface_layer_L2"], "omega L2 (m/s)"),
            "divergence_history": (["divergence_relative", "divergence_bulk_relative"], "||div v|| / ||grad v||")}
        for name, (keys, ylabel) in specifications.items():
            fig, ax = plt.subplots(figsize=(8, 4.5))
            for key in keys:
                ax.plot(t, d[key][:], label=key)
            if name == "slope_history":
                ax.axhline(.1, color="grey", ls="--", label="conservative policy 0.1")
                ax.axhline(.3, color="red", ls=":", label="hard policy 0.3")
            if name == "divergence_history":
                ax.axhline(.05, color="red", ls=":", label="particle caution 5%")
            ax.set(xlabel="t (s)", ylabel=ylabel, title=f"{c['mesh']}, actual alpha={c['alpha_deg']}° run")
            ax.legend(fontsize=8)
            ax.grid(True)
            fig.tight_layout()
            save_figure(fig, output, name)
            plt.close(fig)
        fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
        coords = h["probes/coordinates"][:]
        for ax, component in zip(axes, ("u", "w")):
            values = h[f"probes/{component}"][:]
            for i, point in enumerate(coords.T):
                ax.plot(t, values[:, i], label=f"({point[0]:g},{point[1]:g})")
            ax.set(ylabel=f"{component} (m/s)")
            ax.legend(ncol=4, fontsize=7)
            ax.grid(True)
        axes[-1].set_xlabel("t (s)")
        fig.tight_layout()
        save_figure(fig, output, "velocity_probes")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 4))
        for i, px in enumerate(h["probes/eta_coordinates"]):
            ax.plot(t, h["probes/eta"][:, i], label=f"x={px:g}")
        ax.set(xlabel="t (s)", ylabel="eta (m)")
        ax.legend()
        ax.grid(True)
        fig.tight_layout()
        save_figure(fig, output, "eta_probes")
        plt.close(fig)
        signal = d["modal_eta"][:]
        modal = modal_events(t, signal, c["a"], c["d"], c["g"])
        with (output/"modal_summary.json").open("w") as handle:
            json.dump(modal, handle, indent=2, allow_nan=False)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, signal, label="computed first antisymmetric projection")
        ax.axhline(0, color="grey", lw=.7)
        for e in modal["extrema"]:
            ax.plot(e["time"], e["amplitude"], "ko", ms=3)
        for zero in modal["zero_crossings_s"]:
            ax.axvline(zero, color="grey", ls=":", lw=.6)
        ax.set(xlabel="t (s)", ylabel="modal eta (m)", title=f"Measured period: {modal['estimated_period_s']} s")
        ax.legend()
        ax.grid(True)
        fig.tight_layout()
        save_figure(fig, output, "modal_signal")
        plt.close(fig)
        geometry = FEMGeometry.from_hdf5(h)
        depths = np.unique(np.r_[np.linspace(-c["d"], -1.5, 35), np.linspace(-1.5, 0, 61)])
        maps = []
        for z in depths:
            qp, weights = horizontal_quadrature(geometry, z)
            maps.append((geometry.interpolation(qp), weights))
        fig, ax = plt.subplots(figsize=(6, 6))
        depth_rows = []
        for when in selected[1:]:
            i = np.argmin(abs(t-when))
            u, w = h["fields/u"][i], h["fields/w"][i]
            norms = [float(np.sqrt(weights@(evaluate_map(u, mapping)**2+evaluate_map(w, mapping)**2))) for mapping, weights in maps]
            # Zero is not on a logarithmic axis. Do not draw a spurious
            # long line toward log(0); retain the exact zero in raw data.
            positive = np.asarray(norms) > 0
            ax.semilogx(np.asarray(norms)[positive], depths[positive], label=f"t={t[i]:g}")
            depth_rows.extend({"time": float(t[i]), "z": float(z), "velocity_x_L2": n} for z, n in zip(depths, norms))
        ax.set(xlabel="||v(.,z)|| L2(dx) (m^(3/2)/s)", ylabel="z (m)")
        ax.text(.03, .02, "No-slip bottom: zero norm is off the log axis", transform=ax.transAxes, fontsize=8)
        ax.legend()
        ax.grid(True)
        fig.tight_layout()
        save_figure(fig, output, "depth_penetration")
        plt.close(fig)
        write_tables(output/"raw", "depth_penetration", depth_rows)
    return modal


def comparison_plot(output, name="comparison_full"):
    output = Path(output)
    rows = json.loads((output/"raw"/(name+".json")).read_text())["rows"]
    plt = pyplot()
    fig, axes = plt.subplots(3, 2, figsize=(10, 10))
    for ax, keys in zip(axes.ravel(), (("eta_full_difference", "eta_bulk_difference", "eta_contact_difference"),
                                       ("velocity_difference",), ("omega_full_difference", "omega_wall_left_difference", "omega_surface_layer_difference"),
                                       ("divergence_relative_medium", "divergence_relative_fine"),
                                       ("max_slope_medium", "max_slope_fine"),
                                       ("divergence_bulk_relative_medium", "divergence_bulk_relative_fine"))):
        for key in keys:
            ax.plot([r["time"] for r in rows], [r[key] for r in rows], label=key)
        ax.set_xlabel("t (s)")
        ax.legend(fontsize=7)
        ax.grid(True)
    axes[0, 0].set_ylabel("eta L2 difference (m^(3/2))")
    axes[0, 1].set_ylabel("velocity L2 difference (m²/s)")
    axes[1, 0].set_ylabel("omega L2 difference (m/s)")
    axes[1, 1].set_ylabel("relative strong divergence")
    axes[2, 0].set_ylabel("global max |eta_x|")
    axes[2, 0].axhline(.1, color="grey", ls="--", lw=.8, label="conservative policy 0.1")
    axes[2, 0].axhline(.3, color="red", ls=":", lw=.8, label="hard policy 0.3")
    axes[2, 0].legend(fontsize=7)
    axes[2, 1].set_ylabel("bulk relative strong divergence")
    fig.tight_layout()
    save_figure(fig, output, "medium_vs_fine_short" if name == "comparison_short" else "medium_vs_fine")
    plt.close(fig)
