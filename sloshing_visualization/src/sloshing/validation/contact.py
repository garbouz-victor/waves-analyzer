"""Bulk/contact comparison of unsmoothed piecewise-quadratic surface fields."""

from dataclasses import asdict
import gc
import json
from pathlib import Path

import h5py
import numpy as np

from ..config import SimulationConfig
from ..fem_spaces import FEMSystem
from ..postprocess.vorticity import VorticityProjector
from .artifacts import pyplot, save_figure, write_tables
from .physical import run_history
from .regions import VorticityRegions


def surface_value(x, eta, query):
    query = np.asarray(query)
    edges = x[::2]
    if np.any(query < edges[0]) or np.any(query > edges[-1]):
        raise ValueError("Surface query outside the physical domain")
    i = np.searchsorted(edges, query, side="right")-1
    i = np.where(query == edges[-1], len(edges)-2, i)
    s = (query-edges[i])/(edges[i+1]-edges[i])
    return (2*s*s-3*s+1)*eta[2*i]+4*s*(1-s)*eta[2*i+1]+(2*s*s-s)*eta[2*i+2]


def surface_difference(xa, ea, xb, eb, left=-1., right=1.):
    # Split at BOTH sets of P2 element endpoints. A squared difference then has
    # degree <=4 on each interval, so three Gauss points integrate it exactly.
    interior = np.r_[xa[::2], xb[::2]]
    edges = np.unique(np.r_[left, interior[(interior > left) & (interior < right)], right])
    nodes, weights = np.polynomial.legendre.leggauss(3)
    half = np.diff(edges)/2
    points = (edges[1:]+edges[:-1])[:, None]/2+half[:, None]*nodes
    error = surface_value(xa, ea, points)-surface_value(xb, eb, points)
    return float(np.sqrt(np.sum(half[:, None]*weights*error**2)))


def surface_regions(xa, ea, xb, eb):
    return {"full": surface_difference(xa, ea, xb, eb),
            "bulk": surface_difference(xa, ea, xb, eb, -.9, .9),
            "intermediate": surface_difference(xa, ea, xb, eb, -.98, .98),
            "contact": float(np.hypot(surface_difference(xa, ea, xb, eb, -1, -.98),
                                       surface_difference(xa, ea, xb, eb, .98, 1)))}


def _compatible_prefix(path, config):
    if not path.exists():
        return False
    with h5py.File(path) as h:
        stored = json.loads(h.attrs["config_json"])
        same = all(stored[k] == value for k, value in asdict(config).items()
                   if k not in ("t_end", "snapshot_dt"))
        return (h.attrs["status"] == "complete" and same
                and bool(np.any(np.isclose(h["times"][:], config.t_end, atol=1e-14, rtol=0))))


def contact_mesh_study(output, time_output=Path("validation_results/time")):
    output = Path(output)
    profiles, rows = {}, []
    for mesh in ("coarse", "medium", "fine"):
        config = SimulationConfig(alpha_deg=.2, nu=.01, mesh=mesh, t_end=.1,
                                  dt=.00125, snapshot_dt=.025, integrator="sdirk2")
        path = output/"runs"/f"{mesh}.h5"
        existing = Path(time_output)/"runs"/"nu0.01_sdirk2_dt0.00125.h5"
        if mesh == "medium" and _compatible_prefix(existing, config):
            path = existing
            print("reuse medium prefix at t=0.1", path, flush=True)
        else:
            f = FEMSystem(config)
            projector = VorticityProjector(f)
            regions = VorticityRegions(f, projector)
            run_history(f, config, path, projector, regions)
            del f, projector, regions
            gc.collect()
        with h5py.File(path) as h:
            i = int(np.flatnonzero(np.isclose(h["times"][:], .1, atol=1e-14, rtol=0))[0])
            x, eta = h["surface_x"][:], h["eta"][i]
            profiles[mesh] = (x, eta)
            row = {"mesh": mesh, "nx": config.resolution[0], "nz": config.resolution[1],
                   "time": .1, "dt": config.dt, "integrator": config.integrator,
                   "source": str(path), "eta_antisymmetry_max": float(np.max(abs(eta+eta[::-1])))}
            for key in ("max_slope", "divergence_l2", "weak_divergence_l2", "volume_error", "contact_error",
                        "energy_balance_relative", "omega_L2", "omega_wall_left_L2", "omega_wall_right_L2",
                        "omega_max", "total_energy"):
                row[key] = float(h[f"diagnostics/{key}"][i])
            row["no_slip_max"] = max(float(h[f"diagnostics/{wall}_{c}_max"][i])
                                     for wall in ("left", "right", "bottom") for c in ("u", "w"))
            rows.append(row)
        print("contact", row, flush=True)
    for row in rows:
        differences = surface_regions(*profiles[row["mesh"]], *profiles["fine"])
        row.update({key+"_eta_difference_to_fine": value for key, value in differences.items()})
    consecutive = []
    for a, b in (("coarse", "medium"), ("medium", "fine")):
        consecutive.append({"mesh": a, "reference": b, **surface_regions(*profiles[a], *profiles[b])})
    # Existing medium time-refinement snapshots quantify temporal contamination
    # at the exact time used in this spatial study.
    time_check = {}
    time_coarse = Path(time_output)/"runs"/"nu0.01_sdirk2_dt0.0025.h5"
    if time_coarse.exists():
        with h5py.File(time_coarse) as h:
            i = int(np.flatnonzero(np.isclose(h["times"][:], .1, atol=1e-14, rtol=0))[0])
            time_check = surface_regions(h["surface_x"][:], h["eta"][i], *profiles["medium"])
    write_tables(output, "contact_mesh", rows, {"reference_is_fine_not_exact": True,
                 "consecutive_differences": consecutive, "medium_dt_vs_dt_half_at_t_0_1": time_check,
                 "regions": {"bulk": "|x|<=0.9", "intermediate": "|x|<=0.98", "contact": "0.98<|x|<=1"}})
    print("consecutive", consecutive, "time check", time_check, flush=True)
    plt = pyplot()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for mesh, (x, eta) in profiles.items():
        queries = np.unique(np.r_[np.linspace(-1, 1, 501), x])
        axes[0].plot(queries, surface_value(x, eta, queries), label=mesh)
        axes[1].plot(queries, surface_value(x, eta, queries)-surface_value(*profiles["fine"], queries), label=mesh)
        near = np.unique(np.r_[np.linspace(.95, 1, 201), x[x >= .95]])
        axes[2].plot(near, surface_value(x, eta, near), label=mesh)
    axes[0].set(xlabel="x (m)", ylabel="eta (m)", title="Unsmoothed P2 surface, t=0.1 s")
    axes[1].set(xlabel="x (m)", ylabel="eta - eta_fine (m)", title="Difference to fine (not exact)")
    axes[2].set(xlabel="x (m)", ylabel="eta (m)", title="Right contact region")
    axes[2].plot(1, profiles["fine"][1][-1], "ko", ms=4)
    for ax in axes:
        ax.legend()
        ax.grid(True)
    fig.tight_layout()
    save_figure(fig, output, "contact_profiles")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for region in ("full", "bulk", "intermediate", "contact"):
        axes[0].loglog([r["nx"] for r in rows[:-1]], [r[region+"_eta_difference_to_fine"] for r in rows[:-1]],
                       "o-", label=region)
    axes[0].set(xlabel="Nx", ylabel="L2 difference to fine (m^(3/2))")
    axes[0].legend()
    axes[1].plot([r["nx"] for r in rows], [r["max_slope"] for r in rows], "o-")
    axes[1].set(xlabel="Nx", ylabel="max |eta_x|")
    for ax in axes:
        ax.grid(True)
    fig.tight_layout()
    save_figure(fig, output, "contact_sensitivity")
    plt.close(fig)
    return rows
