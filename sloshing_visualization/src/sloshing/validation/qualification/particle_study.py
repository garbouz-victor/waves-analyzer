"""Mesh, snapshot and ODE-step trajectory qualification, with invalidity checks."""

from dataclasses import replace
import json
from pathlib import Path
import time

import h5py
import numpy as np

from ...config import SimulationConfig
from ..artifacts import pyplot, save_figure, write_tables
from ..particles import FEMVelocityHistory, advect, compare_paths, particle_seeds
from .dataset import dataset_fingerprint, run_dataset, validate_cache


def _aggregate(rows):
    result = {}
    for group in ("bulk", "near_wall"):
        selected = [r for r in rows if r["group"] == group]
        result[group+"_separation_max_m"] = max(r["separation_max_m"] or 0. for r in selected)
        final = [r["separation_final_m"] for r in selected if r["separation_final_m"] is not None]
        result[group+"_separation_final_m"] = max(final) if final else None
        result[group+"_invalid_count"] = sum(r["invalid_a"] is not None or r["invalid_b"] is not None for r in selected)
    return result


def wall_limit(history, mesh):
    rows = []
    epsilon = np.array([.02, .01, .005, .001, .0001, .00001, 0.])
    for t in (.1, .5, 1., 2., 5.):
        if t > history.times[-1]+1e-12:
            continue
        for z in (-.05, -.2, -1.):
            for side in (-1, 1):
                points = np.vstack((side*(1-epsilon), np.full(len(epsilon), z)))
                speed = np.linalg.norm(history.velocity(t, points), axis=0)
                rows.extend({"mesh": mesh, "time": t, "z": z, "side": side,
                             "epsilon": float(e), "speed_m_per_s": float(v)} for e, v in zip(epsilon, speed))
    return rows


def particle_study(output):
    output = Path(output)
    start = time.perf_counter()
    medium, fine = output/"runs/medium-full.h5", output/"runs/fine-full.h5"
    config = SimulationConfig(**validate_cache(medium))
    validate_cache(fine)
    dense_path = output/"runs/sampling-dense.h5"
    dense_config = replace(config, t_end=.5, snapshot_dt=.0025)
    run_dataset(dense_config, dense_path)
    # Verify that this is sampling refinement, not a change to the underlying
    # PDE integrator, initial data, geometry or velocity states.
    same_field_max = 0.
    with h5py.File(medium) as a, h5py.File(dense_path) as b:
        for j, t in enumerate(b["times"]):
            if j % 2:
                continue
            i = int(round(t/config.snapshot_dt))
            for key in ("u", "w", "eta"):
                av, bv = a[f"fields/{key}"][i], b[f"fields/{key}"][j]
                np.testing.assert_allclose(av, bv, atol=2e-13, rtol=1e-9)
                same_field_max = max(same_field_max, float(np.max(abs(av-bv))))
    seeds, groups = particle_seeds()
    histories = [FEMVelocityHistory.open(path) for path in (medium, fine, dense_path)]
    try:
        hm, hf, hd = histories
        print("PARTICLES medium, 40 deterministic seeds", flush=True)
        pm = advect(hm, seeds)
        print("PARTICLES fine", flush=True)
        pf = advect(hf, seeds, times=pm["times"])
        mesh_rows = compare_paths(pm, pf, seeds, groups)
        short_times = pm["times"][pm["times"] <= .5+1e-12]
        short = advect(hm, seeds, times=short_times)
        print("PARTICLES dense snapshot temporal interpolation", flush=True)
        dense = advect(hd, seeds, times=short_times)
        sampling_rows = compare_paths(short, dense, seeds, groups)
        # Independent integration-step check, same piecewise-linear saved field.
        finer_rk = advect(hm, seeds, times=short_times, max_step=.00125)
        rk_rows = compare_paths(short, finer_rk, seeds, groups)
        wall_rows = wall_limit(hm, "medium")+wall_limit(hf, "fine")
    finally:
        for h in histories:
            h.close()
    summary = {"sources": {"medium": dataset_fingerprint(medium), "fine": dataset_fingerprint(fine)},
               "sampling_source": dataset_fingerprint(dense_path),
               "mesh": _aggregate(mesh_rows), "snapshot": _aggregate(sampling_rows),
               "rk_step": _aggregate(rk_rows), "bulk_seeds": int(np.sum(groups == "bulk")),
               "near_wall_seeds": int(np.sum(groups == "near_wall")), "integration": "RK4, max_step=0.0025 s",
               "spatial_interpolation": "actual element P2, no smoothing",
               "temporal_interpolation": "linear between saved FEM states; RK stops at snapshot knots",
               "first_invalid_time_resolution_s": .0025, "sampling_test_end_s": .5,
               "same_underlying_snapshot_field_max_difference": same_field_max,
               "wall_velocity_max_m_per_s": max(r["speed_m_per_s"] for r in wall_rows if r["epsilon"] == 0.),
               "runtime_seconds": time.perf_counter()-start}
    for group in ("bulk", "near_wall"):
        selected = groups == group
        displacement = np.linalg.norm(pf["paths"][:, selected]-seeds[selected], axis=2)
        peak = float(np.nanmax(displacement))
        summary["mesh"][group+"_reference_displacement_max_m"] = peak
        summary["mesh"][group+"_relative_to_peak_displacement"] = summary["mesh"][group+"_separation_max_m"]/max(peak, 1e-30)
    write_tables(output, "particle_mesh", mesh_rows, summary["mesh"])
    write_tables(output, "particle_sampling", sampling_rows, summary["snapshot"])
    write_tables(output, "particle_rk_step", rk_rows, summary["rk_step"])
    wall_compact = []
    for mesh, t, z, side in sorted({(r["mesh"], r["time"], r["z"], r["side"]) for r in wall_rows}):
        group = [r for r in wall_rows if (r["mesh"], r["time"], r["z"], r["side"]) == (mesh, t, z, side)]
        speeds = {r["epsilon"]: r["speed_m_per_s"] for r in group}
        wall_compact.append({"mesh": mesh, "time": t, "z": z, "side": side,
                             "wall_speed": speeds[0.], "speed_at_epsilon_0_02": speeds[.02],
                             "speed_at_epsilon_1e_5": speeds[.00001],
                             "ratio_smallest_to_largest_epsilon": speeds[.00001]/max(speeds[.02], 1e-30)})
    write_tables(output, "wall_velocity_limit", wall_compact, {"wall_velocity_max_m_per_s": summary["wall_velocity_max_m_per_s"]})
    write_tables(output/"raw", "wall_velocity_limit", wall_rows)
    with (output/"particle_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    np.savez_compressed(output/"particle_paths.npz", times=pm["times"], seeds=seeds, groups=groups,
                        medium=pm["paths"], fine=pf["paths"], first_invalid_medium=pm["first_invalid_time"],
                        first_invalid_fine=pf["first_invalid_time"], short_times=short_times,
                        sampling_coarse=short["paths"], sampling_dense=dense["paths"])
    plt = pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for group in ("bulk", "near_wall"):
        selected = groups == group
        mesh_error = np.linalg.norm(pm["paths"][:, selected]-pf["paths"][:, selected], axis=2)
        sampling_error = np.linalg.norm(short["paths"][:, selected]-dense["paths"][:, selected], axis=2)
        axes[0].plot(pm["times"], np.nanmax(mesh_error, axis=1), label=group)
        axes[1].plot(short_times, np.nanmax(sampling_error, axis=1), label=group)
    for ax, title in zip(axes, ("Medium vs fine", "Snapshots 0.005 vs 0.0025 s")):
        ax.set(xlabel="t (s)", ylabel="max particle separation (m)", title=title)
        ax.legend()
        ax.grid(True)
    fig.tight_layout()
    save_figure(fig, output, "particle_qualification")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4))
    for mesh in ("medium", "fine"):
        for t in (.5, 1., 5.):
            rows = [r for r in wall_rows if r["mesh"] == mesh and r["time"] == t and r["side"] == -1 and r["z"] == -.05 and r["epsilon"] > 0]
            ax.loglog([r["epsilon"] for r in rows], [r["speed_m_per_s"] for r in rows], "o-", ms=3, label=f"{mesh}, t={t:g}")
    ax.set(xlabel="distance from left wall (m)", ylabel="|v| (m/s)", title="No fitted power law; exact wall value recorded separately")
    ax.legend(fontsize=8)
    ax.grid(True)
    fig.tight_layout()
    save_figure(fig, output, "wall_velocity_limit")
    plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)
    return summary
