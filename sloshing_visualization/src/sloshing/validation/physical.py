"""Production-geometry temporal refinement; no animation and no smoothing."""

from dataclasses import replace
import gc
import json
from pathlib import Path
import time

import h5py
import numpy as np

from ..config import SimulationConfig
from ..diagnostics import Diagnostics
from ..fem_spaces import FEMSystem
from ..postprocess.vorticity import VorticityProjector
from ..time_integrator import INTEGRATORS
from .artifacts import pyplot, save_figure, write_tables
from .regions import VorticityRegions


PROBES = np.array([[0, -.75, .75, 0, 0, -.95, .95], [-.1, -.2, -.2, -1, -3, -.05, -.05]])
OMEGA_PROBES = np.array([[-.95, .95, -.995, .995], [-.05, -.05, -.005, -.005]])


def run_history(fem, config, path, projector, regions):
    """Save FEM coefficients at snapshot_dt; step diagnostics are sampled at dt."""
    path = Path(path)
    if path.exists():
        try:
            with h5py.File(path) as h:
                if h.attrs.get("status") == "complete" and h.attrs["config_json"] == config.to_json():
                    print("reuse", path, flush=True)
                    return
        except OSError as exc:
            raise RuntimeError(f"Unreadable validation cache: {path}; choose another output directory") from exc
        raise RuntimeError(f"Incomplete or incompatible validation cache: {path}; choose another output directory")
    fem.config = config
    start = time.perf_counter()
    it = INTEGRATORS[config.integrator](fem)
    state = it.initial_state()
    diagnostics = Diagnostics(fem, state)
    vp = fem.scalar.probes(PROBES)
    op = projector.basis.probes(OMEGA_PROBES)
    histories = {key: [] for key in ("times", "velocity", "eta", "omega", "probes_u", "probes_w")}
    rows, steps = [], []
    vertices = fem.mesh.p[:, fem.mesh.t]
    edges = [vertices[:, i]-vertices[:, j] for i, j in ((0, 1), (1, 2), (2, 0))]
    h_min = min(float(np.min(np.linalg.norm(e, axis=0))) for e in edges)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "x") as h:
        h.attrs.update(status="running", config_json=config.to_json(), h_min=h_min,
                       rough_stiffness=config.nu*config.dt/h_min**2,
                       schema="validation-state-history-1.5")
        h.flush()
        try:
            for step in range(config.nsteps+1):
                if step:
                    state = it.advance(state)
                full = fem.expand_velocity(state.velocity)
                omega = projector.evaluate(full)
                norms = regions.norms(omega)
                probes_u, probes_w = [vp@full[ids] for ids in fem.component_dofs]
                omega_probes = op@omega
                metrics = {"time": state.time, "omega_L2": norms["full"],
                           "omega_wall_left_L2": norms["wall_left"],
                           "omega_wall_right_L2": norms["wall_right"],
                           "omega_surface_layer_L2": norms["surface_layer"],
                           "omega_max": float(np.max(abs(omega))),
                           "omega_wall_symmetry": abs(norms["wall_left"]-norms["wall_right"]),
                           **{f"omega_probe_{i}": float(v) for i, v in enumerate(omega_probes)},
                           **{f"probe_u_{i}": float(v) for i, v in enumerate(probes_u)},
                           **{f"probe_w_{i}": float(v) for i, v in enumerate(probes_w)}}
                steps.append(metrics)
                if step % config.save_every == 0 or step == config.nsteps:
                    row = diagnostics.measure(state)
                    diagnostics.validate(row)
                    rows.append({**row, **metrics})
                    for key, value in (("times", state.time), ("velocity", state.velocity),
                                       ("eta", state.eta), ("omega", omega),
                                       ("probes_u", probes_u), ("probes_w", probes_w)):
                        histories[key].append(value)
                if step and step % max(1, config.nsteps//5) == 0:
                    print(f"{path.stem}: t={state.time:.3f}, omegaL2={norms['full']:.5g}", flush=True)
            for key, values in histories.items():
                h.create_dataset(key, data=np.asarray(values), compression="gzip", compression_opts=1)
            for key in rows[0]:
                h.create_dataset(f"diagnostics/{key}", data=[r[key] for r in rows])
            for key in steps[0]:
                h.create_dataset(f"steps/{key}", data=[r[key] for r in steps])
            h.create_dataset("surface_x", data=fem.surface_x)
            h.attrs["status"] = "complete"
        except BaseException as exc:
            h.attrs["status"] = "failed"
            h.attrs["failure"] = repr(exc)
            raise
    write_tables(path.parent, path.stem+"_diagnostics", rows,
                 {"config": json.loads(config.to_json()), "h_min": h_min,
                  "rough_nu_dt_over_hmin_squared": config.nu*config.dt/h_min**2,
                  "runtime_seconds": time.perf_counter()-start})
    write_tables(path.parent, path.stem+"_steps", steps)
    print("completed", path.name, f"{time.perf_counter()-start:.1f}s", flush=True)


def compare_histories(fem, regions, path, reference):
    with h5py.File(path) as a, h5py.File(reference) as b:
        np.testing.assert_allclose(a["times"][:], b["times"][:], atol=1e-14)
        series = []
        reference_norms = {key: 0. for key in ("eta", "velocity", *regions.matrices)}
        for i, t in enumerate(a["times"]):
            dv, de, dw = a["velocity"][i]-b["velocity"][i], a["eta"][i]-b["eta"][i], a["omega"][i]-b["omega"][i]
            vn, en = b["velocity"][i], b["eta"][i]
            wn = regions.norms(b["omega"][i])
            row = {"time": float(t), "eta": float(np.sqrt(de@(fem.S@de))),
                   "velocity": float(np.sqrt(dv@(fem.M@dv))), **regions.norms(dw),
                   "probes_max": float(max(np.max(abs(a["probes_u"][i]-b["probes_u"][i])),
                                            np.max(abs(a["probes_w"][i]-b["probes_w"][i]))))}
            row["energy_difference"] = abs(float(a["diagnostics/total_energy"][i]-b["diagnostics/total_energy"][i]))
            row["slope_difference"] = abs(float(a["diagnostics/max_slope"][i]-b["diagnostics/max_slope"][i]))
            series.append(row)
            for key, value in {"eta": np.sqrt(en@(fem.S@en)), "velocity": np.sqrt(vn@(fem.M@vn)), **wn}.items():
                reference_norms[key] = max(reference_norms[key], float(value))
        result = {}
        for key in series[0]:
            if key == "time":
                continue
            result[key+"_difference_final"] = series[-1][key]
            result[key+"_difference_max"] = max(row[key] for row in series)
            if key in reference_norms:
                result[key+"_relative_max"] = result[key+"_difference_max"]/max(reference_norms[key], 1e-30)
        # Observe errors at INTERNAL dt, before snapshots can alias a sign flip.
        ta, tb = a["steps/time"][:], b["steps/time"][:]
        e = a["steps/omega_probe_2"][:]-np.interp(ta, tb, b["steps/omega_probe_2"][:])
        early = e[(ta > 0) & (ta <= .05+1e-12)]
        result["early_corner_error_max"] = float(np.max(abs(early)))
        result["early_corner_error_sign_changes"] = int(np.sum(early[1:]*early[:-1] < 0))
        result["early_corner_error_alternating_amplitude"] = float(abs(np.mean(early*(-1.)**np.arange(len(early)))))
        result["rough_stiffness"] = float(a.attrs["rough_stiffness"])
    return result


def _compute_viscosity_history(nu, output):
    base = SimulationConfig(alpha_deg=.2, nu=nu, mesh="medium", t_end=.5, snapshot_dt=.025)
    fem = FEMSystem(base)
    projector = VorticityProjector(fem)
    regions = VorticityRegions(fem, projector)
    for method in INTEGRATORS:
        for dt in (.005, .0025, .00125):
            config = replace(base, dt=dt, integrator=method)
            path = Path(output)/"runs"/f"nu{nu:g}_{method}_dt{dt:g}.h5"
            run_history(fem, config, path, projector, regions)


def physical_time_study(output, workers=1):
    output = Path(output)
    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        from multiprocessing import get_context
        with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn")) as pool:
            futures = [pool.submit(_compute_viscosity_history, nu, str(output)) for nu in (.01, .1, 1.)]
            for future in futures:
                future.result()
    rows, cross, physics = [], [], []
    for nu in (.01, .1, 1.):
        base = SimulationConfig(alpha_deg=.2, nu=nu, mesh="medium", t_end=.5, snapshot_dt=.025)
        fem = FEMSystem(base)
        projector = VorticityProjector(fem)
        regions = VorticityRegions(fem, projector)
        paths = {}
        for method in INTEGRATORS:
            for dt in (.005, .0025, .00125):
                config = replace(base, dt=dt, integrator=method)
                path = output/"runs"/f"nu{nu:g}_{method}_dt{dt:g}.h5"
                paths[method, dt] = path
                run_history(fem, config, path, projector, regions)
                with h5py.File(path) as h:
                    summary = {"nu": nu, "integrator": method, "dt": dt}
                    for key in ("volume_error", "contact_error", "energy_balance_relative", "max_step_energy_increase",
                                "divergence_l2", "weak_divergence_l2", "divergence_relative", "omega_wall_symmetry",
                                "max_slope", "omega_max", "rk_energy_correction"):
                        summary[key+"_max"] = float(np.max(abs(h[f"diagnostics/{key}"][:])))
                    summary["no_slip_max"] = max(float(np.max(h[f"diagnostics/{wall}_{c}_max"][:]))
                                                  for wall in ("left", "right", "bottom") for c in ("u", "w"))
                    physics.append(summary)
        for method in INTEGRATORS:
            for dt in (.005, .0025):
                row = {"nu": nu, "integrator": method, "dt": dt, "reference_dt": dt/2,
                       **compare_histories(fem, regions, paths[method, dt], paths[method, dt/2])}
                rows.append(row)
                print("refinement", row, flush=True)
        for dt in (.005, .0025, .00125):
            cross.append({"nu": nu, "dt": dt, **compare_histories(fem, regions, paths["midpoint", dt], paths["sdirk2", dt])})
        # Persist each viscosity before moving on; completed runs are resumable.
        write_tables(output, "physical_refinement", rows)
        write_tables(output, "physical_cross_integrator", cross)
        write_tables(output, "physical_diagnostics", physics)
        plt = pyplot()
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))
        for method in INTEGRATORS:
            for dt in (.005, .0025, .00125):
                with h5py.File(paths[method, dt]) as h:
                    label = f"{method}, dt={dt:g}"
                    for ax, key in zip(axes.ravel(), ("omega_wall_left_L2", "total_energy", "max_slope", "divergence_l2")):
                        ax.plot(h["times"][:], h[f"diagnostics/{key}"][:], label=label)
                        ax.set(xlabel="t (s)", ylabel=key)
        axes[0, 0].legend(fontsize=7)
        fig.suptitle(f"Physical temporal validation: nu={nu:g}, alpha=0.2°, medium")
        fig.tight_layout()
        save_figure(fig, output, f"physical_nu{nu:g}")
        plt.close(fig)
        del fem, projector, regions
        gc.collect()
    from .temporal_plots import internal_refinement, temporal_error_plots
    internal_refinement(output)
    temporal_error_plots(output)
    return rows
