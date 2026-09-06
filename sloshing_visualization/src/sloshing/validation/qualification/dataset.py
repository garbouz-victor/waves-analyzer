"""Streaming, reusable FEM dataset; continuation only from a COMPLETE prefix."""

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import time

import h5py
import numpy as np

from ...diagnostics import Diagnostics
from ...solver import SloshingSolver
from ...postprocess.vorticity import VorticityProjector
from ...time_integrator import State
from ..physical import PROBES
from ..regions import VorticityRegions
from ..contact import surface_value
from .metrics import DerivativeMetrics, modal_projection, slope_regions


SCHEMA = "animation-qualification-1.6"
_VERIFIED_PAYLOADS = set()
RESTART = ("dissipated_energy", "step_energy_residual", "max_step_energy_residual",
           "max_step_energy_increase", "work_energy", "rk_energy_correction")


def dataset_fingerprint(path):
    path = Path(path)
    config = validate_cache(path)
    stat = path.stat()
    return {"path": str(path.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "config_sha256": hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()}


def validate_cache(path, config=None, prefix=False):
    try:
        with h5py.File(path, "r") as h:
            if h.attrs.get("status") != "complete" or h.attrs.get("schema") != SCHEMA:
                raise RuntimeError(f"Unfinished/incompatible dataset: {path}")
            stored = json.loads(h.attrs["config_json"])
            if config is not None:
                for k, v in asdict(config).items():
                    if prefix and k == "t_end":
                        continue
                    if stored[k] != v:
                        raise RuntimeError(f"Cache configuration mismatch {k}: {path}")
            times = h["times"][:]
            count = round(stored["t_end"]/stored["snapshot_dt"])+1
            if len(times) != count or not np.allclose(times, np.arange(count)*stored["snapshot_dt"], atol=1e-13, rtol=0):
                raise RuntimeError(f"Incomplete dataset time grid: {path}")
            for key in ("u", "w", "q", "eta", "omega"):
                if h[f"fields/{key}"].shape[0] != count:
                    raise RuntimeError(f"Incomplete field {key}: {path}")
                if not np.all(np.isfinite(h[f"fields/{key}"][-1])):
                    raise RuntimeError(f"Non-finite endpoint field {key}: {path}")
            if h.attrs.get("snapshots_written") != count:
                raise RuntimeError(f"Snapshot completion marker mismatch: {path}")
            for group in ("diagnostics", "restart"):
                for key, ds in h[group].items():
                    if ds.shape != (count,) or not np.all(np.isfinite(ds[:])):
                        raise RuntimeError(f"Incomplete/non-finite {group}/{key}: {path}")
            if not np.array_equal(h["restart/step"][:], np.rint(times/stored["dt"]).astype(int)):
                raise RuntimeError(f"Restart step/time mismatch: {path}")
            if prefix and times[-1] >= config.t_end:
                raise RuntimeError("Prefix must end before target")
            # Read every compressed payload at least once per file revision.
            # An intact endpoint alone does not certify earlier saved states.
            stat = Path(path).stat()
            revision = (str(Path(path).resolve()), stat.st_size, stat.st_mtime_ns)
            if revision not in _VERIFIED_PAYLOADS:
                for group in ("fields", "probes"):
                    for key, ds in h[group].items():
                        if group == "probes" and key in ("coordinates", "eta_coordinates"):
                            continue
                        if ds.shape[0] != count:
                            raise RuntimeError(f"Incomplete {group}/{key}: {path}")
                        for i in range(count):
                            if not np.all(np.isfinite(ds[i])):
                                raise RuntimeError(f"Non-finite {group}/{key} at snapshot {i}: {path}")
                _VERIFIED_PAYLOADS.add(revision)
            return stored
    except OSError as exc:
        raise RuntimeError(f"Unreadable dataset {path}; choose a new path; never silently overwrite") from exc


def _append(h, key, value, index):
    data = np.asarray(value)
    if key not in h:
        opts = {}
        if data.ndim:
            opts = dict(chunks=(1,)+data.shape, compression="gzip", compression_opts=1, shuffle=True,
                        fletcher32=True)
        h.create_dataset(key, shape=(0,)+data.shape, maxshape=(None,)+data.shape, dtype=data.dtype, **opts)
    ds = h[key]
    ds.resize(index+1, axis=0)
    ds[index] = data


def _metadata(h, fem, projector, config):
    h.attrs.update(schema=SCHEMA, status="running", config_json=config.to_json(),
                   created_utc=datetime.now(timezone.utc).isoformat(),
                   model="linearized Navier-Stokes; fixed domain; strict no-slip; sigma=0; pinned contacts",
                   pressure="instantaneous q; traction-fixed gauge",
                   dependencies_json=json.dumps({n: importlib.metadata.version(n) for n in ("numpy", "scipy", "scikit-fem", "h5py")}))
    for key, values in {"coordinates": fem.mesh.p, "triangles": fem.mesh.t,
                        "velocity_dof_coordinates": fem.scalar.doflocs,
                        "velocity_element_dofs": fem.scalar.element_dofs,
                        "pressure_element_dofs": fem.pressure.element_dofs,
                        "omega_element_dofs": projector.basis.element_dofs,
                        "surface_x": fem.surface_x}.items():
        h.create_dataset(f"mesh/{key}", data=values)
    h.create_dataset("probes/coordinates", data=PROBES)
    h.create_dataset("probes/eta_coordinates", data=[-.75, -.5, 0, .5, .75])


def run_dataset(config, path, prefix=None):
    path = Path(path)
    if path.exists():
        validate_cache(path, config)
        print("REUSE complete dataset", path, flush=True)
        return path
    if config.nsteps % config.save_every:
        raise ValueError("Qualification end time must be on the snapshot grid")
    if prefix is not None:
        validate_cache(prefix, config, prefix=True)
    start = time.perf_counter()
    print(f"ASSEMBLE {config.mesh}, t_end={config.t_end:g}, dt={config.dt:g}", flush=True)
    solver = SloshingSolver(config)
    f, it = solver.fem, solver.integrator
    print(f"FACTORIZED {config.mesh}; velocity unknowns={len(f.free)}, "
          f"pressure unknowns={f.pressure.N}; elapsed={time.perf_counter()-start:.1f}s", flush=True)
    projector = VorticityProjector(f)
    regions = VorticityRegions(f, projector)
    derivatives = DerivativeMetrics(f, regions)
    print(f"REGIONAL QUADRATURE ready; elapsed={time.perf_counter()-start:.1f}s", flush=True)
    vp = f.scalar.probes(PROBES)
    initial = it.initial_state()
    diagnostics = Diagnostics(f, initial)
    state, count = initial, 0
    path.parent.mkdir(parents=True, exist_ok=True)
    if prefix is not None:
        # Exclusive creation closes the race between the initial cache check
        # and long FEM assembly: even another concurrent run cannot be overwritten.
        with Path(prefix).open("rb") as source, path.open("xb") as target:
            shutil.copyfileobj(source, target, length=8*1024*1024)
    with h5py.File(path, "r+" if prefix is not None else "x") as h:
        if prefix is None:
            _metadata(h, f, projector, config)
        else:
            h.attrs.update(status="running", config_json=config.to_json(), continued_from=str(prefix),
                           continued_utc=datetime.now(timezone.utc).isoformat())
            count = len(h["times"])
            full = np.zeros(f.velocity.N)
            for key, ids in zip(("u", "w"), f.component_dofs):
                full[ids] = h[f"fields/{key}"][-1]
            state = State(int(h["restart/step"][-1]), float(h["times"][-1]), full[f.free], h["fields/eta"][-1],
                          **{k: float(h[f"restart/{k}"][-1]) for k in RESTART})
        h.flush()
        try:
            last_progress = time.perf_counter()
            first = state.step if prefix is None else state.step+1
            for step in range(first, config.nsteps+1):
                if step:
                    state = it.advance(state)
                if step % config.save_every:
                    continue
                row = diagnostics.measure(state)
                diagnostics.validate(row)
                full = f.expand_velocity(state.velocity)
                omega = projector.evaluate(full)
                norms = regions.norms(omega)
                row.update(slope_regions(f.surface_x, state.eta))
                row.update(derivatives.evaluate(full))
                row.update(omega_L2=norms["full"], omega_wall_left_L2=norms["wall_left"],
                           omega_wall_right_L2=norms["wall_right"], omega_surface_layer_L2=norms["surface_layer"],
                           omega_max=float(np.max(abs(omega))),
                           omega_wall_symmetry=abs(norms["wall_left"]-norms["wall_right"]),
                           eta_antisymmetry=float(np.max(abs(state.eta+state.eta[::-1]))),
                           modal_eta=modal_projection(f.surface_x, state.eta, config.a))
                pressure = it.last_stage_pressure if config.integrator == "sdirk2" and step else it.instantaneous_pressure(state)
                if config.integrator == "sdirk2":
                    # Initial q needs recovery; subsequent stiffly accurate q is
                    # already computed. Do not retain a second large LU in RAM.
                    it._pressure_lu = None
                arrays = {"times": state.time, "fields/u": full[f.component_dofs[0]],
                          "fields/w": full[f.component_dofs[1]], "fields/q": pressure,
                          "fields/eta": state.eta, "fields/omega": omega,
                          "probes/u": vp@full[f.component_dofs[0]], "probes/w": vp@full[f.component_dofs[1]],
                          "probes/eta": surface_value(f.surface_x, state.eta, h["probes/eta_coordinates"][:]),
                          "restart/step": state.step,
                          **{f"restart/{k}": getattr(state, k) for k in RESTART},
                          **{f"diagnostics/{k}": v for k, v in row.items()}}
                if not all(np.all(np.isfinite(v)) for v in arrays.values()):
                    raise FloatingPointError("Non-finite qualification field")
                for key, values in arrays.items():
                    _append(h, key, values, count)
                count += 1
                h.attrs["snapshots_written"] = count
                if step % max(config.save_every, round(.1/config.dt)) == 0 or time.perf_counter()-last_progress >= 30:
                    # Flush a streaming batch, not dozens of HDF5 metadata
                    # chunks at every frame. A running/failed file is NEVER reused.
                    h.flush()
                    elapsed = time.perf_counter()-start
                    print(f"{config.mesh} t={state.time:.3f}/{config.t_end:g}; snapshots={count}; "
                          f"slope={row['slope_global']:.5g}; div_rel={row['divergence_relative']:.4g}; "
                          f"E={row['total_energy']:.7g}; elapsed={elapsed:.1f}s", flush=True)
                    last_progress = time.perf_counter()
            h.flush()  # Publish all arrays before publishing the complete marker.
            h.attrs.update(status="complete", runtime_this_invocation_seconds=time.perf_counter()-start)
            h.flush()
        except BaseException as exc:
            h.attrs.update(status="failed", failure=repr(exc))
            h.flush()
            raise
    validate_cache(path, config)
    print("COMPLETE", path, f"{time.perf_counter()-start:.1f}s", flush=True)
    return path


def dataset_summary(path):
    validate_cache(path)
    with h5py.File(path) as h:
        times = h["times"][:]
        diag = h["diagnostics"]
        out = {"path": str(path), "complete": True, "snapshots": len(times),
               "parameters": json.loads(h.attrs["config_json"]),
               "runtime_this_invocation_seconds": float(h.attrs["runtime_this_invocation_seconds"])}
        for key in diag:
            v = diag[key][:]
            i = int(np.argmax(abs(v)))
            out[key+"_max"] = float(abs(v[i]))
            out[key+"_final"] = float(v[-1])
            out[key+"_time_at_max"] = float(times[i])
        for cutoff in (.1, .3):
            hit = np.flatnonzero(diag["slope_global"][:] > cutoff)
            out[f"first_slope_above_{cutoff:g}"] = float(times[hit[0]]) if len(hit) else None
        active = diag["velocity_gradient_l2"][:] > .01*np.max(diag["velocity_gradient_l2"][:])
        out["relative_div_above_5_percent_active_fraction"] = float(np.mean(diag["divergence_relative"][:][active] > .05)) if np.any(active) else 0.
        slope = out["slope_global_max"]
        angle = np.deg2rad(out["parameters"]["alpha_deg"])
        out["extrapolated_alpha_for_slope_0_1_deg"] = float(np.rad2deg(np.arctan(np.tan(angle)*.1/slope))) if slope else None
        out["extrapolated_alpha_for_slope_0_05_deg"] = float(np.rad2deg(np.arctan(np.tan(angle)*.05/slope))) if slope else None
        out["slope_policy"] = "failed" if slope > .3 else "conditional" if slope > .1 else "qualified"
        return out
