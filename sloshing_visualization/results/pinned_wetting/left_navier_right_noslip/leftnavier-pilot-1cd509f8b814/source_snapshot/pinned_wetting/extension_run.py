"""Accepted-state Navier execution, extending the existing PW1 checkpoint schema."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np

from .coating import Coating
from .navier import LABEL, NavierConfig, NavierDiagnostics, NavierFEM, NavierIntegrator, NavierState
from .mission import read_json, write_json, sha256


def numerical_identity(root):
    names = ["config.py", "mesh.py", "fem_spaces.py", "problems.py", "time_integrator.py", "diagnostics.py",
             "pinned_wetting/coating.py", "pinned_wetting/navier.py", "pinned_wetting/extension_run.py"]
    src = root / "sloshing_visualization/src/sloshing"
    hashes = {name: sha256(src/name) for name in names}
    return hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest(), hashes


def physical_identity(c):
    return {"a_m": c.a, "d_m": c.d, "g_m_s2": c.g, "nu_m2_s": c.nu, "alpha_deg": c.alpha_deg,
            "slip_length_m": c.slip_length_m, "sigma_N_m": 0., "bottom": "no_slip", "sides": "impermeable_Navier",
            "geometry": "linearized_free_surface", "layer": "irreversible_lower_connected_zero_volume_one_way"}


def execute_case(root, run_dir, config, identity, mission_state, heartbeat, stop_after=None,
                 *, fem_factory=NavierFEM, diagnostics_factory=NavierDiagnostics,
                 contact_method="actual_free_endpoint_trace_Navier", initialize_native=None,
                 record_sides=("L", "R")):
    path = run_dir / "native.h5"
    if path.exists():
        with h5py.File(path, "r") as old:
            if json.loads(old.attrs["identity"]) != identity:
                raise ValueError("Native source/config identity mismatch")
            if old.attrs.get("completed", False):
                return read_json(run_dir / "summary.json")
    f = fem_factory(config)
    integrator = NavierIntegrator(f)
    initial = integrator.initial_state()
    diagnostics = diagnostics_factory(f, initial)
    resumed = path.exists()
    with h5py.File(path, "r+" if resumed else "x") as h:
        if resumed:
            keys = sorted(h["accepted"], key=int)
            if keys != [str(i) for i in range(len(keys))] or not all(h["accepted"][k].attrs.get("accepted", False) for k in keys):
                raise ValueError("Unpublished or non-contiguous checkpoint; retain this run for diagnosis")
            last = h["accepted"][keys[-1]]
            scalars = json.loads(last.attrs["integrator_scalars"])
            s = NavierState(velocity=last["velocity"][:], eta=last["eta"][:], **scalars)
            coat = Coating.from_dict(json.loads(last.attrs["coating"]))
            if (coat.run_id != identity["run_id"] or coat.time != s.time or coat.accepted_step != s.step or
                    not np.array_equal(coat.H, last["H"][:]) or
                    not np.array_equal(last["wet_intervals"][:], [[-config.d, v] for v in coat.H])):
                raise ValueError("Checkpoint coating/solver mismatch")
            running = initial.eta[[0, -1]].copy()
            for k in keys:
                g = h["accepted"][k]
                running = np.maximum(running, g["eta"][:][[0, -1]])
                memory = json.loads(g.attrs["coating"])
                if (not np.array_equal(g["H"][:], running) or memory["H"] != running.tolist() or
                        memory["markers"][:2] != Coating.initial(identity["run_id"], -config.d, *initial.eta[[0,-1]]).to_dict()["markers"]):
                    raise ValueError("Checkpoint wet history/P1/P2 corrupted")
            recent = [h["accepted"][k]["eta"][:][[0, -1]] for k in keys[-2:]]
        else:
            h.attrs.update(identity=json.dumps(identity, sort_keys=True), config=config.to_json(),
                           physical_model_label=identity.get("physical_model_label", LABEL), contact_method=contact_method,
                           wetting_sampling="every_accepted_step_endpoint", layer_model="L_macro union W_i=[-d,H_i]")
            mesh = h.create_group("mesh")
            for key, value in {"points": f.mesh.p, "triangles": f.mesh.t, "surface_x": f.surface_x,
                               "free_velocity_dofs": f.free}.items():
                mesh.create_dataset(key, data=value)
            h.create_group("accepted")
            if initialize_native is not None:
                initialize_native(h, f)
            s = initial
            coat = Coating.initial(identity["run_id"], -config.d, *initial.eta[[0, -1]])
            row = diagnostics.measure(s)
            row["momentum_relative_residual"] = 0.
            publish(h, s, coat, row)
            recent = [s.eta[[0, -1]].copy()]
        heartbeat(s, str(path))
        last_notice = time.monotonic()
        while s.step < config.nsteps:
            if stop_after is not None and s.step >= stop_after:
                h.flush(); heartbeat(s, str(path))
                return {"completed": False, "accepted_step": s.step, "time_s": s.time, "native": str(path)}
            candidate = integrator.advance(s)
            row = diagnostics.measure(candidate)
            row["momentum_relative_residual"] = integrator.momentum_residual
            diagnostics.validate(row)
            if row["momentum_relative_residual"] > 1e-9:
                raise ValueError("Navier variational residual exceeds tolerance")
            current = candidate.eta[[0, -1]].copy()
            coat.accept(candidate.step, candidate.time, current)
            if len(recent) == 2:
                for side_index, side in enumerate(("L", "R")):
                    if side not in record_sides:
                        continue
                    older, peak = recent[0][side_index], recent[1][side_index]
                    prior_marker = max(m.height_m for m in coat.markers if m.side == side)
                    if (peak > older and peak > current[side_index] and peak > prior_marker + .0001 and
                            peak >= coat.H[side_index] - 1e-13):
                        n = sum(m.side == side for m in coat.markers)
                        marker_id = f"P{n+2}" if side == "L" else f"R{n+2}"
                        coat.record_peak(marker_id, side, float(peak), s.time, f"{coat.run_id}:{s.step}")
            publish(h, candidate, coat, row, integrator.stages)
            s = candidate
            recent = (recent + [current])[-2:]
            if time.monotonic() - last_notice > 3 or s.step == config.nsteps:
                h.flush(); heartbeat(s, str(path))
                print(json.dumps({"run_id": identity["run_id"], "accepted_step": s.step,
                                  "t_s": s.time, "T_s": config.t_end, "R_m": current.tolist(), "H_m": coat.H}), flush=True)
                last_notice = time.monotonic()
        rows = [json.loads(h["accepted"][str(i)].attrs["diagnostics"]) for i in range(s.step+1)]
        summary = {"completed": True, "run_id": identity["run_id"], "native": str(path), "time_s": s.time,
            "accepted_states": s.step+1, "max_eta_over_a": max(r["max_eta_over_a"] for r in rows),
            "max_slope": max(r["max_slope"] for r in rows), "max_eta_over_b": max(r["max_eta_over_b"] for r in rows),
            "max_convective_indicator": max(r["convective_indicator"] for r in rows),
            "max_kinematic_indicator": max(r["kinematic_indicator"] for r in rows),
            "max_wall_traction_L2": max(r["wall_traction_L2"] for r in rows),
            "max_momentum_relative_residual": max(r["momentum_relative_residual"] for r in rows),
            "max_weak_divergence": max(r["weak_divergence_l2"] for r in rows),
            "max_mass_relative_error": max(r["volume_error"] for r in rows)/(2*config.a*config.d),
            "max_energy_relative_error": max(abs(r["energy_balance_relative"]) for r in rows),
            "H_final_m": coat.H, "markers": coat.to_dict()["markers"],
            "final_wall_loss_per_density": s.wall_dissipated_energy,
            "final_bulk_loss_per_density": s.bulk_dissipated_energy,
            "final_numerical_remainder_per_density": s.rk_energy_correction}
        for key in ("right_endpoint_error", "right_u_max", "right_w_max", "left_u_max", "bottom_u_max", "bottom_w_max"):
            if key in rows[0]:
                summary["max_"+key] = max(abs(r[key]) for r in rows)
        h.attrs["completed"] = True
        h.flush()
    write_json(run_dir / "summary.json", summary)
    return summary


def publish(h, s, coating, row, stages=None):
    g = h["accepted"].create_group(str(s.step))
    for name, value in (("velocity", s.velocity), ("eta", s.eta), ("H", coating.H),
                        ("wet_intervals", [[coating.bottom, v] for v in coating.H])):
        g.create_dataset(name, data=value, compression="lzf")
    if stages is not None:
        g.create_dataset("stage1_velocity", data=stages[0][0], compression="lzf")
        g.create_dataset("stage1_pressure", data=stages[0][-1], compression="lzf")
        g.create_dataset("stage2_pressure", data=stages[1][-1], compression="lzf")
    scalars = {name: getattr(s, name) for name in s.__dataclass_fields__ if name not in ("velocity", "eta")}
    g.attrs.update(time_s=s.time, state_id=f"{coating.run_id}:{s.step}", coating=json.dumps(coating.to_dict()),
                   diagnostics=json.dumps(row), integrator_scalars=json.dumps(scalars), accepted=True)
