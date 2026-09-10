"""Read-only reuse of the pinned FEM; native retained state at EVERY step."""
from dataclasses import asdict
import json
from pathlib import Path
import time

import h5py
import numpy as np

from ..config import SimulationConfig
from ..diagnostics import Diagnostics
from ..solver import SloshingSolver
from ..time_integrator import State
from .coating import Coating

LABEL = "REFERENCE / NOT TARGET"


def run_reference(run_dir, config, identity, heartbeat, stop_after=None):
    """Append-only accepted groups; each group is also a restart checkpoint.

    Incomplete groups are never used. A publication error fails the run; this
    small adapter intentionally does not try to repair damaged HDF5 journals.
    """
    path = Path(run_dir) / "native.h5"
    solver = SloshingSolver(config)
    f = solver.fem
    initial = solver.integrator.initial_state()
    diagnostic = Diagnostics(f, initial)
    resumed = path.exists()
    with h5py.File(path, "r+" if resumed else "x") as h:
        if resumed:
            if json.loads(h.attrs["identity"]) != identity:
                raise ValueError("Resume numerical/physical identity mismatch")
            groups = h["accepted"]
            keys = sorted(groups.keys(), key=int)
            if keys != [str(i) for i in range(len(keys))]:
                raise ValueError("Non-contiguous accepted history")
            last = groups[keys[-1]]
            if not all(g.attrs.get("accepted", False) for g in groups.values()):
                raise ValueError("Unpublished native group; preserve run and use a new ID")
            scalars = json.loads(last.attrs["integrator_scalars"])
            state = State(velocity=last["velocity"][:], eta=last["eta"][:], **scalars)
            coating = Coating.from_dict(json.loads(last.attrs["coating"]))
            if coating.accepted_step != state.step or coating.run_id != identity["run_id"]:
                raise ValueError("Checkpoint retained state mismatch")
        else:
            h.attrs.update(identity=json.dumps(identity, sort_keys=True),
                           config=config.to_json(), physical_model_label=LABEL,
                           contact_method="actual_fixed_endpoint_trace",
                           wetting_sampling="every_accepted_step_endpoint",
                           layer_model="L_macro union W_i=[-d,H_i], zero-volume leading order")
            mesh = h.create_group("mesh")
            for key, value in {"points": f.mesh.p, "triangles": f.mesh.t,
                               "surface_x": f.surface_x, "free_velocity_dofs": f.free}.items():
                mesh.create_dataset(key, data=value)
            h.create_group("accepted")
            state = initial
            coating = Coating.initial(identity["run_id"], -config.d, *state.eta[[0, -1]])
            _publish(h, state, coating, diagnostic.measure(state))
        saved = time.monotonic()
        heartbeat(state, str(path))
        while state.step < config.nsteps:
            if stop_after is not None and state.step >= stop_after:
                return {"completed": False, "step": state.step, "time_s": state.time}
            candidate = solver.integrator.advance(state)
            row = diagnostic.measure(candidate)
            diagnostic.validate(row)
            coating.accept(candidate.step, candidate.time, candidate.eta[[0, -1]])
            _publish(h, candidate, coating, row)
            state = candidate
            if time.monotonic() - saved >= 3 or state.step == config.nsteps:
                h.flush()
                heartbeat(state, str(path))
                print(json.dumps({"run_id": identity["run_id"], "accepted_step": state.step,
                                  "time_s": state.time, "t_end_s": config.t_end}), flush=True)
                saved = time.monotonic()
        h.attrs["complete_reference"] = True
        h.flush()
    return {"completed": True, "step": state.step, "time_s": state.time}


def _publish(h, state, coating, row):
    group = h["accepted"].create_group(str(state.step))
    group.create_dataset("velocity", data=state.velocity, compression="lzf")
    group.create_dataset("eta", data=state.eta)
    group.create_dataset("H", data=coating.H)
    group.create_dataset("wet_intervals", data=[[-hconfig(h).d, v] for v in coating.H])
    group.attrs["coating"] = json.dumps(coating.to_dict())
    group.attrs["diagnostics"] = json.dumps(row)
    scalar_names = [name for name in state.__dataclass_fields__ if name not in ("velocity", "eta")]
    group.attrs["integrator_scalars"] = json.dumps({name: getattr(state, name) for name in scalar_names})
    group.attrs["time_s"] = state.time
    group.attrs["state_id"] = f"{coating.run_id}:{state.step}"
    group.attrs["accepted"] = True


def hconfig(h):
    return SimulationConfig(**json.loads(h.attrs["config"]))


def load_native(path):
    with h5py.File(path, "r") as h:
        groups = [h["accepted"][k] for k in sorted(h["accepted"], key=int)]
        return {"identity": json.loads(h.attrs["identity"]),
                "config": json.loads(h.attrs["config"]), "x": h["mesh/surface_x"][:],
                "contact_method": h.attrs["contact_method"],
                "times": np.array([g.attrs["time_s"] for g in groups]),
                "eta": np.array([g["eta"][:] for g in groups]),
                "H": np.array([g["H"][:] for g in groups]),
                "intervals": np.array([g["wet_intervals"][:] for g in groups]),
                "ids": [g.attrs["state_id"] for g in groups],
                "accepted": [bool(g.attrs["accepted"]) for g in groups],
                "coatings": [json.loads(g.attrs["coating"]) for g in groups],
                "diagnostics": [json.loads(g.attrs["diagnostics"]) for g in groups],
                "initial_velocity_max": float(np.max(abs(groups[0]["velocity"][:])))}
