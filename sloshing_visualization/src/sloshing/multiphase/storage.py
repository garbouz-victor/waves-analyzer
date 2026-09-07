"""Rank-local exact coefficient checkpoints; no silent failed-cache reuse."""
import json
from pathlib import Path
import h5py
import numpy as np

SCHEMA = "step3-chns-checkpoint-v1"


def write_checkpoint(folder, solver, diagnostics, status="running"):
    path=Path(folder)
    if status not in ("running","complete","failed"):
        raise ValueError("Invalid checkpoint status")
    if solver.comm.rank == 0:
        path.mkdir(parents=True,exist_ok=True)
    solver.comm.barrier()
    target=path/f"rank_{solver.comm.rank:04d}.h5"
    temporary=target.with_suffix(".partial.h5")
    with h5py.File(temporary,"w") as h:
        h.attrs.update(schema=SCHEMA,status=status,config_fingerprint=solver.config.fingerprint(),
                       ranks=solver.comm.size,step=solver.step_number,time=solver.time,
                       dt=solver.config.dt,config_json=json.dumps(solver.config.as_dict()),
                       energy_json=json.dumps({"initial":diagnostics.initial,"previous":diagnostics.previous,
                                               "cumulative":diagnostics.cumulative}))
        for name in ("state","old","older"):
            h.create_dataset(name,data=getattr(solver,name).x.array)
        fields=h.create_group("fields")
        for index,name in enumerate(("velocity","pressure_pi","phi","mu_ch")):
            field=solver.state.sub(index).collapse()
            group=fields.create_group(name)
            group.create_dataset("coefficients",data=field.x.array)
            group.create_dataset("dof_coordinates",data=field.function_space.tabulate_dof_coordinates())
        fields.attrs["mechanical_pressure_recovery"]="p=pi+phi*(mu_ch-rho_prime*Psi)-rho_bar*Psi"
        fields.attrs["units"]="velocity m/s; pi,mu_ch Pa; phi dimensionless; energies J/m; powers W/m"
        h.create_dataset("geometry",data=solver.mesh.geometry.x)
        h.create_dataset("geometry_dofmap",data=solver.mesh.geometry.dofmap)
    temporary.replace(target)
    solver.comm.barrier()
    if solver.comm.rank == 0:
        (path/"status.json").write_text(json.dumps({"schema":SCHEMA,"status":status,
            "ranks":solver.comm.size,"step":solver.step_number,"time":solver.time,
            "config_fingerprint":solver.config.fingerprint()},indent=2)+"\n")


def load_checkpoint(folder, solver, diagnostics, allow_running=False):
    folder=Path(folder)
    manifest=json.loads((folder/"status.json").read_text())
    allowed=("complete","running") if allow_running else ("complete",)
    if manifest["status"] not in allowed:
        raise ValueError("Checkpoint not complete; explicit running restart required; failed is invalid")
    with h5py.File(folder/f"rank_{solver.comm.rank:04d}.h5","r") as h:
        if h.attrs["schema"]!=SCHEMA or h.attrs["config_fingerprint"]!=solver.config.fingerprint():
            raise ValueError("Checkpoint/config mismatch")
        if h.attrs["ranks"]!=solver.comm.size:
            raise ValueError("Restart requires the original MPI rank count and mesh partition")
        if h.attrs["step"]!=manifest["step"] or h.attrs["status"]!=manifest["status"]:
            raise ValueError("Incomplete rank checkpoint transaction")
        if not np.array_equal(h["geometry"][:],solver.mesh.geometry.x):
            raise ValueError("Mesh partition differs")
        for name in ("state","old","older"):
            getattr(solver,name).x.array[:]=h[name][:]
            getattr(solver,name).x.scatter_forward()
        solver.step_number,solver.time=int(h.attrs["step"]),float(h.attrs["time"])
        energy=json.loads(h.attrs["energy_json"])
        diagnostics.initial,diagnostics.previous=energy["initial"],energy["previous"]
        diagnostics.cumulative=energy["cumulative"]
    solver._check_material()
