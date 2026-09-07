"""Measured histories and explicit status for each independent benchmark."""
import csv
import json
import time
import warnings
from pathlib import Path
import numpy as np
from . import __doc__ as _benchmark_scope
from ..solver import CHNSSolver
from ..diagnostics import Diagnostics, interface_resolution
from ..storage import write_checkpoint


def run_case(config, phase, output, velocity=None, observer=None):
    started=time.perf_counter()
    output=Path(output)
    # Creation, not overwriting a successful/failed scientific result.
    from mpi4py import MPI
    comm=MPI.COMM_WORLD
    creation_error=None
    if comm.rank==0:
        try:
            output.mkdir(parents=True,exist_ok=True)
            # Atomic exclusive creation. Only rank zero accesses this status
            # before broadcasting; other ranks must not race a 'running' file.
            with (output/"summary.json").open("x") as stream:
                json.dump({"status":"running","config":config.as_dict(),
                           "fingerprint":config.fingerprint()},stream,indent=2)
        except Exception as error:
            creation_error=f"Cannot create run {output}: {error}; existing cache is not overwritten"
    creation_error=comm.bcast(creation_error,root=0)
    if creation_error:
        raise FileExistsError(creation_error)
    solver=None
    try:
        solver=CHNSSolver(config)
        solver.initialize(phase,velocity)
        diagnostics=Diagnostics(solver)
        rows=[diagnostics.measure()]
        observed=[]
        if observer:
            observed.append(observer(solver))
        resolution=interface_resolution(solver)
        if not resolution["qualified_resolution"] and comm.rank==0:
            warnings.warn(f"Under-resolved diffuse interface (not qualified): {resolution}",RuntimeWarning)
        nsteps=int(round(config.t_end/config.dt))
        if not np.isclose(nsteps*config.dt,config.t_end):
            raise ValueError("t_end must be a multiple of fixed dt")
        budget_warning=False
        for step in range(nsteps):
            log=solver.advance()
            rows.append(diagnostics.measure())
            if abs(rows[-1]["energy_budget_relative"])>.01 and not budget_warning:
                if comm.rank==0:
                    warnings.warn("Time-integrated dissipation is NOT qualified: "
                        f"relative budget defect={rows[-1]['energy_budget_relative']:.6g}. "
                        "Resolve the startup/time quadrature; do not interpret this integral as physical heat.",
                        RuntimeWarning)
                budget_warning=True
            if observer:
                observed.append(observer(solver))
            if comm.rank==0 and (step==0 or (step+1)%max(1,nsteps//10)==0):
                print(f"{config.name}: {step+1}/{nsteps} t={solver.time:.6g} "
                      f"SNES={log['newton_iterations']} residual={log['residual']:.3e} "
                      f"Umax={rows[-1]['speed_max_dof_sample']:.3e} "
                      f"E={rows[-1]['E_total']:.10g}",flush=True)
            if (step+1)%max(1,nsteps//10)==0:
                write_checkpoint(output/"checkpoint",solver,diagnostics,"running")
        write_checkpoint(output/"checkpoint",solver,diagnostics,"complete")
        result={"status":"complete","fingerprint":config.fingerprint(),"config":config.as_dict(),
            "status_note":"complete means solver finished, NOT model/benchmark qualification",
            "energy_accounting":"trapezoidal step-end power integral; requires independent time convergence",
            "runtime_s":time.perf_counter()-started,"mesh":solver.mesh_report,"dimensionless":config.groups(),
            "initial_interface_resolution":resolution,"final_interface_resolution":interface_resolution(solver),
            "max_speed_m_per_s":max(r["speed_max_dof_sample"] for r in rows),
            "max_mass_error_relative":max(abs(r["mass_error_relative_to_domain"]) for r in rows),
            "max_energy_increase_step":float(max(np.diff([r["E_total"] for r in rows]))),
            "max_energy_budget_relative":max(abs(r["energy_budget_relative"]) for r in rows),
            "max_impermeability_L2":max(r["wall_normal_velocity_L2"] for r in rows),
            "initial":rows[0],"final":rows[-1]}
        if comm.rank==0:
            for name,data in (("history",rows),("nonlinear",solver.logs)):
                with (output/(name+".csv")).open("w",newline="") as f:
                    writer=csv.DictWriter(f,fieldnames=data[0].keys())
                    writer.writeheader();writer.writerows(data)
            (output/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
            if observer:
                (output/"observations.json").write_text(json.dumps(observed,indent=2,allow_nan=False)+"\n")
        return solver,rows,result
    except Exception as error:
        if comm.rank==0:
            (output/"summary.json").write_text(json.dumps({"status":"failed","config":config.as_dict(),
                "fingerprint":config.fingerprint(),"error":str(error),"runtime_s":time.perf_counter()-started},indent=2)+"\n")
        raise
