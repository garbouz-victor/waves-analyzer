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
from ..energy_validation import validate_energy
from ..provenance import run_provenance
from ..validation_policy import EnergyPolicy


def run_case(config, phase, output, velocity=None, observer=None):
    started=time.perf_counter()
    output=Path(output)
    historical=Path(__file__).resolve().parents[4]/"validation_results/step3"
    if output.resolve()==historical.resolve() or historical.resolve() in output.resolve().parents:
        raise ValueError("Historical STEP 3 result tree is read-only; choose validation_results/step3a1")
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
    history_stream=None
    try:
        solver=CHNSSolver(config)
        provenance=run_provenance(solver)
        solver.run_metadata=provenance
        if comm.rank==0:
            (output/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
        solver.initialize(phase,velocity)
        diagnostics=Diagnostics(solver)
        rows=[diagnostics.measure()]
        observed=[]
        if observer:
            observed.append(observer(solver))
        resolution=interface_resolution(solver)
        resolutions=[resolution]
        def record_resolution(row,value):
            row["interface_cells_min"]=value.get("cells_across_transition_min",0.)
            row["interface_cells_p05"]=value.get("cells_across_transition_p05",0.)
            row["interface_cells_p50"]=value.get("cells_across_transition_p50",0.)
            row["interface_cells_p95"]=value.get("cells_across_transition_p95",0.)
            row["interface_active_cells"]=value["active_cell_count"]
            worst=value.get("worst_cell",{})
            row["interface_worst_x"]=worst.get("centroid",[0.,0.])[0]
            row["interface_worst_z"]=worst.get("centroid",[0.,0.])[1]
        record_resolution(rows[0],resolution)
        if comm.rank==0:
            history_stream=(output/"history.csv").open("x",newline="")
            history_writer=csv.DictWriter(history_stream,fieldnames=rows[0].keys())
            history_writer.writeheader();history_writer.writerow(rows[0]);history_stream.flush()
        if not resolution["qualified_resolution"] and comm.rank==0:
            warnings.warn(f"Under-resolved diffuse interface (not qualified): {resolution}",RuntimeWarning)
        nsteps=solver.schedule.nsteps
        budget_warning=False
        for step in range(nsteps):
            log=solver.advance()
            rows.append(diagnostics.measure())
            current_resolution=interface_resolution(solver)
            resolutions.append(current_resolution)
            record_resolution(rows[-1],current_resolution)
            if comm.rank==0:
                history_writer.writerow(rows[-1]);history_stream.flush()
            if abs(rows[-1]["energy_budget_relative_to_initial_scale"])>EnergyPolicy().budget_relative_tolerance and not budget_warning:
                if comm.rank==0:
                    warnings.warn("Time-integrated dissipation is NOT qualified: "
                        f"relative budget defect={rows[-1]['energy_budget_relative_to_initial_scale']:.6g}. "
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
        energy_gate=validate_energy(rows)
        result={"status":"complete","qualification_status":"not_assessed",
            "energy_validation":energy_gate,"provenance":provenance,
            "fingerprint":config.fingerprint(),"config":config.as_dict(),
            "status_note":"complete means solver finished, NOT model/benchmark qualification",
            "energy_accounting":"trapezoidal step-end power integral; requires independent time convergence",
            "runtime_s":time.perf_counter()-started,"mesh":solver.mesh_report,"dimensionless":config.groups(),
            "initial_interface_resolution":resolution,"final_interface_resolution":interface_resolution(solver),
            "interface_resolved_all_times":all(r["qualified_resolution"] for r in resolutions),
            "minimum_interface_cells_all_times":min(r.get("cells_across_transition_min",0.) for r in resolutions),
            "worst_interface_resolution":min(resolutions,key=lambda r:r.get("cells_across_transition_min",0.)),
            "max_speed_m_per_s":max(r["speed_max_dof_sample"] for r in rows),
            "max_mass_error_relative":max(abs(r["mass_error_relative_to_domain"]) for r in rows),
            "max_energy_increase_step":float(max(np.diff([r["E_total"] for r in rows]))),
            "max_energy_budget_relative":max(abs(r["energy_budget_relative"]) for r in rows),
            "max_impermeability_L2":max(r["wall_normal_velocity_L2"] for r in rows),
            "initial":rows[0],"final":rows[-1]}
        if comm.rank==0:
            for name,data in (("nonlinear",solver.logs),):
                with (output/(name+".csv")).open("w",newline="") as f:
                    writer=csv.DictWriter(f,fieldnames=data[0].keys())
                    writer.writeheader();writer.writerows(data)
            (output/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
            if observer:
                (output/"observations.json").write_text(json.dumps(observed,indent=2,allow_nan=False)+"\n")
        return solver,rows,result
    except BaseException as error:
        if comm.rank==0:
            (output/"summary.json").write_text(json.dumps({"status":"failed","config":config.as_dict(),
                "fingerprint":config.fingerprint(),"error":str(error),"runtime_s":time.perf_counter()-started},indent=2)+"\n")
        raise
    finally:
        if history_stream is not None:
            history_stream.close()
