"""Measured STEP 3A.2 sequence; each stage has explicit scientific prerequisites."""
import csv
import hashlib
import json
from pathlib import Path
import tarfile
import time
import numpy as np
from ..config import ModelConfig
from ..solver import CHNSSolver
from ..initialization import sessile_drop
from ..equilibrium import (solve_equilibrium, save_prepared, load_prepared, EquilibriumPolicy, EquilibriumFailure,
    first_variations, zero_mass_perturbations, _integral, _measures)
from ..equilibrium_diagnostics import state_report, transient_term_residuals, PreservationDiagnostics
from ..diagnostics import Diagnostics, interface_resolution
from ..provenance import run_provenance
from ..storage import write_checkpoint
from ..energy_validation import validate_energy, positive_energy_scale
from .contact_angle import observe


def write_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")


def reserve(output,kind,config):
    output=Path(output)
    historical=Path(__file__).resolve().parents[4]/"validation_results"
    for name in ("step3","step3a1"):
        if (historical/name).resolve() in (output.resolve(),*output.resolve().parents):
            raise ValueError("Historical result trees are read-only")
    output.mkdir(parents=True,exist_ok=True)
    with (output/"summary.json").open("x") as stream:
        json.dump({"status":"running","kind":kind,"config":config.as_dict()},stream,indent=2)
    source=Path(__file__).resolve().parents[1]
    with tarfile.open(output/"multiphase_source.tar.gz","x:gz") as archive:
        for path in sorted(source.rglob("*.py")):
            archive.add(path,arcname=str(path.relative_to(source)))
    return output


def prepare(config,output):
    started=time.perf_counter()
    output=reserve(output,"prepare_equilibrium",config)
    try:
        solver=CHNSSolver(config)
        if solver.comm.size!=1:
            raise ValueError("Scientific preparation currently requires serial partition")
        solver.initialize(sessile_drop(config.epsilon,config.refinement_circle_radius,60.,config.z_min))
        analytic=state_report(solver)
        write_json(output/"analytic60.json",analytic)
        def monitor(row):
            print("equilibrium "+json.dumps(row),flush=True)
            with (output/"snes_history.jsonl").open("a") as stream:
                stream.write(json.dumps(row)+"\n")
        equilibrium=solve_equilibrium(solver,mu_initial=config.sigma/(2*config.refinement_circle_radius),
                                      history_callback=monitor)
        solver.initialize_prepared_equilibrium(equilibrium)
        prepared=state_report(solver,equilibrium.metadata["mu_star"])
        variations=first_variations(solver,solver.state.sub(2).collapse())
        terms=transient_term_residuals(solver,equilibrium.metadata["mu_star"])
        policy=EquilibriumPolicy()
        gates={"mass_constraint":equilibrium.metadata["mass_residual_relative_to_domain"]<=policy.mass_domain_tolerance,
            "stationary_riesz":prepared["stationarity"]["normalized"]<=policy.normalized_riesz_tolerance,
            "constant_mu":prepared["mu_coefficient_spread"]<=policy.mu_spread_tolerance,
            "initial_CH_near_zero":prepared["CH_dissipation"]<=policy.initial_ch_tolerance,
            "wall_residual_materially_lower":prepared["wall_residual"]["L2"]<=policy.wall_improvement_ratio*analytic["wall_residual"]["L2"],
            "independent_first_variations":all(r["normalized_abs_derivative"]<=policy.first_variation_tolerance for r in variations),
            "certified_resolution":prepared["resolution"]["qualified_resolution"]}
        gates={k:bool(v) for k,v in gates.items()}
        save_prepared(output,equilibrium)
        result={"status":"complete","qualification_status":"passed" if all(gates.values()) else "failed",
            "gates":gates,"config":config.as_dict(),"provenance":run_provenance(solver),
            "mesh":solver.mesh_report,"runtime_s":time.perf_counter()-started,
            "equilibrium_fingerprint":equilibrium.metadata["fingerprint"],"analytic60":analytic,
            "prepared":prepared,"first_variations":variations,"transient_terms_at_prepared_state":terms,
            "free_energy_difference_prepared_minus_analytic":prepared["E_total"]-analytic["E_total"],
            "historical_interpretation":"geometrically compatible, variationally unprepared; not an equilibrium benchmark"}
        write_json(output/"summary.json",result)
        print(json.dumps({"equilibrium_gates":gates,"runtime_s":result["runtime_s"],
            "mass_target":equilibrium.metadata["target_mass"],"mu_star":equilibrium.metadata["mu_star"],
            "riesz":prepared["stationarity"],"certified_cells":prepared["resolution"].get("cells_across_transition_certified_min")}),flush=True)
        return result
    except BaseException as error:
        if isinstance(error,EquilibriumFailure):
            import h5py
            write_json(output/"failure_diagnostics.json",error.diagnostics)
            with h5py.File(output/"incomplete_candidate.h5","x") as h:
                h.attrs["schema"]="step3a2-incomplete-stationary-candidate-diagnostic-v1"
                h.attrs["status"]="incomplete"
                h.create_dataset("phi_coefficients",data=error.coefficients)
        write_json(output/"summary.json",{"status":"failed","qualification_status":"failed",
            "scientific_verdict":"MODEL NOT YET VALIDATED","stage":"PREPARED EQUILIBRIUM NOT OBTAINED",
            "config":config.as_dict(),"error":str(error),"runtime_s":time.perf_counter()-started})
        raise


def transient(equilibrium_folder,output,scheme="be",dt=1e-5,t_end=1e-3,perturb=False):
    started=time.perf_counter()
    equilibrium_folder=Path(equilibrium_folder)
    foundation=json.loads((equilibrium_folder/"summary.json").read_text())
    required=("mass_constraint","stationary_riesz","constant_mu","initial_CH_near_zero",
              "wall_residual_materially_lower","independent_first_variations","certified_resolution")
    if foundation.get("qualification_status")!="passed" or not all(foundation.get("gates",{}).get(k,False) for k in required):
        raise ValueError("STOP: equilibrium foundation has not passed")
    meta=json.loads((equilibrium_folder/"prepared.json").read_text())
    if not perturb and scheme=="bdf2":
        be=json.loads((equilibrium_folder.parent/"preservation_be/summary.json").read_text())
        if be.get("qualification_status")!="passed" or be.get("equilibrium_source_fingerprint")!=meta["fingerprint"]:
            raise ValueError("STOP: BE must preserve this equilibrium before BDF2")
    c=ModelConfig(**meta["config"]).changed(time_scheme=scheme,dt=dt,t_end=t_end,
        startup_dt=None,startup_t_end=0.,startup_stages=())
    output=reserve(output,"perturbation" if perturb else "preservation",c)
    rows=[];resolutions=[];geometry=[]
    stream=None
    try:
        solver=CHNSSolver(c)
        equilibrium=load_prepared(equilibrium_folder,solver,meta["target_mass"])
        solver.initialize_prepared_equilibrium(equilibrium)
        phi_eq=solver.state.sub(2).collapse()
        perturbation=None
        if perturb:
            # Both preservation records must refer to THIS prepared state.
            root=equilibrium_folder.parent
            for name in ("preservation_be","preservation_bdf2"):
                record=json.loads((root/name/"summary.json").read_text())
                if record.get("qualification_status")!="passed" or record.get("equilibrium_source_fingerprint")!=meta["fingerprint"]:
                    raise ValueError("STOP: both schemes must preserve this equilibrium before perturbation")
            from dolfinx import fem
            psi,description=zero_mass_perturbations(solver,phi_eq)["smooth_bulk"]
            perturbed=fem.Function(phi_eq.function_space)
            delta=1e-3
            perturbed.x.array[:]=phi_eq.x.array+delta*psi.x.array
            perturbed.x.scatter_forward()
            # Normal consistent weak mu projection; no equilibrium constant substituted.
            solver.initialize(perturbed)
            dx,_=_measures(solver)
            perturbation={**description,"delta":delta,"mass_change":_integral(solver,(perturbed-phi_eq)*dx),
                "psi_coefficient_sha256":hashlib.sha256(psi.x.array.tobytes()).hexdigest(),
                "initial_phi_coefficient_sha256":hashlib.sha256(perturbed.x.array.tobytes()).hexdigest(),
                "wall_trace_max_abs_psi":float(np.max(abs(psi.x.array[np.isclose(
                    psi.function_space.tabulate_dof_coordinates()[:,1],c.z_min,rtol=0,atol=1e-12)]))),
                "definition":"phi_eq + delta*psi_h; bulk compact support away from wall/contact line, plus roundoff-level mean subtraction"}
            write_json(output/"perturbation.json",perturbation)
        initial_report=state_report(solver,geometry=False)
        write_json(output/"initial_state.json",initial_report)
        if not perturb and initial_report["CH_dissipation"]>EquilibriumPolicy().initial_ch_tolerance:
            raise RuntimeError("STOP: prepared transient representation has nonzero initial CH power; no dt reduction")
        solver.run_metadata={**run_provenance(solver),"equilibrium_source_fingerprint":meta["fingerprint"]}
        write_json(output/"provenance.json",solver.run_metadata)
        diagnostics=Diagnostics(solver)
        drift=PreservationDiagnostics(solver,phi_eq)
        stream=(output/"history.csv").open("x",newline="")
        writer=None
        def sample():
            nonlocal writer
            row=diagnostics.measure();row.update(drift.measure())
            res=interface_resolution(solver)
            row["cells_across_transition_certified_min"]=res.get("cells_across_transition_certified_min",0.)
            row["cells_across_transition_directional_min"]=res.get("cells_across_transition_directional_min",0.)
            area=(c.x_max-c.x_min)*(c.z_max-c.z_min)
            row["mass_error_to_target_relative_to_domain"]=abs(row["phase_mass"]-meta["target_mass"])/area
            # Geometry at every accepted state; preserves all fitted windows/crossings.
            obs=observe(solver)
            if geometry:
                before=sorted(p["coordinate"] for p in geometry[0]["crossings"])
                after=sorted(p["coordinate"] for p in obs["crossings"])
                if len(before)!=2 or len(after)!=2:
                    raise RuntimeError("Contact topology changed in small-equilibrium transient")
                row["contact_displacement_max"]=float(max(abs(np.array(after)-before)))
                row["angle_drift_max"]=max(abs(a["theta_deg"]-b["theta_deg"]) for a,b in zip(obs["fits"],geometry[0]["fits"]))
            else:
                row["contact_displacement_max"]=row["angle_drift_max"]=0.
            geometry.append(obs);resolutions.append(res);rows.append(row)
            if writer is None:
                writer=csv.DictWriter(stream,fieldnames=row.keys());writer.writeheader()
            writer.writerow(row);stream.flush()
            with (output/"resolutions.jsonl").open("a") as f:
                f.write(json.dumps(res)+"\n")
            return row
        first=sample()
        for _ in range(solver.schedule.nsteps):
            log=solver.advance()
            row=sample()
            with (output/"nonlinear.jsonl").open("a") as f:
                f.write(json.dumps(log)+"\n")
            if solver.step_number%10==0 or solver.step_number==1:
                print(f"{output.name}: {solver.step_number}/{solver.schedule.nsteps} "
                    f"SNES={log['newton_iterations']} U={row['speed_max_dof_sample']:.3e} "
                    f"D_CH={row['CH_dissipation']:.3e} phi_drift={row['relative_phi_L2_drift']:.3e}",flush=True)
        write_checkpoint(output/"checkpoint",solver,diagnostics,"complete")
        write_json(output/"geometry.json",geometry)
        scale=positive_energy_scale(first)
        metrics={"max_speed":max(r["speed_max_dof_sample"] for r in rows),
            "relative_phi_L2_drift":max(r["relative_phi_L2_drift"] for r in rows),
            "phi_L2_drift":max(r["phi_L2_drift"] for r in rows),
            "phi_H1_drift":max(r["phi_H1_drift"] for r in rows),
            "mass_error":max(r["mass_error_to_target_relative_to_domain"] for r in rows),
            "max_D_CH":max(r["CH_dissipation"] for r in rows),
            "max_D_visc":max(r["viscous_dissipation"] for r in rows),
            "max_D_slip":max(r["slip_dissipation"] for r in rows),
            "max_kinetic_energy":max(r["E_kin"] for r in rows),
            "energy_drift":max(abs(r["E_total"]-first["E_total"]) for r in rows),
            "relative_energy_drift":max(abs(r["E_total"]-first["E_total"])/scale for r in rows),
            "absolute_budget_defect":max(abs(r["energy_budget_defect"]) for r in rows),
            "relative_budget_defect":max(abs(r["energy_budget_defect"])/scale for r in rows),
            "contact_displacement_max":max(r["contact_displacement_max"] for r in rows),
            "angle_drift_max":max(r["angle_drift_max"] for r in rows)}
        gates={"mass":metrics["mass_error"]<=1e-10,
            "certified_resolution":all(r["qualified_resolution"] for r in resolutions)}
        energy=validate_energy(rows)
        if perturb:
            gates.update(no_energy_growth=energy["energy_growth_ok"],energy_closure=energy["energy_budget_ok"])
            change=abs(energy["energy_change_final_J_per_m"])
            metrics["final_budget_relative_to_actual_change"]=(abs(energy["final_budget_defect_J_per_m"])/change if change else None)
            gates["nonzero_final_budget_over_actual_change"]=(change>0 and metrics["final_budget_relative_to_actual_change"]<=.05)
        else:
            gates.update(speed=metrics["max_speed"]<=1e-8,phase_drift=metrics["relative_phi_L2_drift"]<=1e-8,
                energy_drift=metrics["relative_energy_drift"]<=1e-8,
                absolute_budget=metrics["absolute_budget_defect"]<=1e-10,
                relative_budget=metrics["relative_budget_defect"]<=1e-8,
                CH_power=metrics["max_D_CH"]<=1e-16)
        result={"status":"complete","qualification_status":"passed" if all(gates.values()) else "failed",
            "gates":gates,"kind":"perturbation" if perturb else "preservation","scheme":scheme,"dt":dt,"t_end":t_end,
            "config":c.as_dict(),"provenance":solver.run_metadata,"equilibrium_source_fingerprint":meta["fingerprint"],
            "metrics":metrics,"initial":first,"final":rows[-1],"energy_validation":energy,
            "symmetric_energy_ratio":{"diagnostic_only":True,"E_floor":1e-10},
            "accepted_steps":solver.step_number,"runtime_s":time.perf_counter()-started,
            "qualification_scope":"single-history gates; perturbation also requires the independent temporal-series gate",
            "minimum_certified_cells":min(r["cells_across_transition_certified_min"] for r in rows),
            "perturbation":perturbation}
        if not perturb and not all(gates.values()):
            solver.initialize_prepared_equilibrium(equilibrium)
            result["failure_term_audit"]=transient_term_residuals(solver,meta["mu_star"])
        write_json(output/"summary.json",result)
        print(json.dumps({"output":str(output),"gates":gates,"metrics":metrics,"runtime_s":result["runtime_s"]}),flush=True)
        return result
    except BaseException as error:
        write_json(output/"summary.json",{"status":"failed","qualification_status":"failed","config":c.as_dict(),
            "error":str(error),"accepted_samples":len(rows),"runtime_s":time.perf_counter()-started})
        raise
    finally:
        if stream:
            stream.close()
