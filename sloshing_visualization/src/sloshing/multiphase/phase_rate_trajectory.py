"""Validation-only multi-step nonlinear isolated CH, using retained phase rates.

The physical solver, its residual, Newton defaults and diagnostic mathematics
are unchanged. This driver supplies transactional history and immutable clocks.
"""
import copy
from pathlib import Path
import time
import numpy as np
from .linearized_ch import array_hash
from .nonlinear_accuracy import NonlinearControls, actual_controls, apply_controls
from .phase_rate_history import (archive, load_archive, atomic_json, read_json, file_hash,
    Journal, WriterLock, recover_prefix, validate_complete, VERSION)
from .performance import Performance


class ScientificFailure(RuntimeError):
    pass


class CostStop(RuntimeError):
    pass


def finite_tree(value):
    if isinstance(value,dict): return all(finite_tree(v) for v in value.values())
    if isinstance(value,(list,tuple,np.ndarray)): return all(finite_tree(v) for v in value)
    if isinstance(value,(float,int,np.number)): return bool(np.isfinite(value))
    return True


def rejected_metadata(value):
    """Nonfinite diagnostics stay explicit; never serialize invalid JSON numbers."""
    if isinstance(value,dict): return {k:rejected_metadata(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [rejected_metadata(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value): return "nonfinite:"+str(value)
    return value


def physical_gates(row, resolution, work, policy, target_mass, area, crossings):
    gates={"finite":finite_tree((row,resolution,work)),
        "mass_initial":abs(row["mass_error_relative_to_domain"])<=policy["nonlinear_mass_domain"],
        "mass_target":abs(row["phase_mass"]-target_mass)/area<=policy["nonlinear_mass_domain"],
        "certified_resolution":resolution.get("cells_across_transition_certified_min",0.)>=policy["certified_cells"],
        "topology":row["wall_crossings"]==crossings,
        "energy_non_growth":row["delta_E_interval"]<=policy["energy_growth_J_per_m"]}
    if work is not None:
        gates.update(weak_work=abs(work["weak_work_defect"])<=policy["weak_work_J_per_m"],
            decomposition=abs(work["decomposition_roundoff"])<=policy["weak_work_J_per_m"],
            work_identity=abs(work["weak_work_identity_roundoff"])<=policy["weak_work_J_per_m"])
    return gates


def diagnostic_state(diag):
    return copy.deepcopy({k:getattr(diag,k) for k in ("initial","previous","cumulative")})


class RateTrajectory:
    def __init__(self, engine, schedule, folder, identity, policy, target_mass, *,
                 expected_crossings=2, resume=False, source_guard=lambda:None, expanded_audit=True):
        from dolfinx import fem
        from .diagnostics import Diagnostics
        from .be_work_diagnostics import BEWorkDiagnostics
        from .ch_validation_diagnostics import WallTrace
        self.engine=engine; self.s=engine.solver; self.schedule=schedule; self.folder=Path(folder)
        self.identity=copy.deepcopy(identity); self.policy=policy; self.target_mass=target_mass
        self.expected_crossings=expected_crossings; self.source_guard=source_guard
        self.expanded_audit=expanded_audit; self.lock=WriterLock(folder)
        self.closed=False; self.source_guard()
        if self.s.comm.size!=1: raise ValueError("Scientific checkpoint layout is serial only")
        self.observer=copy.copy(self.s)
        self.observer.state=fem.Function(self.s.space); self.observer.older=fem.Function(self.s.space)
        self.diag=Diagnostics(self.observer); self.work=BEWorkDiagnostics(self.observer)
        self.wall=WallTrace(self.observer)
        self.work_sums={k:0. for k in ("BE_endpoint_integral","B_BE","weak_work_defect","Delta_E")}
        self.area=engine.metric.area; self.controls=NonlinearControls()
        apply_controls(engine.problem.solver,self.controls,"step3a7_trajectory_")
        self.trace=[]
        def monitor(snes,i,norm):
            ksp=snes.getKSP()
            self.trace.append({"iteration":i,"residual":float(norm),
                "KSP_iterations":int(ksp.getIterationNumber()) if i else None,
                "KSP_reason":int(ksp.getConvergedReason()) if i else None})
        engine.problem.solver.setMonitor(monitor)
        if resume:
            if (self.folder/"COMPLETE.json").exists(): raise ValueError("Complete level cannot be continued")
            latest=read_json(self.folder/"LATEST.json")
            arrays,meta=load_archive(self.folder/latest["checkpoint"],identity)
            if meta["semantic_sha256"]!=latest["checkpoint_sha256"] or not meta["accepted"]:
                raise ValueError("Invalid accepted checkpoint pointer")
            self.journal=recover_prefix(self.folder,meta)
            self.restore(arrays,meta)
        else:
            if (self.folder/"scalar_history.jsonl").exists(): raise ValueError("Fresh trajectory would overwrite history")
            self.journal=Journal(self.folder/"scalar_history.jsonl")
            self.initialize()
        self.timing_journal=Journal(self.folder/"timing_history.jsonl")

    def close(self):
        if not self.closed:
            self.engine.problem.solver.cancelMonitor(); self.lock.close(); self.closed=True

    def update_observer(self, state, older, time_value, step, dt):
        self.observer.state.x.array[:]=state; self.observer.older.x.array[:]=older
        self.observer.state.x.scatter_forward(); self.observer.older.x.scatter_forward()
        self.observer.time=time_value; self.observer.step_number=step; self.observer.current_dt=dt

    def observe(self, diag):
        from .diagnostics import interface_resolution
        started=time.perf_counter(); row=diag.measure(); diag_s=time.perf_counter()-started
        started=time.perf_counter(); resolution=interface_resolution(self.observer); resolution_s=time.perf_counter()-started
        started=time.perf_counter(); crossings=self.wall.crossings(); topology_s=time.perf_counter()-started
        row.update(wall_crossings=len(crossings),wall_crossing_positions=crossings,
            cells_across_transition_certified_min=resolution.get("cells_across_transition_certified_min",0.),
            target_mass_error_relative_to_domain=(row["phase_mass"]-self.target_mass)/self.area)
        # Keep exact worst-cell witness; omit potentially huge fallback-index lists.
        compact={k:v for k,v in resolution.items() if k!="diameter_fallback_local_cells_by_rank"}
        return row,compact,{"diagnostics_s":diag_s,"resolution_s":resolution_s,"topology_s":topology_s}

    def initialize(self):
        if self.s.step_number!=0 or self.s.time!=0.: raise ValueError("Fresh level must start at t=0")
        self.update_observer(self.s.state.x.array,self.s.state.x.array,0.,0,None)
        row,res,timing=self.observe(self.diag)
        gates=physical_gates(row,res,None,self.policy,self.target_mass,self.area,self.expected_crossings)
        if not all(gates.values()): raise ScientificFailure("Initial physical gates: "+str(gates))
        self.journal.append({**row,"step":0,"dt":None,"block":None,"solver":None,"work":None,
            "physical_gates":gates,"resolution":res,"timing":timing,"benchmark":"nonlinear isolated CH; u=0 by construction"})
        self.save_field(0); self.checkpoint()

    def physical_arrays(self):
        e=self.engine; s=self.s
        return {**{k:getattr(s,k).x.array.copy() for k in ("state","old","older")},
            "phi":s.state.x.array[e.physical_phi_map].copy(),"mu":s.state.x.array[e.physical_mu_map].copy(),
            "phase_rate":(e.last_accepted_rate if e.last_accepted_rate is not None else e.initial_semidiscrete_rate).copy(),
            "physical_phi_map":e.physical_phi_map.copy(),"physical_mu_map":e.physical_mu_map.copy(),
            "rate_map":e.rate_map.copy(),"mu_map":e.mu_map.copy()}

    def save_field(self, step):
        a=self.physical_arrays()
        return archive(self.folder/"fields"/f"{step:06d}",{k:a[k] for k in ("phi","mu","phase_rate")},
            {"identity":self.identity,"step":step,"time":self.s.time,"dt":self.s.current_dt,
             "rate_role":"retained accepted BE unknown" if step else "semidiscrete initial guess, not accepted BE rate"})

    def checkpoint(self):
        self.source_guard(); s=self.s; step=s.step_number
        path=self.folder/"checkpoints"/f"{step:06d}"
        if path.exists():
            arrays,record=load_archive(path,self.identity)
            if not np.array_equal(arrays["state"],s.state.x.array): raise ValueError("Checkpoint collision")
        else:
            record=archive(path,self.physical_arrays(),{"identity":self.identity,"accepted":True,
                "accepted_steps":step,"time":s.time,"current_dt":s.current_dt,"current_phase":s.current_phase,
                "has_accepted_rate":self.engine.last_accepted_rate is not None,
                "diagnostics":diagnostic_state(self.diag),"work_sums":copy.deepcopy(self.work_sums),
                "journal_prefix":self.journal.prefix()})
        atomic_json(self.folder/"LATEST.json",{"checkpoint":str(path.relative_to(self.folder)),
            "checkpoint_sha256":record["semantic_sha256"]})

    def restore(self, arrays, meta):
        s,e=self.s,self.engine
        for key in ("physical_phi_map","physical_mu_map","rate_map","mu_map"):
            if not np.array_equal(arrays[key],getattr(e,key)): raise ValueError("Checkpoint DOF layout mismatch")
        for key in ("state","old","older"):
            getattr(s,key).x.array[:]=arrays[key]; getattr(s,key).x.scatter_forward()
        s.time=meta["time"]; s.step_number=meta["accepted_steps"]
        s.current_dt=meta["current_dt"]; s.current_phase=meta["current_phase"]
        if s.time!=(self.schedule.intervals[s.step_number-1].end if s.step_number else 0.):
            raise ValueError("Checkpoint clock/schedule mismatch")
        e.last_accepted_rate=arrays["phase_rate"].copy() if meta["has_accepted_rate"] else None
        e.state.x.array[e.rate_map]=arrays["phase_rate"]
        e.state.x.array[e.mu_map]=arrays["mu"]
        e.state.x.scatter_forward()
        e.history=[r["solver"] for r in self.journal.rows[1:]]
        for key,value in meta["diagnostics"].items(): setattr(self.diag,key,copy.deepcopy(value))
        self.work_sums=copy.deepcopy(meta["work_sums"])
        self.update_observer(arrays["state"],arrays["older"],s.time,s.step_number,s.current_dt)

    def advance(self):
        s,e=self.s,self.engine; interval=self.schedule.intervals[s.step_number]
        if s.time!=interval.start: raise ValueError("Schedule advancement/clock mismatch")
        if actual_controls(e.problem.solver)!=self.controls.record(): raise ScientificFailure("Changed Newton controls")
        old_mu=s.state.x.array[e.physical_mu_map].copy()
        guess=(e.last_accepted_rate if e.last_accepted_rate is not None else e.initial_semidiscrete_rate).copy()
        guess_stats={"policy":"previous accepted rate, unscaled" if s.step_number else "semidiscrete rate",
            "sha256":array_hash(guess),"L2":e.metric.norm(guess),"max_abs":float(np.max(abs(guess)))}
        self.trace=[]; e.performance=Performance(); started=time.perf_counter(); pending={}
        def validate(candidate, arrays, dt, endpoint):
            self.update_observer(candidate.x.array,s.state.x.array,endpoint,interval.step,dt)
            trial=copy.copy(self.diag)
            for key,value in diagnostic_state(self.diag).items(): setattr(trial,key,value)
            row,res,timing=self.observe(trial)
            before=time.perf_counter(); work=self.work.measure(dt,self.diag.previous,row)
            timing["BE_work_s"]=time.perf_counter()-before
            gates=physical_gates(row,res,work,self.policy,self.target_mass,self.area,self.expected_crossings)
            gates["physical_arrays_finite"]=bool(np.isfinite(candidate.x.array).all())
            pending.update(row=row,resolution=res,timing=timing,work=work,gates=gates,diagnostics=trial)
            if not all(gates.values()): raise ScientificFailure("Candidate physical gates: "+str(gates))
            if actual_controls(e.problem.solver)!=self.controls.record(): raise ScientificFailure("Changed actual controls")
            audit_started=time.perf_counter()
            if self.expanded_audit and interval.step in self.schedule.audit_indices:
                from .benchmarks.accuracy_ch import residual_audit
                pending["expanded_audit"]=residual_audit(self.observer,{**arrays,"mu_old":old_mu},e.metric,e.space)
            timing["expanded_residual_audit_s"]=time.perf_counter()-audit_started
        try:
            result=e.step(interval.dt,interval.end,candidate_validator=validate)
        except Exception as exc:
            # The existing solver publishes only after our callback succeeds.
            self.update_observer(s.state.x.array,s.older.x.array,s.time,s.step_number,s.current_dt)
            error={"execution_status":"failed","accepted":False,"failed_interval":interval.step,
                "last_accepted_step":s.step_number,"dt":interval.dt,"error":str(exc),
                "SNES_reason":int(e.problem.solver.getConvergedReason()),"trace":self.trace,
                "actual_controls":actual_controls(e.problem.solver),"physical_gates":pending.get("gates"),
                "identity":self.identity,"attempt_wall_s":time.perf_counter()-started}
            error=rejected_metadata(error)
            atomic_json(self.folder/"FAILED.json",error)
            # These are explicitly different rejected Function/iterate states.
            rejected={"last_residual_Function":e.state.x.array.copy(),"PETSc_iterate":e.problem.x.array.copy(),
                "last_accepted_physical_state":s.state.x.array.copy(),"phi_old":e.phi_old.x.array.copy()}
            archive(self.folder/"rejected_attempt",rejected,error)
            raise ScientificFailure(str(exc)) from exc
        self.diag=pending["diagnostics"]
        for key in self.work_sums: self.work_sums[key]+=pending["work"][key]
        timing=pending["timing"]
        timing["SNES_s"]=sum(r["wall_s"] for r in e.performance.records["nonlinear_SNES"])
        initial=self.trace[0]["residual"] if self.trace else None
        result.update(trace=copy.deepcopy(self.trace),actual_controls=actual_controls(e.problem.solver),
            initial_residual=initial,effective_target=max(self.controls.atol,self.controls.rtol*initial) if initial is not None else None,
            guess=guess_stats,retained_rate_sha256=array_hash(e.last_accepted_rate))
        row={**pending["row"],"step":interval.step,"dt":interval.dt,"block":interval.block,
            "DeltaF_from_start":pending["row"]["E_total"]-self.diag.initial["E_total"],
            "parent_step":interval.parent_step,"clock_minus_nominal_dt":interval.end-interval.start-interval.dt,
            "solver":result,"work":pending["work"],"work_sums":copy.deepcopy(self.work_sums),
            "physical_gates":pending["gates"],"resolution":pending["resolution"],"timing":timing}
        io_started=time.perf_counter()
        self.journal.append(row)
        if "expanded_audit" in pending:
            atomic_json(self.folder/"residual_audits"/f"{interval.step:06d}"/"norms.json",pending["expanded_audit"])
        if interval.step in self.schedule.field_indices: self.save_field(interval.step)
        if interval.step in self.schedule.checkpoint_indices: self.checkpoint()
        timing["archival_s"]=time.perf_counter()-io_started
        timing["total_s"]=time.perf_counter()-started
        timing["scope"]="entire accepted interval including journal/fields/checkpoint; excludes this timing-record write"
        self.timing_journal.append({"step":interval.step,"timing":copy.deepcopy(timing)})
        return row

    def run_until(self, count, *, wall_remaining_s=float("inf")):
        if not self.s.step_number<=count<=self.schedule.nsteps: raise ValueError("Invalid immutable prefix target")
        started=time.perf_counter()
        while self.s.step_number<count:
            if time.perf_counter()-started>=wall_remaining_s:
                self.checkpoint(); raise CostStop("Frozen wall budget exhausted; incomplete trajectory")
            row=self.advance()
            if row["step"]%10==0 or row["step"]==count:
                print(f"step={row['step']}/{self.schedule.nsteps} t={row['time']:.10g} "
                      f"Newton={row['solver']['snes_iterations']} D={row['CH_dissipation']:.9g} "
                      f"mass={row['mass_error_relative_to_domain']:.3g} "
                      f"cells={row['cells_across_transition_certified_min']:.9g}",flush=True)
        self.checkpoint()
        if count==self.schedule.nsteps:
            self.source_guard()
            marker={"identity":self.identity,"accepted_steps":count,
                "expected_steps":self.schedule.nsteps,"final_time":self.s.time,
                "final_state_sha256":array_hash(self.s.state.x.array),
                "scalar_history_sha256":file_hash(self.journal.path)}
            validate_complete(self.folder,self.identity,self.schedule,marker=marker)
            atomic_json(self.folder/"COMPLETE.json",marker)
        return self.status()

    def status(self):
        rows=self.journal.rows; last=rows[-1]
        minimum=min(rows,key=lambda r:r["cells_across_transition_certified_min"])
        return {"execution_status":"failed" if (self.folder/"FAILED.json").exists() else "complete" if (self.folder/"COMPLETE.json").exists() else "incomplete",
            "expected_steps":self.schedule.nsteps,"accepted_steps":self.s.step_number,
            "solver_status":"converged" if not (self.folder/"FAILED.json").exists() else "failed",
            "physical_status":"passed" if all(all(r["physical_gates"].values()) for r in rows) else "failed",
            "max_mass_domain":max(abs(r["mass_error_relative_to_domain"]) for r in rows),
            "min_certified_cells":minimum["cells_across_transition_certified_min"],
            "minimum_resolution_time":minimum["time"],"minimum_resolution_witness":minimum["resolution"]["worst_cell"],
            "final_F":last["E_total"],"DeltaF":last["E_total"]-rows[0]["E_total"],
            "integrated_D_CH":last["cumulative_CH_dissipation"],"final_budget_defect":last["energy_budget_defect"],
            "complete_marker":str(self.folder/"COMPLETE.json") if (self.folder/"COMPLETE.json").exists() else None}
