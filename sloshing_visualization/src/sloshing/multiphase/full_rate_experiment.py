"""Independent one-interval bridge experiments, never a temporal refinement."""
import copy
import time
import numpy as np
from .full_phase_rate import CHNSPhaseRateBE
from .phase_rate import PhaseMetric,snes_options
from .nonlinear_accuracy import (apply_controls,actual_controls,synchronized_solution,
    check_candidate_material,publish_candidate)
from .full_rate_diagnostics import FullObserver,assemble_vector
from .phase_rate_trajectory import physical_gates,diagnostic_state,ScientificFailure
from .performance import Performance
from .linearized_ch import array_hash


class AbsoluteBEReference:
    """Private candidate using the UNMODIFIED historical full weak_form.residual.

    This is solely an accuracy comparison adapter, not a production alternative.
    CHNSSolver.advance and the historical solver remain unmodified.
    """
    def __init__(self,solver):
        import ufl
        from dolfinx import fem
        from dolfinx.fem.petsc import NonlinearProblem
        from petsc4py import PETSc
        from .weak_form import residual
        self.solver=solver; self.space=solver.space
        self.maps=[np.asarray(self.space.sub(i).collapse()[1]) for i in range(4)]
        self.u_map,self.pi_map,self.rate_map,self.mu_map=self.maps
        self.physical_phi_map=self.rate_map; self.physical_mu_map=self.mu_map
        self.state=fem.Function(self.space); self.state.x.array[:]=solver.state.x.array
        self.state.x.scatter_forward()
        self.phi_old=solver.state.sub(2).collapse()
        self.metric=PhaseMetric(solver,self.phi_old.function_space)
        self.coefficients=[fem.Constant(solver.mesh,PETSc.ScalarType(v)) for v in (1.,-1.,0.)]
        self.F=residual(self.state,solver.old,solver.older,self.coefficients,solver.config,solver.tags)
        self.J=ufl.derivative(self.F,self.state,ufl.TrialFunction(self.space))
        self.problem=NonlinearProblem(self.F,self.state,J=self.J,bcs=solver.bcs,
            petsc_options_prefix="step3a8_absolute_reference_",petsc_options=snes_options(solver.config))
        self.performance=Performance()

    def step(self,dt,*,candidate_validator=None):
        from dolfinx import fem
        s=self.solver
        if s.step_number!=0: raise ValueError("Absolute reference is one independent interval only")
        for coefficient,value in zip(self.coefficients,(1/dt,-1/dt,0.)): coefficient.value=value
        with self.performance.measure("nonlinear_SNES",dt=dt): self.problem.solve()
        synchronized_solution(self.problem,self.state)
        candidate=fem.Function(s.space); candidate.x.array[:]=self.state.x.array; candidate.x.scatter_forward()
        check_candidate_material(s,candidate)
        arrays={"physical_state":candidate.x.array.copy(),"mixed_Function":self.state.x.array.copy(),
            "PETSc_solution":self.problem.x.array.copy(),"phi_old":self.phi_old.x.array.copy(),"dt":np.array(dt),
            "u_new":candidate.x.array[self.u_map].copy(),"pi_new":candidate.x.array[self.pi_map].copy(),
            "phi_new":candidate.x.array[self.rate_map].copy(),"mu_new":candidate.x.array[self.mu_map].copy()}
        if candidate_validator is not None: candidate_validator(candidate,arrays,dt,dt)
        publish_candidate(s,candidate,dt,dt,"absolute_full_BE_reference")
        self.last_snapshot=arrays


class OneStepAudit:
    def __init__(self,engine,inherited,full_policy,target_mass,*,crossings=2):
        from dolfinx import fem
        from .diagnostics import Diagnostics
        from .be_work_diagnostics import BEWorkDiagnostics
        from .ch_validation_diagnostics import WallTrace
        self.engine=engine; self.s=engine.solver; self.policy=inherited; self.full_policy=full_policy
        self.target_mass=target_mass; self.crossings=crossings; self.timing={}
        self.observer=copy.copy(self.s)
        self.observer.state=fem.Function(self.s.space); self.observer.older=fem.Function(self.s.space)
        self.observer.state.x.array[:]=self.s.state.x.array
        self.observer.older.x.array[:]=self.s.state.x.array
        self.observer.state.x.scatter_forward(); self.observer.older.x.scatter_forward()
        self.diag=Diagnostics(self.observer); self.work=BEWorkDiagnostics(self.observer)
        self.wall=WallTrace(self.observer); self.full=FullObserver(self.observer,engine.metric)
        self.initial,self.initial_resolution=self.sample(self.diag)
        initial_gates=physical_gates(self.initial,self.initial_resolution,None,inherited,target_mass,engine.metric.area,crossings)
        if not all(initial_gates.values()) or not all(self.initial["full_physical_checks"].values()):
            raise ScientificFailure("Initial physical failure: "+str(initial_gates)+str(self.initial))
        self.accepted=False

    def sample(self,diag):
        from .diagnostics import interface_resolution
        started=time.perf_counter(); row=diag.measure(); self.timing["diagnostics_s"]=time.perf_counter()-started
        started=time.perf_counter(); resolution=interface_resolution(self.observer)
        self.timing["resolution_s"]=time.perf_counter()-started
        row["wall_crossings"]=len(self.wall.crossings())
        row["cells_across_transition_certified_min"]=resolution.get("cells_across_transition_certified_min",0.)
        started=time.perf_counter(); row.update(self.full.measure(row,self.full_policy,expanded=True))
        self.timing["full_constraint_audit_s"]=time.perf_counter()-started
        return row,resolution

    def validate(self,candidate,arrays,dt,endpoint):
        self.observer.state.x.array[:]=candidate.x.array; self.observer.state.x.scatter_forward()
        self.observer.time=endpoint; self.observer.current_dt=dt; self.observer.step_number=1
        trial=copy.copy(self.diag)
        for k,v in diagnostic_state(self.diag).items(): setattr(trial,k,v)
        row,resolution=self.sample(trial)
        started=time.perf_counter(); work=self.work.measure(dt,self.initial,row)
        self.timing["BE_work_s"]=time.perf_counter()-started
        gates=physical_gates(row,resolution,work,self.policy,self.target_mass,self.engine.metric.area,self.crossings)
        self.result={"initial":self.initial,"physical":row,"resolution":resolution,"work":work,
            "physical_gates":gates,"full_physical_checks":row["full_physical_checks"],"timing":self.timing.copy()}
        if not all(gates.values()) or not all(row["full_physical_checks"].values()):
            raise ScientificFailure("One-step physical failure: "+str(self.result))
        self.diag=trial; self.accepted=True


def solve_one(engine,dt,controls,audit,prefix):
    snes=engine.problem.solver; apply_controls(snes,controls,prefix)
    trace=[]
    def monitor(snes,i,norm):
        ksp=snes.getKSP()
        trace.append({"iteration":i,"residual":float(norm),
            "KSP_iterations":int(ksp.getIterationNumber()) if i else None,
            "KSP_reason":int(ksp.getConvergedReason()) if i else None})
    snes.setMonitor(monitor); engine.performance=Performance()
    before=actual_controls(snes); guess=engine.state.x.array.copy()
    physical_guess=guess.copy()
    if isinstance(engine,CHNSPhaseRateBE):
        physical_guess[engine.rate_map]=engine.phi_old.x.array+dt*guess[engine.rate_map]
    started=time.perf_counter(); error=None
    def validate(candidate,arrays,dt,endpoint):
        if actual_controls(snes)!=before: raise ScientificFailure("Solver controls changed before publication")
        audit.validate(candidate,arrays,dt,endpoint)
    try:
        engine.step(dt,candidate_validator=validate)
        if actual_controls(snes)!=before: raise ScientificFailure("Solver controls changed")
    except Exception as exc: error=str(exc)
    finally: snes.cancelMonitor()
    result={"execution_status":"complete" if error is None else "failed", "error":error,
        "accepted":error is None and audit.accepted,"nonlinear_solver_status":"converged" if snes.getConvergedReason()>0 else "diverged",
        "SNES_reason":int(snes.getConvergedReason()),"iterations":int(snes.getIterationNumber()),
        "initial_residual":trace[0]["residual"] if trace else None,"final_residual":float(snes.getFunctionNorm()),
        "controls_before":before,"controls_after":actual_controls(snes),"trace":trace,
        "dt":dt,"solve_wall_s":time.perf_counter()-started,"physical_initial_guess_sha256":array_hash(physical_guess),
        "initial_unknown_sha256":array_hash(guess),"physical_old_sha256":array_hash(engine.phi_old.x.array),
        "initial_physical":audit.initial,
        "candidate_audit":getattr(audit,"result",None),"performance":engine.performance.records}
    result["effective_target"]=max(controls.atol,controls.rtol*result["initial_residual"]) if trace else None
    if result["accepted"]:
        assembled=assemble_vector(engine.F)
        result["residual_blocks"]={k:float(np.linalg.norm(assembled[m])) for k,m in zip(("momentum","continuity","phase","chemical"),engine.maps)}
        result["phase_Riesz"]=engine.metric.riesz_norm(assembled[engine.rate_map])
        result["chemical_Riesz"]=engine.metric.riesz_norm(assembled[engine.mu_map])
        result["raw_residual_note"]="Unconstrained weak rows; momentum includes essential-BC reactions, not the SNES norm"
    return result


def compare_states(a,b,metric,observer,absolute_row,rate_row,policy):
    area=metric.area
    errors={"phi_relative":metric.norm(b["phi_new"]-a["phi_new"])/metric.norm(a["phi_new"]),
        "mu_relative":metric.norm(b["mu_new"]-a["mu_new"])/metric.norm(a["mu_new"]),
        "u_scaled":observer.velocity_metric.norm(b["u_new"]-a["u_new"])/(observer.U_ref*np.sqrt(area)),
        "pi_scaled":observer.pressure_metric.norm(b["pi_new"]-a["pi_new"])/(observer.P_ref*np.sqrt(area)),
        "D_relative":abs(rate_row["CH_dissipation"]/absolute_row["CH_dissipation"]-1),
        "energy_absolute":abs(rate_row["E_total"]-absolute_row["E_total"]),
        "mass_absolute":abs(rate_row["phase_mass"]-absolute_row["phase_mass"])}
    limits={k:(policy["moderate_u_pi_scaled"] if k.endswith("scaled") else
        policy["moderate_energy_mass_absolute"] if k.endswith("absolute") else policy["moderate_phi_mu_D_relative"]) for k in errors}
    gates={k:bool(np.isfinite(v) and v<=limits[k]) for k,v in errors.items()}
    return {"errors":errors,"limits":limits,"gates":gates,"passed":all(gates.values())}
