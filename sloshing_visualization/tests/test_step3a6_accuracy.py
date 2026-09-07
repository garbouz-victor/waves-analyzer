"""Small accuracy/control/transaction checks; never production 96x48 cases."""
from types import SimpleNamespace
import copy
import numpy as np
import pytest
from sloshing.multiphase.nonlinear_accuracy import (NonlinearControls, actual_controls,
    apply_controls, finite_arrays, rate_arrays, write_snapshot, read_snapshot)


class FakeSNES:
    def __init__(self):
        self.controls = dict(atol=1e-10, rtol=1e-9, stol=0., max_it=30)
    def setTolerances(self, **kwargs): self.controls.update(kwargs)
    def setOptionsPrefix(self, prefix): self.prefix=prefix
    def getTolerances(self): return tuple(self.controls[k] for k in ("rtol","atol","stol","max_it"))
    def getType(self): return "newtonls"
    def getLineSearch(self): return SimpleNamespace(getType=lambda:"bt")
    def getKSP(self): return SimpleNamespace(getType=lambda:"preonly",getPC=lambda:SimpleNamespace(
        getType=lambda:"lu",getFactorSolverType=lambda:"mumps"))
    def getConvergedReason(self): return 2
    def getIterationNumber(self): return 2
    def getFunctionNorm(self): return 2e-14


def test_controls_are_separate_from_physics_and_cannot_leak():
    from sloshing.multiphase.config import ModelConfig
    config=ModelConfig.from_json("configs/step3/benchmark_flat.json")
    before=config.fingerprint()
    moderate,tiny=FakeSNES(),FakeSNES()
    apply_controls(moderate,NonlinearControls(3e-13,3e-13),"moderate_")
    apply_controls(tiny,NonlinearControls(),"tiny_")
    assert actual_controls(moderate)["atol"]==3e-13
    assert actual_controls(tiny)==NonlinearControls().record()
    assert config.fingerprint()==before


def test_NewtonTrace_target_uses_actual_not_config_tolerances():
    from sloshing.multiphase.benchmarks.phase_rate_ch import NewtonTrace
    snes=FakeSNES(); apply_controls(snes,NonlinearControls(3e-13,3e-13),"trace_")
    trace=NewtonTrace.__new__(NewtonTrace)
    trace.engine=SimpleNamespace(problem=SimpleNamespace(solver=snes),solver=SimpleNamespace(config=None))
    trace.interval=1; trace.rows=[{"interval":1,"residual":.5}]
    result=trace.summary()
    assert result["effective_target"]==3e-13
    assert result["rtol_target"]==1.5e-13
    assert result["actual_controls"]==actual_controls(snes)


def test_retained_rate_snapshot_roundtrip_and_independent_arrays(tmp_path):
    old=np.array([1.,-.9]); rate=np.array([.123456789,-.987654321]); mu=np.array([.1,.2])
    arrays=rate_arrays(old,rate,mu,1e-12,rate_map=np.array([0,2]),mu_map=np.array([1,3]))
    original=rate.copy(); rate[:]=0.; old[:]=0.; mu[:]=0.
    assert np.array_equal(arrays["phase_rate"],original)
    assert not np.array_equal(arrays["phase_rate"],arrays["phase_rate_recovered_from_rounded_phi"])
    for i,a in enumerate(arrays.values()):
        assert all(not np.shares_memory(a,b) for b in list(arrays.values())[i+1:])
    record=write_snapshot(tmp_path,arrays,{"accepted":True,"execution_status":"complete"})
    restored,loaded=read_snapshot(tmp_path)
    assert record==loaded
    assert all(np.array_equal(v,restored[k]) for k,v in arrays.items())


@pytest.mark.parametrize("bad",[float("nan"),float("inf"),-float("inf")])
def test_finite_and_snapshot_reject_nonfinite(bad):
    assert not finite_arrays([1.,bad])
    with pytest.raises(ValueError): rate_arrays([1.],[bad],[.1],1e-12)


def test_cross_PASS_never_overrides_physical_FAIL_and_P1_cross_FAIL_can_proceed():
    from sloshing.multiphase.benchmarks.accuracy_ch import comparison_gates,may_continue_accuracy
    checks=comparison_gates({"phi_relative":0.,"mu_relative":0.,"D_relative":0.,
        "energy_absolute":0.,"mass_absolute":0.},{"mass":False},{"mass":True})
    assert all(checks["cross_formulation_checks"].values()) and not checks["passed"]
    status={"nonlinear_solver_status":"converged","physical_checks":"passed","cross":"failed"}
    assert may_continue_accuracy(status)
    status["physical_checks"]="failed"
    assert not may_continue_accuracy(status)


def test_strict_tiny_FEM_moderate_comparison_and_matched_guess(tmp_path):
    from test_step3a5_phase_rate import tiny_solver
    from sloshing.multiphase.isolated_ch import IsolatedCH
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate,PhaseMetric
    from sloshing.multiphase.benchmarks.accuracy_ch import residual_audit
    fields={}
    for name in ("absolute","semidiscrete","zero"):
        s=tiny_solver()
        initial=lambda x:.2+.1*np.cos(2*x[0])*np.cos(3*x[1])
        engine=IsolatedCH(s,initial) if name=="absolute" else IsolatedCHPhaseRate(s,initial,initial_guess=name)
        mu_old=s.state.sub(3).collapse().x.array.copy()
        apply_controls(engine.problem.solver,NonlinearControls(3e-13,3e-13),"tiny_accuracy_"+name+"_")
        # Coarse CI mesh h=.25 has far larger mass rows than production h~.01.
        # Use a toy moderate dt above its absolute cancellation floor. The
        # production dt=1e-6 and target tiny dt are NOT changed by this test.
        engine.step(1e-3)
        assert engine.problem.solver.getConvergedReason()>0
        fields[name]=(s.state.sub(2).collapse().x.array.copy(),s.state.sub(3).collapse().x.array.copy())
        metric=PhaseMetric(s,s.state.sub(2).collapse().function_space)
        if name!="absolute":
            folder=tmp_path/name; folder.mkdir()
            arrays={**engine.last_snapshot,"mu_old":mu_old}
            write_snapshot(folder,arrays,{"accepted":True,"execution_status":"complete"})
            restored,_=read_snapshot(folder)
            assert np.array_equal(restored["phase_rate"],engine.state.x.array[engine.rate_map])
            audit=residual_audit(s,restored,metric,engine.state.function_space)
            assert audit["retained_rate_UFL"]["algebraic_norm"]<3e-13
            assert audit["chemical_rate_expression"]["algebraic_norm"]<3e-13
    for name in ("semidiscrete","zero"):
        for i in (0,1):
            assert metric.norm(fields[name][i]-fields["absolute"][i])/metric.norm(fields["absolute"][i])<1e-10


def test_converged_SNES_material_failure_preserves_every_accepted_field(monkeypatch):
    from test_step3a5_phase_rate import tiny_solver
    s=tiny_solver()
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    from sloshing.multiphase.diagnostics import Diagnostics
    import sloshing.multiphase.nonlinear_accuracy as helpers
    engine=IsolatedCHPhaseRate(s,lambda x:np.ones(x.shape[1]),initial_guess="zero")
    engine.step(1e-12)  # Exercise preservation of nonempty accepted history/rate.
    diag=Diagnostics(s); diag.measure()
    before={k:getattr(s,k).x.array.copy() for k in ("state","old","older")}
    attrs={k:getattr(s,k) for k in ("time","step_number","current_dt","current_phase")}
    cumulative=copy.deepcopy(diag.cumulative); history=copy.deepcopy(engine.history)
    previous=engine.last_accepted_rate
    snapshot={k:v.copy() for k,v in engine.last_snapshot.items()}
    def fail(*args): raise RuntimeError("injected material candidate failure")
    monkeypatch.setattr(helpers,"check_candidate_material",fail)
    with pytest.raises(RuntimeError,match="material candidate"):
        engine.step(1e-12)
    assert engine.problem.solver.getConvergedReason()>0
    assert all(np.array_equal(getattr(s,k).x.array,v) for k,v in before.items())
    assert all(getattr(s,k)==v for k,v in attrs.items())
    assert engine.last_accepted_rate is previous
    assert engine.history==history and diag.cumulative==cumulative
    assert all(np.array_equal(engine.last_snapshot[k],v) for k,v in snapshot.items())


def test_physical_validator_failure_does_not_publish():
    from test_step3a5_phase_rate import tiny_solver
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    s=tiny_solver(); engine=IsolatedCHPhaseRate(s,lambda x:np.ones(x.shape[1]),initial_guess="zero")
    before=s.state.x.array.copy()
    def fail(*args): raise RuntimeError("physical gate")
    with pytest.raises(RuntimeError,match="physical gate"):
        engine.step(1e-12,candidate_validator=fail)
    assert s.step_number==0 and s.time==0. and not engine.history
    assert np.array_equal(s.state.x.array,before)
