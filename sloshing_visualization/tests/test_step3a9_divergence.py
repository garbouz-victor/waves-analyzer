import copy
import numpy as np
import pytest
from sloshing.multiphase.divergence_policy import (AUDIT,spatial_decision,production_policy,
    projection_checks,divergence_gates,bridge_divergence)


def rows():
    return [dict(mesh_shape=shape,h=h,strong_divergence_L2=h**2,
        weak_continuity_Riesz=1e-20,projection_agrees=True,pythagorean_pass=True,
        physical_pass=True,SNES_reason=2,KSP_success=True)
        for shape,h in zip(AUDIT["meshes"],[.0125,.01,1/120])]


def test_predeclared_selection_and_no_override():
    result=spatial_decision(rows()); assert result["passed"]
    assert result["p_div"]==pytest.approx(2.)
    for key,value in [("weak_continuity_Riesz",2e-12),("projection_agrees",False),
            ("pythagorean_pass",False),("physical_pass",False),("SNES_reason",-3),("KSP_success",False)]:
        r=rows();r[1][key]=value
        assert not spatial_decision(r)["passed"]
    r=rows();r[2]["strong_divergence_L2"]=1.
    assert not spatial_decision(r)["passed"]
    r=rows()
    for v in r:v["strong_divergence_L2"]=v["h"]**.3
    assert not spatial_decision(r)["passed"]
    with pytest.raises(ValueError):production_policy(spatial_decision(r))
    assert not spatial_decision(rows()[:2])["passed"]


def test_new_policy_retires_cap_not_continuity_or_cost():
    p=production_policy(spatial_decision(rows()))
    assert "strong_divergence_max" not in p
    assert not p["historical_metadata"]["active_hard_gate"]
    assert p["historical_metadata"]["STEP3A8_envelope"]==3.6880704671415807e-6
    assert p["weak_continuity"]==1e-12
    assert (p["preliminary_wall_s"],p["full_L0_wall_s"],p["total_wall_s"])==(1800,14400,19800)


def test_production_finite_weak_fraction_and_same_root():
    r=dict(strong_divergence_L2=1e-4,projected_divergence_L2=1e-20,
        orthogonal_divergence_L2=1e-4,weak_continuity_Riesz=1e-20,projection_fraction=1e-16,
        projection_agrees=True,pythagorean_pass=True)
    assert all(divergence_gates(r,production=True).values())
    for key,value in [("weak_continuity_Riesz",2e-12),("projection_fraction",2e-6),
            ("strong_divergence_L2",np.nan),("strong_divergence_L2",np.inf)]:
        bad=dict(r,**{key:value})
        assert not all(divergence_gates(bad,production=True).values())
    assert bridge_divergence(4e-6,4.0001e-6)["passed"]
    assert not bridge_divergence(4e-6,4.001e-6)["passed"]
    check=projection_checks(1.,.6,.6,.8)
    assert check["projection_agrees"] and check["pythagorean_pass"]
    import json
    json.dumps(check,allow_nan=False)


def tiny_solver():
    pytest.importorskip("dolfinx")
    from test_step3a8_full_phase_rate import tiny_solver as create
    return create()


def test_independent_pressure_projection_on_nontrivial_velocity():
    import math
    s=tiny_solver()
    from sloshing.multiphase.full_rate_diagnostics import FullObserver,assemble_vector
    from sloshing.multiphase.phase_rate import PhaseMetric
    from sloshing.multiphase.divergence_audit import DivergenceAudit
    s.state.sub(0).interpolate(lambda x:np.array([np.sin(2*x[0])*x[1],np.cos(3*x[1])*x[0]]))
    s.state.x.scatter_forward()
    observer=FullObserver(s,PhaseMetric(s,s.space.sub(2).collapse()[0]))
    weak=observer.pressure_metric.riesz_norm(assemble_vector(observer.continuity))
    audit=DivergenceAudit(s); r=audit.measure(weak)
    assert r["projection_agrees"] and r["pythagorean_pass"],r
    assert r["projected_divergence_L2"]==pytest.approx(weak,rel=1e-12)
    a,stats=audit.cells()
    assert stats["sum_divergence_squared"]==pytest.approx(r["strong_divergence_L2"]**2,rel=1e-12)
    assert not np.shares_memory(audit.M.data,observer.pressure_metric.M.data)


def test_Taylor_Hood_weak_zero_not_strong_zero():
    from sloshing.multiphase.divergence_audit import kernel_demonstration
    s=tiny_solver(); before=s.state.x.array.copy()
    r=kernel_demonstration(s)
    assert r["passed"],r
    assert np.array_equal(s.state.x.array,before)


def test_new_policy_candidate_weak_failure_keeps_accepted_state():
    s=tiny_solver()
    from sloshing.multiphase.full_phase_rate import CHNSPhaseRateBE
    from sloshing.multiphase.divergence_observer import DivergenceOneStepAudit
    from sloshing.multiphase.full_rate_experiment import solve_one
    from sloshing.multiphase.full_rate_policy import CONTROLS
    from sloshing.multiphase.phase_rate_history import read_json
    e=CHNSPhaseRateBE(s); p=production_policy(spatial_decision(rows()))
    inherited=read_json("validation_results/step3a4/policy/policy.json")["policy"]
    mass=float(e.metric.mass@s.state.x.array[e.rate_map])
    audit=DivergenceOneStepAudit(e,inherited,p,mass,crossings=0,production=True)
    before=[f.x.array.copy() for f in (s.state,s.old,s.older)]
    history=copy.deepcopy((s.time,s.step_number,audit.diag.initial,audit.diag.previous))
    measure=audit.full.projection.measure
    def bad(weak):
        r=measure(weak);r["weak_continuity_Riesz"]=2e-12;return r
    audit.full.projection.measure=bad
    result=solve_one(e,1e-3,CONTROLS["P2"],audit,"step3a9_badweak_")
    assert result["SNES_reason"]>0 and not result["accepted"]
    for f,a in zip((s.state,s.old,s.older),before):assert np.array_equal(f.x.array,a)
    assert (s.time,s.step_number,audit.diag.initial,audit.diag.previous)==history
    assert not result["candidate_audit"]["full_physical_checks"]["weak_continuity"]
