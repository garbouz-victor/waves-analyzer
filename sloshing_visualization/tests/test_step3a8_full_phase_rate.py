import numpy as np
import pytest
from sloshing.multiphase.full_phase_rate import physical_coefficients,transformed_direction,require_scope
from sloshing.multiphase.full_rate_policy import POLICY,CONTROLS


def tiny_solver(theta=90.):
    pytest.importorskip("dolfinx")
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    c=ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4,nz=4,epsilon=1.,theta_equilibrium_deg=theta,g=0.,a_x=0.,quadrature_degree=12)
    s=CHNSSolver(c)
    # Toy only: not the immutable production stiff bump.
    s.initialize(lambda x:.2+1e-4*np.cos(2*x[0])*np.cos(3*x[1]))
    return s


def test_transformation_is_columns_and_physical_copy():
    old=np.arange(8,dtype=float); y=np.ones(8); m=np.array([4,5]); dt=.01
    physical=physical_coefficients(old,y,m,dt)
    assert np.array_equal(physical[m],old[m]+dt*y[m])
    assert np.array_equal(transformed_direction(y,m,dt)[m],dt*y[m])
    physical[:]=0; assert np.all(y==1)
    with pytest.raises(ValueError): physical_coefficients(old,y*np.nan,m,dt)


def test_scope_rejects_unmatched_force_or_BDF2():
    from sloshing.multiphase.config import ModelConfig
    c=ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(quadrature_degree=12,g=0,a_x=0)
    require_scope(c)
    for change in ({"rho_gas":.5},{"g":1.},{"a_x":1.}):
        with pytest.raises(ValueError): require_scope(c.changed(**change))
    with pytest.raises(ValueError): require_scope(c,"bdf2")


def test_arbitrary_full_residual_Jacobian_BC_bridge():
    from sloshing.multiphase.full_rate_bridge import audit
    result=audit(tiny_solver(theta=60.))
    assert result["passed"],result


def one(formulation,dt=1e-3):
    from sloshing.multiphase.full_phase_rate import CHNSPhaseRateBE
    from sloshing.multiphase.full_rate_experiment import AbsoluteBEReference,OneStepAudit,solve_one
    from sloshing.multiphase.phase_rate_history import read_json
    s=tiny_solver(); e=CHNSPhaseRateBE(s) if formulation=="rate" else AbsoluteBEReference(s)
    p=read_json("validation_results/step3a4/policy/policy.json")["policy"]
    mass=float(e.metric.mass@s.state.x.array[e.rate_map])
    audit=OneStepAudit(e,p,POLICY,mass,crossings=0)
    result=solve_one(e,dt,CONTROLS["P2" if dt==1e-3 else "P0"],audit,"step3a8_test_"+formulation+"_")
    assert result["accepted"],result
    return e,audit,result


def test_moderate_full_physical_equivalence():
    from sloshing.multiphase.full_rate_experiment import compare_states
    a,aa,ar=one("absolute"); b,ba,br=one("rate")
    result=compare_states(a.last_snapshot,b.last_snapshot,b.metric,ba.full,
        ar["candidate_audit"]["physical"],br["candidate_audit"]["physical"],POLICY)
    assert result["passed"],result


def test_tiny_coupled_rate_retained_archive_and_equilibrium(tmp_path):
    from sloshing.multiphase.phase_rate_history import archive,load_archive
    from sloshing.multiphase.full_rate_diagnostics import assemble_vector
    e,a,r=one("rate",dt=1e-14)
    assert r["phase_Riesz"]<1e-8
    archive(tmp_path/"accepted",e.last_snapshot,{"accepted":True})
    saved,_=load_archive(tmp_path/"accepted")
    assert np.array_equal(saved["phase_rate"],e.last_accepted_rate)
    assert np.array_equal(saved["phi_new"],saved["phi_old"]+float(saved["dt"])*saved["phase_rate"])
    assert not np.array_equal(saved["phase_rate"],saved["phase_rate_recovered_from_rounded_phi"])
    assert np.linalg.norm(assemble_vector(e.F)[e.rate_map])<1e-10
    from sloshing.multiphase.full_phase_rate import CHNSPhaseRateBE
    s=tiny_solver(); eq=CHNSPhaseRateBE(s,lambda x:np.zeros(x.shape[1]),initial_guess="zero")
    eq.step(1e-12)
    assert max(abs(eq.last_accepted_rate))<1e-12
    assert max(abs(s.state.x.array[eq.rate_map]))<1e-15
