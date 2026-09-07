"""Pure change-of-variables tests and tiny pinned-FEM checks, never 96x48 CI."""
from decimal import Decimal, localcontext
import numpy as np
import pytest
from sloshing.multiphase.phase_rate import reconstruct_phase


def test_scalar_moderate_equivalence_and_jacobian_chain_rule():
    old, dt, stiffness = 1., .125, 3.
    new = old/(1+dt*stiffness)
    rate = -stiffness*old/(1+dt*stiffness)
    assert old+dt*rate == pytest.approx(new, abs=1e-16)
    assert (new-old)/dt+stiffness*new == pytest.approx(0., abs=1e-15)
    assert rate+stiffness*(old+dt*rate) == pytest.approx(0., abs=1e-15)
    assert (1/dt+stiffness)*dt == 1+dt*stiffness


def test_tiny_absolute_subtraction_loses_digits_but_rate_retains_root():
    old, dt = 1., 1e-12
    with localcontext() as context:
        context.prec = 70
        exact_dt = Decimal.from_float(dt)
        exact_rate = -Decimal(1)/(1+exact_dt)
        rate = float(exact_rate)
    new = old+dt*rate
    assert abs((new-old)/dt+new) > 1e-6
    assert abs(rate+new) < 1e-15
    assert abs(Decimal.from_float(rate)-exact_rate) < Decimal("1e-16")


def test_reconstruction_is_exact_coefficient_expression_without_aliases():
    old = np.array([1., -.2, .3]); rate = np.array([-3., 2., 1.])
    saved = old.copy()
    result = reconstruct_phase(old, rate, .125)
    assert np.array_equal(result, old+.125*rate)
    result[:] = 0.
    assert np.array_equal(old, saved)
    with pytest.raises(ValueError):
        reconstruct_phase(old, rate, 0.)


def tiny_solver():
    pytest.importorskip("dolfinx")
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=60., g=0., a_x=0., quadrature_degree=12)
    return CHNSSolver(c)


def test_tiny_FE_residual_and_jacobian_equivalence_and_FD():
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    from sloshing.multiphase.isolated_ch import IsolatedCH
    s = tiny_solver()
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    rate = IsolatedCHPhaseRate(s, lambda x: .2+.1*np.cos(2*x[0])*np.cos(3*x[1]))
    absolute = IsolatedCH(s, s.state.sub(2).collapse())
    dt = .125
    rate.dt.value = absolute.dt.value = dt
    rate.state.sub(0).interpolate(lambda x: .1*np.sin(2*x[0]))
    rate.state.sub(1).interpolate(lambda x: .3*np.cos(x[1]))
    phi = reconstruct_phase(rate.phi_old.x.array, rate.state.sub(0).collapse().x.array, dt)
    _, am0 = absolute.state.function_space.sub(0).collapse()
    _, am1 = absolute.state.function_space.sub(1).collapse()
    absolute.state.x.array[am0] = phi
    absolute.state.x.array[am1] = rate.state.sub(1).collapse().x.array
    rate.state.x.scatter_forward(); absolute.state.x.scatter_forward()
    def vector(form):
        value = fp.assemble_vector(fem.form(form)); out = value.array.copy(); value.destroy(); return out
    old_indices = np.r_[am0, am1]
    new_indices = np.r_[rate.rate_map, rate.mu_map]
    a, b = vector(absolute.F)[old_indices], vector(rate.F)[new_indices]
    assert np.linalg.norm(a-b)/np.linalg.norm(a) < 1e-12
    matrices = []
    for form, state, indices in ((absolute.F, absolute.state, old_indices), (rate.F, rate.state, new_indices)):
        J = ufl.derivative(form, state, ufl.TrialFunction(state.function_space))
        matrix = fp.assemble_matrix(fem.form(J)); matrix.assemble()
        # Tiny test only, no production densification.
        from scipy.sparse import csr_matrix
        ptr, col, val = matrix.getValuesCSR()
        dense = csr_matrix((val, col, ptr), shape=matrix.getSize()).toarray()
        matrices.append(dense[np.ix_(indices, indices)])
        matrix.destroy()
    expected = matrices[0]*np.r_[np.full(len(am0), dt), np.ones(len(am1))][None, :]
    assert np.linalg.norm(matrices[1]-expected)/np.linalg.norm(expected) < 1e-12
    original = rate.state.x.array.copy()
    direction = np.random.default_rng(3505).normal(size=len(original))
    h = 1e-5
    rate.state.x.array[:] = original+h*direction
    plus = vector(rate.F)[new_indices]
    rate.state.x.array[:] = original-h*direction
    minus = vector(rate.F)[new_indices]
    rate.state.x.array[:] = original
    actual = matrices[1] @ direction[new_indices]
    assert np.linalg.norm((plus-minus)/(2*h)-actual)/np.linalg.norm(actual) < 1e-6


def test_prepared_uniform_rate_equilibrium_is_preserved():
    s = tiny_solver()
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    # phi=+1 is a known exact state: wall derivative also vanishes at +1.
    s.initialize(lambda x: np.ones(x.shape[1]))
    s.state.sub(3).interpolate(lambda x: np.zeros(x.shape[1]))
    s.state.x.scatter_forward()
    engine = IsolatedCHPhaseRate(s, initial_guess="zero")
    old = s.state.x.array.copy()
    for dt in (1e-10, 2e-10, 1e-9):
        result = engine.step(dt)
        assert result["reason"] > 0
        assert result["rate_L2"] < 1e-9
    assert max(abs(s.state.x.array-old)) < 1e-9


def test_tiny_nonlinear_rate_convergence_and_diagnostic_trace(tmp_path):
    s = tiny_solver()
    from petsc4py import PETSc
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    from sloshing.multiphase.benchmarks.phase_rate_ch import NewtonTrace
    PETSc.Log.begin()
    engine = IsolatedCHPhaseRate(s, lambda x: 1.+1e-3*np.cos(2*x[0])*np.cos(3*x[1]))
    trace = NewtonTrace(engine, tmp_path)
    result = engine.step(4.657657028534679e-12)
    assert result["reason"] > 0
    assert trace.summary()["accepted"]
    assert (tmp_path/"line_search.txt").exists()
    assert (tmp_path/"snes_iterations.jsonl").exists()


@pytest.mark.parametrize("raises", [True, False])
def test_rate_failed_SNES_does_not_publish_or_advance(raises):
    from types import SimpleNamespace
    from sloshing.multiphase.phase_rate import IsolatedCHPhaseRate
    calls = []
    def fail():
        calls.append(True)
        if raises:
            raise RuntimeError("simulated SNES failure")
    def vector(a):
        return SimpleNamespace(array=np.array(a, dtype=float), scatter_forward=lambda: None)
    old = np.array([.3, .4, .5, .6])
    observer = SimpleNamespace(state=SimpleNamespace(x=vector(old)), time=0., step_number=0)
    engine = IsolatedCHPhaseRate.__new__(IsolatedCHPhaseRate)
    engine.solver, engine.performance = observer, None
    engine.physical_phi_map = engine.rate_map = np.array([0, 1])
    engine.physical_mu_map = engine.mu_map = np.array([2, 3])
    engine.phi_old = SimpleNamespace(x=vector(old[:2]))
    engine.state = SimpleNamespace(x=vector([0., 0., .5, .6]))
    engine.last_accepted_rate = None
    engine.dt = SimpleNamespace(value=1.)
    engine.problem = SimpleNamespace(solve=fail, solver=SimpleNamespace(getConvergedReason=lambda: -6))
    with pytest.raises(RuntimeError, match="SNES fail"):
        engine.step(1e-12)
    assert calls == [True]
    assert np.array_equal(observer.state.x.array, old)
    assert observer.time == 0. and observer.step_number == 0


def test_failed_moderate_gate_forbids_linear_and_downstream_work(tmp_path, monkeypatch):
    from sloshing.multiphase.benchmarks import phase_rate_ch as stages
    monkeypatch.setattr(stages, "RESULTS", tmp_path)
    folder = tmp_path/"moderate_dt"; folder.mkdir()
    stages.write_json(folder/"status.json", {"qualification_status": "failed"})
    with pytest.raises(ValueError, match="prerequisite not passed"):
        stages.linear_tiny()  # must stop before importing PETSc/assembling anything
    assert not (tmp_path/"tiny_step_linear").exists()


def test_unscaled_phase_rate_guess_across_dt_changes_is_pure_kinematics():
    old, rate = np.array([.7, -.1]), np.array([2., -3.])
    first = reconstruct_phase(old, rate, 1e-4)
    second = reconstruct_phase(first, rate, 2e-4)
    assert np.array_equal(second, first+2e-4*rate)
    assert np.allclose(second, old+3e-4*rate, atol=2e-16, rtol=0.)
