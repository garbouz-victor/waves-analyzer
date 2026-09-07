"""Pure dense/sparse comparisons plus explicitly tiny FEniCSx operator tests."""
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from sloshing.multiphase.linearized_ch import LinearizedCH


def synthetic(n=8):
    rng = np.random.default_rng(33)
    a = rng.normal(size=(n, n))
    M = a.T @ a+np.eye(n)
    p = np.eye(n)-np.ones((n, n))/n
    a = rng.normal(size=(n, n))
    K = p @ (a.T @ a+np.eye(n)) @ p
    a = rng.normal(size=(n, n))
    H = a.T @ a+np.eye(n)
    return LinearizedCH(M, K, H, mobility=.7)


def test_mass_projection_symmetry_and_conservation():
    op = synthetic()
    q = op.project(np.arange(op.size, dtype=float))
    assert abs(op.mass @ q) < 1e-12
    assert op.project(q) == pytest.approx(q, abs=1e-14)
    checks = op.matrix_checks()
    assert checks["max_relative_mass_defect"] < 1e-11
    assert checks["max_projection_relative_defect"] < 1e-11
    assert max(checks["frobenius_symmetry_defects"].values()) < 1e-11
    assert max(checks["bilinear_symmetry_defects"].values()) < 1e-11
    assert checks["max_relative_energy_identity_defect"] < 1e-11


def test_dense_tiny_operator_and_coupled_BE():
    op = synthetic()
    M, K, H = (a.toarray() for a in (op.M0, op.K, op.H))
    L = -.7*np.linalg.solve(M, K @ np.linalg.solve(M, H))
    q = op.project(np.arange(op.size, dtype=float))
    assert op.raw(q) == pytest.approx(L @ q, rel=1e-12, abs=1e-12)
    assert op.be_step(q, .01) == pytest.approx(np.linalg.solve(np.eye(op.size)-.01*L, q), abs=1e-12)
    assert np.linalg.eigvals(L).real.max() < 1e-10


def test_changed_operator_fingerprint_and_invalid_inputs():
    a, b = synthetic(), synthetic()
    assert a.fingerprint == b.fingerprint
    b = LinearizedCH(a.M0, a.K, 2*a.H, a.mobility)
    assert a.fingerprint != b.fingerprint
    with pytest.raises(ValueError, match="Positive"):
        a.be_step(np.ones(a.size), 0.)


def tiny_fem_operator():
    pytest.importorskip("dolfinx")
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.equilibrium import solve_equilibrium
    from sloshing.multiphase.linearized_ch import assemble_linearized_ch
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=90., g=0., a_x=0., quadrature_degree=12)
    s = CHNSSolver(c)
    s.initialize(lambda x: np.ones(x.shape[1]))
    eq = solve_equilibrium(s)
    s.initialize_prepared_equilibrium(eq)
    phi = s.state.sub(2).collapse()
    return s, phi, assemble_linearized_ch(s, phi)


def test_tiny_fem_hessian_linearization_and_stable_spectrum():
    from sloshing.multiphase.linearized_ch import NonlinearCHRate
    s, phi, op = tiny_fem_operator()
    checks = op.matrix_checks()
    assert checks["max_relative_mass_defect"] < 1e-11
    assert checks["frobenius_symmetry_defects"]["H"] < 1e-11
    rng = np.random.default_rng(332)
    psi = op.project(rng.normal(size=op.size)); psi /= max(abs(psi))
    rows = NonlinearCHRate(s, phi.function_space, op).finite_difference_checks(phi.x.array, psi)
    assert rows[-1]["central_rate_relative_L2_error"] < 1e-6
    assert rows[-1]["hessian_relative_M_dual_error"] < 1e-6
    assert rows[-1]["central_rate_relative_L2_error"] < rows[0]["central_rate_relative_L2_error"]
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    values = np.linalg.eigvals(dense)
    assert values.real.max() < 1e-8*max(abs(values))
    assert max(abs(values.imag)) < 1e-8*max(abs(values))


def test_pinned_MUMPS_resolvent_matches_dense():
    pytest.importorskip("petsc4py")
    op = synthetic()
    op.resolvent_backend = "petsc_mumps"
    q = op.project(np.arange(op.size, dtype=float))
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    for dt in (1e-10, .001, .1):
        assert op.be_step(q, dt) == pytest.approx(np.linalg.solve(np.eye(op.size)-dt*dense, q), abs=1e-12)
    op.clear_be_cache()


def test_sparse_shift_inverse_full_operator_residual():
    pytest.importorskip("petsc4py")
    from sloshing.multiphase.linearized_ch import shifted_modes
    op = synthetic(16)
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    exact = np.linalg.eigvals(dense)
    rows, modes = shifted_modes(op, -.1, count=4)
    for row in rows:
        assert row["relative_full_operator_residual"] < 1e-8
        assert min(abs(exact-row["real"])) < 1e-8
        assert row["relative_mass"] < 1e-12
