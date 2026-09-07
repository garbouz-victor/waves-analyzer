import numpy as np
import pytest
from scipy.linalg import expm
from sloshing.multiphase.linearized_ch import LinearizedCH
from sloshing.multiphase.ch_krylov import MassArnoldi, compare_references


def known_system(n=10):
    rng = np.random.default_rng(3303)
    M = np.diag(np.linspace(1., 2., n))
    mass = M @ np.ones(n)
    # Euclidean nullspace of sqrt(M)*1 gives an M-orthonormal mass-zero basis.
    from scipy.linalg import null_space
    V = np.diag(1/np.sqrt(np.diag(M))) @ null_space(np.sqrt(np.diag(M))[None, :])
    K = M @ V @ np.diag(np.geomspace(1., 100., n-1)) @ V.T @ M
    op = LinearizedCH(M, K, M)
    q = op.project(rng.normal(size=n))
    return op, q, V


def test_M_orthogonality_mass_and_known_Ritz_values():
    op, q, _ = known_system()
    arnoldi = MassArnoldi(op, q, 9)
    arnoldi.extend(9)
    checks = arnoldi.diagnostics()
    assert checks["orthogonality_defect"] < 1e-12
    assert checks["max_checked_Arnoldi_relative_residual"] < 1e-12
    assert checks["maximum_basis_mass"] < 1e-12
    values = sorted(r["real"] for r in arnoldi.ritz())
    assert values == pytest.approx(sorted(-np.geomspace(1., 100., 9)), rel=1e-10)


def test_exponential_reference_and_increasing_dimension():
    op, q, _ = known_system()
    arnoldi = MassArnoldi(op, q, 9)
    times = np.array([0., .0001, .001, .01, .1])
    refs = []
    for m in (3, 6, 9):
        arnoldi.extend(m)
        refs.append(arnoldi.reference(times))
    a, b = compare_references(op, refs[0], refs[2]), compare_references(op, refs[1], refs[2])
    assert b["max_relative_L2_difference"] < a["max_relative_L2_difference"]
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    exact = np.array([expm(t*dense) @ q for t in times])
    assert refs[-1]["states"] == pytest.approx(exact, abs=1e-12)
    assert max(refs[-1]["projected_energy_identity_relative_defect"]) < 1e-11
    with pytest.raises(ValueError, match="identical"):
        compare_references(op, refs[0], {**refs[1], "times": times+1})


def test_zero_mass_seed_required_and_fingerprint_deterministic():
    op, q, _ = known_system()
    with pytest.raises(ValueError, match="Nonzero"):
        MassArnoldi(op, np.zeros(op.size))
    a, b = MassArnoldi(op, q, 6), MassArnoldi(op, q, 6)
    a.extend(6); b.extend(6)
    assert a.fingerprint() == b.fingerprint()


def test_restarted_SLEPc_M_exponential_matches_dense():
    pytest.importorskip("slepc4py")
    from sloshing.multiphase.ch_krylov import RestartedExponential
    op, q, _ = known_system(12)
    reference = RestartedExponential(op, dimension=6)
    try:
        result, log = reference.action(q, .2)
        dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
        exact = expm(.2*dense) @ q
        assert op.norm(result-exact)/op.norm(exact) < 1e-9
        assert abs(op.mass @ result) < 1e-12
        assert log["reason"] > 0
    finally:
        reference.close()


def test_fixed_pole_rational_M_space_against_dense_exponential():
    from sloshing.multiphase.ch_krylov import RationalMassKrylov
    op, q, _ = known_system(12)
    rational = RationalMassKrylov(op, q, (.001, .01, .1), max_dimension=11)
    rational.extend(11)
    assert rational.diagnostics()["orthogonality_defect"] < 1e-12
    times = np.array([0., .01, .1, 1.])
    reference = rational.reference(times)
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    exact = np.array([expm(t*dense) @ q for t in times])
    assert reference["states"] == pytest.approx(exact, abs=1e-11)
    assert max(r["relative_Ritz_residual"] for r in rational.ritz()) < 1e-10
    op.clear_be_cache()


def test_Ritz_gate_uses_actual_operator_not_only_terminal_recurrence():
    op, q, _ = known_system()
    arnoldi = MassArnoldi(op, q, 9); arnoldi.extend(9)
    arnoldi.hessenberg[9, 8] = 0.
    arnoldi.hessenberg[:9, :9] += 1e-4*np.eye(9)
    ritz = arnoldi.ritz()
    assert max(r["Arnoldi_recurrence_residual_estimate"] for r in ritz) == 0.
    assert max(r["relative_Ritz_residual"] for r in ritz) > 1e-6


def test_rational_happy_breakdown_cannot_append_zero_basis_vectors():
    from sloshing.multiphase.ch_krylov import RationalMassKrylov
    op, _, V = known_system()
    rational = RationalMassKrylov(op, V[:, 0], (.01, .1), max_dimension=8)
    rational.extend(4)
    assert rational.breakdown
    assert rational.dimension == 1
    rational.extend(8)
    assert rational.dimension == 1
    assert rational.reference([0., .1])["states"][1] == pytest.approx(np.exp(-.1)*V[:, 0], abs=1e-12)
    op.clear_be_cache()
