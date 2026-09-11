"""Independent algebra checks for the PW2 implementation acceleration.

These explicitly synthetic small-mesh fixtures are not physical trajectories.
The reference operators below are the original vector variational forms, not
production's scalar-block assembler or its operator-action implementation.
"""
from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import bmat, csc_matrix
from scipy.sparse.linalg import spsolve
from skfem import (Basis, BilinearForm, ElementTriP1, ElementTriP2,
                   ElementVector, LinearForm, MeshTri, MeshTri2, asm)
from skfem.helpers import ddot, div, dot, sym_grad

from sloshing.pinned_wetting import free_boundary as production
from sloshing.pinned_wetting import free_boundary_remesh as transfer


@BilinearForm
def _reference_mass(u, v, _):
    return dot(u, v)


@BilinearForm
def _reference_strain(u, v, _):
    return 2*.01*ddot(sym_grad(u), sym_grad(v))


@BilinearForm
def _reference_divergence(u, q, _):
    return q*div(u)


@BilinearForm
def _reference_left_robin(u, v, _):
    return .01/.50*u[1]*v[1]


@BilinearForm
def _reference_geometric_mass(u, v, w):
    return .5*w.dilation*dot(u, v)


@LinearForm
def _reference_body_gravity(v, _):
    return -9.81*v[1]


@LinearForm
def _reference_hydrostatic_surface(v, w):
    return -9.81*w.x[1]*dot(v, w.n)


@LinearForm
def _reference_momentum_x(v, _):
    return v[0]


@LinearForm
def _reference_momentum_z(v, _):
    return v[1]


def _synthetic_mesh(curved):
    linear = MeshTri.init_tensor(np.linspace(-1., 1., 5), np.linspace(-10., 0., 5))
    linear = linear.with_boundaries({"left": lambda x: x[0] == -1.,
        "right": lambda x: x[0] == 1., "bottom": lambda x: x[1] == -10.,
        "surface": lambda x: x[1] == 0.})
    points = linear.p.copy()
    points[1] += np.tan(np.deg2rad(2.))*points[0]*(points[1]+10.)/10.
    mesh = replace(MeshTri2.from_mesh(replace(linear, doflocs=points)),
                   _boundaries=linear.boundaries)
    if curved:
        points = mesh.p.copy()
        height_fraction = (points[1]+10.)/(10.+np.tan(np.deg2rad(2.))*points[0])
        points[1] += .03*np.sin(np.pi*(points[0]+1.)/2.)*height_fraction**2
        mesh = replace(mesh, doflocs=points)
    return mesh


def _reference_operators(mesh, controls, advecting_velocity):
    velocity = Basis(mesh, ElementVector(ElementTriP2()), intorder=controls.intorder)
    pressure = Basis(mesh, ElementTriP1(), intorder=controls.intorder)
    field = velocity.interpolate(advecting_velocity)
    mass = asm(_reference_mass, velocity).tocsr()
    geometric = (asm(_reference_geometric_mass, velocity, dilation=div(field)).tocsr()
                 if controls.skew_divergence else mass*0.)
    body = asm(_reference_body_gravity, velocity)
    # FacetBasis is safe on this deliberately mildly curved tiny fixture. It
    # provides a separate route from production's edge-polynomial integration.
    surface = asm(_reference_hydrostatic_surface, velocity.boundary("surface"))
    return {"M": mass, "Kb": asm(_reference_strain, velocity).tocsr(),
            "Kl": asm(_reference_left_robin, velocity.boundary("left")).tocsr(),
            "B": asm(_reference_divergence, velocity, pressure).tocsr(),
            "Kg": geometric, "f_body": body,
            "f": surface if controls.hydrostatic_split else body}


def _reference_mixed_solution(op, matrices, v0, dt):
    ids = op.free
    mass = matrices["M"][ids][:, ids]
    viscous = (matrices["Kb"]+matrices["Kl"]+matrices["Kg"])[ids][:, ids]
    divergence = matrices["B"][:, ids]
    operator = bmat([[mass/dt+.5*viscous, -divergence.T],
                     [.5*divergence, None]], format="csc")
    rhs = np.r_[(mass/dt-.5*viscous)@v0[ids]+matrices["f"][ids],
                -.5*divergence@v0[ids]]
    solution = spsolve(operator, rhs)
    velocity = np.zeros_like(v0)
    velocity[ids] = solution[:len(ids)]
    return velocity, solution[len(ids):], operator, rhs


@pytest.mark.parametrize("curved", [False, True])
@pytest.mark.parametrize("hydrostatic", [False, True])
@pytest.mark.parametrize("skew", [False, True])
def test_accelerated_full_operators_equal_original_vector_forms(curved, hydrostatic, skew):
    mesh = _synthetic_mesh(curved)
    controls = production.Controls(nx=4, nz=4, intorder=8,
        hydrostatic_split=hydrostatic, skew_divergence=skew)
    n = Basis(mesh, ElementVector(ElementTriP2()), intorder=8).N
    advecting = np.random.default_rng(8132).normal(scale=.02, size=n)
    op = production.Operators(mesh, controls, advecting)
    reference = _reference_operators(mesh, controls, advecting)
    # Every full matrix entry counts, including rows hidden by the boundary
    # constraints. An accidental RIGHT Robin term must not become invisible.
    for name in ("M", "Kb", "Kl", "B", "Kg"):
        np.testing.assert_allclose(getattr(op, name).toarray(), reference[name].toarray(),
                                   rtol=4e-12, atol=3e-13, err_msg=name)
    for name in ("f_body", "f"):
        np.testing.assert_allclose(getattr(op, name), reference[name],
                                   rtol=4e-12, atol=3e-12, err_msg=name)
    np.testing.assert_array_equal(op.Kl[op.wall["right"]].toarray(), 0.)
    np.testing.assert_allclose(op.Kb.toarray(), op.Kb.T.toarray(), rtol=0., atol=3e-13)
    assert float(advecting@(op.Kb@advecting)) >= 0.
    assert float(advecting@(op.Kl@advecting)) >= 0.


@pytest.mark.parametrize("curved", [False, True])
def test_accelerated_mixed_solve_equals_independent_vector_direct(curved):
    mesh = _synthetic_mesh(curved)
    controls = production.Controls(nx=4, nz=4, dt=.005,
        hydrostatic_split=True, skew_divergence=True)
    n = Basis(mesh, ElementVector(ElementTriP2()), intorder=8).N
    advecting = np.random.default_rng(317).normal(scale=.01, size=n)
    op = production.Operators(mesh, controls, advecting)
    v0 = np.random.default_rng(123).normal(scale=.02, size=n)
    v0[op.fixed] = 0.
    reference = _reference_operators(mesh, controls, advecting)
    expected_v, expected_p, operator, rhs = _reference_mixed_solution(op, reference, v0, controls.dt)
    actual_v, actual_p = op.solve(v0, controls.dt)
    np.testing.assert_allclose(actual_v, expected_v, rtol=2e-10, atol=2e-11)
    np.testing.assert_allclose(actual_p, expected_p, rtol=2e-10, atol=2e-10)
    residual = operator@np.r_[actual_v[op.free], actual_p]-rhs
    assert np.max(abs(residual[:len(op.free)])) < 1e-9
    assert np.max(abs(residual[len(op.free):])) < 1e-11
    assert np.max(abs(actual_v[op.fixed])) == 0.


def _synthetic_saddle(variation):
    """Nonsymmetric saddle convention exactly matching the material step."""
    hessian = np.diag([2., 3., 5.])+variation*np.diag([.03, .02, .04])
    divergence = np.array([[1., 0., .5], [0., 1., -.25]])
    divergence += variation*np.array([[0., .002, 0.], [.003, 0., 0.]])
    return csc_matrix(np.block([[hessian, -divergence.T],
                                [.5*divergence, np.zeros((2, 2))]]))


def test_lagged_lu_solves_current_saddle_not_stale_operator():
    solver = production.PicardLinearSolver()
    first, current = _synthetic_saddle(0.), _synthetic_saddle(1.)
    first_rhs = np.array([.01, -.02, .03, -.01, .02])
    current_rhs = first_rhs+np.array([.001, .002, -.001, .003, -.002])
    solver.solve(first, first_rhs, 3)
    initial_factor = solver.factor
    expected = spsolve(current, current_rhs)
    assert np.max(abs(first@expected-current_rhs)) > 1e-5
    actual = solver.solve(current, current_rhs, 3)
    np.testing.assert_allclose(actual, expected, rtol=0., atol=1e-13)
    assert np.max(abs(current@actual-current_rhs)) < 1e-13
    assert solver.factor is initial_factor
    assert solver.statistics["factorizations"] == 1
    assert solver.statistics["direct_fallbacks"] == 0
    assert 0 < solver.statistics["krylov_iterations"] <= 120
    assert solver.statistics["maximum_true_block_residual"] <= 1e-13


@pytest.mark.parametrize("failure", ["divergence_only", "momentum_only", "stale_operator",
                                      "nonfinite", "iteration_limit", "exception"])
def test_gmres_false_success_or_failure_uses_current_direct_fallback(monkeypatch, failure):
    solver = production.PicardLinearSolver()
    first, current = _synthetic_saddle(0.), _synthetic_saddle(1.)
    rhs = np.array([.01, -.02, .03, -.01, .02])
    solver.solve(first, rhs, 3)
    expected = spsolve(current, rhs)
    candidate = expected.copy()
    if failure == "divergence_only":
        # An error invisible to momentum rows: H*dv-B.T*dp=0, B*dv != 0.
        dp = np.array([1e-5, -2e-5])
        matrix = current.toarray()
        dv = np.linalg.solve(matrix[:3, :3], -matrix[:3, 3:]@dp)
        candidate += np.r_[dv, dp]
        residue = current@candidate-rhs
        assert np.max(abs(residue[:3])) < 1e-13
        assert np.max(abs(residue[3:])) > 1e-7
    elif failure == "momentum_only":
        candidate[-1] += 1e-5
        assert np.max(abs((current@candidate-rhs)[3:])) < 1e-13
    elif failure == "stale_operator":
        candidate = spsolve(first, rhs)
    elif failure == "nonfinite":
        candidate[:] = np.nan

    calls = []
    def invalid_gmres(matrix, supplied_rhs, **kwargs):
        calls.append(1)
        np.testing.assert_array_equal(matrix.toarray(), current.toarray())
        np.testing.assert_array_equal(supplied_rhs, rhs)
        assert kwargs["callback_type"] == "pr_norm"
        assert 0 < kwargs["restart"]*kwargs["maxiter"] <= 120
        kwargs["callback"](0.)  # Report false preconditioned convergence.
        if failure == "exception":
            raise RuntimeError("intentional tiny Krylov failure")
        return candidate, 1 if failure == "iteration_limit" else 0

    monkeypatch.setattr(production, "gmres", invalid_gmres)
    actual = solver.solve(current, rhs, 3)
    assert len(calls) == 1
    np.testing.assert_allclose(actual, expected, rtol=0., atol=1e-14)
    assert np.max(abs(current@actual-rhs)) < 1e-13
    assert solver.statistics["direct_fallbacks"] == 1
    assert solver.statistics["factorizations"] == 2


@pytest.mark.parametrize("failure", ["nonfinite", "bad_residual"])
def test_direct_fallback_cannot_publish_a_bad_linear_solution(monkeypatch, failure):
    class InvalidFactor:
        def solve(self, rhs):
            return np.full_like(rhs, np.nan) if failure == "nonfinite" else np.zeros_like(rhs)
    monkeypatch.setattr(production, "splu", lambda *args, **kwargs: InvalidFactor())
    solver = production.PicardLinearSolver()
    with pytest.raises(ValueError, match="current mixed linear system residual"):
        solver.solve(_synthetic_saddle(0.), np.ones(5)*.01, 3)


def test_curved_operator_picard_helper_matches_independent_vector_direct():
    controls = production.Controls(nx=4, nz=4, dt=.005,
        hydrostatic_split=True, skew_divergence=True)
    solver = production.PicardLinearSolver()
    rng = np.random.default_rng(346)
    first_mesh = _synthetic_mesh(True)
    n = Basis(first_mesh, ElementVector(ElementTriP2()), intorder=8).N
    advecting = rng.normal(scale=.01, size=n)
    first = production.Operators(first_mesh, controls, advecting)
    v0 = rng.normal(scale=.01, size=n)
    v0[first.fixed] = 0.
    first.solve(v0, controls.dt, linear_solver=solver)
    points = first_mesh.p.copy()
    points[1] += 1e-5*np.sin(np.pi*(points[0]+1.)/2.)*((points[1]+10.)/10.)**2
    current_mesh = replace(first_mesh, doflocs=points)
    current = production.Operators(current_mesh, controls, advecting+.001)
    reference = _reference_operators(current_mesh, controls, advecting+.001)
    expected_v, expected_p, operator, rhs = _reference_mixed_solution(current, reference, v0, controls.dt)
    actual_v, actual_p = current.solve(v0, controls.dt, linear_solver=solver)
    np.testing.assert_allclose(actual_v, expected_v, rtol=2e-10, atol=2e-11)
    np.testing.assert_allclose(actual_p, expected_p, rtol=2e-10, atol=2e-10)
    residual = operator@np.r_[actual_v[current.free], actual_p]-rhs
    assert np.max(abs(residual[:len(current.free)])) < 1e-9
    assert np.max(abs(residual[len(current.free):])) < 1e-11
    assert np.max(abs(actual_v[current.fixed])) == 0.
    assert solver.statistics["krylov_iterations"] > 0


def _projection_fixture(epsilon=None):
    mass = np.diag([2., 3., 4., 5.])
    divergence = np.array([[1., 0., 0., 0.]])
    momenta = (np.array([[0., 1., 0., 0.], [0., 0., 1., 0.]]) if epsilon is None else
               np.array([[0., 1., 0., 0.], [0., 1., epsilon, 0.]]))
    known_v = np.array([0., .02, -.03, .04])
    known_p = np.array([.015])
    known_lambda = np.array([.035, -.025])
    rhs = mass@known_v+divergence.T@known_p+momenta.T@known_lambda
    target = momenta@known_v
    constraints = np.vstack((divergence, momenta))
    full = np.block([[mass, constraints.T],
                     [constraints, np.zeros((len(constraints), len(constraints)))]])
    full_rhs = np.r_[rhs, np.zeros(len(divergence)), target]
    return mass, divergence, momenta, rhs, target, full, full_rhs, known_v


@pytest.mark.parametrize("epsilon", [None, 1e-3, 1e-5])
def test_transfer_schur_equals_full_dense_system_including_near_dependence(epsilon):
    mass, divergence, momenta, rhs, target, full, full_rhs, known_v = _projection_fixture(epsilon)
    solution, statistics = transfer.solve_transfer(csc_matrix(mass), csc_matrix(divergence),
                                                   momenta, rhs, target)
    reference = np.linalg.solve(full, full_rhs)
    np.testing.assert_allclose(solution[:len(mass)], known_v, rtol=0., atol=3e-10)
    np.testing.assert_allclose(solution[:len(mass)], reference[:len(mass)], rtol=0., atol=3e-10)
    assert np.max(abs(full@solution-full_rhs)) <= 1e-12
    assert statistics["full_mixed_residual"] <= 1e-12
    if epsilon is None:
        np.testing.assert_allclose(solution, reference, rtol=0., atol=1e-13)
        assert not statistics["full_direct_fallback"]
    else:
        assert statistics["schur_condition_number"] > 1e5


def test_curved_transfer_schur_matches_independent_full_vector_projection():
    mesh = _synthetic_mesh(True)
    controls = production.Controls(nx=4, nz=4, intorder=8)
    basis = Basis(mesh, ElementVector(ElementTriP2()), intorder=8)
    op = production.Operators(mesh, controls)
    reference = _reference_operators(mesh, controls, np.zeros(basis.N))
    ids = op.free
    mass = reference["M"][ids][:, ids]
    divergence = reference["B"][:, ids]
    momenta = np.vstack((asm(_reference_momentum_x, basis),
                         asm(_reference_momentum_z, basis)))[:, ids]
    rhs = np.random.default_rng(413).normal(scale=.001, size=len(ids))
    target = np.zeros(2)
    constraints = np.vstack((divergence.toarray(), momenta))
    full = np.block([[mass.toarray(), constraints.T],
                     [constraints, np.zeros((len(constraints), len(constraints)))]])
    full_rhs = np.r_[rhs, np.zeros(len(constraints))]
    expected = np.linalg.solve(full, full_rhs)
    actual, statistics = transfer.solve_transfer(mass, divergence, momenta, rhs, target)
    np.testing.assert_allclose(actual[:len(ids)], expected[:len(ids)], rtol=1e-9, atol=2e-10)
    assert np.max(abs(full@actual-full_rhs)) < 1e-12
    assert np.max(abs(divergence@actual[:len(ids)])) < 1e-12
    assert np.max(abs(momenta@actual[:len(ids)]-target)) < 1e-12
    assert statistics["full_mixed_residual"] < 1e-12


def test_transfer_exactly_dependent_consistent_constraints_keep_all_original_rows():
    mass, divergence, momenta, rhs, target, full, full_rhs, known_v = _projection_fixture(0.)
    assert np.linalg.matrix_rank(full) < len(full)
    # Multipliers need not be unique, but constrained minimum-energy velocity
    # is unique. A rank-revealing small solve must satisfy BOTH original rows.
    actual, statistics = transfer.solve_transfer(csc_matrix(mass), csc_matrix(divergence),
                                                 momenta, rhs, target)
    np.testing.assert_allclose(actual[:len(mass)], known_v, rtol=0., atol=1e-13)
    assert np.max(abs(full@actual-full_rhs)) < 1e-12
    assert np.max(abs(momenta@actual[:len(mass)]-target)) < 1e-12
    assert statistics["schur_numerical_rank"] == 1
    assert not statistics["full_direct_fallback"]


def test_transfer_dependent_inconsistent_momenta_are_never_silently_dropped():
    mass, divergence, momenta, rhs, target, _, _, _ = _projection_fixture(0.)
    target[1] += 1e-3  # Identical rows cannot have different required integrals.
    with pytest.raises(ValueError, match="conservative transfer mixed solve failed"):
        transfer.solve_transfer(csc_matrix(mass), csc_matrix(divergence), momenta, rhs, target)


@pytest.mark.parametrize("corruption", ["wrong_schur_solution", "nonfinite_local_solve"])
def test_transfer_false_schur_success_is_rechecked_against_full_system(monkeypatch, corruption):
    mass, divergence, momenta, rhs, target, full, full_rhs, _ = _projection_fixture()
    expected = np.linalg.solve(full, full_rhs)
    if corruption == "wrong_schur_solution":
        real_solve = np.linalg.solve
        def bad_small_solve(matrix, vector):
            result = real_solve(matrix, vector)
            if np.shape(matrix) == (2, 2):
                result = result+.05  # A claimed 2x2 answer, not a valid full solve.
            return result
        monkeypatch.setattr(transfer.np.linalg, "solve", bad_small_solve)
    else:
        class NonfiniteFactor:
            def solve(self, loads):
                return np.full_like(loads, np.nan)
        monkeypatch.setattr(transfer, "splu", lambda *args, **kwargs: NonfiniteFactor())
    actual, statistics = transfer.solve_transfer(csc_matrix(mass), csc_matrix(divergence),
                                                 momenta, rhs, target)
    assert statistics["full_direct_fallback"]
    np.testing.assert_allclose(actual, expected, rtol=0., atol=1e-13)
    assert np.max(abs(full@actual-full_rhs)) < 1e-12


def test_transfer_bad_full_direct_fallback_is_rejected(monkeypatch):
    mass, divergence, momenta, rhs, target, full, _, _ = _projection_fixture()
    def unavailable_local(*args, **kwargs):
        raise RuntimeError("intentional local factorization failure")
    monkeypatch.setattr(transfer, "splu", unavailable_local)
    monkeypatch.setattr(transfer, "spsolve", lambda *args, **kwargs: np.zeros(len(full)))
    with pytest.raises(ValueError, match="conservative transfer mixed solve failed"):
        transfer.solve_transfer(csc_matrix(mass), csc_matrix(divergence), momenta, rhs, target)
