"""Tiny independent tests of an OPTIONAL, non-physical mesh objective.

Synthetic fixtures here are not solver trajectories or resolution evidence.
No test moves a physical interface as an optimization choice: the full current
P2 boundary is fixed during every ordinary optimizer-coordinate evaluation.
"""
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from skfem import Basis, ElementTriP2, MeshTri, MeshTri2

from sloshing.pinned_wetting import free_boundary_mesh_quality as quality


def _orientation(mesh):
    vertices = mesh.p[:, mesh.t]
    first = vertices[:, 1]-vertices[:, 0]
    second = vertices[:, 2]-vertices[:, 0]
    return np.sign(first[0]*second[1]-first[1]*second[0])


def _mesh_pair():
    linear = MeshTri.init_tensor(np.array([-1., -.15, .8, 1.]),
                                np.array([-10., -1., -.1, 0.]))
    linear = linear.with_boundaries({"left": lambda x: x[0] == -1.,
        "right": lambda x: x[0] == 1., "bottom": lambda x: x[1] == -10.,
        "surface": lambda x: x[1] == 0.})
    vertices = linear.p.copy()
    slope = np.tan(np.deg2rad(2.))
    vertices[1] += slope*vertices[0]*(vertices[1]+10.)/10.
    reference = replace(MeshTri2.from_mesh(replace(linear, doflocs=vertices)),
                        _boundaries=linear.boundaries)
    points = reference.p.copy()
    height_fraction = (points[1]+10.)/(10.+slope*points[0])
    points[1] += .015*np.sin(np.pi*(points[0]+1.)/2.)*height_fraction**2
    return reference, replace(reference, doflocs=points)


def _single_triangle():
    linear = MeshTri(np.array([[0., 1., 0.], [0., 0., 1.]]),
                     np.array([[0], [1], [2]]))
    mesh = MeshTri2.from_mesh(linear)
    # Isolate the geometry guard: no physical wall claim is made for this
    # reference triangle, and all node coordinates are test variables.
    return replace(mesh, _boundaries={name: np.array([], dtype=int)
                    for name in ("left", "right", "bottom", "surface")})


@pytest.mark.parametrize("angle", [0., .61, -1.9])
def test_mesh_density_rotation_is_stress_free(angle):
    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
    value, stress, determinant, inverse_transpose = quality.density_and_stress(rotation)
    assert abs(value) < 2e-15
    np.testing.assert_allclose(stress, 0., rtol=0., atol=2e-15)
    assert abs(determinant-1.) < 2e-15
    np.testing.assert_allclose(inverse_transpose, rotation, rtol=0., atol=2e-15)


@pytest.mark.parametrize("matrix", [np.array([[1.2, .15], [-.2, .9]]),
                                    np.array([[.3, .07], [.04, 1.1]])])
def test_mesh_density_stress_matches_independent_directional_derivative(matrix):
    direction = np.array([[.3, -.2], [.4, -.1]])
    _, stress, _, _ = quality.density_and_stress(matrix)
    step = 1e-6
    plus = quality.density_and_stress(matrix+step*direction)[0]
    minus = quality.density_and_stress(matrix-step*direction)[0]
    numerical = (plus-minus)/(2*step)
    predicted = float(np.sum(stress*direction))
    np.testing.assert_allclose(numerical, predicted, rtol=2e-8, atol=2e-9)


def test_mesh_density_penalizes_isotropic_collapse_and_rejects_wrong_orientation():
    values = [quality.density_and_stress(scale*np.eye(2))[0]
              for scale in (1., .1, .01, .001)]
    assert np.all(np.diff(values) > 0.) and values[-1] > 100.
    for matrix in (np.diag([-1., 1.]), np.diag([0., 1.]), np.full((2, 2), np.nan)):
        with pytest.raises(ValueError, match="positive relative Jacobian"):
            quality.density_and_stress(matrix)


def test_projected_scaled_nodal_gradient_matches_finite_differences():
    reference, current = _mesh_pair()
    objective = quality.MeshObjective(reference, current, _orientation(reference))
    rng = np.random.default_rng(1287)
    variables = rng.normal(scale=1e-3, size=int(objective.allowed.sum()))
    frozen_scale = objective.scale.copy()
    value, gradient = objective.value_gradient(variables)
    assert np.isfinite(value) and np.isfinite(gradient).all()
    for _ in range(4):
        direction = rng.normal(size=len(variables))
        direction /= np.linalg.norm(direction)
        step = 1e-6
        plus = objective.value_gradient(variables+step*direction)[0]
        minus = objective.value_gradient(variables-step*direction)[0]
        numerical = (plus-minus)/(2*step)
        predicted = float(gradient@direction)
        np.testing.assert_allclose(numerical, predicted, rtol=2e-5, atol=2e-7)
    np.testing.assert_array_equal(objective.scale, frozen_scale)


def test_optimizer_coordinates_keep_full_surface_bottom_and_wall_normals_bitwise_fixed():
    reference, current = _mesh_pair()
    objective = quality.MeshObjective(reference, current, _orientation(reference))
    variables = np.linspace(-.002, .003, int(objective.allowed.sum()))
    points = objective.coordinates(variables)
    basis = Basis(reference, ElementTriP2())
    surface = basis.get_dofs("surface").all()
    bottom = basis.get_dofs("bottom").all()
    left, right = (basis.get_dofs(name).all() for name in ("left", "right"))
    np.testing.assert_array_equal(points[:, surface], current.p[:, surface])
    np.testing.assert_array_equal(points[:, bottom], current.p[:, bottom])
    np.testing.assert_array_equal(points[0, np.union1d(left, right)], current.p[0, np.union1d(left, right)])
    np.testing.assert_array_equal(points[~objective.allowed], current.p[~objective.allowed])
    right_interior = np.setdiff1d(right, np.union1d(surface, bottom))
    assert len(right_interior) > 0
    assert np.any(points[1, right_interior] != current.p[1, right_interior])
    # Tangential MESH redistribution is allowed; it creates no physical slip
    # velocity and does not change either endpoint or the continuous curve.


def test_exact_vertex_negative_rejected_despite_positive_interior_gauss_jacobians():
    reference = _single_triangle()
    objective = quality.MeshObjective(reference, reference, _orientation(reference))
    r, s = reference.p
    epsilon = -1e-3
    points = np.array([r, epsilon*s+r*s+.5*s*s])
    quadrature = Basis(reference, ElementTriP2(), intorder=4).X
    # Direct analytic determinant of X=(r, eps*s+r*s+s²/2), not production J.
    assert np.all(epsilon+quadrature.sum(axis=0) > 0.)
    assert epsilon < 0.  # det J at the closed-triangle vertex (0,0).
    variables = ((points-reference.p)/objective.scale)[objective.allowed]
    value, gradient = objective.value_gradient(variables)
    assert np.isinf(value) and np.isfinite(gradient).all()
    assert objective.invalid_evaluations == 1
    np.testing.assert_array_equal(objective.best, reference.p)


def test_invalid_initial_current_geometry_cannot_be_returned_as_best():
    reference = _single_triangle()
    r, s = reference.p
    invalid = replace(reference, doflocs=np.array([r, -.001*s+r*s+.5*s*s]))
    with pytest.raises(ValueError):
        quality.MeshObjective(reference, invalid, _orientation(reference))
    with pytest.raises(ValueError):
        quality.optimize_mesh(reference, invalid, _orientation(reference), max_iterations=1)


def test_nonaffine_reference_is_not_silently_used_as_affine():
    reference = _single_triangle()
    r, s = reference.p
    curved = replace(reference, doflocs=np.array([r, s+.05*r*s]))
    with pytest.raises(ValueError, match="affine"):
        quality.MeshObjective(curved, curved, _orientation(reference))


def test_optimizer_nonconvergence_does_not_invent_improvement_or_success(monkeypatch):
    reference, current = _mesh_pair()
    def stopped(fun, initial, max_iterations):
        value, _ = fun(initial)
        return SimpleNamespace(x=initial.copy(), fun=value, success=False,
                               nit=0, message="intentional bounded test stop")
    monkeypatch.setattr(quality, "feasible_lbfgs", stopped)
    selected, info = quality.optimize_mesh(reference, current, _orientation(reference), max_iterations=1)
    np.testing.assert_array_equal(selected.p, current.p)
    assert info["optimizer_success"] is False
    assert info["initial_min_relative_jacobian"] == info["selected_min_relative_jacobian"]
    assert info["physical_constitutive_law"] is False
    assert info["full_free_boundary_displacement_m"] == 0.


def _edge_minimum_pair(matrix, epsilon=.04, location=.25):
    """Analytic P2 map with det DQ=(r-location)²+eps+s*(1+2r).

    Its unique minimum on the closed triangle lies inside s=0, not at a
    vertex or quadrature point.  A nonorthogonal affine reference checks the
    physical-to-reference gradient transformation independently.
    """
    standard = _single_triangle()
    r, s = standard.p
    reference = replace(standard, doflocs=matrix@standard.p)
    q = np.array([r+.5*r*r+s,
                  s+r*s+(1-location*location-epsilon)*r+(1+location)*r*r+s*s])
    return reference, replace(reference, doflocs=matrix@q), r


@pytest.mark.parametrize("matrix", [np.eye(2), np.array([[1.3, .4], [-.2, .9]]),
                                    np.array([[-1.3, .4], [.2, .9]])])
def test_exact_nonvertex_minimum_witness_matches_analytic_quadratic_map(matrix):
    epsilon, location = .04, .25
    reference, current, _ = _edge_minimum_pair(matrix, epsilon, location)
    minima, witnesses = quality.minimum_witnesses(current, _orientation(reference))
    np.testing.assert_allclose(minima, abs(np.linalg.det(matrix))*epsilon,
                               rtol=3e-13, atol=3e-15)
    np.testing.assert_allclose(witnesses[:, 0], [location, 0.], rtol=0., atol=2e-14)
    # All three vertex values exceed the edge minimum: a vertex-only barrier
    # cannot represent this event (the actual saved failure is also nonvertex).
    assert epsilon < min(location**2+epsilon,
                         (1-location)**2+epsilon, location**2+epsilon+1.)


@pytest.mark.parametrize("matrix", [np.eye(2), np.array([[1.3, .4], [-.2, .9]])])
def test_nonvertex_minimum_barrier_value_and_envelope_gradient_are_independent(matrix):
    epsilon = .04
    reference, current, r = _edge_minimum_pair(matrix, epsilon)
    objective = quality.MeshObjective(reference, current, _orientation(reference))
    initial = np.zeros(int(objective.allowed.sum()))
    total, gradient = objective.value_gradient(initial)
    objective.minimum_weight = 0.
    without, without_gradient = objective.value_gradient(initial)
    objective.minimum_weight = .25
    np.testing.assert_allclose(total-without, .25*(epsilon-1-np.log(epsilon)),
                               rtol=1e-13, atol=2e-14)
    # Vary epsilon analytically: Q_y changes by -r*d_eps, so the exact
    # minimum increases by d_eps while its reference witness stays at .25.
    displacement = matrix@np.array([np.zeros_like(r), -r])
    direction = (displacement/objective.scale)[objective.allowed]
    np.testing.assert_allclose((gradient-without_gradient)@direction,
                               .25*(1-1/epsilon), rtol=1e-12, atol=1e-12)
    # General nodal perturbations move the active witness.  Envelope
    # differentiation must still agree without differentiating that minimizer.
    rng = np.random.default_rng(48021)
    for _ in range(3):
        direction = rng.normal(size=len(initial))
        direction /= np.linalg.norm(direction)
        step = 1e-7
        def barrier(variables):
            objective.minimum_weight = .25
            full = objective.value_gradient(variables)[0]
            objective.minimum_weight = 0.
            base = objective.value_gradient(variables)[0]
            objective.minimum_weight = .25
            return full-base
        numerical = (barrier(initial+step*direction)-barrier(initial-step*direction))/(2*step)
        predicted = (gradient-without_gradient)@direction
        np.testing.assert_allclose(numerical, predicted, rtol=2e-6, atol=2e-7)


def test_feasible_armijo_recovers_after_infeasible_first_trial():
    calls = []
    def objective(point):
        x = float(point[0]); calls.append(x)
        if x <= 0:
            return np.inf, np.zeros_like(point)
        return (x-.1)**2, np.array([2*(x-.1)])
    result = quality.feasible_lbfgs(objective, np.array([.5]), max_iterations=8)
    assert calls[1] < 0.  # first descent proposal is deliberately infeasible
    assert len(calls) >= 3 and calls[2] > 0.
    assert result.success
    assert result.fun < 1e-20
    np.testing.assert_allclose(result.x, [.1], rtol=0., atol=1e-10)


def test_feasible_armijo_stall_cannot_claim_unchanged_objective_convergence():
    calls = []
    def stalled(point):
        calls.append(point.copy())
        return 1., np.ones_like(point)
    initial = np.array([0., 0.])
    result = quality.feasible_lbfgs(stalled, initial, max_iterations=3)
    assert not result.success
    assert "line search exhausted" in result.message
    np.testing.assert_array_equal(result.x, initial)
    assert result.fun == 1.
    assert len(calls) == 61  # bounded initial evaluation + 60 trial evaluations


@pytest.mark.parametrize("value,gradient", [(np.inf, np.array([0.])),
                                           (1., np.array([np.nan]))])
def test_feasible_armijo_nonfinite_initial_evaluation_cannot_succeed(value, gradient):
    result = quality.feasible_lbfgs(lambda point: (value, gradient.copy()),
                                    np.array([0.]), max_iterations=2)
    assert not result.success
