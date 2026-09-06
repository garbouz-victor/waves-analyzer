import numpy as np

from sloshing.config import SimulationConfig
from sloshing.fem_spaces import FEMSystem
from sloshing.validation.particles import FEMVelocityHistory, advect
from sloshing.validation.qualification.fem_evaluation import FEMGeometry


def analytic_history(spacing, upward=False):
    # Independent analytic field ONLY for interpolation/RK tests, never a
    # production dataset. P2 represents a constant velocity exactly in space.
    f = FEMSystem(SimulationConfig(nx=4, nz=6))
    geom = FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs)
    times = np.linspace(0, 1, round(1/spacing)+1)
    velocity = np.exp(-times[:, None])*np.ones((len(times), f.scalar.N))
    zero = np.zeros_like(velocity)
    return FEMVelocityHistory(geom, times, zero if upward else .1*velocity,
                              velocity if upward else zero,
                              np.zeros((len(times), len(f.surface_x))), f.surface_x)


def test_linear_snapshot_interpolation_converges_second_order():
    errors = []
    for spacing in (.1, .05, .025):
        result = advect(analytic_history(spacing), np.array([[0, -.3]]), np.linspace(0, 1, 11), max_step=.01)
        exact_x = .1*(1-np.exp(-1))
        errors.append(abs(result["paths"][-1, 0, 0]-exact_x))
        assert not np.any(np.isfinite(result["first_invalid_time"]))
    assert 3.9 < errors[0]/errors[1] < 4.1
    assert 3.9 < errors[1]/errors[2] < 4.1


def test_invalid_surface_crossing_is_stopped_not_reflected():
    result = advect(analytic_history(.025, upward=True), np.array([[0, -.1]]), max_step=.005)
    assert .1 < result["first_invalid_time"][0] < .12
    assert np.all(np.isnan(result["paths"][-1]))


def test_particle_rk4_has_fourth_order_for_exact_fem_linear_flow():
    f = FEMSystem(SimulationConfig(nx=4, nz=6))
    geom = FEMGeometry(f.mesh.p, f.mesh.t, f.scalar.element_dofs)
    u = np.tile(f.scalar.doflocs[0], (2, 1))
    history = FEMVelocityHistory(geom, [0., 1.], u, np.zeros_like(u),
                                 np.zeros((2, len(f.surface_x))), f.surface_x)
    errors = []
    for step in (.2, .1, .05):
        result = advect(history, np.array([[.1, -.3]]), max_step=step)
        errors.append(abs(result["paths"][-1, 0, 0]-.1*np.e))
    assert 13 < errors[0]/errors[1] < 17
    assert 13 < errors[1]/errors[2] < 17
