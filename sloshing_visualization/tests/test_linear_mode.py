from dataclasses import replace

import numpy as np

from sloshing import SimulationConfig, SloshingSolver


def test_release_direction_is_computed(short_run):
    solver, snapshots = short_run
    state = snapshots[1][0]
    points = np.array([[0, -.75, .75], [-.1, -.2, -.2]])
    velocity = solver.fem.velocity_at(solver.fem.expand_velocity(state.velocity), points)
    assert velocity[0, 0] < 0   # Surface flow from right (high) to left (low).
    assert velocity[1, 1] > 0   # Left up.
    assert velocity[1, 2] < 0   # Right down.


def test_linear_amplitude_and_sign_response():
    c = SimulationConfig(nx=8, nz=16, t_end=.05)
    small = SloshingSolver(c)
    alpha2 = np.rad2deg(np.arctan(-2*c.slope))
    large = SloshingSolver(replace(c, alpha_deg=alpha2))
    s1 = list(small.snapshots())[-1][0]
    s2 = list(large.snapshots())[-1][0]
    np.testing.assert_allclose(s2.velocity, -2*s1.velocity, rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(s2.eta, -2*s1.eta, rtol=1e-8, atol=1e-12)


def test_constant_surface_fixes_pressure_gauge_without_spurious_motion():
    solver = SloshingSolver(SimulationConfig(nx=8, nz=16, t_end=.025))
    eta = .02 * np.ones_like(solver.fem.surface_x)
    state = list(solver.snapshots(eta_initial=eta))[-1][0]
    assert np.max(abs(state.velocity)) < 1e-12
    q = solver.integrator.instantaneous_pressure(state)
    np.testing.assert_allclose(q, solver.config.g * .02, atol=1e-10)


def test_second_order_time_convergence():
    solutions = []
    for dt in (.01, .005, .0025):
        solver = SloshingSolver(SimulationConfig(nx=8, nz=16, dt=dt, t_end=.1, snapshot_dt=.1))
        state = list(solver.snapshots())[-1][0]
        solutions.append(state)
    f = solver.fem
    def difference(a, b):
        return np.sqrt(2 * sum(f.energy(a.velocity-b.velocity, a.eta-b.eta)))
    ratio = difference(solutions[0], solutions[1]) / difference(solutions[1], solutions[2])
    assert 3.5 < ratio < 4.5
