import numpy as np
import pytest

from sloshing import SimulationConfig, SloshingSolver
from sloshing.diagnostics import Diagnostics, PhysicsViolation


def test_energy_budget_on_every_internal_step(short_run):
    solver, _ = short_run
    f = solver.fem
    state = solver.integrator.initial_state()
    initial = sum(f.energy(state.velocity, state.eta))
    assert initial == pytest.approx(solver.config.g * solver.config.slope**2 * solver.config.a**3 / 3)
    for _ in range(solver.config.nsteps):
        next_state = solver.integrator.advance(state)
        vmid = f.expand_velocity((state.velocity + next_state.velocity) / 2)
        grad = f.velocity.interpolate(vmid).grad
        strain = (grad + grad.swapaxes(0, 1)) / 2
        dissipation = 2 * solver.config.nu * np.sum(strain**2 * f.velocity.dx[None, None])
        e0 = sum(f.energy(state.velocity, state.eta))
        e1 = sum(f.energy(next_state.velocity, next_state.eta))
        assert e1 <= e0 + 1e-14
        assert abs(e1 - e0 + solver.config.dt * dissipation) < 1e-12 * initial
        state = next_state


def test_zero_angle_stays_at_rest():
    solver = SloshingSolver(SimulationConfig(nx=8, nz=16, alpha_deg=0, t_end=.05))
    for state, row in solver.snapshots():
        assert np.count_nonzero(state.velocity) == 0
        assert np.count_nonzero(state.eta) == 0
        assert row["total_energy"] == 0


@pytest.mark.parametrize("key", ["left_u_max", "right_w_max", "bottom_u_max", "volume_error",
                                 "weak_divergence_l2", "contact_error", "energy_balance_relative",
                                 "max_step_energy_increase"])
def test_invalid_diagnostics_fail_loudly(short_run, key):
    solver, snapshots = short_run
    diagnostic = Diagnostics(solver.fem, snapshots[0][0])
    row = dict(snapshots[0][1])
    row[key] = 1.0
    with pytest.raises(PhysicsViolation):
        diagnostic.validate(row)
