import numpy as np
import pytest

from sloshing import SimulationConfig, SloshingSolver


def test_all_six_wall_residuals_and_independent_boundary_probes(short_run):
    solver, snapshots = short_run
    f, c = solver.fem, solver.config
    z = np.linspace(-c.d, 0, 301)
    x = np.linspace(-c.a, c.a, 101)
    points = np.concatenate([np.vstack((-c.a * np.ones_like(z), z)),
                             np.vstack((c.a * np.ones_like(z), z)),
                             np.vstack((x, -c.d * np.ones_like(x)))], axis=1)
    for state, row in snapshots:
        for wall in ("left", "right", "bottom"):
            for component in ("u", "w"):
                assert row[f"{wall}_{component}_max"] == 0.0
        full = f.expand_velocity(state.velocity)
        assert np.max(np.abs(f.velocity_at(full, points))) < 1e-13


@pytest.mark.parametrize("nu", [.001, .01, .1, 1.])
def test_every_viscosity_preset_obeys_constraints(nu):
    solver = SloshingSolver(SimulationConfig(nx=8, nz=16, nu=nu, t_end=.025))
    snapshots = list(solver.snapshots())
    assert snapshots[-1][1]["total_energy"] < snapshots[0][1]["total_energy"]
    assert snapshots[-1][1]["contact_error"] == 0
