import numpy as np
import pytest

from sloshing.diagnostics import Diagnostics, SLOPE_WARNING, surface_slopes


def test_contact_points_are_exactly_pinned(short_run):
    solver, snapshots = short_run
    expected = solver.config.slope * np.array([-solver.config.a, solver.config.a])
    for state, row in snapshots:
        np.testing.assert_array_equal(state.eta[[0, -1]], expected)
        assert row["contact_error"] == 0


def test_exact_quadratic_surface_slope():
    x = np.array([-1., -.75, -.5, .25, 1.])
    eta = x**2 + 3*x
    left, right = surface_slopes(x, eta)
    np.testing.assert_allclose(left, 2*x[:-2:2] + 3)
    np.testing.assert_allclose(right, 2*x[2::2] + 3)


def test_warning_uses_dimensional_eta_slope_without_extra_amplitude_factor(short_run):
    solver, snapshots = short_run
    diagnostic = Diagnostics(solver.fem, snapshots[0][0])
    row = dict(snapshots[0][1], max_slope=.4)
    with pytest.warns(RuntimeWarning, match=SLOPE_WARNING):
        diagnostic.validate(row)
