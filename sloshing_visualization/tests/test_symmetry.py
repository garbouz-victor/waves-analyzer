import numpy as np
import pytest

from sloshing import SimulationConfig, SloshingSolver
from sloshing.diagnostics import quadratic_edge_max


@pytest.mark.parametrize("nx", [8, 9])
def test_mesh_and_release_preserve_reflection_symmetry(nx):
    solver = SloshingSolver(SimulationConfig(nx=nx, nz=16, t_end=.05))
    f = solver.fem
    state = list(solver.snapshots())[-1][0]
    np.testing.assert_allclose(state.eta, -state.eta[::-1], atol=1e-12)
    points = np.array([[-.9, -.7, -.3, -.1], [-.05, -.2, -1., -3.]])
    reflected = points * np.array([[-1], [1]])
    full = f.expand_velocity(state.velocity)
    left = f.velocity_at(full, points)
    right = f.velocity_at(full, reflected)
    np.testing.assert_allclose(left[0], right[0], atol=1e-12)
    np.testing.assert_allclose(left[1], -right[1], atol=1e-12)


def test_wall_max_includes_quadratic_extrema_between_dofs():
    # p(s)=1-(s-.25)^2; its maximum is not a nodal value.
    left, middle, right = (np.array([1-(s-.25)**2]) for s in (0., .5, 1.))
    assert quadratic_edge_max(left, middle, right) == pytest.approx(1.)
