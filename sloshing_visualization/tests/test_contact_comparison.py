import numpy as np

from sloshing.validation.contact import surface_difference, surface_regions, surface_value


def test_p2_evaluation_and_cross_mesh_quadrature():
    xa = np.array([-1., -.75, -.5, 0, .5, .75, 1.])
    xb = np.array([-1., -.5, 0, .5, 1.])
    ea, eb = 1+xa**2, 1+xb**2
    q = np.linspace(-1, 1, 101)
    np.testing.assert_allclose(surface_value(xa, ea, q), 1+q*q, atol=1e-14)
    assert surface_difference(xa, ea, xb, eb) < 1e-14
    np.testing.assert_allclose(surface_difference(xa, ea, xb, 0*eb)**2, 56/15, atol=1e-14)
    norms = surface_regions(xa, ea, xb, 0*eb)
    np.testing.assert_allclose(norms["full"]**2, norms["intermediate"]**2+norms["contact"]**2, atol=1e-14)
