import numpy as np
import pytest
import sympy as sp

from sloshing.validation.stiffness import amplification, scalar_step, scalar_time_rows


@pytest.mark.parametrize("lam", [-1., -100., -10000.])
@pytest.mark.parametrize("method", ["midpoint", "sdirk2"])
def test_actual_scalar_stages_match_analytic_stability_function(lam, method):
    assert scalar_step(method, 1., lam, .01) == pytest.approx(amplification(method, lam*.01))


def test_second_order_scalar_decay():
    for row in scalar_time_rows():
        if "order" in row:
            assert 1.95 < row["order"] < 2.1


def test_l_stability_is_not_a_stability():
    # R_mid(z)=(1+z/2)/(1-z/2) -> -1. R_sdirk(z)=(1+(1-2g)z)/(1-gz)^2 -> 0.
    assert abs(amplification("sdirk2", -1e8)) < 5e-8
    assert amplification("midpoint", -1e8) < -.9999999
    for z in -np.logspace(-4, 8, 30):
        assert abs(amplification("sdirk2", z)) <= 1


def test_a_stability_on_imaginary_boundary_and_infinity():
    y = sp.symbols("y", real=True)
    z = sp.symbols("z")
    gamma = 1-1/sp.sqrt(2)
    # |denominator(iy)|²-|numerator(iy)|² = gamma^4*y^4 >= 0.
    difference = (1+gamma**2*y**2)**2-(1+(1-2*gamma)**2*y**2)
    assert sp.simplify(difference-gamma**4*y**4) == 0
    assert gamma > 0  # The only pole 1/gamma is strictly in the right half-plane.
    rational = (1+(1-2*gamma)*z)/(1-gamma*z)**2
    assert sp.limit(rational, z, sp.oo) == 0
    # Analyticity in Re(z)<=0 and the maximum-modulus principle give A-stability;
    # the zero limit at infinity then supplies L-stability.
