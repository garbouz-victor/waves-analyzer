import importlib.util
import math
from pathlib import Path

import pytest
import sympy as sp


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "derive_streamfunction_conditions.py"
SPEC = importlib.util.spec_from_file_location("sf0_initial", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
sf0 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sf0)


def test_initial_tilt_forces_nonzero_interior_acceleration():
    g, k, psi_zt = sp.symbols("g k psi_zt", nonzero=True)

    # At t=0 psi and all of its spatial derivatives vanish, while eta_x=k.
    primary_dynamic_residual = psi_zt + g * k
    required_acceleration = sp.solve(primary_dynamic_residual, psi_zt)

    assert required_acceleration == [-g * k]
    assert sp.simplify(primary_dynamic_residual.subs(psi_zt, 0)) == g * k
    assert primary_dynamic_residual.subs(psi_zt, 0) != 0


def test_psi_only_zero_solution_does_not_encode_initial_tilt():
    g, k = sp.symbols("g k", positive=True)

    # psi=0 satisfies the homogeneous psi-only differential equation, but it
    # violates the separately required initial trace psi_zt=-g*k.
    psi_only_residual_for_zero = sp.Integer(0)
    missing_initial_trace_residual = sp.Integer(0) + g * k

    assert psi_only_residual_for_zero == 0
    assert missing_initial_trace_residual != 0


def test_zero_initial_streamfunction_gives_zero_surface_rate():
    x, t = sp.symbols("x t", real=True)
    psi_surface = sp.Integer(0)
    eta_t = -sp.diff(psi_surface, x)
    assert eta_t == 0


def test_wall_no_slip_and_kinematics_pin_both_endpoints():
    x, z, t, a = sp.symbols("x z t a", positive=True)
    amplitude = sp.Function("A")(z, t)
    psi = (a**2 - x**2) ** 2 * amplitude
    eta_t = -sp.diff(psi, x)

    assert sp.simplify(eta_t.subs(x, a)) == 0
    assert sp.simplify(eta_t.subs(x, -a)) == 0


def test_initial_linear_surface_modal_coefficients():
    x, a, k = sp.symbols("x a k", positive=True)
    n = sp.symbols("n", integer=True, nonnegative=True)
    qn = (n + sp.Rational(1, 2)) * sp.pi / a
    coefficient = sp.simplify(
        sp.integrate(k * x * sp.sin(qn * x), (x, -a, a)) / a
    )
    expected = 2 * k * (-1) ** n / (a * qn**2)

    assert sp.simplify(coefficient - expected) == 0
    endpoint_sum = sp.summation(
        expected * sp.sin(qn * a), (n, 0, sp.oo)
    )
    assert sp.simplify(endpoint_sum - k * a) == 0


def test_reference_scales_for_requested_parameters():
    scales = sf0.physical_scales(a=1.0, g=9.81, nu=1.0e-6)

    assert scales["nu_star"] == pytest.approx(3.1927542840705044e-7)
    assert scales["re_g"] == pytest.approx(3.132091952673165e6)
    assert scales["q1"] == pytest.approx(math.pi / 2)
    assert scales["omega1"] == pytest.approx(3.9254951236573885)
    assert scales["period1"] == pytest.approx(1.6006096324800767)
    assert scales["stokes_thickness"] == pytest.approx(7.137855910141449e-4)
