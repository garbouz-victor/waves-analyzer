import importlib.util
from pathlib import Path

import sympy as sp


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "derive_streamfunction_conditions.py"
SPEC = importlib.util.spec_from_file_location("sf0_derivation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
sf0 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sf0)


def test_all_scripted_sign_regressions_simplify_to_zero():
    checks = sf0.symbolic_checks()
    assert checks
    assert {name: sp.simplify(value) for name, value in checks.items()} == {
        name: 0 for name in checks
    }


def test_componentwise_curl_is_negative_streamfunction_bulk_residual():
    x, z, t = sp.symbols("x z t", real=True)
    rho, nu = sp.symbols("rho nu", positive=True)
    psi = sp.Function("psi")(x, z, t)
    pressure = sp.Function("p")(x, z, t)
    u, w = sp.diff(psi, z), -sp.diff(psi, x)
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, z, 2)

    horizontal = sp.diff(u, t) + sp.diff(pressure, x) / rho - nu * lap(u)
    vertical = sp.diff(w, t) + sp.diff(pressure, z) / rho - nu * lap(w)
    primitive_curl = sp.diff(vertical, x) - sp.diff(horizontal, z)
    bulk = sp.diff(lap(psi), t) - nu * lap(lap(psi))

    assert sp.simplify(primitive_curl + bulk) == 0


def test_free_surface_elimination_has_three_xxz_terms_and_positive_solved_viscosity():
    x, z, t = sp.symbols("x z t", real=True)
    nu, g = sp.symbols("nu g", positive=True)
    psi = sp.Function("psi")(x, z, t)
    eta = sp.Function("eta")(x, t)

    px_over_rho = g * sp.diff(eta, x) - 2 * nu * sp.diff(psi, x, 2, z)
    solved_horizontal = -px_over_rho + nu * (
        sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3)
    )
    expected = -g * sp.diff(eta, x) + nu * (
        3 * sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3)
    )

    assert sp.simplify(solved_horizontal - expected) == 0


def test_hydrostatic_linearization_gives_requested_normal_stress_sign():
    atmospheric, rho, g, eta, pressure, mu, w_z = sp.symbols(
        "P_atm rho g eta p mu w_z", real=True
    )
    base_pressure_at_surface = atmospheric - rho * g * eta
    traction_balance_residual = (
        -(base_pressure_at_surface + pressure) + 2 * mu * w_z + atmospheric
    )
    expected_residual = -(pressure - 2 * mu * w_z - rho * g * eta)

    assert sp.simplify(traction_balance_residual - expected_residual) == 0


def test_psi_only_condition_retains_positive_solved_viscous_term():
    x, z, t = sp.symbols("x z t", real=True)
    nu, g = sp.symbols("nu g", positive=True)
    psi = sp.Function("psi")(x, z, t)
    eta = sp.Function("eta")(x, t)
    primary = (
        sp.diff(psi, z, t)
        - nu * (3 * sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3))
        + g * sp.diff(eta, x)
    )
    derived = sp.diff(primary, t).subs(
        sp.diff(eta, x, t), -sp.diff(psi, x, 2)
    )
    solved_form = sp.diff(psi, z, t, 2) - (
        g * sp.diff(psi, x, 2)
        + nu
        * (3 * sp.diff(psi, x, 2, z, t) + sp.diff(psi, z, 3, t))
    )

    assert sp.simplify(derived - solved_form) == 0


def test_reconstructed_pressure_curl_equals_negative_bulk_residual():
    x, z, t = sp.symbols("x z t", real=True)
    nu = sp.symbols("nu", positive=True)
    psi = sp.Function("psi")(x, z, t)
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, z, 2)
    px = -sp.diff(psi, z, t) + nu * (
        sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3)
    )
    pz = sp.diff(psi, x, t) - nu * (
        sp.diff(psi, x, 3) + sp.diff(psi, x, z, 2)
    )
    pressure_curl = sp.diff(px, z) - sp.diff(pz, x)
    bulk = sp.diff(lap(psi), t) - nu * lap(lap(psi))

    assert sp.simplify(pressure_curl + bulk) == 0


def test_proposed_mms_forcings_close_every_forced_equation():
    mms = sf0.manufactured_residuals()
    x, z, t, nu, g = sp.symbols("x z t nu g", real=True)
    psi, eta = mms["psi"], mms["eta"]
    lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, z, 2)

    assert sp.simplify(sp.diff(lap(psi), t) - nu * lap(lap(psi)) - mms["bulk_forcing"]) == 0
    assert sp.simplify(sp.diff(eta, t) + sp.diff(psi, x).subs(z, 0) - mms["kinematic_forcing"]) == 0
    assert sp.simplify((sp.diff(psi, z, 2) - sp.diff(psi, x, 2)).subs(z, 0) - mms["shear_forcing"]) == 0
    assert sp.simplify(
        (
            sp.diff(psi, z, t)
            - nu * (3 * sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3))
            + g * sp.diff(eta, x)
        ).subs(z, 0)
        - mms["dynamic_forcing"]
    ) == 0
    for name in (
        "left_wall_psi",
        "left_wall_psix",
        "right_wall_psi",
        "right_wall_psix",
    ):
        assert mms[name] == 0
