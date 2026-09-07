#!/usr/bin/env python3
"""Reproducible SF-0 sign checks for the linear streamfunction formulation.

This script is deliberately symbolic and cheap.  It is a regression guard, not
a proof and not a numerical sloshing solver.
"""

from __future__ import annotations

import math

import sympy as sp


def laplacian(expression: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> sp.Expr:
    return sp.diff(expression, x, 2) + sp.diff(expression, z, 2)


def symbolic_checks() -> dict[str, sp.Expr]:
    """Return expressions that must all simplify to zero."""

    x, z, t = sp.symbols("x z t", real=True)
    rho, nu, g = sp.symbols("rho nu g", positive=True)
    psi = sp.Function("psi")(x, z, t)
    eta = sp.Function("eta")(x, t)
    pressure = sp.Function("p")(x, z, t)

    u = sp.diff(psi, z)
    w = -sp.diff(psi, x)
    omega = sp.diff(w, x) - sp.diff(u, z)
    delta_psi = laplacian(psi, x, z)
    bulk_residual = sp.diff(delta_psi, t) - nu * laplacian(delta_psi, x, z)

    # Primitive momentum is written with every term on the left.  Taking
    # d_x(vertical residual) - d_z(horizontal residual) is an independent
    # componentwise route to the vorticity equation.
    horizontal_residual = (
        sp.diff(u, t)
        + sp.diff(pressure, x) / rho
        - nu * laplacian(u, x, z)
    )
    vertical_residual = (
        sp.diff(w, t)
        + sp.diff(pressure, z) / rho
        - nu * laplacian(w, x, z)
    )
    primitive_curl = sp.diff(vertical_residual, x) - sp.diff(
        horizontal_residual, z
    )

    shear_from_velocity = sp.diff(u, z) + sp.diff(w, x)
    shear_from_psi = sp.diff(psi, z, 2) - sp.diff(psi, x, 2)

    # p/rho - 2 nu w_z = g eta must be identical to
    # p/rho + 2 nu psi_xz = g eta.
    normal_primitive = pressure / rho - 2 * nu * sp.diff(w, z) - g * eta
    normal_streamfunction = (
        pressure / rho + 2 * nu * sp.diff(psi, x, z) - g * eta
    )

    # Differentiate normal stress in x, solve for p_x/rho, and substitute it
    # in horizontal momentum.
    px_over_rho = g * sp.diff(eta, x) - 2 * nu * sp.diff(psi, x, 2, z)
    eliminated_dynamic = (
        sp.diff(psi, z, t)
        + px_over_rho
        - nu * (sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3))
    )
    expected_dynamic = (
        sp.diff(psi, z, t)
        - nu * (3 * sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3))
        + g * sp.diff(eta, x)
    )

    dynamic_time_derivative = sp.diff(expected_dynamic, t)
    psi_only_from_primary = dynamic_time_derivative.subs(
        sp.diff(eta, x, t), -sp.diff(psi, x, 2)
    )
    expected_psi_only = (
        sp.diff(psi, z, t, 2)
        - g * sp.diff(psi, x, 2)
        - nu
        * (
            3 * sp.diff(psi, x, 2, z, t)
            + sp.diff(psi, z, 3, t)
        )
    )

    px_reconstructed = -sp.diff(psi, z, t) + nu * (
        sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3)
    )
    pz_reconstructed = sp.diff(psi, x, t) - nu * (
        sp.diff(psi, x, 3) + sp.diff(psi, x, z, 2)
    )
    pressure_curl = sp.diff(px_reconstructed, z) - sp.diff(
        pz_reconstructed, x
    )

    exx = sp.diff(u, x)
    ezz = sp.diff(w, z)
    exz = (sp.diff(u, z) + sp.diff(w, x)) / 2
    twice_e_contract_e = 2 * (exx**2 + ezz**2 + 2 * exz**2)
    streamfunction_dissipation_density = (
        4 * sp.diff(psi, x, z) ** 2
        + (sp.diff(psi, z, 2) - sp.diff(psi, x, 2)) ** 2
    )

    return {
        "incompressibility": sp.diff(u, x) + sp.diff(w, z),
        "vorticity_definition": omega + delta_psi,
        "componentwise_vorticity_equation": primitive_curl + bulk_residual,
        "tangential_stress": shear_from_velocity - shear_from_psi,
        "normal_stress_streamfunction_sign": (
            normal_primitive - normal_streamfunction
        ),
        "pressure_eliminated_surface_condition": (
            eliminated_dynamic - expected_dynamic
        ),
        "psi_only_surface_condition": psi_only_from_primary - expected_psi_only,
        "pressure_gradient_compatibility": pressure_curl + bulk_residual,
        "strain_dissipation_identity": (
            twice_e_contract_e - streamfunction_dissipation_density
        ),
    }


def manufactured_residuals() -> dict[str, sp.Expr]:
    """Return proposed nondimensional forced-MMS residual expressions."""

    x, z, t, nu, g = sp.symbols("x z t nu g", real=True)
    psi = sp.exp(t + z) * (1 - x**2) ** 2 * (1 + z)
    eta = sp.exp(t) * (x + x**3 / 5)
    delta_psi = laplacian(psi, x, z)

    return {
        "psi": psi,
        "eta": eta,
        "bulk_forcing": sp.simplify(
            sp.diff(delta_psi, t) - nu * laplacian(delta_psi, x, z)
        ),
        "kinematic_forcing": sp.simplify(
            sp.diff(eta, t) + sp.diff(psi, x).subs(z, 0)
        ),
        "shear_forcing": sp.simplify(
            (sp.diff(psi, z, 2) - sp.diff(psi, x, 2)).subs(z, 0)
        ),
        "dynamic_forcing": sp.simplify(
            (
                sp.diff(psi, z, t)
                - nu
                * (3 * sp.diff(psi, x, 2, z) + sp.diff(psi, z, 3))
                + g * sp.diff(eta, x)
            ).subs(z, 0)
        ),
        "left_wall_psi": sp.simplify(psi.subs(x, -1)),
        "left_wall_psix": sp.simplify(sp.diff(psi, x).subs(x, -1)),
        "right_wall_psi": sp.simplify(psi.subs(x, 1)),
        "right_wall_psix": sp.simplify(sp.diff(psi, x).subs(x, 1)),
    }


def physical_scales(
    *, a: float = 1.0, g: float = 9.81, nu: float = 1.0e-6
) -> dict[str, float]:
    q1 = math.pi / (2 * a)
    omega1 = math.sqrt(g * q1)
    return {
        "nu_star": nu / (a * math.sqrt(g * a)),
        "re_g": math.sqrt(g * a) * a / nu,
        "q1": q1,
        "omega1": omega1,
        "period1": 2 * math.pi / omega1,
        "stokes_thickness": math.sqrt(2 * nu / omega1),
    }


def main() -> None:
    checks = symbolic_checks()
    failed = {
        name: sp.simplify(expression)
        for name, expression in checks.items()
        if sp.simplify(expression) != 0
    }
    if failed:
        for name, expression in failed.items():
            print(f"FAIL {name}: {expression}")
        raise SystemExit(1)

    for name in checks:
        print(f"PASS {name}")

    scales = physical_scales()
    for name, value in scales.items():
        print(f"{name} = {value:.16g}")


if __name__ == "__main__":
    main()
