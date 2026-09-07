import sympy as sp


def test_streamfunction_velocity_strain_and_dissipation_identity():
    x, z, t, mu = sp.symbols("x z t mu", real=True)
    psi = sp.exp(t + z) * (1 - x**2) ** 2 * (1 + z)
    u, w = sp.diff(psi, z), -sp.diff(psi, x)

    exx = sp.diff(u, x)
    ezz = sp.diff(w, z)
    exz = (sp.diff(u, z) + sp.diff(w, x)) / 2
    primitive_density = 2 * mu * (exx**2 + ezz**2 + 2 * exz**2)
    streamfunction_density = mu * (
        4 * sp.diff(psi, x, z) ** 2
        + (sp.diff(psi, z, 2) - sp.diff(psi, x, 2)) ** 2
    )

    assert sp.simplify(sp.diff(u, x) + sp.diff(w, z)) == 0
    assert sp.simplify(exx - sp.diff(psi, x, z)) == 0
    assert sp.simplify(ezz + sp.diff(psi, x, z)) == 0
    assert sp.simplify(primitive_density - streamfunction_density) == 0


def test_free_surface_boundary_work_converts_to_potential_energy_rate():
    rho, g, eta, eta_t = sp.symbols("rho g eta eta_t", real=True)
    w = eta_t
    tangential_traction = sp.Integer(0)
    normal_traction = -rho * g * eta
    u = sp.symbols("u", real=True)

    surface_power_density = u * tangential_traction + w * normal_traction
    potential_energy_rate_density = rho * g * eta * eta_t

    assert sp.simplify(surface_power_density + potential_energy_rate_density) == 0


def test_manufactured_dissipation_is_nonnegative_and_nontrivial():
    x, z = sp.symbols("x z", real=True)
    psi = sp.exp(z) * (1 - x**2) ** 2
    density = sp.expand(
        4 * sp.diff(psi, x, z) ** 2
        + (sp.diff(psi, z, 2) - sp.diff(psi, x, 2)) ** 2
    )
    total = sp.integrate(density, (x, -1, 1), (z, -sp.oo, 0))

    assert sp.simplify(total) > 0


def test_normal_stress_sign_gives_negative_gravitational_surface_work():
    rho, nu, g = sp.symbols("rho nu g", positive=True)
    p, w_z, eta = sp.symbols("p w_z eta", real=True)
    normal_condition = sp.Eq(p - 2 * rho * nu * w_z, rho * g * eta)
    fluid_normal_traction = -p + 2 * rho * nu * w_z

    solved_traction = sp.simplify(
        fluid_normal_traction.subs(p, sp.solve(normal_condition, p)[0])
    )
    assert solved_traction == -rho * g * eta
