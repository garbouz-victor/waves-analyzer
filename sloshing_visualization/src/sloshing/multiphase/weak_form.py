"""Monolithic conservative-CH / AGG-skew-momentum residual, dimensional SI."""
import ufl
from .free_energy import bulk_derivative, lambda_from_sigma, wall_derivative
from .material import density, density_derivative, viscosity
from .mesh import WALL_IDS


def residual(state, old, older, coefficients, config, tags):
    u, pressure_pi, phi, mu = ufl.split(state)
    u0, _, phi0, _ = ufl.split(old)
    u1, _, phi1, _ = ufl.split(older)
    v, q, s, chi = ufl.TestFunctions(state.function_space)
    a0, a1, a2 = coefficients
    rho, rho0, rho1 = [density(p, config) for p in (phi, phi0, phi1)]
    eta = viscosity(phi, config)
    drho = density_derivative(config)
    dx = ufl.Measure("dx", domain=state.function_space.mesh,
                     metadata={"quadrature_degree": config.quadrature_degree})
    ds = ufl.Measure("ds", domain=state.function_space.mesh, subdomain_data=tags,
                     metadata={"quadrature_degree": config.quadrature_degree})
    x = ufl.SpatialCoordinate(state.function_space.mesh)
    potential = config.g*x[1]-config.a_x*x[0]
    J = -drho*config.mobility*ufl.grad(mu)
    flux = rho*u+J
    time_momentum = (a0*rho*u+a1*rho0*u0+a2*rho1*u1
                     - 0.5*(a0*rho+a1*rho0+a2*rho1)*u)
    inertia = (ufl.inner(time_momentum, v)
               + 0.5*(ufl.inner(ufl.dot(ufl.grad(u), flux), v)
                      - ufl.inner(ufl.dot(ufl.grad(v), flux), u))) * dx
    momentum = (inertia + (2*eta*ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v)))
                - pressure_pi*ufl.div(v) + phi*ufl.dot(ufl.grad(mu), v))*dx)
    for wall, marker in WALL_IDS.items():
        if wall not in config.no_slip_walls and wall not in config.free_slip_walls:
            n = ufl.FacetNormal(state.function_space.mesh)
            ut = u-ufl.dot(u, n)*n
            vt = v-ufl.dot(v, n)*n
            momentum += eta/config.slip_length_m*ufl.inner(ut, vt)*ds(marker)
    continuity = q*ufl.div(u)*dx
    phase = ((a0*phi+a1*phi0+a2*phi1)*s - phi*ufl.dot(u, ufl.grad(s))
             + config.mobility*ufl.dot(ufl.grad(mu), ufl.grad(s))) * dx
    chemical = ((mu-drho*potential-bulk_derivative(phi, config.sigma, config.epsilon))*chi
                 -lambda_from_sigma(config.sigma)*config.epsilon*
                 ufl.dot(ufl.grad(phi), ufl.grad(chi))) * dx
    for wall in config.wetting_walls:
        chemical -= wall_derivative(phi, config.sigma, config.theta_equilibrium_deg)*chi*ds(WALL_IDS[wall])
    return momentum+continuity+phase+chemical
