"""Pure constitutive/operator checks; no FEniCSx dependency in ordinary CI."""
from pathlib import Path
import runpy
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.free_energy import (lambda_from_sigma, bulk_energy,
    wall_energy, wall_derivative, transition_width)
from sloshing.multiphase.material import density, viscosity, assert_admissible
from sloshing.multiphase.initialization import flat_interface, effective_acceleration
from sloshing.multiphase.time_integrator import derivative_coefficients
from sloshing.multiphase.film import classify_film, nusselt_flux

CONFIG = Path(__file__).resolve().parents[1]/"configs/step3/benchmark_flat.json"


def test_step3_tanh_normalization_symbolic_and_quadrature():
    s, e, lam = sp.symbols("s e lam", positive=True)
    phi = sp.tanh(s/(sp.sqrt(2)*e))
    residual = lam/e*(phi**3-phi)-lam*e*sp.diff(phi, s, 2)
    assert sp.simplify(residual) == 0
    # Substitute y=tanh(s/(sqrt(2)e)), ds=sqrt(2)e/(1-y²)dy.
    y = sp.symbols("y")
    sigma_exact = sp.integrate(lam/sp.sqrt(2)*(1-y*y), (y, -1, 1))
    assert sp.simplify(sigma_exact-2*sp.sqrt(2)*lam/3) == 0
    sigma, eps = 0.073, 0.017
    energy = quad(lambda x: 2*bulk_energy(np.tanh(x/(np.sqrt(2)*eps)), sigma, eps),
                  -20*eps, 20*eps, epsabs=1e-12)[0]
    assert energy == pytest.approx(sigma, rel=2e-10)


def test_step3_mass_identity_and_momentum_flux_orientation():
    x, z, t, rb, rd, mobility = sp.symbols("x z t rb rd mobility")
    phi = x*x*z*sp.exp(-t)
    mu = sp.sin(x)*z*z
    u, w = x*x, -2*x*z  # divergence-free, independent of phi/mu
    rho = rb+rd*phi
    jx, jz = -rd*mobility*sp.diff(mu, x), -rd*mobility*sp.diff(mu, z)
    ch = sp.diff(phi, t)+sp.diff(phi*u, x)+sp.diff(phi*w, z)-mobility*(sp.diff(mu,x,2)+sp.diff(mu,z,2))
    mass = sp.diff(rho,t)+sp.diff(rho*u+jx,x)+sp.diff(rho*w+jz,z)
    assert sp.simplify(mass-rd*ch) == 0
    # div(u tensor flux) - (flux.grad)u = u div(flux), component by component.
    flux = [rho*u+jx, rho*w+jz]
    for component in (u, w):
        lhs = sum(sp.diff(component*f, y)-f*sp.diff(component,y)
                  for f,y in zip(flux,(x,z)))
        rhs = component*sum(sp.diff(f,y) for f,y in zip(flux,(x,z)))
        assert sp.simplify(lhs-rhs) == 0


def test_step3_pressure_transform_and_gravity_not_double_counted():
    x,z,g,ax,rb,rd = sp.symbols("x z g ax rb rd")
    phi, mu0, p = x*x+z, sp.sin(x*z), x*z*z
    psi = g*z-ax*x
    rho, mu = rb+rd*phi, mu0+rd*psi
    pi = p-phi*mu0+rb*psi
    for y in (x,z):
        original = -sp.diff(p,y)+mu0*sp.diff(phi,y)-rho*sp.diff(psi,y)
        transformed = -sp.diff(pi,y)-phi*sp.diff(mu,y)
        assert sp.simplify(original-transformed) == 0


@pytest.mark.parametrize("theta", [60,90,120])
def test_step3_young_sign_and_planar_wall_condition(theta):
    sigma, epsilon = .1,.02
    assert wall_energy(-1,sigma,theta)-wall_energy(1,sigma,theta) == pytest.approx(
        sigma*np.cos(np.deg2rad(theta)), abs=1e-15)
    phi=np.array([-.8,0.,.8])
    # Bottom outward normal -ez, gradient into liquid gives partial_n phi
    # = cos(theta)*(1-phi²)/(sqrt(2)epsilon).
    dn = np.cos(np.deg2rad(theta))*(1-phi*phi)/(np.sqrt(2)*epsilon)
    assert np.max(np.abs(lambda_from_sigma(sigma)*epsilon*dn+
                         wall_derivative(phi,sigma,theta))) < 1e-15


def test_step3_parameters_material_overshoot_and_groups():
    c=ModelConfig.from_json(CONFIG)
    assert density(1,c) == c.rho_liquid and viscosity(-1,c) == c.mu_gas
    assert c.groups()["capillary_length_infinite"]
    assert c.groups()["Bo"] == 0
    with pytest.raises(ValueError):
        c.changed(slip_length_m=0)
    with pytest.raises(RuntimeError, match="no clipping"):
        assert_admissible(np.array([-1.1]), c.changed(rho_liquid=100))
    assert c.fingerprint() != c.changed(mobility=2).fingerprint()


def test_step3_positive_tilt_and_effective_gravity():
    points=np.array([[-.5,.5],[0.,0.]])
    assert flat_interface(.02,alpha_deg=2)(points)[1] > 0
    assert flat_interface(.02,alpha_deg=2)(points)[0] < 0
    assert effective_acceleration(2,9.81)/9.81 == pytest.approx(np.tan(np.deg2rad(2)))


def test_step3_bdf2_derivative_polynomial_and_be_startup():
    dt=.02
    b=np.array(derivative_coefficients(2,dt))
    times=np.array([1.,1.-dt,1.-2*dt])
    assert b@times**2 == pytest.approx(2.,abs=1e-13)
    assert derivative_coefficients(1,dt)[2] == 0


def test_step3_film_classification_not_adsorption():
    assert classify_film(.01,.01,.001)=="unresolved"
    assert classify_film(.03,.01,.001)=="marginal"
    assert classify_film(.04,.01,.001)=="resolved_candidate"
    assert classify_film(.04,.001,.01)=="marginal"
    assert nusselt_flux(.01,9.81,.01)==pytest.approx(9.81*.01**3/(3*.01))
    assert nusselt_flux(.01,9.81,.01,.001) > nusselt_flux(.01,9.81,.01)
    assert transition_width(.01)/.01 == pytest.approx(4.164065537,rel=1e-8)


@pytest.mark.parametrize("benchmark", ["flat","contact","operators","startup-energy"])
def test_step3_analysis_mode_never_falls_through_to_solver(benchmark):
    cli=runpy.run_path(str(CONFIG.parents[2]/"scripts/step3_benchmarks.py"))
    with pytest.raises(SystemExit) as exc:
        cli["parse_arguments"](["--benchmark",benchmark,"--analyze-only"])
    assert exc.value.code == 2


def test_step3_analysis_uses_saved_parameters_and_rejects_unimplemented_study():
    cli=runpy.run_path(str(CONFIG.parents[2]/"scripts/step3_benchmarks.py"))
    parse=cli["parse_arguments"]
    assert parse(["--benchmark","laplace","--analyze-only"]).analyze_only
    for args in (["--benchmark","laplace","--analyze-only","--dt",".001"],
                 ["--benchmark","contact","--study"],
                 ["--estimate-only","--analyze-only"]):
        with pytest.raises(SystemExit) as exc:
            parse(args)
        assert exc.value.code == 2
