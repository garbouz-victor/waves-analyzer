import numpy as np
import pytest
from sloshing.multiphase.be_energy_work import bulk_remainder, wall_remainder
from sloshing.multiphase.free_energy import bulk_energy, bulk_derivative, wall_energy, wall_derivative


def test_quadratic_stiff_first_interval_large_trapezoid_small_BE_work():
    rate, dt, old = 10000., .01, 1.
    new = old/(1+dt*rate)
    change = (new*new-old*old)/2
    B = (new-old)**2/2
    old_power, new_power = rate*old*old, rate*new*new
    weak = change+dt*new_power+B
    continuum = change+dt/2*(old_power+new_power)
    gap = dt/2*(old_power-new_power)
    assert old_power > 1e4*new_power
    assert abs(weak) < 1e-15
    assert continuum > 40
    assert continuum == pytest.approx(weak-B+gap)


@pytest.mark.parametrize("old,new", [(-.1, .1), (.8, .9), (-1., .9)])
def test_quartic_wall_cubic_exact_work_and_negative_pieces(old, new):
    sigma, epsilon, theta = .1, .025, 60.
    B = bulk_remainder(old, new, sigma, epsilon)
    assert bulk_energy(new, sigma, epsilon)-bulk_energy(old, sigma, epsilon)+B == pytest.approx(
        bulk_derivative(new, sigma, epsilon)*(new-old), abs=1e-14)
    W = wall_remainder(old, new, sigma, theta)
    assert wall_energy(new, sigma, theta)-wall_energy(old, sigma, theta)+W == pytest.approx(
        wall_derivative(new, sigma, theta)*(new-old), abs=1e-14)
    if (old, new) == (-.1, .1):
        assert B < 0


def test_tiny_isolated_CH_mode_decay_and_first_interval_work():
    pytest.importorskip("dolfinx")
    from scipy.linalg import eig
    from dolfinx import fem
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.equilibrium import solve_equilibrium
    from sloshing.multiphase.linearized_ch import assemble_linearized_ch
    from sloshing.multiphase.isolated_ch import IsolatedCH
    from sloshing.multiphase.diagnostics import Diagnostics
    from sloshing.multiphase.be_work_diagnostics import BEWorkDiagnostics
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=90., g=0., a_x=0., quadrature_degree=12)
    s = CHNSSolver(c); s.initialize(lambda x: np.ones(x.shape[1]))
    prepared = solve_equilibrium(s); s.initialize_prepared_equilibrium(prepared)
    phi = s.state.sub(2).collapse()
    op = assemble_linearized_ch(s, phi)
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    values, vectors = eig(dense)
    index = np.argsort(values.real)[-3]
    rate = values[index].real
    q = op.project(vectors[:, index].real); q /= max(abs(q))
    endpoints = []
    for steps in (10, 20, 40):
        trial = fem.Function(phi.function_space)
        trial.x.array[:] = phi.x.array+1e-5*q
        engine = IsolatedCH(s, trial)
        diagnostics, work = Diagnostics(s), BEWorkDiagnostics(s)
        previous = diagnostics.measure()
        dt = -1/rate/steps
        for i in range(steps):
            engine.step(dt, (i+1)*dt)
            row = diagnostics.measure()
            audit = work.measure(dt, previous, row)
            assert abs(audit["weak_work_defect"]) < 1e-12
            assert abs(audit["decomposition_roundoff"]) < 1e-12
            assert abs(row["mass_error_relative_to_domain"]) < 1e-10
            previous = row
        delta = s.state.sub(2).collapse().x.array-phi.x.array
        amplitude = op.inner(q, delta).real/op.inner(q, q).real/1e-5
        endpoints.append(amplitude)
        assert amplitude == pytest.approx((1-rate*dt)**(-steps), abs=1e-4)
    errors = abs(np.asarray(endpoints)-np.exp(-1.))
    assert errors[0] > errors[1] > errors[2]
    assert np.log2(errors[1]/errors[2]) > .8


def test_tiny_FEM_stiff_first_interval_keeps_t0_power_and_small_weak_work():
    pytest.importorskip("dolfinx")
    from scipy.linalg import eig
    from dolfinx import fem
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.equilibrium import solve_equilibrium
    from sloshing.multiphase.linearized_ch import assemble_linearized_ch
    from sloshing.multiphase.isolated_ch import IsolatedCH
    from sloshing.multiphase.diagnostics import Diagnostics
    from sloshing.multiphase.be_work_diagnostics import BEWorkDiagnostics
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=90., g=0., a_x=0., quadrature_degree=12)
    s = CHNSSolver(c); s.initialize(lambda x: np.ones(x.shape[1]))
    eq = solve_equilibrium(s); s.initialize_prepared_equilibrium(eq)
    phi = s.state.sub(2).collapse(); op = assemble_linearized_ch(s, phi)
    dense = np.column_stack([op.raw(e) for e in np.eye(op.size)])
    values, modes = eig(dense)
    j = np.argmin(values.real)
    q = op.project(modes[:, j].real); q /= max(abs(q))
    initial = fem.Function(phi.function_space); initial.x.array[:] = phi.x.array+1e-5*q
    engine = IsolatedCH(s, initial)
    diagnostics, work = Diagnostics(s), BEWorkDiagnostics(s)
    old = diagnostics.measure()
    dt = 100/abs(values[j].real)
    engine.step(dt, dt)
    new = diagnostics.measure(); audit = work.measure(dt, old, new)
    assert audit["D_old"] > 1000*audit["D_new"]
    assert abs(audit["weak_work_defect"]) < 1e-12
    assert abs(audit["decomposition_roundoff"]) < 1e-12
    assert abs(audit["continuum_local_defect"]) > 10*abs(audit["Delta_E"])
    assert new["cumulative_CH_dissipation"] == pytest.approx(audit["trapezoid_integral"])
