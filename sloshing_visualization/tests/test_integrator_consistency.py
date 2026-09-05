from dataclasses import replace

import numpy as np
import pytest

from sloshing import SimulationConfig, SloshingSolver
from sloshing.fem_spaces import FEMSystem
from sloshing.time_integrator import GAMMA, INTEGRATORS
from sloshing.validation.manufactured import ManufacturedProblem
from sloshing.validation.spatial import manufactured_run


@pytest.mark.parametrize("method", INTEGRATORS)
def test_source_stage_times_and_forced_energy_identity(method):
    c = SimulationConfig(a=1, d=1, nx=8, nz=8, x_grading=0, z_grading=0, t_end=.025)
    fem = FEMSystem(c)
    integrator = INTEGRATORS[method](fem, ManufacturedProblem(c))
    state = integrator.initial_state()
    called = []
    original = integrator.loads
    def record(t):
        called.append(t)
        return original(t)
    integrator.loads = record
    final = integrator.advance(state)
    expected = [c.dt/2] if method == "midpoint" else [GAMMA*c.dt, c.dt]
    np.testing.assert_allclose(called, expected)
    assert abs(final.step_energy_residual) < 1e-12
    assert np.max(abs(fem.B@final.velocity)) < 1e-12


def test_sdirk_stages_satisfy_monolithic_dae():
    c = SimulationConfig(nx=8, nz=16, integrator="sdirk2")
    solver = SloshingSolver(c)
    f, it = solver.fem, solver.integrator
    initial = it.initial_state()
    basev, basee = initial.velocity, initial.eta
    for time in (GAMMA*c.dt, c.dt):
        v, eta, kv, ke, load, q = it._stage(basev, basee, time)
        np.testing.assert_allclose(f.M@kv+f.K@v-f.B.T@q+c.g*f.C@eta, load, atol=1e-11)
        np.testing.assert_allclose(f.B@v, 0, atol=1e-12)
        np.testing.assert_allclose((eta-basee)/(GAMMA*c.dt), f.R@v, atol=1e-12)
        basev = initial.velocity+c.dt*(1-GAMMA)*kv
        basee = initial.eta+c.dt*(1-GAMMA)*ke


@pytest.mark.parametrize("nu", [.01, .1, 1.])
def test_sdirk_production_invariants_and_rk_budget(nu):
    c = SimulationConfig(nx=12, nz=24, nu=nu, alpha_deg=.2, t_end=.1, integrator="sdirk2")
    rows = [r for _, r in SloshingSolver(c).snapshots()]
    assert rows[-1]["total_energy"] < rows[0]["total_energy"]
    for row in rows:
        assert abs(row["energy_balance_relative"]) < 1e-9
        assert row["contact_error"] == 0
        assert row["volume_error"] < 1e-12


def test_sdirk_manufactured_smoke():
    row = manufactured_run(8, integrator_class=INTEGRATORS["sdirk2"])
    assert row["velocity_l2"] < .003
    assert row["pressure_l2"] < .006
    assert row["eta_l2"] < .004


@pytest.mark.parametrize("method", INTEGRATORS)
def test_mms_temporal_second_order_on_fixed_mesh(method):
    c = SimulationConfig(a=1, d=1, nx=8, nz=8, x_grading=0, z_grading=0,
                         t_end=.1, snapshot_dt=.1, nu=.1)
    states = []
    for dt in (.01, .005, .0025):
        f = FEMSystem(replace(c, dt=dt))
        it = INTEGRATORS[method](f, ManufacturedProblem(f.config))
        state = it.initial_state()
        for _ in range(f.config.nsteps):
            state = it.advance(state)
        states.append(state)
    def diff(a, b):
        return np.sqrt(2*sum(f.energy(a.velocity-b.velocity, a.eta-b.eta)))
    ratio = diff(states[0], states[1])/diff(states[1], states[2])
    assert 3.5 < ratio < 4.5
