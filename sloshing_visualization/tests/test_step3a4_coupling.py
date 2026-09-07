"""Tiny CH-vs-CHNS smoke evidence is NOT a production coupling qualification."""
import numpy as np
import pytest


def test_tiny_matched_density_coupling_decreases_with_amplitude():
    pytest.importorskip("dolfinx")
    from dolfinx import fem
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.isolated_ch import IsolatedCH
    from sloshing.multiphase.linearized_ch import assemble_linearized_ch
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=90., g=0., a_x=0., quadrature_degree=12,
        dt=1e-6, t_end=3e-6, time_scheme="be")
    observer, full = CHNSSolver(c), CHNSSolver(c)
    full.initialize(lambda x: np.ones(x.shape[1]))
    eq = full.state.sub(2).collapse(); op = assemble_linearized_ch(full, eq)
    psi = fem.Function(eq.function_space)
    psi.interpolate(lambda x: np.cos(3*x[0])*np.cos(2*x[1]))
    q = op.project(psi.x.array)
    discrepancies = []
    for amplitude in (1e-4, 5e-5):
        initial = fem.Function(eq.function_space)
        initial.x.array[:] = 1+amplitude*q
        engine = IsolatedCH(observer, initial)
        full.initialize(initial)
        full.step_number, full.time = 0, 0.
        for i in range(3):
            engine.step(c.dt, (i+1)*c.dt)
            full.advance()
        a = observer.state.sub(2).collapse().x.array
        b = full.state.sub(2).collapse().x.array
        error = op.norm(a-b)
        discrepancies.append(error)
        assert error/max(op.norm(a-1), 1e-30) < 1e-3
        assert abs(op.mass @ (a-initial.x.array))/op.area < 1e-10
        assert abs(op.mass @ (b-initial.x.array))/op.area < 1e-10
    assert discrepancies[1] <= max(discrepancies[0], 1e-12)
