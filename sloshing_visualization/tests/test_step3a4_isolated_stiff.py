"""Tiny nonlinear benchmark tests, never production stiff-bump evidence."""
import numpy as np
import pytest


def test_failed_SNES_is_not_retried_or_published_as_accepted():
    from types import SimpleNamespace
    from sloshing.multiphase.isolated_ch import IsolatedCH
    from sloshing.multiphase.performance import Performance
    calls = []
    def fail():
        calls.append(True)
        raise RuntimeError("SNES failure")
    engine = IsolatedCH.__new__(IsolatedCH)
    engine.solver = SimpleNamespace(step_number=0, time=0.)
    engine.dt = SimpleNamespace(value=1.)
    engine.problem = SimpleNamespace(solve=fail)
    engine.performance = Performance()
    with pytest.raises(RuntimeError, match="SNES failure"):
        engine.step(1e-12, 1e-12)
    assert calls == [True]
    assert engine.solver.step_number == 0 and engine.solver.time == 0.
    assert engine.performance.snapshot()["categories"]["nonlinear_SNES"]["count"] == 1


def test_chemical_and_phase_form_equivalence_with_zero_velocity():
    pytest.importorskip("dolfinx")
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.isolated_ch import IsolatedCH
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=60., g=0., a_x=0., quadrature_degree=12)
    full = CHNSSolver(c)
    initial = lambda x: .3+.2*np.cos(3*x[0])*np.cos(2*x[1])
    iso = IsolatedCH(full, initial)
    assert iso.state.sub(1).collapse().x.array == pytest.approx(
        full.state.sub(3).collapse().x.array, abs=0, rel=0)
    changed = lambda x: .31+.21*np.cos(3*x[0])*np.cos(2*x[1])
    mu = lambda x: .2+np.sin(2*x[0])*np.cos(x[1])
    for state, phase_index, mu_index in ((full.state, 2, 3), (iso.state, 0, 1)):
        state.sub(phase_index).interpolate(changed)
        state.sub(mu_index).interpolate(mu)
        state.x.scatter_forward()
    iso.dt.value = 1/16
    for constant, value in zip(full.coefficients, (16., -16., 0.)):
        constant.value = value
    a = fp.assemble_vector(fem.form(full.F))
    b = fp.assemble_vector(fem.form(iso.F))
    for fi, ii in ((2, 0), (3, 1)):
        _, fm = full.space.sub(fi).collapse()
        _, im = iso.state.function_space.sub(ii).collapse()
        af, bi = a.array[fm], b.array[im]
        assert np.linalg.norm(af-bi)/max(np.linalg.norm(af), np.linalg.norm(bi), 1e-30) <= 1e-12
    a.destroy(); b.destroy()


def test_nonlinear_isolated_stiff_tiny_three_level_convergence():
    pytest.importorskip("dolfinx")
    from scipy.linalg import eig
    from dolfinx import fem
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.equilibrium import solve_equilibrium
    from sloshing.multiphase.linearized_ch import assemble_linearized_ch
    from sloshing.multiphase.isolated_ch import IsolatedCH
    from sloshing.multiphase.diagnostics import Diagnostics
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(
        nx=4, nz=4, theta_equilibrium_deg=90., g=0., a_x=0., quadrature_degree=12)
    s = CHNSSolver(c); s.initialize(lambda x: np.ones(x.shape[1]))
    eq = solve_equilibrium(s); s.initialize_prepared_equilibrium(eq)
    phi = s.state.sub(2).collapse(); op = assemble_linearized_ch(s, phi)
    values, vectors = eig(np.column_stack([op.raw(e) for e in np.eye(op.size)]))
    order = np.argsort(values.real)
    q = op.project(vectors[:, order[0]].real+vectors[:, order[len(order)//2]].real)
    q /= max(abs(q))
    horizon = 1/abs(values[order[0]].real)
    ends, powers, defects = [], [], []
    for steps in (10, 20, 40):
        initial = fem.Function(phi.function_space)
        initial.x.array[:] = phi.x.array+1e-4*q
        engine = IsolatedCH(s, initial)
        diagnostic = Diagnostics(s); previous = diagnostic.measure()
        original = engine.state.x.array.copy()
        for j in range(steps):
            engine.step(horizon/steps, (j+1)*horizon/steps)
            row = diagnostic.measure()
            assert abs(row["mass_error_relative_to_domain"]) <= 1e-10
            assert row["E_total"]-previous["E_total"] <= 1e-14
            assert np.array_equal(engine.old.x.array, engine.state.x.array)
            previous = row
        assert not np.array_equal(original, engine.state.x.array)
        ends.append(s.state.sub(2).collapse().x.array.copy())
        powers.append(row["cumulative_CH_dissipation"])
        defects.append(abs(row["energy_budget_defect"]))
    differences = [op.norm(b-a) for a, b in zip(ends, ends[1:])]
    assert differences[0] > differences[1] > 0
    assert np.log2(differences[0]/differences[1]) > .8
    assert defects[0] > defects[1] > defects[2]
    changes = abs(np.diff(powers))
    assert changes[0] > changes[1] and changes[1]/powers[-1] <= .05


def test_cached_wall_trace_equals_uncached_and_contour_roots():
    pytest.importorskip("dolfinx")
    from sloshing.multiphase.config import ModelConfig
    from sloshing.multiphase.solver import CHNSSolver
    from sloshing.multiphase.ch_validation_diagnostics import WallTrace
    from sloshing.multiphase.interface import contour_segments, connected_polylines
    from sloshing.multiphase.contact_line import all_wall_crossings
    c = ModelConfig.from_json("configs/step3/benchmark_flat.json").changed(nx=4, nz=4)
    s = CHNSSolver(c)
    cache = WallTrace(s)
    center = .5*(c.x_min+c.x_max)
    for fraction in (.23, .16):
        radius = fraction*(c.x_max-c.x_min)
        s.state.sub(2).interpolate(lambda x: (x[0]-center)**2-radius**2)
        s.state.x.scatter_forward()
        cached = cache.crossings()
        assert cached == WallTrace(s).crossings()
        assert cached == pytest.approx([center-radius, center+radius], abs=1e-12)
        direct = all_wall_crossings(connected_polylines(contour_segments(s)), c.z_min, axis=1)
        assert cached == pytest.approx(sorted(r["coordinate"] for r in direct), abs=1e-12)
