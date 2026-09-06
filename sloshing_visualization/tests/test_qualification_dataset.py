from dataclasses import replace

import h5py
import numpy as np
import pytest

from sloshing.config import SimulationConfig
from sloshing.solver import SloshingSolver
from sloshing.fem_spaces import FEMSystem
from sloshing.postprocess.vorticity import VorticityProjector
from sloshing.validation.regions import VorticityRegions
from sloshing.validation.qualification.dataset import run_dataset, validate_cache
from sloshing.validation.qualification.metrics import DerivativeMetrics, pinned_discrete_equilibrium, slope_regions


def test_exact_regional_slopes_include_cut_inside_p2_edge():
    x = np.linspace(-1, 1, 9)
    slopes = slope_regions(x, x*x)
    np.testing.assert_allclose([slopes[k] for k in ("slope_global", "slope_bulk", "slope_intermediate", "slope_contact")],
                               [2, 1.8, 1.96, 2], atol=1e-14)


def test_regional_derivatives_against_analytic_polynomial_integral():
    f = FEMSystem(SimulationConfig(a=1, d=1, nx=4, nz=6))
    p = VorticityProjector(f)
    metrics = DerivativeMetrics(f, VorticityRegions(f, p))
    v = np.zeros(f.velocity.N)
    x, z = f.scalar.doflocs
    v[f.component_dofs[0]], v[f.component_dofs[1]] = x*x, z*z
    result = metrics.evaluate(v)
    exact = np.sqrt(4*(2*.9**3/3+1.8/3))
    np.testing.assert_allclose(result["divergence_bulk_l2"], exact, atol=2e-12)
    np.testing.assert_allclose(result["gradient_bulk_l2"], exact, atol=2e-12)


def test_static_surface_interpretation_matches_actual_trace_nullspace():
    f = FEMSystem(SimulationConfig(nx=8, nz=16, alpha_deg=.02))
    eta, info = pinned_discrete_equilibrium(f.surface_x, f.initial_surface()[[0, -1]], f.config.g)
    np.testing.assert_allclose(f.C@eta, 0, atol=1e-18)
    assert abs(f.surface_weights@eta) < 1e-18
    np.testing.assert_allclose(info["potential_energy"], .5*f.config.g*eta@(f.S@eta), atol=1e-20)


def test_complete_prefix_continuation_matches_uninterrupted(tmp_path):
    c = SimulationConfig(nx=4, nz=6, alpha_deg=.02, integrator="sdirk2",
                         dt=.00125, snapshot_dt=.005, t_end=.01)
    prefix = run_dataset(c, tmp_path/"prefix.h5")
    full = replace(c, t_end=.02)
    continued = run_dataset(full, tmp_path/"continued.h5", prefix)
    direct = run_dataset(full, tmp_path/"direct.h5")
    with h5py.File(continued) as a, h5py.File(direct) as b:
        for group in ("fields", "diagnostics", "restart"):
            for key in a[group]:
                np.testing.assert_allclose(a[f"{group}/{key}"][:], b[f"{group}/{key}"][:], atol=1e-16, rtol=1e-11)
        assert a["fields/q"].shape[0] == 5 and "visualization" not in a
    validate_cache(prefix, c)  # Original complete prefix was not changed.
    assert run_dataset(full, continued) == continued
    with h5py.File(continued, "r+") as h:
        original_value = h["fields/u"][1, 1]
        h["fields/u"][1, 1] = np.nan
    with pytest.raises(RuntimeError, match="Non-finite fields/u at snapshot 1"):
        validate_cache(continued)
    with h5py.File(continued, "r+") as h:
        h["fields/u"][1, 1] = original_value
    with h5py.File(continued, "r+") as h:
        h.attrs["status"] = "failed"
    with pytest.raises(RuntimeError, match="Unfinished"):
        run_dataset(full, continued)


@pytest.mark.parametrize("nu", [.01, 1.])
def test_stiffly_accurate_stage_pressure_is_instantaneous_pressure(nu):
    solver = SloshingSolver(SimulationConfig(nx=8, nz=16, alpha_deg=.02, nu=nu, integrator="sdirk2"))
    it = solver.integrator
    state = it.initial_state()
    for _ in range(5):
        state = it.advance(state)
        np.testing.assert_allclose(it.last_stage_pressure, it.instantaneous_pressure(state), atol=3e-13, rtol=2e-9)


def test_prefix_target_cannot_be_overwritten_by_concurrent_creation(tmp_path, monkeypatch):
    import sloshing.validation.qualification.dataset as module
    c = SimulationConfig(nx=4, nz=6, alpha_deg=.02, integrator="sdirk2",
                         dt=.00125, snapshot_dt=.005, t_end=.005)
    prefix = run_dataset(c, tmp_path/"prefix.h5")
    target = tmp_path/"another_owner.h5"
    original = module.SloshingSolver
    def concurrent_creation(config):
        solver = original(config)
        target.write_bytes(b"another owner created this after the initial cache check")
        return solver
    monkeypatch.setattr(module, "SloshingSolver", concurrent_creation)
    with pytest.raises(FileExistsError):
        run_dataset(replace(c, t_end=.01), target, prefix)
    assert target.read_bytes().startswith(b"another owner")
