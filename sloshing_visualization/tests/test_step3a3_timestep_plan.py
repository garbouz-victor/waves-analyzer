import copy
import json
import numpy as np
import pytest
from scipy.integrate import quad
from sloshing.multiphase.ch_timestep_planner import (ReducedCHReference, make_plan,
    validate_plan, FrozenSchedule, PlannerPolicy)


@pytest.fixture(scope="module")
def planned():
    A = -np.diag([1., 100., 10000.])
    model = ReducedCHReference(A, np.eye(3), -A, np.ones(3))
    hashes = {k: k+"-fixture" for k in ("operator", "equilibrium", "perturbation", "Krylov")}
    plan = make_plan(model, .1, 10000., hashes)
    return model, hashes, plan


def test_stiff_plan_starts_small_grows_and_passes(planned):
    model, hashes, plan = planned
    assert plan["qualification_status"] == "linear_reduced_pass"
    assert plan["dt0"] <= .05/10000
    assert max(b["dt"] for b in plan["blocks"]) > 2*plan["dt0"]
    errors = plan["predicted_errors"]
    assert errors["max_checkpoint_relative_L2_error"] <= 1e-3
    assert errors["max_quadratic_energy_relative_error"] <= 1e-3
    assert errors["integrated_D2_relative_error"] <= 1e-2
    assert not plan["nonlinear_execution_authorized"]
    integral = quad(lambda t: model.power(model.state(t)), 0., .1, epsabs=1e-12)[0]
    assert model.integrated_power(.1) == pytest.approx(integral, rel=1e-11)


def test_deterministic_JSON_and_fingerprint_rejection(planned):
    model, hashes, plan = planned
    again = make_plan(model, .1, 10000., hashes)
    assert json.dumps(plan, sort_keys=True) == json.dumps(again, sort_keys=True)
    validate_plan(plan, hashes)
    for key in ("operator", "perturbation"):
        with pytest.raises(ValueError, match="fingerprint"):
            validate_plan(plan, {**hashes, key: "changed"})
    changed = copy.deepcopy(plan); changed["blocks"][0]["dt"] *= 2
    with pytest.raises(ValueError, match="integrity"):
        validate_plan(changed, hashes)


def test_observed_nonlinear_state_cannot_change_schedule(planned):
    _, hashes, plan = planned
    frozen = FrozenSchedule(plan, hashes, analysis_only=True)
    original = [frozen.step(i) for i in range(frozen.nsteps)]
    observed_field = np.zeros(100)
    observed_field[:] = 100*np.arange(100)
    # Even modifying the caller's plan object cannot alter the frozen schedule.
    saved = plan["blocks"][0]["dt"]
    plan["blocks"][0]["dt"] = 999.
    assert [frozen.step(i) for i in range(frozen.nsteps)] == original
    plan["blocks"][0]["dt"] = saved
    half = FrozenSchedule(plan, hashes, factor=2, analysis_only=True)
    assert half.nsteps == 2*frozen.nsteps
    assert half.step(0)["dt"] == original[0]["dt"]/2
    assert half.step(half.nsteps-1)["time"] == pytest.approx(.1)


def test_failed_candidate_never_authorizes_nonlinear(planned):
    model, hashes, _ = planned
    plan = make_plan(model, .1, 10000., hashes, PlannerPolicy(max_candidate_refinements=0))
    assert plan["qualification_status"] == "failed"
    assert not plan["nonlinear_execution_authorized"]


def test_reduced_pass_alone_cannot_authorize_nonlinear(planned):
    _, hashes, plan = planned
    with pytest.raises(ValueError, match="Unqualified"):
        FrozenSchedule(plan, hashes)


def test_sparse_schedule_matches_reduced_on_known_tangent_system():
    from scipy.linalg import null_space
    from sloshing.multiphase.linearized_ch import LinearizedCH
    from sloshing.multiphase.ch_timestep_planner import verify_full_sparse_schedule, verify_schedule
    basis = null_space(np.ones((1, 4)))
    rate = -np.diag([1., 100., 10000.])
    op = LinearizedCH(np.eye(4), basis @ (-rate) @ basis.T, np.eye(4))
    reference = ReducedCHReference(rate, np.eye(3), -rate, np.ones(3))
    blocks = [{"block": 0, "t_start": 0., "t_end": .01, "dt": .0001, "steps": 100}]
    reduced = verify_schedule(reference, blocks)
    full = verify_full_sparse_schedule(op, reference, basis, blocks)
    for key in ("endpoint_relative_L2_error", "max_checkpoint_relative_L2_error",
                "max_quadratic_energy_relative_error", "integrated_D2_relative_error"):
        assert full[key] == pytest.approx(reduced[key], abs=1e-11)
    assert full["max_mass_domain"] < 1e-14


def test_analysis_cannot_reserve_historical_or_external_result_paths(tmp_path, monkeypatch):
    from sloshing.multiphase.benchmarks import spectral_ch
    monkeypatch.setattr(spectral_ch, "RESULTS", tmp_path/"step3a3")
    for name in ("../step3a2", str(tmp_path/"outside"), "."):
        with pytest.raises(ValueError, match="Historical/external"):
            spectral_ch.reserve(name)
    assert not (tmp_path/"step3a2").exists()
