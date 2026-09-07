import copy
import json
import numpy as np
import pytest
from sloshing.multiphase.ch_schedule_optimizer import optimize, OptimizerGuard, legal_blocks
from sloshing.multiphase.ch_timestep_planner import (ReducedCHReference, refine_blocks,
    verify_schedule, validate_plan)
from sloshing.multiphase.ch_validation_policy import linear_gates, LINEAR_LIMITS


def fixture():
    model = ReducedCHReference(-np.diag([1., 100., 10000.]), np.eye(3),
        np.diag([1., 100., 10000.]), np.ones(3))
    blocks = [dict(block=0, t_start=0., t_end=2e-5, dt=5e-7, steps=40),
              dict(block=1, t_start=2e-5, t_end=1e-4, dt=1e-6, steps=80)]
    hashes = {k: k for k in ("operator", "equilibrium", "perturbation", "Krylov")}
    return model, blocks, hashes


def test_blockwise_beats_uniform_with_identical_hard_thresholds():
    model, blocks, hashes = fixture()
    assert not all(linear_gates(verify_schedule(model, blocks)).values())
    assert all(linear_gates(verify_schedule(model, refine_blocks(blocks, 2))).values())
    result = optimize(model, blocks, hashes, reference_qualified=True)
    assert result["total_steps"] == 200 < 240
    assert result["refinements"] == [0, 1]
    assert all(linear_gates(result["predicted_errors"]).values())
    assert result["policy"]["linear_limits"] == LINEAR_LIMITS
    assert result["locally_minimal_under_legal_moves"]
    assert not result["global_optimum_claimed"]
    assert legal_blocks(result["blocks"])
    assert any(r["action"] == "trial_merge" for r in result["optimization_history"])


def test_byte_reproducible_history_and_semantic_fingerprint():
    model, blocks, hashes = fixture()
    a = optimize(model, blocks, hashes, reference_qualified=True)
    b = optimize(model, blocks, hashes, reference_qualified=True)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)
    validate_plan(a, hashes)
    for key in ("operator", "perturbation"):
        with pytest.raises(ValueError, match="fingerprint"):
            validate_plan(a, {**hashes, key: "changed"})
    changed = copy.deepcopy(a); changed["blocks"][0]["dt"] *= 2
    with pytest.raises(ValueError, match="integrity"):
        validate_plan(changed, hashes)


def test_no_nonlinear_observations_or_unqualified_reference_inputs():
    model, blocks, hashes = fixture()
    with pytest.raises(TypeError):
        optimize(model, blocks, hashes, reference_qualified=True, observed_field=np.zeros(3))
    with pytest.raises(ValueError, match="qualified"):
        optimize(model, blocks, hashes, reference_qualified=False)
    with pytest.raises(ValueError, match="qualified"):
        optimize(np.eye(3), blocks, hashes, reference_qualified=True)


def test_impossible_with_allowed_hierarchy_fails_finitely_without_relaxation():
    model, blocks, hashes = fixture()
    result = optimize(model, blocks, hashes, reference_qualified=True, guard=OptimizerGuard(max_splits=0))
    assert result["qualification_status"] == "failed"
    assert len(result["optimization_history"]) == 1
    assert result["policy"]["linear_limits"] == LINEAR_LIMITS
    timed = optimize(model, blocks, hashes, reference_qualified=True, guard=OptimizerGuard(evaluations=0))
    assert timed["qualification_status"] == "failed"
    assert "finite optimizer resource guard" in timed["reason"]


def test_illegal_time_growth_and_boundary_changes_rejected():
    _, blocks, _ = fixture()
    assert legal_blocks(blocks)
    blocks[1]["dt"] *= 2; blocks[1]["steps"] //= 2
    assert not legal_blocks(blocks)


def test_role_aware_frozen_schedule_never_changes_with_nonlinear_observation():
    from sloshing.multiphase.ch_validation_diagnostics import ValidationSchedule, observation_indices
    model, blocks, hashes = fixture()
    plan = optimize(model, blocks, hashes, reference_qualified=True)
    sparse = dict(qualification_status="passed", plan_sha256=plan["plan_sha256"], input_fingerprints=hashes)
    with pytest.raises(ValueError, match="not authorized"):
        ValidationSchedule(plan, hashes, sparse, "isolated_CH", {})
    schedule = ValidationSchedule(plan, hashes, sparse, "isolated_CH", {"isolated_CH_execution_authorized": True})
    initial = [schedule.step(i) for i in range(schedule.nsteps)]
    observed_field = np.zeros(3); observed_field[:] = 100.
    plan["blocks"][0]["dt"] = 999.
    assert [schedule.step(i) for i in range(schedule.nsteps)] == initial
    selected = observation_indices(schedule.descriptor["blocks"], 1e-8, 1e-6)
    assert selected["BE_work"][:20] == list(range(1, 21))
    assert all(t in [r["time"] for r in initial] for t in selected["probe_times"])
    assert selected["probe_times"][-1] <= 1e-5
