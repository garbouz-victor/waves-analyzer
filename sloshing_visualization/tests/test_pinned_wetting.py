"""SYNTHETIC negative fixtures never become scientific outputs."""
from copy import deepcopy
import fcntl
import json

import h5py
import numpy as np
import pytest

from sloshing.config import SimulationConfig
from sloshing.solver import SloshingSolver
from sloshing.pinned_wetting.audit import blocker_check, check_consistency
from sloshing.pinned_wetting.coating import Coating
from sloshing.pinned_wetting.mission import lock_busy
from sloshing.pinned_wetting.reference import run_reference


def fixture():
    """SYNTHETIC source: one qualified inner peak between endpoint samples."""
    times = np.array([0., 1., 2., 5.])
    raw = np.array([[-.03, .03], [.05, -.01], [.01, .02], [-.02, -.01]])
    reach = raw.copy()
    reach[1, 0] = .06
    H = np.maximum.accumulate(reach, axis=0)
    eta = np.column_stack([raw[:, 0], np.zeros(4), raw[:, 1]])
    markers = [{"marker_id": "P1", "side": "L", "height_m": -.03, "created_at_s": 0.},
               {"marker_id": "P2", "side": "R", "height_m": .03, "created_at_s": 0.},
               {"marker_id": "P3", "side": "L", "height_m": .06, "created_at_s": 2.}]
    intervals = np.stack([np.full_like(H, -10.), H], axis=-1)
    frames = [{"state_index": i, "time_s": t, "H": H[i].copy(), "intervals": intervals[i].copy(),
               "markers": deepcopy([m for m in markers if m["created_at_s"] <= t])} for i, t in enumerate(times)]
    return {"times": times, "eta": eta, "raw_R": raw, "reach": reach, "H": H,
        "x": np.array([-1., 0., 1.]), "slope": .03, "mean": 0., "bottom": -10.,
        "native_R": raw.copy(), "native_reach": reach.copy(), "native_H": H.copy(),
        "intervals": intervals, "accepted": [True]*4, "frames": frames, "markers": markers,
        "source_run_id": "SYNTHETIC", "manifest_run_id": "SYNTHETIC", "refinement_run_id": "SYNTHETIC",
        "source_config_hash": "fixture", "manifest_config_hash": "fixture", "refinement_config_hash": "fixture",
        "synthetic": True, "source_kind": "SYNTHETIC", "contact_method": "synthetic_test_only",
        "initial_velocity_max": 0.}


def test_positive_fixture_is_explicitly_synthetic():
    d = fixture()
    assert check_consistency(d, t_end=5.) == []
    assert "reference_not_target" in check_consistency(d, require_target=True, t_end=5.)


def corrupt(d, case):
    if case == "H_decreases":
        d["H"][-1, 0] -= .01
    elif case == "H_without_reach":
        d["H"][-1, 0] += .01
    elif case == "P2_moves":
        d["frames"][-1]["markers"][1]["height_m"] += .01
    elif case == "P3_disappears":
        d["frames"][-1]["markers"].pop()
    elif case == "dry_hole":
        d["intervals"][-1, 0, 0] = -.01
    elif case == "curved_initial":
        d["eta"][0, 1] = .001
    elif case == "extrapolated_contact":
        d["contact_method"] = "fixed_endpoint_near_wall_extrapolation"
    elif case == "mixed_run":
        d["refinement_run_id"] = "OTHER_RUN"
    elif case == "mixed_physics":
        d["refinement_config_hash"] = "OTHER_PHYSICS"
    elif case == "missing_t5":
        d["times"][-1] = 4.9
    elif case == "lost_inner_peak":
        d["reach"][1, 0] = d["raw_R"][1, 0]
    elif case == "render_state_disagrees":
        d["frames"][-1]["H"][0] += .01
    elif case == "synthetic_scientific":
        d["source_kind"] = "real_reference_solver"
    elif case == "rejected_state":
        d["accepted"][-1] = False
    return d


@pytest.mark.parametrize("case, expected", [
    ("H_decreases", "decreasing_memory"), ("H_without_reach", "memory_without_reach"),
    ("P2_moves", "marker_moved"), ("P3_disappears", "marker_disappeared"),
    ("dry_hole", "dry_hole_in_W"), ("curved_initial", "initial_not_straight"),
    ("extrapolated_contact", "invalid_contact_method"), ("mixed_run", "mixed_run_ids"),
    ("mixed_physics", "mixed_physics"), ("missing_t5", "missing_final_solver_time"),
    ("lost_inner_peak", "missed_detector_peak"), ("render_state_disagrees", "render_memory_mismatch"),
    ("synthetic_scientific", "synthetic_claimed_scientific"), ("rejected_state", "rejected_state_in_history")])
def test_required_negative_cases(case, expected):
    assert expected in check_consistency(corrupt(fixture(), case), t_end=5.)


def test_connected_retained_state_and_immutable_ids():
    coat = Coating.initial("SYNTHETIC", -10., -.03, .03)
    coat.accept(1, .1, [.05, -.02], reach=[.06, .01])
    coat.accept(2, .2, [-.02, -.04])
    coat.record_peak("P3", "L", .06, .1, "SYNTHETIC:1")
    restored = Coating.from_dict(json.loads(json.dumps(coat.to_dict())))
    assert restored.to_dict() == coat.to_dict()
    for side, H in zip(("L", "R"), [.06, .03]):
        assert all(restored.liquid_present(side, float(z)) for z in np.linspace(-10, H, 1001))
        assert not restored.liquid_present(side, H+1e-15)
    with pytest.raises(ValueError, match="reused"):
        restored.record_peak("P2", "R", .1, .1, "SYNTHETIC:1")
    with pytest.raises(ValueError):
        restored.accept(3, .3, [np.nan, 0.])
    assert restored.H == [.06, .03]


def test_real_checkpoint_continuation(tmp_path):
    c = SimulationConfig(nx=8, nz=12, dt=.0025, t_end=.025, snapshot_dt=.005, integrator="sdirk2")
    identity = {"run_id": "TEST_REAL_REFERENCE", "source_hash": "test_same_source"}
    uninterrupted, resumed = tmp_path / "uninterrupted", tmp_path / "resumed"
    uninterrupted.mkdir(); resumed.mkdir()
    noop = lambda *args: None
    run_reference(uninterrupted, c, identity, noop)
    partial = run_reference(resumed, c, identity, noop, stop_after=5)
    assert partial == {"completed": False, "step": 5, "time_s": .0125}
    run_reference(resumed, c, identity, noop)
    with h5py.File(uninterrupted / "native.h5") as a, h5py.File(resumed / "native.h5") as b:
        for key in a["accepted"]:
            for field in ("eta", "velocity", "H", "wet_intervals"):
                np.testing.assert_array_equal(a["accepted"][key][field][:], b["accepted"][key][field][:])
            assert a["accepted"][key].attrs["coating"] == b["accepted"][key].attrs["coating"]
    with pytest.raises(ValueError, match="identity"):
        run_reference(resumed, c, {**identity, "source_hash": "changed"}, noop)
    with h5py.File(resumed / "native.h5", "r+") as corrupt_native:
        last = corrupt_native["accepted"]["10"]
        bad = json.loads(last.attrs["coating"])
        bad["time"] += .001
        last.attrs["coating"] = json.dumps(bad)
    with pytest.raises(ValueError, match="Checkpoint"):
        run_reference(resumed, c, identity, noop)


def test_independent_zero_disturbance_equilibrium():
    """Exact gravity/hydrostatic reference: no disturbance implies zero motion."""
    solver = SloshingSolver(SimulationConfig(nx=8, nz=12, alpha_deg=0., t_end=.02, dt=.0025, snapshot_dt=.01))
    state = solver.integrator.initial_state()
    for _ in range(8):
        state = solver.integrator.advance(state)
    assert np.max(abs(state.velocity)) == 0.
    assert np.max(abs(state.eta)) == 0.
    assert state.dissipated_energy == 0.


def test_existing_job_lock_prevents_duplicate(tmp_path):
    path = tmp_path / "job.lock"
    assert not lock_busy(path)
    assert not path.exists()
    with path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert lock_busy(path)
    assert not lock_busy(path)


def test_reproducible_two_path_blocker():
    r = blocker_check()
    assert r["status"] == "NEEDS_MODEL_DECISION"
    sharp = r["sharp_material_no_slip"]
    assert sharp["endpoint_trace_nonzeros"] == [0, 0]
    assert sharp["endpoint_displacement_m"] == [0., 0.]
    assert sharp["interior_surface_change_max_m"] > 0
    assert all(value == 0. for value in r["existing_diffuse_sigma_zero"]["energy_terms_max_abs"].values())
