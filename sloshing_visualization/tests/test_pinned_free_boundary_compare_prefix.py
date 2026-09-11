"""Explicit temporal-prefix diagnostics; tiny real solves are not mission results."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest

from sloshing.pinned_wetting import free_boundary as production
from sloshing.pinned_wetting import free_boundary_run as execution
from sloshing.pinned_wetting import free_boundary_compare as comparison

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def prefix_natives(tmp_path_factory):
    directory = tmp_path_factory.mktemp("tiny_real_prefix_comparisons")
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "scipy", "scikit-fem", "h5py")}
    paths = {}
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(execution, "REL_RUNS", str(directory))
        patch.setattr(execution, "doctor", lambda root: {"versions": versions,
                                                       "execution_HEAD": "prefix_unit_HEAD"})
        for role, nx, nz, dt, count in (("short", 4, 4, .0025, 2),
                                        ("long", 4, 4, .0025, 4),
                                        ("timefine", 4, 4, .00125, 6),
                                        ("spacefine", 6, 6, .00125, 6)):
            result = execution.run_case(ROOT, production.Controls(nx=nx, nz=nz, dt=dt, t_end=.02),
                "prefix-"+role, {}, lambda *args: None, stop_after=count)
            paths[role] = Path(result["native"])
    return paths


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_prefix_is_explicit_and_does_not_mutate_or_truncate_native(prefix_natives):
    first, second = prefix_natives["short"], prefix_natives["long"]
    before = [sha(first), sha(second)]
    with pytest.raises(ValueError, match="refinement_horizons_differ"):
        comparison.compare_runs(first, second, require_target=False)
    result = comparison.compare_runs(first, second, require_target=False, until_time_s=.005)
    assert result["passed"] and result["explicit_common_interval"]
    assert result["comparison_scope"] == "EXPLICIT_ACCEPTED_COMMON_PREFIX_DIAGNOSTIC"
    assert result["prefix_excludes_native_future_states"]
    with h5py.File(second, "r") as h:
        actual_long_end = float(h["states/00000004"].attrs["time_s"])
    assert actual_long_end > .005  # a tiny real solve may accept adaptive halvings
    assert result["full_native_horizons_s"] == [.005, actual_long_end]
    assert result["full_native_state_counts"] == [3, 5]
    assert result["compared_native_state_counts"] == [3, 3]
    assert result["comparison_times"] == 3
    assert result["common_time_policy"] == "union_of_all_accepted_times_in_explicit_prefix"
    assert result["time_end_s"] == .005 and not result["full_5s_interval_compared"]
    assert [sha(first), sha(second)] == before
    assert json.loads(json.dumps(result))["explicit_common_interval"] is True


@pytest.mark.parametrize("end", [.003, .00375, .0075, .03])
def test_prefix_requires_an_actual_accepted_endpoint_in_both_natives(prefix_natives, end):
    with pytest.raises(ValueError, match="match_one_accepted_state"):
        comparison.compare_runs(prefix_natives["short"], prefix_natives["long"],
                                require_target=False, until_time_s=end)


@pytest.mark.parametrize("end", [0., -.1, np.inf, np.nan, "not-seconds"])
def test_prefix_rejects_invalid_end_without_loading_native(end):
    with pytest.raises(ValueError, match="positive_finite_seconds"):
        comparison.compare_runs("missing", "missing", require_target=False, until_time_s=end)


def test_prefix_cannot_be_used_as_five_second_release_or_sampled_substitute(prefix_natives):
    first, second = prefix_natives["short"], prefix_natives["long"]
    with pytest.raises(ValueError, match="partial_prefix_cannot_qualify"):
        comparison.compare_runs(first, second, until_time_s=.005)
    with pytest.raises(ValueError, match="match_one_accepted_state"):
        comparison.compare_runs(first, second, until_time_s=5.)
    with pytest.raises(ValueError, match="requires_union"):
        comparison.compare_runs(first, second, require_target=False, until_time_s=.005, sample_dt=.005)


def test_clock_tolerance_selects_exact_native_endpoint_fields(prefix_natives, monkeypatch):
    path = prefix_natives["short"]
    with h5py.File(path, "r") as h:
        ids = h["topology/interface_nodes"][:]
        actual = h["states/00000002/geometry"][:][:, ids]
    original = comparison.curve_distance; calls = []
    def distance(one, two, **kwargs):
        calls.append((one.copy(), two.copy()))
        return original(one, two, **kwargs)
    monkeypatch.setattr(comparison, "curve_distance", distance)
    requested = .005+5e-13
    result = comparison.compare_runs(path, path, require_target=False, until_time_s=requested)
    assert result["requested_until_time_s"] == requested
    assert result["compared_native_end_times_s"] == [.005, .005]
    np.testing.assert_array_equal(calls[-1][0], actual)
    np.testing.assert_array_equal(calls[-1][1], actual)
    assert not result["prefix_excludes_native_future_states"]


def test_future_marker_catalog_is_not_read_back_into_a_prefix(prefix_natives, tmp_path):
    """Only a labelled negative input copy receives a fabricated future marker."""
    destination = tmp_path/"future-marker-negative"
    shutil.copytree(prefix_natives["long"].parent, destination)
    path = destination/"native.h5"
    with h5py.File(path, "r+") as h:
        group = h["states/00000004"]
        future_source_time = float(h["states/00000003"].attrs["time_s"])
        future_confirmation_time = float(group.attrs["time_s"])
        coating = json.loads(group.attrs["coating"])
        coating["markers"].append({"marker_id": "P3", "side": "L", "height_m": .01,
            "created_at_s": future_confirmation_time, "state_id": str(h.attrs["run_id"])+":3"})
        group.attrs["coating"] = json.dumps(coating)
    before = sha(path)
    result = comparison.compare_runs(prefix_natives["short"], path,
                                     require_target=False, until_time_s=.005)
    assert [m["marker_id"] for m in result["second_markers"]] == ["P1", "P2"]
    assert result["peak_differences"]["only_in_second"] == []
    loaded_full = comparison._load(path, lambda: None)
    assert loaded_full["markers"][-1]["marker_id"] == "P3"
    assert loaded_full["markers"][-1]["source_time_s"] == future_source_time > .005
    assert sha(path) == before


@pytest.mark.parametrize("mutation,expected", [
    ("physics", "different_physics"),
    ("sources", "different_numerical_sources"),
    ("gamma", "different_numerical_formulations"),
    ("not_accepted", "nonaccepted_native_state"),
])
def test_prefix_does_not_bypass_physics_sources_formulation_or_acceptance(
        prefix_natives, tmp_path, mutation, expected):
    destination = tmp_path/mutation
    shutil.copytree(prefix_natives["long"].parent, destination)
    path = destination/"native.h5"
    with h5py.File(path, "r+") as h:
        ident = json.loads(h.attrs["numerical_identity"])
        if mutation == "physics":
            ident["physical_contract"]["boundaries"]["right"]["type"] = "Navier"
        elif mutation == "sources":
            ident["source_hashes"]["pinned_wetting/free_boundary.py"] = "not-the-source"
        elif mutation == "gamma":
            ident["controls"]["grad_div_gamma_m2_s"] = 1.
        else:
            h["states/00000002"].attrs["accepted"] = False
        h.attrs["numerical_identity"] = json.dumps(ident)
    with pytest.raises(ValueError, match=expected):
        comparison.compare_runs(prefix_natives["short"], path,
                                require_target=False, until_time_s=.005)


def test_prefix_keeps_independent_dt_h_guards_and_full_parametric_curve(prefix_natives, monkeypatch):
    original = comparison.curve_distance; node_counts = []
    def distance(one, two, **kwargs):
        node_counts.append((one.shape[1], two.shape[1]))
        return original(one, two, **kwargs)
    monkeypatch.setattr(comparison, "curve_distance", distance)
    temporal = comparison.compare_runs(prefix_natives["short"], prefix_natives["timefine"],
        kind="temporal", require_target=False, until_time_s=.005)
    spatial = comparison.compare_runs(prefix_natives["timefine"], prefix_natives["spacefine"],
        kind="spatial", require_target=False, until_time_s=.005)
    assert temporal["passed"] and spatial["passed"]
    assert temporal["comparison_times"] == spatial["comparison_times"] == 5
    assert node_counts == [(9, 9)]*5+[(9, 13)]*5
    with pytest.raises(ValueError, match="temporal_refinement"):
        comparison.compare_runs(prefix_natives["short"], prefix_natives["long"],
            kind="temporal", require_target=False, until_time_s=.005)
    with pytest.raises(ValueError, match="spatial_refinement"):
        comparison.compare_runs(prefix_natives["short"], prefix_natives["timefine"],
            kind="spatial", require_target=False, until_time_s=.005)
