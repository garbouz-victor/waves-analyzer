"""Real tiny nonlinear solves plus intentional native/operator corruptions.

Synthetic geometry is used only for the explicitly labelled algebra tests.
No fixture is a mission result or a five-second qualification.
"""
from dataclasses import replace
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest
from skfem import MeshTri, MeshTri2, asm

from sloshing.pinned_wetting import free_boundary as production
from sloshing.pinned_wetting import free_boundary_run as execution
from sloshing.pinned_wetting.free_boundary_compare import compare_runs, curve_distance
from sloshing.pinned_wetting.free_boundary_geometry_audit import right_branch_gap
from sloshing.pinned_wetting.free_boundary_verify import (
    MODEL_ID, _actions, _basis, _boundary_facets, _canonical, _determinant, _facet_quadrature, _read_matrix,
    _expected_rezone, _left_height_is_global_record, _scalar_mass, _scatter_scalar, _solid_contact_events,
    _quadratic_curve_intersections, _self_intersects,
    curve_samples, determinant_minima,
    physical_contract, rebuild, verifier_dependencies, verify_native,
)

ROOT = Path(__file__).resolve().parents[2]


def test_full_quadratic_roots_detect_intersection_missed_by_endpoint_chords():
    # Synthetic P(s)=(s,s^2), Q(t)=(1-t,1-3t+3t^2).  Their sole
    # topological join is P(1)=Q(0), but P(.5)=Q(.5) is an extra root.
    geometry = np.array([[0., .5, 1., .5, 0.], [0., .25, 1., .25, 1.]])
    assert not _self_intersects(geometry[:, ::2].T)
    roots = _quadratic_curve_intersections(geometry, np.arange(5))
    assert roots and all(event["edges"] == [0, 1] for event in roots)
    root = next(event for event in roots if event.get("strictly_interior"))
    assert json.loads(json.dumps(roots)) == roots  # real failure reports must serialize
    np.testing.assert_allclose(root["parameters"], [.5, .5], atol=1e-14, rtol=0.)
    np.testing.assert_allclose(root["position_m"], [.5, .25], atol=1e-14, rtol=0.)
    assert root["equation_residual_m"] <= 1e-14 and root["transverse"]


def test_quadratic_root_arbitrarily_close_to_shared_vertex_is_not_deleted_by_radius():
    epsilon = 2.**-30
    # A real second root only ~1e-9 in parameter from the known join.
    geometry = np.array([[0., .5, 1., .5, 0.],
                         [0., 0., 0., .25-.5*epsilon, 1.-epsilon]])
    roots = _quadratic_curve_intersections(geometry, np.arange(5))
    root = next(event for event in roots if event.get("strictly_interior"))
    np.testing.assert_allclose(root["parameters"], [1.-epsilon, epsilon], atol=1e-14, rtol=0.)


@pytest.mark.parametrize("slope", [0., np.tan(np.deg2rad(2.)), -.21])
def test_quadratic_straight_initial_line_has_only_allowed_shared_vertices(slope):
    x = np.linspace(-1., 1., 43)
    assert _quadratic_curve_intersections(np.array([x, slope*x]), np.arange(len(x))) == []


def test_curved_adjacent_edges_with_only_true_shared_endpoint_are_not_intersections():
    # Two restrictions of y=x^2 to [0,1] and [1,2], no overlapping arc.
    geometry = np.array([[0., .5, 1., 1.5, 2.], [0., .25, 1., 2.25, 4.]])
    assert _quadratic_curve_intersections(geometry, np.arange(5)) == []


def test_collinear_nonincident_overlap_is_not_lost_as_zero_resultant():
    x = np.array([0., .5, 1., 1.5, 2., 1.25, .5])
    events = _quadratic_curve_intersections(np.array([x, np.zeros_like(x)]), np.arange(7))
    assert any(event["edges"] == [0, 2] and event["kind"] == "collinear_quadratic_overlap"
               for event in events)


def test_coincident_curved_arc_is_not_lost_as_zero_resultant():
    geometry = np.array([[0., .5, 1., .5, 0.], [0., .25, 1., .25, 0.]])
    events = _quadratic_curve_intersections(geometry, np.arange(5))
    assert any(event["kind"] == "coincident_quadratic_overlap" for event in events)


def test_parallel_nonincident_segments_with_tiny_positive_gap_are_not_collinear_overlap():
    gap = 2.**-48  # below the root residual scale; this is not a physical cutoff
    geometry = np.array([[0., .5, 1., 1., 1., .5, 0.],
                         [0., 0., 0., .5*gap, gap, gap, gap]])
    events = _quadratic_curve_intersections(geometry, np.arange(7))
    assert not any(event["edges"] == [0, 2] for event in events)


def test_nonincident_exact_endpoint_touch_is_not_a_permitted_topological_join():
    geometry = np.array([[0., .5, 1., 1., 1., .5, 0.], [0., 0., 0., .5, 1., .5, 0.]])
    events = _quadratic_curve_intersections(geometry, np.arange(7))
    assert any(event["edges"] == [0, 2] for event in events)


@pytest.mark.parametrize("gap", [0., 2.**-50])
def test_exact_tangency_and_positive_quadratic_gap_are_distinct(gap):
    # Edge0 lies on y=0; edge2 has y=(t-.5)^2+gap.  This very small
    # positive gap must not become a false contact through a residual cutoff.
    geometry = np.array([[0., .5, 1., 1., 1., .5, 0.],
                         [0., 0., 0., .125+.5*gap, .25+gap, gap, .25+gap]])
    events = [event for event in _quadratic_curve_intersections(geometry, np.arange(7))
              if event["edges"] == [0, 2]]
    if gap:
        assert events == []
    else:
        assert events
        assert any(event["kind"] == "exact_quadratic_contact" and not event["transverse"]
                   for event in events)


@pytest.mark.parametrize("strategy", ["harmonic_p2", "straight_interior"])
def test_adaptive_rezone_never_changes_historical_strategy_schedule(strategy):
    controls = {"rezone_strategy": strategy, "rezone_interval_s": .05}
    # A tiny current Jacobian must not invent an adaptive event in old modes.
    assert not _expected_rezone(controls, .021, .02, np.array([1e-8]), np.array([1.]))
    assert _expected_rezone(controls, .07, .02, np.array([1.]), np.array([1.]))


def test_adaptive_rezone_uses_per_cell_fixed_reference_ratios_not_global_min_ratio():
    controls = {"rezone_strategy": "quality_optimized", "rezone_interval_s": .05}
    initial = np.array([10., .1])
    current = np.array([.05, .005])
    assert current.min()/initial.min() > .01  # incorrect global-ratio shortcut
    assert (current/initial).min() < .01
    expected = _expected_rezone(controls, .021, .02, current, initial)
    assert expected  # omitting the actual remap must fail the native schedule check
    assert expected != False
    assert not _expected_rezone(controls, .021, .02, initial, initial)
    # A fabricated remap on an unchanged, healthy early geometry also mismatches.
    assert _expected_rezone(controls, .021, .02, initial, initial) != True


def test_adaptive_rezone_threshold_is_strict_dimensionless_and_requires_enabled_rezoning():
    controls = {"rezone_strategy": "quality_optimized", "rezone_interval_s": .05}
    assert not _expected_rezone(controls, .021, .02, np.array([.01]), np.array([1.]))
    assert _expected_rezone(controls, .021, .02, np.array([.0099]), np.array([1.]))
    for scale in (1e-8, 1e6):
        assert _expected_rezone(controls, .021, .02, scale*np.array([.0099]), scale*np.array([1.]))
    controls["rezone_interval_s"] = 0.
    assert not _expected_rezone(controls, 5., 0., np.array([1e-10]), np.array([1.]))


@pytest.mark.parametrize("curved", [False, True])
@pytest.mark.parametrize("hydrostatic", [False, True])
@pytest.mark.parametrize("skew", [False, True])
def test_tiny_all_row_actions_equal_full_independent_operators(curved, hydrostatic, skew):
    """Synthetic algebra fixture, not a trajectory or physical acceptance."""
    mesh = MeshTri2.from_mesh(MeshTri.init_tensor(np.linspace(-1., 1., 5),
                                                np.linspace(-10., 0., 5)))
    geometry = mesh.p.copy()
    geometry[1] += np.tan(np.deg2rad(2.))*geometry[0]*(geometry[1]+10.)/10.
    boundaries = _boundary_facets(MeshTri2(doflocs=geometry.copy(), t=mesh.t))
    if curved:
        geometry[1] += .02*np.sin(np.pi*(geometry[0]+1.)/2.)*((geometry[1]+10.)/10.)**2
    basis = _basis(geometry, mesh.t, boundaries, 8)
    rng = np.random.default_rng(642)
    # Deliberately nonzero even on constrained rows: equality of an action
    # must not rely on the right/bottom Dirichlet values hiding operator errors.
    velocity = rng.normal(scale=.03, size=basis["velocity"].N)
    acceleration = rng.normal(scale=.06, size=len(velocity))
    pressure = rng.normal(scale=.1, size=basis["pressure"].N)
    matrices = rebuild(geometry, mesh.t, boundaries, 8, hydrostatic_split=hydrostatic,
                       skew_divergence=skew, advecting_velocity=velocity)
    action = _actions(basis, velocity, acceleration, pressure, 8,
                      hydrostatic_split=hydrostatic, skew_divergence=skew)
    expected = (matrices["M"]@acceleration
                +(matrices["bulk"]+matrices["left"]+matrices["geometric"])@velocity
                -matrices["B"].T@pressure-matrices["f"])
    np.testing.assert_allclose(action["momentum"], expected, rtol=0., atol=2e-13)
    np.testing.assert_allclose(action["weak_divergence"], matrices["B"]@velocity,
                               rtol=0., atol=2e-13)
    for key, matrix in (("bulk_power", "bulk"), ("left_power", "left"),
                        ("geometric_power", "geometric")):
        assert abs(action[key]-velocity@(matrices[matrix]@velocity)) < 2e-13
    assert action["bulk_power"] >= 0. and action["left_power"] >= 0.
    assert abs(action["hydrostatic_power"]-velocity@(matrices["f"]-matrices["f_body"])) < 2e-13
    np.testing.assert_array_equal(basis["fixed"], matrices["fixed"])
    assert not {"M", "bulk", "left", "B", "geometric"}.intersection(basis)


@pytest.fixture(scope="module")
def real_native(tmp_path_factory):
    directory = tmp_path_factory.mktemp("real_PW2_native_unit_tests")
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "scipy", "scikit-fem", "h5py")}
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(execution, "REL_RUNS", str(directory))
        patch.setattr(execution, "doctor", lambda root: {"versions": versions, "execution_HEAD": "unit_test_HEAD"})
        result = execution.run_case(ROOT, production.Controls(nx=4, nz=4, dt=.0025, t_end=.005),
                                    "independent-unit", {}, lambda *args: None)
    return Path(result["native"])


@pytest.fixture(scope="module")
def recovery_native(tmp_path_factory):
    directory = tmp_path_factory.mktemp("real_PW2_rezone_unit_tests")
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "scipy", "scikit-fem", "h5py")}
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(execution, "REL_RUNS", str(directory))
        patch.setattr(execution, "doctor", lambda root: {"versions": versions, "execution_HEAD": "recovery_unit_HEAD"})
        controls = production.Controls(nx=4, nz=4, dt=.0025, t_end=.01,
            rezone_interval_s=.005, hydrostatic_split=True, skew_divergence=True)
        result = execution.run_case(ROOT, controls, "recovery-unit", {}, lambda *args: None)
    return Path(result["native"])


def copy_native(real_native, tmp_path):
    directory = tmp_path/"mutated"
    shutil.copytree(real_native.parent, directory)
    return directory/"native.h5"


def replace_matrix(group, matrix):
    matrix = matrix.copy().tocsr()
    matrix.sum_duplicates(); matrix.eliminate_zeros(); matrix.sort_indices()
    for key in ("data", "indices", "indptr"):
        del group[key]
        group.create_dataset(key, data=getattr(matrix, key))
    group.attrs["shape"] = matrix.shape


def alter_identity(path, mutate):
    with h5py.File(path, "r+") as h:
        identity = json.loads(h.attrs["numerical_identity"])
        mutate(identity)
        h.attrs["numerical_identity"] = json.dumps(identity)
        h.attrs["numerical_digest"] = _canonical(identity)
    # Test input mutations, not application writes; use normal serializer APIs.
    with open(path.parent/"identity.json", "w") as stream:
        json.dump(identity, stream)
    resolved = json.loads((path.parent/"resolved_case.json").read_text())
    resolved["physical_contract"] = identity["physical_contract"]
    with open(path.parent/"resolved_case.json", "w") as stream:
        json.dump(resolved, stream)


def test_real_material_native_reconstructed_independently(real_native):
    report = verify_native(real_native)
    assert report["passed"], report
    assert report["accepted_states"] == 3 and report["time_end_s"] == .005
    assert report["actual_5s_trajectory"] is False
    assert report["geometry_resolution_passed"] is None
    assert report["quadratic_curve_configurations_checked"] == 5
    assert report["quadratic_interface_intersection_events"] == []
    assert "not rigorous" in report["quadratic_intersection_qualification"]
    assert report["maxima"]["P2_position_drift_m"] == 0.
    assert report["maxima"]["right_velocity_max_m_s"] == 0.
    assert report["maxima"]["material_GCL_relative"] < 1e-12
    assert report["maxima"]["momentum_residual"] < 1e-8
    assert report["maxima"]["accepted_endpoint_weak_divergence"] >= 0.
    target = verify_native(real_native, require_target=True)
    assert not target["passed"] and "not_real_0_to_5s_trajectory" in target["failures"]


@pytest.mark.parametrize("configuration,expected", [
    ("accepted_endpoint", "full_quadratic_interface_intersection_or_overlap"),
    ("midpoint_collocation", "midpoint_full_quadratic_interface_intersection_or_overlap"),
])
def test_full_quadratic_gate_checks_every_native_geometry_even_if_chords_pass(
        real_native, monkeypatch, configuration, expected):
    """Gate-only negative: no fabricated geometry is labelled solver data."""
    from sloshing.pinned_wetting import free_boundary_verify as module
    with h5py.File(real_native, "r") as h:
        start = h["states/00000000/geometry"][:]
        end = h["states/00000001/geometry"][:]
    wanted = end if configuration == "accepted_endpoint" else .5*(start+end)
    original = module._quadratic_curve_intersections
    fixture = np.array([[0., .5, 1., .5, 0.], [0., .25, 1., .25, 1.]])
    events = original(fixture, np.arange(5))
    assert events and not module._self_intersects(fixture[:, ::2].T)
    calls = []
    def classify(geometry, ids):
        calls.append(np.array(geometry, copy=True))
        return events if np.array_equal(geometry, wanted) else original(geometry, ids)
    monkeypatch.setattr(module, "_quadratic_curve_intersections", classify)
    report = module.verify_native(real_native)
    assert len(calls) == report["quadratic_curve_configurations_checked"] == 5
    assert not report["passed"] and expected in report["failures"], report
    assert all(event["configuration"] == configuration
               for event in report["quadratic_interface_intersection_events"])


@pytest.mark.parametrize("mutation,expected", [
    ("crossing", "midpoint_interface_self_intersection_at_controlled_tessellation"),
    ("wall_contact", "midpoint_unsupported_extra_solid_contact"),
])
def test_midpoint_geometry_gate_cannot_be_skipped_when_endpoint_checks_pass(real_native, monkeypatch,
                                                                         mutation, expected):
    """Diagnostic-negative fixture, not a claimed solver trajectory.

    Keep every accepted native field unchanged. Inject a geometric failure
    ONLY when the auditor asks for the actual first midpoint configuration,
    to prove that passing endpoint geometry cannot bypass the stage gate.
    Exact curve/contact classification is exercised by separate geometry tests.
    """
    from sloshing.pinned_wetting import free_boundary_verify as module
    with h5py.File(real_native, "r") as h:
        midpoint = .5*(h["states/00000000/geometry"][:]+h["states/00000001/geometry"][:])
    injected = []
    if mutation == "crossing":
        original = module.curve_samples
        def curve(geometry, ids, tolerance_m=1e-6):
            if np.array_equal(geometry, midpoint):
                injected.append(True)
                # Two transverse nonadjacent edges: actual classifier returns true.
                return np.array([[-1., 0.], [.6, .4], [-.6, .4], [1., 0.]]), 0.
            return original(geometry, ids, tolerance_m)
        monkeypatch.setattr(module, "curve_samples", curve)
    else:
        original = module._solid_contact_events
        def contacts(geometry, ids):
            if np.array_equal(geometry, midpoint):
                injected.append(True)
                return [{"wall": "right", "edge": 1, "edge_parameter": .5,
                         "position_m": [1., 0.], "kind": "additional_exact_contact",
                         "signed_distance_m": 0., "input_coordinate_8ULP_scale_m": 1e-15}]
            return original(geometry, ids)
        monkeypatch.setattr(module, "_solid_contact_events", contacts)
    report = module.verify_native(real_native)
    assert injected
    assert not report["passed"] and expected in report["failures"], report
    if mutation == "wall_contact":
        assert any(event.get("configuration") == "midpoint_collocation"
                   for event in report["additional_solid_contact_events"])


@pytest.mark.parametrize("mutation,expected", [
    ("right_freed", "right_freed_or_wrong_velocity_constraints"),
    ("right_robin", "forbidden_RIGHT_Robin_even_on_constrained_DOFs"),
    ("left_removed", "full_LEFT_Robin_matrix_mismatch"),
    ("left_wrong_sign", "full_LEFT_Robin_matrix_mismatch"),
    ("normal_components_swapped", "wrong_normal_tangential_components"),
    ("right_endpoint_moves", "limit:P2_position_drift_m"),
    ("right_velocity_moves", "limit:right_velocity_max_m_s"),
    ("wrong_ALE_convection", "limit:material_kinematic_defect_m"),
    ("broken_GCL", "limit:material_GCL_relative"),
    ("remesh_volume_leak", "limit:volume_relative"),
    ("decreased_H", "limit:retained_running_max_error_m"),
    ("right_H_wrong", "limit:right_H_error_m"),
    ("dry_hole", "retained_state_not_single_lower_connected_intervals"),
    ("restart_history_lost", "coating_history_lost_on_restart_or_rejected_trial"),
    ("fake_5s_from_1s", "saved_dt_not_physical_step_interval"),
    ("not_solver_fields", "limit:momentum_residual"),
    ("linear_reference", "not_new_full_nonlinear_model"),
    ("synthetic_source", "synthetic_or_non_solver_source"),
    ("thin_branch_removed", "missing_or_reordered_interface_branch"),
    ("P2_label_moved", "P1_P2_not_immutable_initial_coordinates"),
    ("rejected_trial_promoted", "rejected_or_wrong_step_in_accepted_history"),
    ("fitted_energy_correction", "limit:saved_work_discrepancy"),
])
def test_native_mutations_rejected(real_native, tmp_path, mutation, expected):
    path = copy_native(real_native, tmp_path)
    with h5py.File(path, "r+") as h:
        topology = h["topology"]
        initial = topology["initial_geometry"][:]
        triangles = topology["triangles"][:]
        mesh = MeshTri2(doflocs=initial, t=triangles)
        op = rebuild(initial, triangles, _boundary_facets(mesh))
        final = h["states/00000002"]
        ids = topology["interface_nodes"][:]
        if mutation == "right_freed":
            fixed = topology["fixed_velocity_dofs"][:]
            freed = np.intersect1d(op["wall"]["right"], op["components"][1])[-1]
            del topology["fixed_velocity_dofs"]
            topology.create_dataset("fixed_velocity_dofs", data=fixed[fixed != freed])
        elif mutation == "right_robin":
            scalar = .02*asm(_scalar_mass, op["scalar"].boundary("right", intorder=8))
            extra = _scatter_scalar(scalar, op["components"], op["velocity"].N, 1)
            # Keep only actual right trace DOFs: exactly zero reduced operator.
            support = np.zeros(op["velocity"].N); support[op["wall"]["right"]] = 1.
            extra = extra.multiply(support[:, None]).multiply(support[None, :]).tocsr()
            assert extra[op["free"]][:, op["free"]].nnz == 0 or np.max(abs(extra[op["free"]][:, op["free"]].data)) == 0.
            replace_matrix(h["initial_left_wall_operator"], _read_matrix(h["initial_left_wall_operator"])+extra)
        elif mutation in ("left_removed", "left_wrong_sign"):
            replace_matrix(h["initial_left_wall_operator"], _read_matrix(h["initial_left_wall_operator"])*(0. if mutation == "left_removed" else -1.))
        elif mutation == "normal_components_swapped":
            topology["component_dofs"][:] = topology["component_dofs"][:][::-1]
        elif mutation == "right_endpoint_moves":
            final["geometry"][1, ids[-1]] -= 1e-4
        elif mutation == "right_velocity_moves":
            final["velocity"][op["components"][1][ids[-1]]] = .01
        elif mutation in ("wrong_ALE_convection", "broken_GCL"):
            boundary_nodes = np.unique(np.concatenate([op["scalar"].get_dofs(side).all()
                for side in ("left", "right", "bottom", "surface")]))
            interior = np.setdiff1d(np.arange(initial.shape[1]), boundary_nodes)[0]
            final["geometry"][0, interior] += 1e-3
        elif mutation == "remesh_volume_leak":
            final["geometry"][1, ids[0]] += .01
        elif mutation in ("decreased_H", "right_H_wrong", "dry_hole", "restart_history_lost", "P2_label_moved"):
            coating = json.loads(final.attrs["coating"])
            if mutation == "decreased_H": coating["H"][0] -= .001
            elif mutation == "right_H_wrong": coating["H"][1] += .001
            elif mutation == "dry_hole": coating["dry_intervals"] = [[-1., -.5]]
            elif mutation == "restart_history_lost": coating["accepted_step"] = 0
            else: coating["markers"][1]["height_m"] += .001
            final.attrs["coating"] = json.dumps(coating)
        elif mutation == "fake_5s_from_1s":
            final.attrs["time_s"] = 5.
            h.attrs["last_accepted_time_s"] = 5.
        elif mutation == "not_solver_fields":
            final["pressure_midpoint"][:] = final["pressure_midpoint"][:]+.1
        elif mutation == "linear_reference":
            h.attrs["model_id"] = "PW1_LEFT_NAVIER_RIGHT_NOSLIP_V1"
        elif mutation == "synthetic_source":
            h.attrs["synthetic"] = True
        elif mutation == "thin_branch_removed":
            del topology["interface_nodes"]
            topology.create_dataset("interface_nodes", data=np.r_[ids[:-3], ids[-1]])
        elif mutation == "rejected_trial_promoted":
            final.attrs["accepted"] = False
        elif mutation == "fitted_energy_correction":
            row = json.loads(final.attrs["diagnostics"])
            row["geometry_potential_work"] += 1e-3
            final.attrs["diagnostics"] = json.dumps(row)
    report = verify_native(path)
    assert not report["passed"] and expected in report["failures"], report


@pytest.mark.parametrize("mutation", ["b", "sigma", "swapped_sides", "snapshot_different_physics"])
def test_locked_contract_cannot_be_relabelled(real_native, tmp_path, mutation):
    path = copy_native(real_native, tmp_path)
    def corrupt(identity):
        physical = identity["physical_contract"]
        if mutation == "b": physical["boundaries"]["left"]["slip_length_m"] = 1.
        elif mutation == "sigma": physical["material"]["surface_tension_N_m"] = .001
        elif mutation == "swapped_sides":
            physical["boundaries"]["left"], physical["boundaries"]["right"] = physical["boundaries"]["right"], physical["boundaries"]["left"]
        else: physical["initial"]["alpha_deg"] = 1.
    alter_identity(path, corrupt)
    report = verify_native(path)
    assert not report["passed"] and "locked_physical_contract_mismatch" in report["failures"], report


def test_contract_equals_repository_current_user_contract():
    actual = json.loads((ROOT/"mission/pinned_wetting/solve_to_animation_v2/CONTRACT.json").read_text())
    assert actual["physical_contract"] == physical_contract()


def test_compatible_nonlinear_checkpoint_survives_execution_HEAD_change(tmp_path, monkeypatch):
    """The commit is provenance; numerical bytes and physical contract are fixed."""
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "scipy", "scikit-fem", "h5py")}
    head = ["before_documentation_commit"]
    monkeypatch.setattr(execution, "REL_RUNS", str(tmp_path/"runs"))
    monkeypatch.setattr(execution, "doctor", lambda root: {"versions": versions, "execution_HEAD": head[0]})
    controls = production.Controls(nx=4, nz=4, dt=.0025, t_end=.005)
    first = execution.run_case(ROOT, controls, "restart-unit", {}, lambda *args: None, stop_after=1)
    first_identity = json.loads((Path(first["native"]).parent/"identity.json").read_text())
    head[0] = "after_documentation_commit"
    resumed = execution.run_case(ROOT, controls, "restart-unit", {}, lambda *args: None)
    assert resumed["run_id"] == first["run_id"] and resumed["completed"]
    native = Path(resumed["native"])
    assert json.loads((native.parent/"identity.json").read_text()) == first_identity
    provenance = json.loads((native.parent/"provenance.json").read_text())
    assert [item["HEAD"] for item in provenance["execution_HEAD_history"]] == [
        "before_documentation_commit", "after_documentation_commit"]
    report = verify_native(native)
    assert report["passed"], report


def test_hydrostatic_skew_and_rezoning_reconstructed_independently(recovery_native):
    report = verify_native(recovery_native)
    assert report["passed"], report
    assert report["pressure_variable"] == "pi=p+gz"
    assert report["independently_checked_rezones"] == 1
    assert report["maxima"]["rezone_projection_residual"] < 1e-9
    assert report["maxima"]["rezone_surface_geometry_change_m"] == 0.
    assert report["maxima"]["rezone_momentum_error"] < 1e-10
    assert report["maxima"]["rezone_inverse_map_error_m"] < 1e-10
    assert report["maxima"]["split_energy_relative"] < 1e-6


@pytest.mark.parametrize("mutation,expected", [
    ("wrong_pressure_variable", "pressure_variable_does_not_match_equations"),
    ("remap_surface_cut", "limit:rezone_surface_geometry_change_m"),
    ("remap_source_point", "limit:rezone_inverse_map_error_m"),
    ("remap_projection_multiplier", "limit:rezone_projection_residual"),
    ("remap_right_velocity", "limit:rezone_constrained_velocity_m_s"),
    ("remap_energy_fitted", "limit:rezone_saved_energy_work_discrepancy"),
    ("missing_remap_evidence", "rezone_schedule_or_native_transfer_evidence_missing"),
    ("wrong_geometric_work", "limit:saved_work_discrepancy"),
    ("wrong_surface_gravity_work", "limit:saved_work_discrepancy"),
])
def test_recovery_variant_native_mutations(recovery_native, tmp_path, mutation, expected):
    path = copy_native(recovery_native, tmp_path)
    with h5py.File(path, "r+") as h:
        groups = [group for group in h["states"].values() if "rezoning_before_step" in group]
        group = groups[0]
        transfer = group["rezoning_before_step"]
        if mutation == "wrong_pressure_variable":
            h.attrs["pressure_variable"] = "gauge_p"
        elif mutation == "remap_surface_cut":
            node = h["topology/interface_nodes"][1]
            transfer["geometry"][1, node] += 1e-4
        elif mutation == "remap_source_point":
            transfer["source_reference_coordinates"][0, 0] += 1e-5
        elif mutation == "remap_projection_multiplier":
            transfer["projection_multipliers"][0] += 1.
        elif mutation == "remap_right_velocity":
            right_endpoint = h["topology/interface_nodes"][-1]
            node = h["topology/component_dofs"][1, right_endpoint]
            transfer["velocity"][node] += .01
        elif mutation == "remap_energy_fitted":
            transfer.attrs["energy_work"] += 1e-3
        elif mutation == "missing_remap_evidence":
            del group["rezoning_before_step"]
        else:
            row = json.loads(group.attrs["diagnostics"])
            name = "geometric_divergence_work" if mutation == "wrong_geometric_work" else "hydrostatic_split_work"
            row[name] += 1e-3
            group.attrs["diagnostics"] = json.dumps(row)
    report = verify_native(path)
    assert not report["passed"] and expected in report["failures"], report


def test_synthetic_exact_curved_Jacobian_and_zero_failure_evidence():
    mesh = MeshTri2.from_mesh(MeshTri(np.array([[0., 1., 0.], [0., 0., 1.]]), np.array([[0], [1], [2]])))
    original = _determinant(mesh, np.array([[0.], [0.]]))[:, 0]
    sign = np.sign(original)
    assert determinant_minima(mesh, sign)[0] == pytest.approx(abs(original[0]))
    coordinates = mesh.p.copy()
    coordinates[:, 3] = coordinates[:, 0]
    curved = replace(mesh, doflocs=coordinates)
    # A zero Jacobian is diagnostic data, not an unclassified library exception.
    assert np.isfinite(determinant_minima(curved, sign)).all()
    assert determinant_minima(curved, sign)[0] <= 0.


def test_synthetic_parametric_overhang_is_not_sorted_or_cut():
    geometry = np.array([[-1., 0., .8, 1.05, 1.], [0., -.01, -.02, .01, .03]])
    curve, error = curve_samples(geometry, np.arange(5), tolerance_m=1e-7)
    assert np.array_equal(curve[0], geometry[:, 0]) and np.array_equal(curve[-1], geometry[:, -1])
    assert np.max(curve[:, 0]) > 1.
    assert np.any(np.diff(curve[:, 0]) < 0.)
    assert error <= 1e-7


def test_synthetic_curve_distance_detects_steep_branch_and_bounds_interior_maximum():
    first = np.array([[-1., 0., .9, .99, 1.], [0., -.01, -.015, .02, .035]])
    second = first.copy(); second[0, 3] -= .01
    difference = curve_distance(first, second, maximum_spacing_m=.0005)
    assert difference["sampled_m"] > .001
    assert difference["lower_bound_m"] <= difference["sampled_m"] <= difference["upper_bound_m"]
    assert difference["upper_bound_m"]-difference["lower_bound_m"] < .00026
    same = curve_distance(first, first, maximum_spacing_m=.001)
    assert same["sampled_m"] < 1e-12


def test_full_history_comparison_uses_native_energy_and_complete_curve(real_native):
    result = compare_runs(real_native, real_native, require_target=False)
    assert result["passed"] and result["comparison_times"] == 3
    assert result["common_time_policy"] == "union_of_all_accepted_times"
    assert result["maxima"]["surface_sampled_m"] < 1e-12
    assert result["maxima"]["surface_upper_bound_m"] < .00051
    assert result["maxima"]["P2_drift_m"] == 0.
    assert all(value == 0. for value in result["energy_observable_differences_per_density"].values())
    with pytest.raises(ValueError, match="not_0_to_5s"):
        compare_runs(real_native, real_native)
    with pytest.raises(ValueError, match="temporal_refinement"):
        compare_runs(real_native, real_native, require_target=False, kind="temporal")
    with pytest.raises(ValueError, match="spatial_refinement"):
        compare_runs(real_native, real_native, require_target=False, kind="spatial")


def test_synthetic_normal_gap_is_a_diagnostic_not_a_film_cutoff():
    vertices = np.array([[-1., 1., 1.], [0., np.tan(np.deg2rad(2.)), -10.]])
    mesh = MeshTri2.from_mesh(MeshTri(vertices, np.array([[0], [1], [2]])))
    facet = int(np.flatnonzero(np.all(mesh.facets == np.array([[0], [1]]), axis=0))[0])
    from skfem import Basis, ElementTriP2
    scalar = Basis(mesh, ElementTriP2())
    midpoint = int(scalar.dofs.facet_dofs[0, facet])
    geometry = mesh.p.copy(); geometry[:, midpoint] = [.475, -.1]
    boundaries = {"surface": np.array([facet], dtype=int)}
    report = right_branch_gap(geometry, mesh.t, boundaries, np.array([0, midpoint, 1]))
    assert report["samples"] > 0
    witness = report["minimum_gap_over_cell_width_witness"]
    assert witness["right_wall_intersection_m"][0] == pytest.approx(1.)
    assert witness["gap_normal_m"] > 0. and witness["adjacent_element_normal_width_m"] > 0.
    assert report["molecular_retained_film_thickness_computed"] is False
    assert "passed" not in report


def test_synthetic_curved_boundary_normals_satisfy_Green_identity_without_inverse():
    from skfem import Basis, ElementTriP2
    vertices = np.array([[-1., 1., 1.], [0., np.tan(np.deg2rad(2.)), -10.]])
    mesh = MeshTri2.from_mesh(MeshTri(vertices, np.array([[0], [1], [2]])))
    scalar = Basis(mesh, ElementTriP2(), intorder=8)
    facet = int(np.flatnonzero(np.all(mesh.facets == np.array([[0], [1]]), axis=0))[0])
    coordinates = mesh.p.copy()
    coordinates[:, int(scalar.dofs.facet_dofs[0, facet])] = [.475, -.1]
    mesh = replace(mesh, doflocs=coordinates)
    scalar = Basis(mesh, ElementTriP2(), intorder=8)
    normal_integral = np.zeros(2); moment_integral = np.zeros(2)
    for q in _facet_quadrature(mesh, scalar, mesh.boundary_facets(), 8):
        normal_integral += np.sum(q["normal"]*q["measure"], axis=1)
        moment_integral += np.sum(q["normal"]*q["coordinates"][1]*q["measure"], axis=1)
    assert np.max(abs(normal_integral)) < 1e-12
    assert np.max(abs(moment_integral-np.array([0., np.sum(scalar.dx)]))) < 1e-12


def test_actual_tiny_dt_and_h_comparisons_use_same_physics(tmp_path, monkeypatch):
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "scipy", "scikit-fem", "h5py")}
    monkeypatch.setattr(execution, "REL_RUNS", str(tmp_path/"runs"))
    monkeypatch.setattr(execution, "doctor", lambda root: {"versions": versions, "execution_HEAD": "comparison_unit_HEAD"})
    cases = []
    for role, nx, nz, dt in (("timebase", 4, 4, .0025), ("timefine", 4, 4, .00125), ("spacefine", 8, 8, .00125)):
        result = execution.run_case(ROOT, production.Controls(nx=nx, nz=nz, dt=dt, t_end=.005),
                                    role, {}, lambda *args: None)
        cases.append(result["native"])
    temporal = compare_runs(cases[0], cases[1], kind="temporal", require_target=False)
    spatial = compare_runs(cases[1], cases[2], kind="spatial", require_target=False)
    assert temporal["passed"] and temporal["actual_triangles"] == [32, 32]
    assert spatial["passed"] and spatial["actual_triangles"] == [32, 128]
    assert temporal["comparison_times"] == spatial["comparison_times"] == 5
    assert temporal["maxima"]["P2_drift_m"] == spatial["maxima"]["P2_drift_m"] == 0.


def test_synthetic_record_comparison_reports_absence_and_small_excess():
    from sloshing.pinned_wetting.free_boundary_compare import _marker_differences
    def marker(name, height, source, confirmed):
        return {"marker_id": name, "side": "R" if name == "P2" else "L",
                "height_m": height, "source_time_s": source, "created_at_s": confirmed}
    first = [marker("P1", -.035, 0., 0.), marker("P2", .035, 0., 0.),
             marker("P3", .024, .78, .79), marker("P4", .024002, 2.4, 2.41)]
    second = [marker("P1", -.035, 0., 0.), marker("P2", .035, 0., 0.),
              marker("P3", .024005, .785, .79)]
    result = _marker_differences(first, second)
    assert result["only_in_first"] == ["P4"] and result["only_in_second"] == []
    assert len(result["matched_chronological_records"]) == 1
    difference = result["matched_chronological_records"][0]
    assert difference["height_difference_m"] == pytest.approx(5e-6)
    assert difference["source_time_difference_s"] == pytest.approx(.005)
    assert difference["confirmation_time_difference_s"] == 0.
    second.append(marker("P4", .024015, 2.42, 2.425))
    result = _marker_differences(first, second)
    assert result["matched_chronological_records"][1]["record_excess_not_larger_than_height_difference"]


def test_independent_cache_identity_includes_geometry_helpers(real_native):
    report = verify_native(real_native)
    assert report["passed"] and report["verifier_dependencies"] == verifier_dependencies()
    assert set(report["verifier_dependencies"]) == {
        "free_boundary_verify.py", "free_boundary_geometry_audit.py"}
    assert report["verifier_sha256"] == report["verifier_dependencies"]["free_boundary_verify.py"]


def test_actual_PW2_refinement_CLI_survives_real_documentation_commit(tmp_path):
    """Real tiny CLI/restart in a separate git repo; no parent git mutation."""
    root = tmp_path/"repo"
    mission = root/"mission/pinned_wetting"
    current = mission/"solve_to_animation_v2"
    current.mkdir(parents=True)
    for name in ("CONTRACT.json", "MODEL_NOTES.md"):
        shutil.copyfile(ROOT/"mission/pinned_wetting/solve_to_animation_v2"/name, current/name)
    shutil.copyfile(ROOT/"mission/pinned_wetting/CONFIG.yaml", mission/"CONFIG.yaml")
    (mission/"state.json").write_text(json.dumps({"run_ids": [], "active_jobs": [],
        "budget": {"actual_local_job_wall_s": 0., "auxiliary_wall_time_upper_bound_s": 0.}}))
    scripts = root/"sloshing_visualization/scripts"
    scripts.mkdir(parents=True)
    shutil.copyfile(ROOT/"sloshing_visualization/scripts/pinned_wetting_mission.py", scripts/"pinned_wetting_mission.py")
    (root/"sloshing_visualization/src").symlink_to(ROOT/"sloshing_visualization/src", target_is_directory=True)
    def git(*args):
        return subprocess.check_output(["git", "-c", "user.name=PW2 numerical test",
            "-c", "user.email=pw2@example.invalid", *args], cwd=root, text=True,
            stderr=subprocess.STDOUT).strip()
    git("init", "-q")
    (root/"note.txt").write_text("fixture creation\n")
    git("add", "note.txt"); git("commit", "-qm", "fixture creation")
    initial_HEAD = git("rev-parse", "HEAD")
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    def cli(command, *extra):
        process = subprocess.run([sys.executable, str(scripts/"pinned_wetting_mission.py"), command,
            "--free-boundary", "--mesh-level", "0", "--dt-scale", "0.5", "--horizon", ".01",
            "--case-role", "refinementcheck", "--hydrostatic-split", "--skew-divergence",
            "--rezone-interval", ".005", "--rezone-strategy", "straight_interior", *extra],
            cwd=root, env=env, capture_output=True, text=True, timeout=90)
        assert process.returncode == 0, process.stdout+process.stderr
    cli("run", "--stop-after", "2")
    directory = next((root/execution.REL_RUNS).iterdir())
    identity_bytes = (directory/"identity.json").read_bytes()
    with h5py.File(directory/"native.h5", "r") as h:
        original_geometry = h["states/00000002/geometry"][:]
        original_velocity = h["states/00000002/velocity"][:]
        original_coating = h["states/00000002"].attrs["coating"]
        assert h.attrs["last_accepted_step"] == 2
    (root/"note.txt").write_text("documentation-only commit; numerical sources unchanged\n")
    git("add", "note.txt"); git("commit", "-qm", "documentation only")
    updated_HEAD = git("rev-parse", "HEAD")
    assert updated_HEAD != initial_HEAD
    cli("resume", "--run-id", directory.name, "--stop-after", "4")
    with h5py.File(directory/"native.h5", "r") as h:
        np.testing.assert_array_equal(h["states/00000002/geometry"][:], original_geometry)
        np.testing.assert_array_equal(h["states/00000002/velocity"][:], original_velocity)
        assert h["states/00000002"].attrs["coating"] == original_coating
        assert h.attrs["last_accepted_step"] == 4
    cli("run")
    assert len(list((root/execution.REL_RUNS).iterdir())) == 1
    assert (directory/"identity.json").read_bytes() == identity_bytes
    provenance = json.loads((directory/"provenance.json").read_text())
    assert [entry["HEAD"] for entry in provenance["execution_HEAD_history"]] == [initial_HEAD, updated_HEAD]
    with h5py.File(directory/"native.h5", "r") as h:
        assert h.attrs["last_accepted_step"] == 8 and h.attrs["last_accepted_time_s"] == .01
        assert json.loads(h.attrs["numerical_identity"])["controls"]["dt"] == .00125
        assert "rezoning_before_step" in h["states/00000005"]
    before = (directory/"native.h5").read_bytes()
    cli("run")
    assert (directory/"native.h5").read_bytes() == before
    report = verify_native(directory/"native.h5")
    assert report["passed"] and report["independently_checked_rezones"] == 1, report


def test_synthetic_unmarked_plateau_is_in_global_record_history():
    heights = np.array([-.035, .03, .03, .02, .025, .022])
    # The later local maximum beats P1 but not the earlier unmarked plateau.
    assert heights[3] < heights[4] > heights[5]
    assert not _left_height_is_global_record(heights, 4)
    assert _left_height_is_global_record(heights, 1)
    heights[4] = np.nextafter(.03, 0.)
    assert not _left_height_is_global_record(heights, 4)
    heights[4] = .031
    assert _left_height_is_global_record(heights, 4)


def test_native_marker_below_earlier_unmarked_plateau_is_rejected(tmp_path, monkeypatch):
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "scipy", "scikit-fem", "h5py")}
    monkeypatch.setattr(execution, "REL_RUNS", str(tmp_path/"runs"))
    monkeypatch.setattr(execution, "doctor", lambda root: {"versions": versions, "execution_HEAD": "plateau_negative_fixture"})
    result = execution.run_case(ROOT, production.Controls(nx=4, nz=4, dt=.0025, t_end=.0125),
                                "plateau-negative-fixture", {}, lambda *args: None)
    path = Path(result["native"])
    # Intentional native corruption, not a physical solution or a solver input.
    with h5py.File(path, "r+") as h:
        ids = h["topology/interface_nodes"][:]
        heights = [h["states/00000000/geometry"][1, ids[0]], .03, .03, .02, .025, .022]
        for index, height in enumerate(heights):
            state = h["states"][f"{index:08d}"]
            state["geometry"][1, ids[0]] = height
            coating = json.loads(state.attrs["coating"])
            coating["H"][0] = float(max(heights[:index+1]))
            if index == 5:
                coating["markers"].append({"marker_id": "P3", "side": "L", "height_m": heights[4],
                    "created_at_s": float(state.attrs["time_s"]), "state_id": str(h.attrs["run_id"])+":4"})
            state.attrs["coating"] = json.dumps(coating)
    report = verify_native(path)
    assert not report["passed"] and "marker_not_global_LEFT_record_at_source" in report["failures"], report


def test_exact_quadratic_interior_wall_touch_not_at_a_native_node_is_detected():
    # Edge 1 has x(s)=1-(s-1/4)^2.  Every one of its three native nodes
    # is strictly inside, but it touches RIGHT at s=1/4 without crossing.
    geometry = np.array([[-1., 0., .9375, .9375, .4375, .7, 1.],
                         [-.035, -.02, -.01, .0, .015, .025, .035]])
    events = _solid_contact_events(geometry, np.arange(7))
    assert any(event["wall"] == "right" and event["edge"] == 1
               and event["edge_parameter"] == .25 and event["kind"] == "additional_exact_contact"
               for event in events)


@pytest.mark.parametrize("wall", ["left", "bottom"])
def test_exact_quadratic_other_solid_contacts_are_detected(wall):
    if wall == "left":
        geometry = np.array([[-1., -.95, -.9, -1., -.9, 0., 1.],
                             [-.035, -.02, -.01, .0, .015, .025, .035]])
    else:
        geometry = np.array([[-1., -.5, 0., .5, 1.], [-.035, -10., -.035, 0., .035]])
    assert any(event["wall"] == wall and event["kind"] == "additional_exact_contact"
               for event in _solid_contact_events(geometry, np.arange(geometry.shape[1])))


def test_solid_contact_check_excludes_only_known_roots_and_has_no_physical_gap_cutoff():
    # Double roots at the two known contacts are not additional contacts.
    geometry = np.array([[-1., -.75, 0., .75, 1.], [-.035, -.02, 0., .02, .035]])
    assert _solid_contact_events(geometry, np.arange(5)) == []
    geometry = np.array([[-1., 0., .9, 1.-1e-10, .9, .95, 1.],
                         [-.035, -.02, -.01, .0, .015, .025, .035]])
    # A 0.1-nanometre gap is not forbidden by an invented physical-film cutoff.
    assert _solid_contact_events(geometry, np.arange(7)) == []
    geometry[0, 3] = np.nextafter(1., 0.)
    events = _solid_contact_events(geometry, np.arange(7))
    assert any(event["kind"] == "precision_unresolved_contact" and event["signed_distance_m"] > 0.
               for event in events)
    assert not any(event["kind"] == "additional_exact_contact" for event in events)
    geometry[0, 3] = 1.+1e-14
    assert any(event["kind"] == "crossing_or_outside"
               for event in _solid_contact_events(geometry, np.arange(7)))


def test_coincident_wall_surface_edge_is_not_exempt_as_known_endpoint():
    geometry = np.array([[-1., -1., -1., 0., 1.], [-.035, -.02, 0., .02, .035]])
    assert any(event["kind"] == "coincident_surface_wall_segment"
               for event in _solid_contact_events(geometry, np.arange(5)))


def _synthetic_curved_locator_fixture():
    from skfem import Basis, ElementTriP2
    mesh = MeshTri2.from_mesh(MeshTri(np.array([[0., 1., 0.], [0., .2, 1.]]),
                                    np.array([[0], [1], [2]])))
    scalar = Basis(mesh, ElementTriP2())
    geometry = mesh.p.copy()
    geometry[:, scalar.element_dofs[3, 0]] = [.5, -.2]
    return replace(mesh, doflocs=geometry), scalar.element_dofs


def _independent_test_lagrange(reference):
    first, second = reference
    zero = 1.-first-second
    return np.array([zero*(2*zero-1), first*(2*first-1), second*(2*second-1),
                     4*zero*first, 4*first*second, 4*second*zero])


def test_pure_inverse_fallback_covers_curved_image_outside_nodal_bounds():
    from sloshing.pinned_wetting.free_boundary_remesh import _bounded_inverse_fallback
    mesh, elements = _synthetic_curved_locator_fixture()
    source = np.array([[.42, .8, .1], [.0001, .05, .8]])
    points = mesh.p[:, elements[:, 0]]@_independent_test_lagrange(source)
    assert points[1, 0] < np.min(mesh.p[1])
    # This synthetic map has determinant 1+1.2*x > 0 on the triangle.
    assert np.min(_determinant(mesh, source)) > 0.
    cells = np.full(3, -1, dtype=int); reference = np.zeros((2, 3)); errors = np.zeros(3)
    _bounded_inverse_fallback(mesh, elements, points, cells, reference, errors)
    assert np.array_equal(cells, np.zeros(3, dtype=int))
    assert np.min(reference) >= -1e-10 and np.max(reference.sum(axis=0)) <= 1.+1e-10
    phi = _independent_test_lagrange(reference)
    mapped = mesh.p[:, elements[:, 0]]@phi
    assert np.max(np.linalg.norm(mapped-points, axis=0)) < 1e-10
    np.testing.assert_allclose(reference, source, rtol=0., atol=1e-10)
    # Physical affine fields reproduce independently of the inverse search.
    A = np.array([[2., -.3], [.7, 1.2]]); constant = np.array([1., -.4])
    nodal = A@mesh.p[:, elements[:, 0]]+constant[:, None]
    assert np.max(abs(nodal@phi-(A@points+constant[:, None]))) < 3e-10


def test_pure_inverse_fallback_does_not_accept_outside_target_inside_Bernstein_bounds():
    from sloshing.pinned_wetting.free_boundary_remesh import _bounded_inverse_fallback
    mesh, elements = _synthetic_curved_locator_fixture()
    x = .42
    boundary_z = -x+1.2*x*x
    points = np.array([[x, x], [boundary_z-.02, -2.]])
    # First point lies above Bernstein's z=-.5 lower bound but outside the
    # actual curved triangle. The second is outside even that enclosure.
    assert points[1, 0] > -.5
    cells = np.full(2, -1, dtype=int); reference = np.zeros((2, 2)); errors = np.zeros(2)
    _bounded_inverse_fallback(mesh, elements, points, cells, reference, errors)
    assert np.array_equal(cells, np.full(2, -1, dtype=int))
    assert np.array_equal(reference, np.zeros((2, 2)))
