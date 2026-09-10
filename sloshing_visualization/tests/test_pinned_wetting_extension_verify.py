"""Small real trajectories and deliberately corrupted copies; no release data."""
from dataclasses import asdict
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest

from sloshing.pinned_wetting import extension_run
from sloshing.pinned_wetting.extension_verify import _review_signoff, _signed_native_audit, comparison_metrics, verify_native
from sloshing.pinned_wetting.mission import sha256, write_json
from sloshing.pinned_wetting.navier import LABEL, NavierConfig, NavierDiagnostics, NavierFEM

ROOT = Path(__file__).resolve().parents[2]


def make_case(directory, *, slip=.5, t_end=.02):
    import hashlib
    directory.mkdir()
    c = NavierConfig(nx=8, nz=12, dt=.005, t_end=t_end, snapshot_dt=.01, integrator="sdirk2", slip_length_m=slip)
    source_hash, files = extension_run.numerical_identity(ROOT)
    physical = extension_run.physical_identity(c)
    identity = {"run_id": directory.name, "source_hash": source_hash, "source_files": files,
                "physical_config_hash": hashlib.sha256(json.dumps(physical, sort_keys=True).encode()).hexdigest(),
                "physical": physical, "physical_model_label": LABEL, "synthetic": False,
                "source_kind": "declared_Navier_linear_FEM", "role": "unit_test_real_trajectory"}
    write_json(directory / "identity.json", identity)
    write_json(directory / "resolved_case.json", {"controls": asdict(c), "physical": physical,
        "physical_model_label": LABEL, "film_thickness_m": None,
        "retained_layer_volume_feedback": "neglected_at_leading_order",
        "extension_contract": {"record_height_tolerance_m": .0001}})
    for name in files:
        dest = directory / "source_snapshot" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "sloshing_visualization/src/sloshing" / name, dest)
    extension_run.execute_case(ROOT, directory, c, identity, {}, lambda *args: None)
    return directory / "native.h5"


def test_independent_positive_native_and_self_comparison(tmp_path):
    native = make_case(tmp_path / "valid")
    report = verify_native(native)
    assert report["passed"], report
    metrics = comparison_metrics(native, native)
    assert all(v is None or v == 0.0 for v in metrics.values())


@pytest.mark.parametrize("corruption, expected", [
    ("omitted_wall_loss", "limit:split_energy_relative"),
    ("fixed_endpoint", "fixed_vertical_endpoint_or_wrong_free_DOFs"),
    ("no_slip_relabel", "not_declared_Navier_source"),
    ("wrong_momentum", "limit:stage_momentum_relative"),
    ("coating_hole", "dry_hole_or_wrong_W"),
    ("memory_decrease", "memory_without_accepted_contact"),
    ("memory_unreached_increase", "memory_without_accepted_contact"),
    ("moving_P2", "marker_not_native_resolved_record"),
    ("curved_initial_surface", "initial_not_straight"),
    ("extrapolated_contact", "contact_not_actual_free_endpoint_trace"),
    ("synthetic_source", "not_declared_Navier_source"),
    ("missing_accepted_final", "native_time_cadence_or_horizon"),
])
def test_corrupt_native_rejected(tmp_path, corruption, expected):
    native = make_case(tmp_path / corruption)
    with h5py.File(native, "r+") as h:
        if corruption == "omitted_wall_loss":
            group = h["accepted/4"]
            values = json.loads(group.attrs["integrator_scalars"])
            values["wall_dissipated_energy"] = 0.0
            group.attrs["integrator_scalars"] = json.dumps(values)
        elif corruption == "fixed_endpoint":
            values = h["mesh/free_velocity_dofs"][:]
            values[-1] = -1
            h["mesh/free_velocity_dofs"][:] = values
        elif corruption == "no_slip_relabel":
            value = json.loads(h.attrs["identity"])
            value["source_kind"] = "real_reference_solver"
            h.attrs["identity"] = json.dumps(value)
        elif corruption == "wrong_momentum":
            values = h["accepted/4/stage2_pressure"][:]
            values[0] += 1.0
            h["accepted/4/stage2_pressure"][:] = values
        elif corruption == "coating_hole":
            values = h["accepted/4/wet_intervals"][:]
            values[0, 0] += 1.0
            h["accepted/4/wet_intervals"][:] = values
        elif corruption.startswith("memory_"):
            values = h["accepted/4/H"][:]
            values[1] += -.001 if corruption == "memory_decrease" else .001
            h["accepted/4/H"][:] = values
        elif corruption == "moving_P2":
            group = h["accepted/4"]
            coating = json.loads(group.attrs["coating"])
            coating["markers"][1]["height_m"] += .001
            group.attrs["coating"] = json.dumps(coating)
        elif corruption == "curved_initial_surface":
            values = h["accepted/0/eta"][:]
            values[len(values) // 2] += .001
            h["accepted/0/eta"][:] = values
        elif corruption == "extrapolated_contact":
            h.attrs["contact_method"] = "near_wall_linear_extrapolation"
        elif corruption == "synthetic_source":
            value = json.loads(h.attrs["identity"])
            value["synthetic"] = True
            h.attrs["identity"] = json.dumps(value)
        elif corruption == "missing_accepted_final":
            del h["accepted/4"]
    report = verify_native(native)
    assert not report["passed"]
    assert expected in report["failures"], report


@pytest.fixture(scope="module")
def tiny_record_native(tmp_path_factory):
    """One genuine inexpensive trajectory; no changes to any mission run."""
    return make_case(tmp_path_factory.mktemp("record_source") / "tiny_record", t_end=1.0)


@pytest.mark.parametrize("corruption,expected", [
    ("disappearing_P3", "marker_not_native_resolved_record"),
    ("frame_sampled_peak", "memory_without_accepted_contact"),
])
def test_real_contact_peak_corruption_rejected(tmp_path, tiny_record_native, corruption, expected):
    directory = tmp_path / corruption
    shutil.copytree(tiny_record_native.parent, directory)
    native = directory / "native.h5"
    with h5py.File(native, "r+") as h:
        groups = [h["accepted"][name] for name in sorted(h["accepted"], key=int)]
        final = json.loads(groups[-1].attrs["coating"])
        assert any(marker["marker_id"] == "P3" for marker in final["markers"])
        if corruption == "disappearing_P3":
            final["markers"] = [marker for marker in final["markers"] if marker["marker_id"] != "P3"]
            groups[-1].attrs["coating"] = json.dumps(final)
        else:
            # Pretend coarse output frames are sufficient to evolve H. Even
            # self-consistent W/checkpoint edits must fail the accepted-step law.
            sampled_H = groups[0]["H"][:]
            missed_peak = 0.0
            for i, group in enumerate(groups):
                if i % 20 == 0:
                    sampled_H = np.maximum(sampled_H, group["eta"][:][[0, -1]])
                missed_peak = max(missed_peak, float(np.max(group["H"][:] - sampled_H)))
                group["H"][:] = sampled_H
                intervals = group["wet_intervals"][:]
                intervals[:, 1] = sampled_H
                group["wet_intervals"][:] = intervals
                coating = json.loads(group.attrs["coating"])
                coating["H"] = sampled_H.tolist()
                group.attrs["coating"] = json.dumps(coating)
            assert missed_peak > 1e-8
    report = verify_native(native)
    assert not report["passed"]
    assert expected in report["failures"], report


def test_wrong_wall_friction_sign_cannot_self_validate(tmp_path, monkeypatch):
    class NegativeFriction(NavierFEM):
        def __init__(self, config):
            super().__init__(config)
            self.K_wall = -self.K_wall
            self.K = self.K_bulk + self.K_wall
    # Deliberately disable production guards to manufacture a physically wrong
    # trajectory; the independent verifier still uses positive Navier friction.
    monkeypatch.setattr(extension_run, "NavierFEM", NegativeFriction)
    monkeypatch.setattr(NavierDiagnostics, "validate", lambda *args: None)
    native = make_case(tmp_path / "wrong_friction_sign")
    report = verify_native(native)
    assert not report["passed"]
    assert "limit:stage_momentum_relative" in report["failures"], report
    assert "limit:split_energy_relative" in report["failures"], report


def test_mixed_slip_is_not_numerical_refinement(tmp_path):
    first = make_case(tmp_path / "b_half", slip=.5)
    second = make_case(tmp_path / "b_one", slip=1.0)
    with pytest.raises(ValueError, match="mixed_physics_or_slip"):
        comparison_metrics(first, second)
    assert comparison_metrics(first, second, same_physics=False)["normalized_macro_Linf"] > 0.0


def test_missing_independent_review_is_never_complete(tmp_path):
    assert _review_signoff(tmp_path, "SYNTHETIC", "digest", set()) == {
        "passed": False, "status": "PENDING", "failures": ["independent_review_missing"], "evidence_paths": []}


def test_review_binds_exact_artifacts_and_native_hash(tmp_path):
    """SYNTHETIC signature fixture tests the gate, not physical acceptance."""
    import xml.etree.ElementTree as ET
    suite = ET.Element("testsuite")
    names = ["test_wrong_wall_friction_sign_cannot_self_validate", "test_mixed_slip_is_not_numerical_refinement",
             "test_corrupt_native_rejected"] + [f"synthetic_{i}" for i in range(5)]
    for name in names:
        ET.SubElement(suite, "testcase", name=name)
    evidence = tmp_path / "tests.xml"
    ET.ElementTree(suite).write(evidence)
    payload = tmp_path / "payload.csv"
    payload.write_text("SYNTHETIC payload\n")
    write_json(tmp_path / "verification.json", {"status": "SYNTHETIC_PRECHECK"})
    review = {"status": "PASS_DECLARED_EXTENSION", "run_id": "SYNTHETIC", "native_sha256": "digest",
              "visual_native_reviewed": True, "artifact_hashes": {"payload.csv": sha256(payload)},
              "test_evidence": {"path": "tests.xml", "sha256": sha256(evidence),
                                "tests": 8, "failures": 0, "errors": 0, "skipped": 0}}
    write_json(tmp_path / "independent_review.json", review)
    assert _review_signoff(tmp_path, "SYNTHETIC", "digest", {"payload.csv"})["passed"]
    assert not _review_signoff(tmp_path, "SYNTHETIC", "changed_native", {"payload.csv"})["passed"]
    payload.write_text("SYNTHETIC changed after visual review\n")
    result = _review_signoff(tmp_path, "SYNTHETIC", "digest", {"payload.csv"})
    assert not result["passed"]
    assert "independent_review_artifact_hash:payload.csv" in result["failures"]


@pytest.mark.parametrize("mutation,expected", [
    ("none", None),
    ("unsigned", "unsigned_native_audit"),
    ("tampered_audit", "signed_native_audit_hash_mismatch"),
    ("tampered_native", "signed_native_data_hash_mismatch"),
    ("wrong_report_scope", "signed_native_report_identity_or_scope"),
])
def test_signed_native_audit_reuse_is_narrow(tmp_path, mutation, expected):
    """Real tiny native/operator audit; SYNTHETIC package/signature gate fixture."""
    native = make_case(tmp_path / "real_tiny")
    checked = verify_native(native)
    assert checked["passed"], checked
    output = tmp_path / "SYNTHETIC_HANDOFF_FIXTURE"
    output.mkdir()
    source = ROOT / "sloshing_visualization/src/sloshing/pinned_wetting/extension_verify.py"
    snapshot = output / "auditor_source.py"
    shutil.copyfile(source, snapshot)
    audit = {"scope": "DECLARED_EXTENSION_FINAL_ACCEPTANCE", "scientific_checks_passed": True,
             "native_verification_mode": "FULL_OPERATOR_RECOMPUTATION",
             "native_reports": {checked["run_id"]: checked},
             "auditor_source_sha256": sha256(snapshot), "auditor_source_snapshot_path": str(snapshot)}
    audit_path = output / "native_audit.json"
    write_json(audit_path, audit)
    manifest = {"run_id": checked["run_id"], "native_data_paths": [str(native)],
                "native_runs": [{"path": str(native), "run_id": checked["run_id"], "sha256": sha256(native)}]}
    review = {"status": "PASS_DECLARED_EXTENSION", "visual_native_reviewed": True,
              "run_id": checked["run_id"], "native_sha256": checked["native_sha256"],
              "native_audit": {"path": "native_audit.json", "sha256": sha256(audit_path),
                               "scope": "FULL_STAGE_OPERATOR_NATIVE_AUDIT", "auditor_source_sha256": sha256(snapshot)}}
    if mutation == "tampered_audit":
        audit["changed_after_signature"] = True
        write_json(audit_path, audit)
    elif mutation == "tampered_native":
        with h5py.File(native, "r+") as h:
            h["accepted/4/H"][0] += .001
    elif mutation == "wrong_report_scope":
        audit["native_reports"][checked["run_id"]]["scope"] = "SOLVER_SELF_DIAGNOSTIC"
        write_json(audit_path, audit)
        review["native_audit"]["sha256"] = sha256(audit_path)
    if mutation != "unsigned":
        write_json(output / "independent_review.json", review)
    if expected is None:
        saved, evidence = _signed_native_audit(output, manifest)
        assert saved[checked["run_id"]] == checked
        assert evidence["all_native_hashes_rechecked"] is True
    else:
        with pytest.raises(ValueError, match=expected):
            _signed_native_audit(output, manifest)
