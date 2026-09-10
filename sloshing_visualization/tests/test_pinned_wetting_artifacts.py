"""Semantic negative tests on copies of a REAL bundle, with fresh checksums.

PW1_AUDIT_BUNDLE selects the existing reference. Source/native files stay read-only.
The absence of the explicit bundle is a visible skip, never an acceptance PASS.
"""
import csv
import json
import os
from pathlib import Path
import shutil

import pytest

from sloshing.pinned_wetting.audit import verify_output
from sloshing.pinned_wetting.mission import sha256, write_json


@pytest.fixture
def bundle(tmp_path):
    source = os.environ.get("PW1_AUDIT_BUNDLE")
    if not source:
        pytest.skip("PW1_AUDIT_BUNDLE not provided; real-artifact negatives NOT RUN")
    shutil.copytree(source, tmp_path / "bundle")
    return tmp_path / "bundle"


def alter_csv(path, mutate):
    with path.open(newline="") as file:
        reader = csv.DictReader(file)
        fields, rows = reader.fieldnames, list(reader)
    mutate(rows)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fields)
        writer.writeheader(); writer.writerows(rows)


def refresh_hashes(bundle):
    path = bundle / "manifest.json"
    manifest = json.loads(path.read_text())
    for item in manifest["artifacts"]:
        item["sha256"] = sha256(bundle / item["path"])
    write_json(path, manifest)


@pytest.mark.parametrize("case, expected", [
    ("track", "marker_tracks_not_native"),
    ("catalog_and_render", "marker_catalog_not_native"),
    ("surface_x", "surface_coordinate_or_state_map"),
    ("surface_state", "surface_coordinate_or_state_map"),
    ("refinement_run", "mixed_run_ids"),
    ("refinement_physics", "mixed_physics"),
    ("frame_time", "frame_time_vs_native"),
    ("frame_axes", "moving_render_axes"),
    ("source_relabel", "manifest_native_identity"),
    ("lost_track", "marker_tracks_not_native"),
])
def test_real_export_semantic_corruption(bundle, case, expected):
    if case in ("track", "lost_track"):
        def change(rows):
            if case == "track":
                rows[-1]["height_m"] = "0.10"
            else:
                rows.pop()
        alter_csv(bundle / "marker_tracks.csv", change)
    elif case == "catalog_and_render":
        alter_csv(bundle / "marker_catalog.csv", lambda rows: rows[1].update(height_m="0.10"))
        p = bundle / "render_evidence.json"
        val = json.loads(p.read_text())
        for f in val["frame_map"]:
            f["markers"][1]["height_m"] = .10
        write_json(p, val)
    elif case.startswith("surface"):
        key, val = ("x_m", "0.42") if case == "surface_x" else ("state_id", "other:0")
        alter_csv(bundle / "surface_history.csv", lambda rows: rows[-1].update({key: val}))
    elif case.startswith("refinement"):
        p = bundle / "refinement.json"
        val = json.loads(p.read_text())
        val["run_id" if case == "refinement_run" else "physical_config_hash"] = "WRONG"
        write_json(p, val)
    elif case.startswith("frame"):
        p = bundle / "render_evidence.json"
        val = json.loads(p.read_text())
        if case == "frame_time":
            val["frame_map"][-1]["time_s"] -= .001
        else:
            val["frame_map"][-1]["ylim_mm"][1] += 10
        write_json(p, val)
    elif case == "source_relabel":
        p = bundle / "manifest.json"
        val = json.loads(p.read_text())
        val["physical_model_label"] = "TARGET"
        val["source_kind"] = "qualified_target_solver"
        write_json(p, val)
    # Recompute checksums deliberately: these must fail semantic cross-checks,
    # not merely the generic tamper detector.
    refresh_hashes(bundle)
    result = verify_output(bundle, reference_only=True)
    assert not result["passed"]
    assert expected in result["failures"], result
