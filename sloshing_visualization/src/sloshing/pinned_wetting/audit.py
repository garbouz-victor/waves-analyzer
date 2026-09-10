"""Read-only cross-checks. A valid reference is explicitly NOT target acceptance."""
import csv
import json
import math
from pathlib import Path
import subprocess

import h5py
import numpy as np

from ..config import SimulationConfig
from ..fem_spaces import FEMSystem
from ..multiphase.free_energy import bulk_energy, bulk_derivative, lambda_from_sigma, wall_energy, wall_derivative
from ..solver import SloshingSolver
from .reference import LABEL, load_native


def blocker_check():
    """Two executable diagnostics, not two failed attempts at a target CFD run."""
    config = SimulationConfig(nx=8, nz=12, t_end=.02, dt=.0025, snapshot_dt=.01)
    solver = SloshingSolver(config)
    f = solver.fem
    s = solver.integrator.initial_state()
    initial = s.eta[[0, -1]].copy()
    for _ in range(config.nsteps):
        s = solver.integrator.advance(s)
    phi = np.linspace(-1, 1, 51)
    energy_terms = {"bulk_energy": bulk_energy(phi, 0., .01),
                    "bulk_derivative": bulk_derivative(phi, 0., .01),
                    "gradient_coefficient": lambda_from_sigma(0.) * .01,
                    "wall_energy": wall_energy(phi, 0., 60.),
                    "wall_derivative": wall_derivative(phi, 0., 60.)}
    terms = {k: float(np.max(abs(v))) for k, v in energy_terms.items()}
    from types import SimpleNamespace
    from ..multiphase.full_phase_rate import require_scope
    rejected_scopes = {}
    for name, rho_gas, g in (("gravity", 1000., 9.81), ("unequal_density", 1., 0.)):
        probe = SimpleNamespace(rho_liquid=1000., rho_gas=rho_gas, g=g, a_x=0.,
                                phase_degree=2, quadrature_degree=12)
        try:
            require_scope(probe)
        except ValueError as error:
            rejected_scopes[name] = str(error)
        else:
            raise AssertionError("Historical CHNS scope changed; reassess target suitability")
    evidence = {"status": "NEEDS_MODEL_DECISION", "blocker_code": "PW1_CONTACT_CLOSURE",
        "sharp_material_no_slip": {"endpoint_trace_nonzeros": [f.R.getrow(i).nnz for i in (0, f.R.shape[0]-1)],
                                   "steps": config.nsteps, "t_end_s": s.time,
                                   "endpoint_displacement_m": (s.eta[[0, -1]]-initial).tolist(),
                                   "interior_surface_change_max_m": float(np.max(abs(s.eta-f.initial_surface())))},
        "existing_diffuse_sigma_zero": {"energy_terms_max_abs": terms,
            "rejected_scope_probes": rejected_scopes,
            "interpretation": "Finite mobility/free energy/contact relaxation needs a separate closure; gravity diffusion alone is not a qualified target law."},
        "missing": "A contact transport/relaxation law and scale, not coating thickness.",
        "no_blanket_impossibility_claim": "No-slip diffuse contact is possible with a specified model.",
        "primary_sources": ["https://arxiv.org/abs/1310.1255", "https://arxiv.org/abs/1507.08945"]}
    if evidence["sharp_material_no_slip"]["endpoint_trace_nonzeros"] != [0, 0] or any(terms.values()):
        raise AssertionError("The diagnosed mathematical obstruction changed; review source")
    return evidence


def check_consistency(d, *, require_target=False, t_end=None):
    """Array audit kernel, shared with destructive-copy negative fixtures.

    d contains native source, export, and render state separately so tampering
    one cannot silently redefine its independent evidence.
    """
    errors = []
    def check(ok, reason):
        if not ok:
            errors.append(reason)
    t, eta, raw, reach, H = (np.asarray(d[k]) for k in ("times", "eta", "raw_R", "reach", "H"))
    check(all(np.isfinite(a).all() for a in (t, eta, raw, reach, H)), "nonfinite")
    check(t[0] == 0 and np.all(np.diff(t) > 0), "time_order")
    if t_end is not None:
        check(abs(t[-1] - t_end) < 1e-10, "missing_final_solver_time")
    check(np.max(abs(eta[0] - (d["mean"] + d["slope"]*d["x"]))) < 1e-12, "initial_not_straight")
    check(np.allclose(raw, d["native_R"], rtol=0, atol=1e-12), "contact_not_native")
    check(np.all(reach >= raw-1e-13), "interval_reach_below_endpoint")
    check(np.allclose(reach, d["native_reach"], rtol=0, atol=1e-12), "missed_detector_peak")
    check(np.allclose(H[0], raw[0], rtol=0, atol=1e-12), "bad_initial_memory")
    check(np.all(np.diff(H, axis=0) >= -1e-13), "decreasing_memory")
    expected = np.maximum.accumulate(np.vstack([raw[0], reach[1:]]), axis=0)
    check(np.allclose(H, expected, rtol=0, atol=1e-12), "memory_without_reach")
    check(np.allclose(H, d["native_H"], rtol=0, atol=1e-12), "memory_not_native")
    check(np.all(np.asarray(d["accepted"])), "rejected_state_in_history")
    check(d["source_run_id"] == d["manifest_run_id"], "mixed_run_ids")
    check(d["source_config_hash"] == d["manifest_config_hash"], "mixed_physics")
    if d["refinement_run_id"] is not None:
        check(d["source_run_id"] == d["refinement_run_id"], "mixed_run_ids")
    if d["refinement_config_hash"] is not None:
        check(d["source_config_hash"] == d["refinement_config_hash"], "mixed_physics")
    check(not d["synthetic"] or d["source_kind"] == "SYNTHETIC", "synthetic_claimed_scientific")
    intervals = np.asarray(d["intervals"])
    check(intervals.shape == (len(t), 2, 2) and
          np.allclose(intervals[:, :, 0], d["bottom"], rtol=0, atol=1e-12) and
          np.allclose(intervals[:, :, 1], H, rtol=0, atol=1e-12), "dry_hole_in_W")
    for fr in d["frames"]:
        i = fr["state_index"]
        check(np.allclose(fr["H"], H[i], rtol=0, atol=1e-12), "render_memory_mismatch")
        check(np.allclose(fr["intervals"], intervals[i], rtol=0, atol=1e-12), "render_coating_mismatch")
        expected_markers = {m["marker_id"]: m for m in d["markers"] if m["created_at_s"] <= fr["time_s"]+1e-12}
        actual = {m["marker_id"]: m for m in fr["markers"]}
        check(set(actual) == set(expected_markers), "marker_disappeared")
        for key in set(actual) & set(expected_markers):
            check(actual[key]["side"] == expected_markers[key]["side"] and
                  abs(actual[key]["height_m"] - expected_markers[key]["height_m"]) < 1e-12, "marker_moved")
    check(d["initial_velocity_max"] == 0., "initial_velocity_not_zero")
    check(d["contact_method"] != "fixed_endpoint_near_wall_extrapolation", "invalid_contact_method")
    if require_target:
        check(d["source_kind"] == "qualified_target_solver", "reference_not_target")
        check(d["contact_method"] != "actual_fixed_endpoint_trace", "target_contact_not_qualified")
    return sorted(set(errors))


def verify_output(output, reference_only=False):
    from .mission import sha256
    output = Path(output)
    if not (output / "manifest.json").exists():
        return {"passed": False, "status": "MISSING_MANIFEST", "path": str(output)}
    m = json.loads((output / "manifest.json").read_text())
    native = Path(m["native_data_paths"][0])
    failures = []
    for item in m["artifacts"]:
        path = output / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            failures.append("artifact_integrity:" + item["path"])
    if not native.is_file() or sha256(native) != m["native_sha256"]:
        failures.append("native_integrity")
        return {"passed": False, "failures": failures, "rendering": {"status": "NOT_RUN"}}
    d = load_native(native)
    c = d["config"]
    if any(m.get(key) != d["identity"].get(key) for key in
           ("source_hash", "physical_config_hash", "physical_model_label", "synthetic")):
        failures.append("manifest_native_identity")
    if m.get("contact_method") != d["contact_method"] or d["identity"].get("physical_model_label") != LABEL:
        failures.append("reference_source_label")
    resolved = json.loads((output / "resolved_case.json").read_text())
    if resolved != json.loads((native.parent / "resolved_case.json").read_text()):
        failures.append("resolved_case_not_native")
    import hashlib
    if hashlib.sha256(json.dumps(resolved, sort_keys=True).encode()).hexdigest() != m["physical_config_hash"]:
        failures.append("resolved_case_hash")
    for name, digest in d["identity"]["source_files"].items():
        if sha256(native.parent / "source_snapshot" / name) != digest:
            failures.append("source_snapshot_integrity")
    def rows(name):
        with (output / name).open(newline="") as file:
            return list(csv.DictReader(file))
    wet = rows("wetting_history.csv")
    surface = rows("surface_history.csv")
    render = json.loads((output / "render_evidence.json").read_text())
    refinement = json.loads((output / "refinement.json").read_text())
    # No fabricated refinement identities: NOT_RUN contains no claimed runs.
    if refinement.get("status") != "NOT_RUN" or refinement.get("runs"):
        failures.append("unverified_refinement_claim")
    n = len(d["times"])
    if len(wet) != n or len(surface) != n*len(d["x"]):
        return {"passed": False, "failures": ["truncated_csv"], "rendering": {"status": "NOT_RUN"}}
    raw = np.array([[float(w[k]) for k in ("R_L_m", "R_R_m")] for w in wet])
    reach = np.array([[float(w[k]) for k in ("reach_L_m", "reach_R_m")] for w in wet])
    H = np.array([[float(w[k]) for k in ("H_L_m", "H_R_m")] for w in wet])
    if [w["state_id"] for w in wet] != d["ids"] or not np.allclose([float(w["time_s"]) for w in wet], d["times"], atol=1e-12, rtol=0):
        failures.append("csv_state_map")
    se = np.array([float(r["z_m"]) for r in surface]).reshape(d["eta"].shape)
    if not np.allclose(se, d["eta"], atol=1e-12, rtol=0):
        failures.append("surface_not_native")
    for i, r in enumerate(surface):
        step, point = divmod(i, len(d["x"]))
        if (int(r["point_index"]) != point or abs(float(r["x_m"])-d["x"][point]) > 1e-12 or
                abs(float(r["time_s"])-d["times"][step]) > 1e-12 or
                r["state_id"] != d["ids"][step] or r["component_id"] != "main"):
            failures.append("surface_coordinate_or_state_map"); break
    id_to_index = {sid: i for i, sid in enumerate(d["ids"])}
    if any(fr["state_id"] not in id_to_index for fr in render["frame_map"]):
        return {"passed": False, "failures": ["unknown_render_state"], "rendering": {"status": "NOT_RUN"}}
    frames = [{"state_index": id_to_index[f["state_id"]], "time_s": f["time_s"], "H": f["H_m"],
               "intervals": f["wet_intervals"], "markers": f["markers"]} for f in render["frame_map"]]
    markers = [{**r, "height_m": float(r["height_m"]), "created_at_s": float(r["created_at_s"])} for r in rows("marker_catalog.csv")]
    if markers != d["coatings"][-1]["markers"]:
        failures.append("marker_catalog_not_native")
    expected_initial = [("P1", "L", float(d["eta"][0,0])), ("P2", "R", float(d["eta"][0,-1]))]
    for key, side, height in expected_initial:
        found = [mr for mr in markers if mr["marker_id"] == key]
        if (len(found) != 1 or found[0]["side"] != side or found[0]["height_m"] != height or
                found[0]["created_at_s"] != 0 or found[0]["state_id"] != d["ids"][0]):
            failures.append("initial_marker_geometry")
    expected_times = np.arange(round(c["t_end"]/.01)+1) * .01
    if len(frames) != len(expected_times):
        failures.append("frame_count_vs_physical_time")
    for number, fr in enumerate(render["frame_map"]):
        state_i = id_to_index[fr["state_id"]]
        if (fr["frame"] != number or abs(fr["time_s"]-d["times"][state_i]) > 1e-12 or
                number >= len(expected_times) or abs(fr["time_s"]-expected_times[number]) > 1e-12):
            failures.append("frame_time_vs_native")
        if (fr["xlim"] != render["frame_map"][0]["xlim"] or
                fr["ylim_mm"] != render["frame_map"][0]["ylim_mm"]):
            failures.append("moving_render_axes")
    tracks = rows("marker_tracks.csv")
    expected_tracks = [(fr["time_s"], mr["marker_id"], mr["side"], mr["height_m"])
                       for fr in frames for mr in markers if mr["created_at_s"] <= fr["time_s"]+1e-12]
    actual_tracks = [(float(tr["time_s"]), tr["marker_id"], tr["side"], float(tr["height_m"])) for tr in tracks]
    if actual_tracks != expected_tracks:
        failures.append("marker_tracks_not_native")
    checks = {"times": d["times"], "eta": se, "raw_R": raw, "reach": reach, "H": H,
        "x": d["x"], "slope": math.tan(math.radians(c["alpha_deg"])), "mean": 0., "bottom": -c["d"],
        "native_R": d["eta"][:, [0, -1]], "native_reach": d["eta"][:, [0, -1]], "native_H": d["H"],
        "intervals": d["intervals"], "accepted": d["accepted"], "frames": frames, "markers": markers,
        "source_run_id": d["identity"]["run_id"], "manifest_run_id": m["run_id"], "refinement_run_id": refinement.get("run_id"),
        "source_config_hash": d["identity"]["physical_config_hash"], "manifest_config_hash": m["physical_config_hash"],
        "refinement_config_hash": refinement.get("physical_config_hash"), "synthetic": m["synthetic"], "source_kind": m["source_kind"],
        "contact_method": d["contact_method"], "initial_velocity_max": d["initial_velocity_max"]}
    failures.extend(check_consistency(checks, require_target=not reference_only, t_end=c["t_end"] if reference_only else 5.))
    if failures:
        return {"passed": False, "scope": "REFERENCE_ONLY" if reference_only else "TARGET_ACCEPTANCE",
                "scientific_status": "TARGET_NOT_VALIDATED", "failures": sorted(set(failures)),
                "rendering": {"status": "NOT_RUN_DUE_TO_INCONSISTENT_DATA"}}
    # Rebuild numerical operators and evaluate native velocities; do not trust diagnostic PASS flags.
    fem = FEMSystem(SimulationConfig(**c))
    weak_max, mass_max, wall_max, energy_defect = 0., 0., 0., 0.
    E0 = .5*c["g"]*float(d["eta"][0] @ (fem.S @ d["eta"][0]))
    with h5py.File(native, "r") as h:
        for i, eta in enumerate(d["eta"]):
            gr = h["accepted"][str(i)]
            v = gr["velocity"][:]
            if not np.isfinite(v).all():
                failures.append("nonfinite_native_velocity"); break
            b = fem.B @ v
            weak_max = max(weak_max, float(np.sqrt(max(0., b @ fem.pressure_mass_lu.solve(b)))))
            mass_max = max(mass_max, abs(float(fem.surface_weights @ (eta-d["eta"][0]))) / (2*c["a"]*c["d"]))
            wall_max = max(wall_max, float(np.max(abs(fem.expand_velocity(v)[fem.fixed]))))
            s = json.loads(gr.attrs["integrator_scalars"])
            E = sum(fem.energy(v, eta))
            energy_defect = max(energy_defect, abs(E-E0+s["dissipated_energy"]+s["rk_energy_correction"]-s["work_energy"]) / max(E0, abs(E-E0)))
            if json.loads(gr.attrs["coating"])["H"] != gr["H"][:].tolist():
                failures.append("checkpoint_coating_mismatch")
    if weak_max > 1e-9 or mass_max > 1e-6 or wall_max > 1e-12 or energy_defect > .05:
        failures.append("reference_physical_budget")
    rendering = {"status": "NOT_RUN"}
    try:
        probe = json.loads(subprocess.check_output(["ffprobe", "-v", "error", "-count_frames", "-select_streams", "v:0", "-show_entries",
            "stream=width,height,avg_frame_rate,nb_read_frames", "-of", "json", str(output / "preview_reference.mp4")], text=True))["streams"][0]
        subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-i", str(output / "preview_reference.mp4"), "-f", "null", "-"], check=True, capture_output=True, timeout=120)
        if (int(probe["nb_read_frames"]) != len(frames) or probe["avg_frame_rate"] != "25/1" or
                (probe["width"], probe["height"]) != (1920, 1080)):
            failures.append("video_cadence")
        rendering = {"status": "PASS_REFERENCE_DECODE", "probe": probe, "evidence_paths": ["preview_reference.mp4", "render_evidence.json"]}
    except (subprocess.SubprocessError, FileNotFoundError) as error:
        failures.append("video_decode:" + str(error))
    if not reference_only:
        for name in ("final_animation.mp4", "wall_detail.mp4"):
            if not (output / name).is_file():
                failures.append("missing_final_artifact:" + name)
        failures.append("target_refinement_NOT_RUN")
    return {"passed": not failures, "scope": "REFERENCE_ONLY" if reference_only else "TARGET_ACCEPTANCE",
        "mission_status": "NEEDS_MODEL_DECISION", "scientific_status": "TARGET_NOT_VALIDATED",
        "failures": sorted(set(failures)), "contract": {"status": "REFERENCE_ONLY", "evidence_paths": [str(native)]},
        "physics": {"status": "PASS_REFERENCE" if not failures else "SEE_FAILURES", "weak_divergence_max": weak_max,
                    "mass_relative_max": mass_max, "wall_velocity_max": wall_max, "energy_budget_relative_max": energy_defect},
        "wetting": {"status": "CHECKED_AGAINST_NATIVE_REFERENCE", "evidence_paths": [str(native), "wetting_history.csv"]},
        "numerical_refinement": {"status": "NOT_RUN", "evidence_paths": ["refinement.json"]},
        "rendering": rendering, "artifact_integrity": {"status": "PASS" if not any("integrity" in x for x in failures) else "FAIL"},
        "independent_review": {"status": "PENDING", "evidence_paths": []}}
