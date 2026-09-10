"""Independent read-only verification of the declared finite-Navier-slip case.

The Navier solver/integrator and their PASS flags are deliberately not imported.
Operators are rebuilt from the historical FEM and scalar boundary mass; stage
equations, energy integrals and contact provenance are checked from native data.
"""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from skfem import asm

from ..config import SimulationConfig
from ..diagnostics import quadratic_edge_max, surface_slopes
from ..fem_spaces import FEMSystem, scalar_mass
from ..validation.contact import surface_value
from .mission import sha256
from .reference import load_native

LABEL = "DECLARED_EXTENSION — конечное скольжение Навье"
CONTACT = "actual_free_endpoint_trace_Navier"
GAMMA = 1.0 - 1.0 / np.sqrt(2.0)
METRIC_NAMES = ("normalized_surface_Linf", "normalized_macro_Linf", "normalized_H_Linf",
                "first_peak_time_fraction", "first_peak_height_normalized")


def _json(path):
    return json.loads(Path(path).read_text())


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _physical(c):
    return {"a_m": c["a"], "d_m": c["d"], "g_m_s2": c["g"], "nu_m2_s": c["nu"],
            "alpha_deg": c["alpha_deg"], "slip_length_m": c["slip_length_m"],
            "sigma_N_m": 0.0, "bottom": "no_slip", "sides": "impermeable_Navier",
            "geometry": "linearized_free_surface",
            "layer": "irreversible_lower_connected_zero_volume_one_way"}


def _operators(c):
    """Rebuild the law without using the production NavierFEM implementation."""
    b = float(c["slip_length_m"])
    if not np.isfinite(b) or b <= 0 or c["integrator"] != "sdirk2":
        raise ValueError("Finite positive physical slip and SDIRK2 required")
    f = FEMSystem(SimulationConfig(**{k: v for k, v in c.items() if k != "slip_length_m"}))
    wall = np.union1d(f.wall_dofs["left"], f.wall_dofs["right"])
    fixed = np.union1d(np.intersect1d(wall, f.component_dofs[0]), f.wall_dofs["bottom"])
    free = np.setdiff1d(np.arange(f.velocity.N), fixed)
    scalar_wall = sum(asm(scalar_mass, f.scalar.boundary(side, intorder=6))
                      for side in ("left", "right")).tocoo()
    vertical = f.component_dofs[1]
    kw_full = coo_matrix((c["nu"] / b * scalar_wall.data,
                         (vertical[scalar_wall.row], vertical[scalar_wall.col])),
                        shape=(f.velocity.N, f.velocity.N)).tocsr()
    op = {"M": f.M_full[free][:, free].tocsr(),
          "bulk": f.K_full[free][:, free].tocsr(),
          "wall": kw_full[free][:, free].tocsr(), "wall_full": kw_full,
          "B": f.B_full[:, free].tocsr(), "R": f.R_full[:, free].tocsr(),
          "free": free, "fixed": fixed, "fem": f}
    op["C"] = (op["R"].T @ f.S).tocsr()
    return op


def _energy(op, velocity, eta):
    f = op["fem"]
    return (.5 * float(velocity @ (op["M"] @ velocity)),
            .5 * f.config.g * float(eta @ (f.S @ eta)))


def _record_markers(run_id, times, eta, H, floor):
    """Derive the catalog at every accepted step from contact samples alone."""
    catalog = [{"marker_id": "P1", "side": "L", "height_m": float(eta[0, 0]),
                "created_at_s": 0.0, "state_id": run_id + ":0"},
               {"marker_id": "P2", "side": "R", "height_m": float(eta[0, -1]),
                "created_at_s": 0.0, "state_id": run_id + ":0"}]
    history = [list(catalog)]
    for i in range(1, len(times)):
        if i >= 2:
            for side, column in (("L", 0), ("R", -1)):
                peak = float(eta[i - 1, column])
                prior = [m for m in catalog if m["side"] == side]
                if (peak > eta[i - 2, column] and peak > eta[i, column]
                        and peak > max(m["height_m"] for m in prior) + floor
                        and peak >= H[i, 0 if side == "L" else 1] - 1e-13):
                    catalog.append({"marker_id": ("P" if side == "L" else "R") + str(len(prior) + 2),
                                    "side": side, "height_m": peak, "created_at_s": float(times[i]),
                                    "state_id": f"{run_id}:{i - 1}"})
        history.append(list(catalog))
    return history


def verify_native(native):
    """Verify one complete native trajectory, including pilots; never writes files."""
    path = Path(native)
    errors = []
    def require(ok, reason):
        if not ok and reason not in errors:
            errors.append(reason)
    report = {"passed": False, "native": str(path), "scope": "DECLARED_NAVIER_NATIVE"}
    try:
        d = load_native(path)
        c, identity = d["config"], d["identity"]
        physical = _physical(c)
        require(identity.get("physical_model_label") == LABEL and identity.get("source_kind") == "declared_Navier_linear_FEM"
                and identity.get("synthetic") is False, "not_declared_Navier_source")
        require(d["contact_method"] == CONTACT, "contact_not_actual_free_endpoint_trace")
        require(identity.get("physical") == physical and identity.get("physical_config_hash") == _hash(physical),
                "physical_identity_mismatch")
        require(_json(path.parent / "identity.json") == identity, "native_identity_file_mismatch")
        resolved = _json(path.parent / "resolved_case.json")
        require(resolved["controls"] == c and resolved["physical"] == physical and
                resolved["physical_model_label"] == LABEL, "resolved_case_mismatch")
        require(resolved.get("film_thickness_m") is None and
                resolved.get("retained_layer_volume_feedback") == "neglected_at_leading_order", "layer_model_changed")
        source = identity["source_files"]
        require(identity.get("source_hash") == _hash(source), "source_aggregate_hash")
        for name, digest in source.items():
            rel = Path(name)
            require(not rel.is_absolute() and ".." not in rel.parts, "invalid_source_path")
            if not rel.is_absolute() and ".." not in rel.parts:
                require(sha256(path.parent / "source_snapshot" / rel) == digest, "source_snapshot_integrity")
        require("pinned_wetting/navier.py" in source and "pinned_wetting/extension_run.py" in source,
                "no_slip_source_relabelled")
        op = _operators(c)
        f, R, dt = op["fem"], op["R"], c["dt"]
        times, eta, H = d["times"], d["eta"], d["H"]
        n = len(times)
        require(n == round(c["t_end"] / dt) + 1 and
                np.allclose(times, np.arange(n) * dt, rtol=0, atol=1e-12), "native_time_cadence_or_horizon")
        require(n > 1 and times[0] == 0, "missing_initial_or_advanced_state")
        require(np.array_equal(d["x"], f.surface_x), "native_surface_mesh_mismatch")
        require(eta.shape == (n, len(f.surface_x)) and H.shape == (n, 2), "native_surface_shape")
        require(np.isfinite(eta).all() and np.isfinite(H).all(), "nonfinite_surface_or_memory")
        require(np.max(abs(eta[0] - f.config.slope * f.surface_x)) < 1e-12,
                "initial_not_straight")
        require(eta[0, -1] > eta[0, 0] and d["initial_velocity_max"] == 0, "initial_geometry_or_velocity")
        require(d["ids"] == [f"{identity['run_id']}:{i}" for i in range(n)] and all(d["accepted"]),
                "accepted_state_identity")
        require([R.getrow(i).nnz for i in (0, R.shape[0] - 1)] == [1, 1], "fixed_vertical_endpoint")
        expected_H = np.maximum.accumulate(eta[:, [0, -1]], axis=0)
        require(np.allclose(H, expected_H, rtol=0, atol=1e-12), "memory_without_accepted_contact")
        expected_intervals = np.stack([np.full_like(H, -c["d"]), H], axis=-1)
        require(np.allclose(d["intervals"], expected_intervals, rtol=0, atol=1e-12), "dry_hole_or_wrong_W")
        floor = resolved["extension_contract"]["record_height_tolerance_m"]
        require(floor == .0001, "undeclared_record_detector_floor")
        catalogs = _record_markers(identity["run_id"], times, eta, H, floor)
        for i, coating in enumerate(d["coatings"]):
            require(coating["run_id"] == identity["run_id"] and coating["accepted_step"] == i
                    and abs(coating["time"] - times[i]) < 1e-12 and coating["bottom"] == -c["d"],
                    "checkpoint_coating_identity")
            require(np.allclose(coating["H"], H[i], rtol=0, atol=1e-12), "checkpoint_coating_H")
            require(coating["markers"] == catalogs[i], "marker_not_native_resolved_record")
        maxima = {key: 0.0 for key in ("stage_momentum_relative", "weak_divergence", "mass_relative",
                  "endpoint_kinematic_absolute", "split_energy_relative", "energy_budget_relative",
                  "max_eta_over_a", "max_slope", "max_eta_over_b", "convective_indicator",
                  "kinematic_indicator", "wall_traction_L2")}
        bulk_loss = wall_loss = numerical = 0.0
        weights = (1.0 - GAMMA, GAMMA)
        surface_basis = f.velocity.boundary("surface", intorder=6)
        scalar_surface = f.scalar.boundary("surface", intorder=6)
        sides = [(f.velocity.boundary(side, intorder=6), sign) for side, sign in (("left", -1), ("right", 1))]
        with h5py.File(path, "r") as h:
            require(bool(h.attrs.get("completed", False)), "native_not_complete")
            require(h.attrs.get("physical_model_label") == LABEL, "native_header_model")
            require(h.attrs.get("wetting_sampling") == "every_accepted_step_endpoint", "unknown_contact_sampling")
            require(sorted(h["accepted"], key=int) == [str(i) for i in range(n)], "noncontiguous_native_groups")
            require(np.array_equal(h["mesh/free_velocity_dofs"][:], op["free"]), "fixed_vertical_endpoint_or_wrong_free_DOFs")
            require(np.array_equal(h["mesh/points"][:], f.mesh.p) and np.array_equal(h["mesh/triangles"][:], f.mesh.t),
                    "native_bulk_mesh_mismatch")
            if errors:
                report["failures"] = errors
                return report
            previous_v = h["accepted/0/velocity"][:]
            E0 = sum(_energy(op, previous_v, eta[0]))
            energy_scale = max(E0, 1e-30)
            for i in range(n):
                group = h["accepted"][str(i)]
                v = group["velocity"][:]
                require(v.shape == (len(op["free"]),) and np.isfinite(v).all(), "invalid_native_velocity")
                scalars = json.loads(group.attrs["integrator_scalars"])
                require(scalars["step"] == i and abs(scalars["time"] - times[i]) < 1e-12, "integrator_state_time")
                require(all(np.isfinite(value) for value in scalars.values()), "nonfinite_integrator_scalars")
                if i:
                    v1, p1, p2 = (group[name][:] for name in ("stage1_velocity", "stage1_pressure", "stage2_pressure"))
                    require(v1.shape == v.shape and p1.shape == p2.shape == (f.pressure.N,)
                            and all(np.isfinite(q).all() for q in (v1, p1, p2)), "invalid_stage_fields")
                    k1 = (v1 - previous_v) / (GAMMA * dt)
                    k2 = (v - previous_v - dt * weights[0] * k1) / (GAMMA * dt)
                    e1 = eta[i - 1] + GAMMA * dt * (R @ v1)
                    predicted = eta[i - 1] + dt * (weights[0] * (R @ v1) + weights[1] * (R @ v))
                    maxima["endpoint_kinematic_absolute"] = max(maxima["endpoint_kinematic_absolute"],
                                                              float(np.max(abs(predicted - eta[i]))))
                    for vv, pp, kk, ee in ((v1, p1, k1, e1), (v, p2, k2, eta[i])):
                        terms = (op["M"] @ kk, op["bulk"] @ vv, op["wall"] @ vv,
                                 c["g"] * (op["C"] @ ee), -(op["B"].T @ pp))
                        residual = np.linalg.norm(sum(terms)) / max(sum(np.linalg.norm(term) for term in terms), 1e-30)
                        maxima["stage_momentum_relative"] = max(maxima["stage_momentum_relative"], float(residual))
                        weak = op["B"] @ vv
                        maxima["weak_divergence"] = max(maxima["weak_divergence"],
                            float(np.sqrt(max(0.0, weak @ f.pressure_mass_lu.solve(weak)))))
                    bulk_loss += dt * (weights[0] * float(v1 @ (op["bulk"] @ v1)) + weights[1] * float(v @ (op["bulk"] @ v)))
                    wall_loss += dt * (weights[0] * float(v1 @ (op["wall"] @ v1)) + weights[1] * float(v @ (op["wall"] @ v)))
                    numerical += dt**2 * GAMMA**2 * (sum(_energy(op, k2, R @ v)) - sum(_energy(op, k1, R @ v1)))
                kinetic, potential = _energy(op, v, eta[i])
                energy = kinetic + potential
                remainder = energy - E0 + bulk_loss + wall_loss + numerical
                maxima["energy_budget_relative"] = max(maxima["energy_budget_relative"], abs(remainder) / max(E0, abs(energy - E0)))
                comparisons = {"bulk_dissipated_energy": bulk_loss, "wall_dissipated_energy": wall_loss,
                               "dissipated_energy": bulk_loss + wall_loss, "rk_energy_correction": numerical, "work_energy": 0.0}
                for key, expected in comparisons.items():
                    maxima["split_energy_relative"] = max(maxima["split_energy_relative"], abs(scalars[key] - expected) / energy_scale)
                diagnostics = d["diagnostics"][i]
                for key, expected in {"kinetic_energy": kinetic, "gravitational_potential_energy": potential,
                         "cumulative_bulk_dissipation": bulk_loss, "cumulative_wall_dissipation": wall_loss,
                         "cumulative_viscous_dissipation": bulk_loss + wall_loss,
                         "cumulative_work": 0.0, "rk_energy_correction": numerical,
                         "energy_balance_residual": remainder}.items():
                    require(abs(diagnostics[key] - expected) <= 1e-8 * energy_scale + 1e-14, "diagnostic_energy_not_native")
                maxima["mass_relative"] = max(maxima["mass_relative"], abs(float(f.surface_weights @ (eta[i] - eta[0]))) / (2*c["a"]*c["d"]))
                amplitude = quadratic_edge_max(eta[i, :-2:2], eta[i, 1::2], eta[i, 2::2])
                slopes = surface_slopes(f.surface_x, eta[i])
                maxima["max_eta_over_a"] = max(maxima["max_eta_over_a"], amplitude / c["a"])
                maxima["max_eta_over_b"] = max(maxima["max_eta_over_b"], amplitude / c["slip_length_m"])
                maxima["max_slope"] = max(maxima["max_slope"], max(float(np.max(abs(s))) for s in slopes))
                full = np.zeros(f.velocity.N)
                full[op["free"]] = v
                field = f.velocity.interpolate(full)
                adv = np.einsum("ij...,j...->i...", field.grad, np.asarray(field))
                max_convective = float(np.max(np.sqrt(np.sum(adv**2, axis=0)))) / (c["g"] * max(abs(f.config.slope), 1e-12))
                trace = surface_basis.interpolate(full)
                scalar = np.zeros(f.scalar.N)
                scalar[f.surface_dofs] = eta[i]
                graph = scalar_surface.interpolate(scalar)
                omitted = np.asarray(graph) * trace.grad[1, 1] - np.asarray(trace)[0] * graph.grad[0]
                max_kinematic = float(np.max(abs(omitted))) / (np.sqrt(c["g"] * c["a"]) * max(abs(f.config.slope), 1e-12))
                maxima["convective_indicator"] = max(maxima["convective_indicator"], max_convective)
                maxima["kinematic_indicator"] = max(maxima["kinematic_indicator"], max_kinematic)
                boundary_residual = 0.0
                for basis, normal in sides:
                    field = basis.interpolate(full)
                    defect = c["nu"] * (normal * (field.grad[1, 0] + field.grad[0, 1]) + np.asarray(field)[1] / c["slip_length_m"])
                    boundary_residual += float(np.sum(defect**2 * basis.dx))
                maxima["wall_traction_L2"] = max(maxima["wall_traction_L2"], float(np.sqrt(boundary_residual)))
                previous_v = v
        thresholds = {"stage_momentum_relative": 1e-9, "weak_divergence": 1e-9,
                      "endpoint_kinematic_absolute": 1e-11, "mass_relative": 1e-6,
                      "split_energy_relative": 1e-8, "energy_budget_relative": .05,
                      "max_eta_over_a": .05, "max_slope": .30}
        for key, bound in thresholds.items():
            require(maxima[key] <= bound, "limit:" + key)
        require(bulk_loss >= 0 and wall_loss > 0, "nonpositive_or_omitted_physical_losses")
        report.update(run_id=identity["run_id"], source_hash=identity["source_hash"],
                      physical_config_hash=identity["physical_config_hash"], physical=physical,
                      native_sha256=sha256(path), time_range=[float(times[0]), float(times[-1])],
                      accepted_states=n, maxima=maxima,
                      final_integrals_per_density={"bulk": bulk_loss, "wall": wall_loss, "signed_RK": numerical},
                      H_final_m=H[-1].tolist(), markers=catalogs[-1],
                      endpoint_trace_rows_nonzeros=[1, 1],
                      applicability_note="Amplitude/slope screens passed only if indicated; omitted-term indicators are not a rigorous physical error bound.")
    except (KeyError, IndexError, ValueError, TypeError, OSError, FloatingPointError) as error:
        errors.append("native_schema_or_read_error:" + str(error))
    report["passed"], report["failures"] = not errors, errors
    return report


def _first_peak(d):
    r = d["eta"][:, 0]
    ids = np.flatnonzero((r[1:-1] > r[:-2]) & (r[1:-1] > r[2:])) + 1
    running = np.maximum.accumulate(r)
    for value in ids:
        i = int(value)
        if r[i] > r[0] + .0001 and r[i] >= running[i + 1] - 1e-13:
            return float(d["times"][i]), float(r[i])
    return None


def comparison_metrics(first, second, *, same_physics=True):
    """Exact spatial P2 sup difference on common accepted physical times."""
    a, b = load_native(first), load_native(second)
    if same_physics and _physical(a["config"]) != _physical(b["config"]):
        raise ValueError("mixed_physics_or_slip_in_refinement")
    ta, tb = np.round(a["times"], 12), np.round(b["times"], 12)
    common, ia, ib = np.intersect1d(ta, tb, return_indices=True)
    if len(common) < 2 or common[0] != 0 or common[-1] != ta[-1] or common[-1] != tb[-1]:
        raise ValueError("refinement_common_horizon")
    edges = np.unique(np.r_[a["x"][::2], b["x"][::2]])
    query = np.column_stack([edges[:-1], .5*(edges[:-1]+edges[1:]), edges[1:]])
    normalization = abs(float(a["eta"][0,-1] - a["eta"][0,0]))
    maximum = 0.0
    for i, j in zip(ia, ib):
        error = surface_value(a["x"], a["eta"][i], query) - surface_value(b["x"], b["eta"][j], query)
        maximum = max(maximum, quadratic_edge_max(error[:, 0], error[:, 1], error[:, 2]))
    pa, pb = _first_peak(a), _first_peak(b)
    result = {"normalized_surface_Linf": maximum / normalization,
              "normalized_macro_Linf": float(np.max(abs(a["eta"][ia][:,[0,-1]] - b["eta"][ib][:,[0,-1]]))) / normalization,
              "normalized_H_Linf": float(np.max(abs(a["H"][ia] - b["H"][ib]))) / normalization,
              "first_peak_time_fraction": None if pa is None or pb is None else abs(pa[0] - pb[0]) / pb[0],
              "first_peak_height_normalized": None if pa is None or pb is None else abs(pa[1] - pb[1]) / normalization}
    if (pa is None) != (pb is None):
        raise ValueError("first_peak_identity_unresolved")
    return result


def _csv(path):
    with Path(path).open(newline="") as file:
        return list(csv.DictReader(file))


def _detected_records(d, indices, floor):
    times, eta = d["times"][indices], d["eta"][indices]
    H = np.maximum.accumulate(eta[:, [0, -1]], axis=0)
    catalog = _record_markers("DETECTOR_SENSITIVITY", times, eta, H, floor)[-1]
    return [{"side": m["side"], "height_m": m["height_m"],
             "peak_time_s": float(times[int(m["state_id"].rsplit(":", 1)[1])]),
             "confirmed_at_s": m["created_at_s"]} for m in catalog if m["marker_id"] not in ("P1", "P2")]


def _channel_reference(c):
    small = {**c, "nx": 8, "nz": 12}
    op = _operators(small)
    f = op["fem"]
    vertical = f.component_dofs[1]
    G = .001
    matrix = (f.K_full + op["wall_full"])[vertical][:, vertical]
    rhs = G * (f.M_full[vertical][:, vertical] @ np.ones(len(vertical)))
    velocity = spsolve(matrix, rhs)
    x = f.velocity.doflocs[0, vertical]
    exact = G / (2 * c["nu"]) * (c["a"]**2 - x*x + 2*c["a"]*c["slip_length_m"])
    return {"independent_channel_relative_error": float(np.max(abs(velocity - exact)) / np.max(abs(exact))),
            "independent_channel_max_absolute_error_m_s": float(np.max(abs(velocity - exact)))}


def _surface_export(rows, d, indices):
    count = len(d["x"])
    if len(rows) != len(indices) * count:
        return False
    for row_number, row in enumerate(rows):
        step, point = divmod(row_number, count)
        i = indices[step]
        if (int(row["point_index"]) != point or row["state_id"] != d["ids"][i] or row["component_id"] != "main"
                or abs(float(row["time_s"]) - d["times"][i]) > 1e-12
                or abs(float(row["x_m"]) - d["x"][point]) > 1e-12
                or abs(float(row["z_m"]) - d["eta"][i,point]) > 1e-12):
            return False
    return True


def _review_signoff(output, run_id, native_digest, required_hashes):
    """Check a separate human/agent-context attestation bound to immutable files."""
    errors = []
    path = output / "independent_review.json"
    if not path.is_file():
        return {"passed": False, "status": "PENDING", "failures": ["independent_review_missing"], "evidence_paths": []}
    try:
        review = _json(path)
        if (review.get("status") != "PASS_DECLARED_EXTENSION" or review.get("run_id") != run_id
                or review.get("native_sha256") != native_digest or review.get("visual_native_reviewed") is not True):
            errors.append("independent_review_scope_or_native_hash")
        signed = review.get("artifact_hashes", {})
        for name in required_hashes:
            if name not in signed or not (output / name).is_file() or sha256(output / name) != signed[name]:
                errors.append("independent_review_artifact_hash:" + name)
        evidence = review["test_evidence"]
        test_path = Path(evidence["path"])
        if not test_path.is_absolute():
            test_path = output / test_path
        if sha256(test_path) != evidence["sha256"]:
            errors.append("negative_test_evidence_hash")
        cases = ET.parse(test_path).getroot().findall(".//testcase")
        failed = sum(case.find("failure") is not None for case in cases)
        errored = sum(case.find("error") is not None for case in cases)
        skipped = sum(case.find("skipped") is not None for case in cases)
        if (len(cases) < 8 or failed or errored or skipped or evidence["tests"] != len(cases)
                or evidence["failures"] != failed or evidence["errors"] != errored or evidence["skipped"] != skipped):
            errors.append("negative_tests_not_all_passed")
        names = {case.attrib.get("name", "").split("[", 1)[0] for case in cases}
        if not {"test_wrong_wall_friction_sign_cannot_self_validate", "test_mixed_slip_is_not_numerical_refinement",
                "test_corrupt_native_rejected"} <= names:
            errors.append("missing_required_extension_negative_test_evidence")
        if not (output / "verification.json").is_file():
            errors.append("missing_verification_artifact")
    except (KeyError, TypeError, ValueError, OSError, ET.ParseError) as error:
        errors.append("independent_review_schema:" + str(error))
    return {"passed": not errors, "status": "PASS" if not errors else "FAIL",
            "failures": errors, "evidence_paths": [str(path)]}


def _report_sections(report, output):
    """OUTPUT_SCHEMA section summaries retain the scope of the actual checks."""
    def paths(names):
        return [str(output / name) for name in names if (output / name).is_file()]
    scientific = report.get("scientific_checks_passed", False)
    native = report.get("native_reports", {})
    native_paths = [value["native"] for value in native.values() if Path(value["native"]).is_file()]
    physics_status = ("PASS" if all(value["passed"] for value in native.values()) else "FAIL") if native else "NOT_RUN"
    report["contract"] = {"status": "PASS_DECLARED_EXTENSION" if scientific else "NOT_PASSED",
                          "evidence_paths": paths(["resolved_case.json", "model_decision.md", "ACCEPTANCE_ADDENDUM.yaml"])}
    report["physics"] = {"status": physics_status,
                         "evidence_paths": native_paths + paths(["boundary_benchmark.json"])}
    report["wetting"] = {"status": "PASS_NATIVE_CROSSCHECK" if scientific else "NOT_PASSED",
                         "evidence_paths": native_paths + paths(["wetting_history.csv", "marker_catalog.csv", "marker_tracks.csv", "detector_sensitivity.json", "record_stability.json"])}
    report["numerical_refinement"] = {"status": "PASS_RECOMPUTED" if scientific else "NOT_PASSED",
                                      "evidence_paths": paths(["refinement.json", "sensitivity.json"])}
    report["rendering"] = {"status": "PASS_FULL_DECODE_AND_STATE_MAP" if scientific else "NOT_PASSED",
                           "evidence_paths": paths(["final_animation.mp4", "wall_detail.mp4", "render_evidence.json", "key_frames.png"])}
    report.setdefault("artifact_integrity", {"status": "NOT_RUN", "evidence_paths": paths(["manifest.json"])})
    report.setdefault("independent_review", {"status": "NOT_RUN", "evidence_paths": []})
    return report


def verify_extension(output):
    """Verify a complete exported 0–5 s package against every native source."""
    output = Path(output)
    errors, natives = [], {}
    def require(ok, reason):
        if not ok and reason not in errors:
            errors.append(reason)
    report = {"passed": False, "scientific_checks_passed": False,
              "scope": "DECLARED_EXTENSION_FINAL_ACCEPTANCE", "output": str(output)}
    try:
        m = _json(output / "manifest.json")
        required = {"index.html", "final_animation.mp4", "wall_detail.mp4", "resolved_case.json", "model_decision.md",
                    "surface_history.csv", "surface_snapshots.csv", "wetting_history.csv", "marker_catalog.csv",
                    "marker_tracks.csv", "runup_summary.json", "energy_history.csv", "refinement.json", "sensitivity.json",
                    "key_frames.png", "FINAL_REPORT.md", "render_evidence.json", "ACCEPTANCE_ADDENDUM.yaml",
                    "boundary_benchmark.json", "case_summaries.json", "detector_sensitivity.json", "record_stability.json"}
        listed = {item["path"] for item in m["artifacts"]}
        require(required <= listed, "missing_mandatory_manifest_artifacts")
        for item in m["artifacts"]:
            relative = Path(item["path"])
            require(not relative.is_absolute() and ".." not in relative.parts, "invalid_artifact_path")
            if not relative.is_absolute() and ".." not in relative.parts:
                require((output / relative).is_file() and sha256(output / relative) == item["sha256"], "artifact_integrity:" + str(relative))
        report["artifact_integrity"] = {"status": "PASS" if not errors else "FAIL",
                                         "evidence_paths": [str(output / "manifest.json")]}
        if errors:
            report["failures"] = errors
            return _report_sections(report, output)
        for item in m["native_runs"]:
            path = Path(item["path"])
            require(sha256(path) == item["sha256"], "native_checksum:" + item["role"])
            checked = verify_native(path)
            require(checked["passed"], "native_failed:" + item["role"])
            require(checked.get("run_id") == item["run_id"], "native_run_id_mismatch")
            if item["run_id"] in natives:
                require(natives[item["run_id"]][0] == path, "duplicate_native_run_id")
            natives[item["run_id"]] = (path, checked)
        if any(reason.startswith("native_checksum:") for reason in errors):
            report["artifact_integrity"]["status"] = "FAIL"
        primary = Path(m["native_data_paths"][0])
        d = load_native(primary)
        identity, c = d["identity"], d["config"]
        require(identity["run_id"] in natives and natives[identity["run_id"]][0] == primary, "primary_not_verified_native")
        require(all(m.get(k) == identity[k] for k in ("run_id", "source_hash", "physical_config_hash", "physical_model_label")),
                "manifest_native_identity")
        require(m.get("schema_version") == "PW1_NAVIER_1.0" and m.get("resolved_case") == "resolved_case.json"
                and m.get("execution_HEAD") == identity.get("execution_HEAD")
                and isinstance(m.get("execution_HEAD"), str) and bool(m["execution_HEAD"]), "manifest_provenance_schema")
        require(m.get("frame_map") == "render_evidence.json" and
                m.get("contact_source_map") == "wetting_history.csv: state_id maps to native accepted group, reach=actual endpoint R",
                "manifest_contact_or_frame_source_map")
        require(m.get("units") == {"length": "m (video mm explicitly)", "time": "physical seconds",
                                   "energy_csv": "J/m", "native_energy": "per density"}
                and isinstance(m.get("limitations"), list) and len(m["limitations"]) >= 5,
                "manifest_units_or_limitations")
        require(m.get("physical_model_label") == LABEL and identity["physical"]["slip_length_m"] == .5,
                "wrong_selected_declared_model")
        require(all(c[k] == value for k, value in
                    {"a": 1.0, "d": 10.0, "g": 9.81, "nu": .01, "alpha_deg": 2.0}.items()),
                "unapproved_change_to_frozen_case")
        require(m["selected_role"] == identity["role"], "selected_run_role_mismatch")
        require(m["time_range"] == [0, 5] and d["times"][0] == 0 and abs(d["times"][-1] - 5) < 1e-12,
                "missing_final_solver_time")
        require(_json(output / "resolved_case.json") == _json(primary.parent / "resolved_case.json"), "resolved_case_not_native")
        if errors:
            report.update(failures=errors, native_reports={k: v[1] for k,v in natives.items()})
            return _report_sections(report, output)
        require(_surface_export(_csv(output / "surface_history.csv"), d, list(range(len(d["times"])))), "surface_export_not_native")
        slices = [int(round(t / c["dt"])) for t in np.arange(0, 5.01, .5)]
        require(_surface_export(_csv(output / "surface_snapshots.csv"), d, slices), "missing_or_incorrect_11_snapshots")
        wet = _csv(output / "wetting_history.csv")
        require(len(wet) == len(d["times"]), "wetting_csv_length")
        for i, row in enumerate(wet):
            if i >= len(d["times"]):
                break
            values = [float(row[k]) for k in ("R_L_m", "R_R_m", "reach_L_m", "reach_R_m", "H_L_m", "H_R_m")]
            expected = np.r_[d["eta"][i,[0,-1]], d["eta"][i,[0,-1]], d["H"][i]]
            require(row["state_id"] == d["ids"][i] and abs(float(row["time_s"]) - d["times"][i]) < 1e-12
                    and np.allclose(values, expected, atol=1e-12, rtol=0), "wetting_not_native_contact_or_memory")
        markers = [{**row, "height_m": float(row["height_m"]), "created_at_s": float(row["created_at_s"])}
                   for row in _csv(output / "marker_catalog.csv")]
        require(markers == d["coatings"][-1]["markers"], "marker_catalog_not_native")
        energy = _csv(output / "energy_history.csv")
        require(len(energy) == len(d["times"]), "energy_csv_length")
        rho = _json(primary.parent / "resolved_case.json")["original_case"]["physics"]["density_kg_m3"]
        require(rho == 1000.0, "unapproved_energy_density_scale")
        columns = {"kinetic_J_per_m": "kinetic_energy", "potential_J_per_m": "gravitational_potential_energy",
                   "physical_dissipation_integral_J_per_m": "cumulative_viscous_dissipation",
                   "external_work_J_per_m": "cumulative_work", "numerical_energy_remainder_J_per_m": "rk_energy_correction",
                   "budget_defect_J_per_m": "energy_balance_residual",
                   "bulk_viscous_dissipation_integral_J_per_m": "cumulative_bulk_dissipation",
                   "wall_friction_dissipation_integral_J_per_m": "cumulative_wall_dissipation"}
        for i, row in enumerate(energy):
            if i >= len(d["times"]):
                break
            require(row["state_id"] == d["ids"][i] and abs(float(row["time_s"]) - d["times"][i]) < 1e-12, "energy_state_map")
            for key, native_key in columns.items():
                require(np.isclose(float(row[key]), rho*d["diagnostics"][i][native_key], atol=1e-10, rtol=1e-10),
                        "energy_export_not_native:" + key)
            require(float(row["capillary_J_per_m"]) == float(row["wall_J_per_m"]) == 0.0, "undeclared_capillary_or_adsorption_energy")
        refinement = _json(output / "refinement.json")
        computed = {}
        for direction in ("time", "space"):
            study = refinement[direction]
            ids = study["run_ids"]
            require(len(ids) == 2 and all(i in natives for i in ids), "unknown_refinement_runs:" + direction)
            a, b = [load_native(natives[i][0]) for i in ids]
            require(a["identity"]["physical"] == b["identity"]["physical"] == identity["physical"], "mixed_slip_or_physics_in_refinement")
            require(study["physical_hashes"] == [a["identity"]["physical_config_hash"], b["identity"]["physical_config_hash"]], "refinement_physical_hashes")
            ac, bc = a["config"], b["config"]
            mesh_keys = ("nx", "nz", "x_grading", "z_grading", "a", "d")
            if direction == "time":
                require(all(ac[k] == bc[k] for k in mesh_keys) and ac["dt"] != bc["dt"], "time_refinement_not_independent")
            else:
                require(ac["dt"] == bc["dt"] and (ac["nx"], ac["nz"]) != (bc["nx"], bc["nz"]), "space_refinement_not_independent")
            metrics = comparison_metrics(natives[ids[0]][0], natives[ids[1]][0])
            computed[direction] = metrics
            for key, value in metrics.items():
                saved = study["metrics"][key]
                require((value is None and saved is None) or
                        (value is not None and saved is not None and np.isclose(value, saved, rtol=2e-6, atol=1e-9)),
                        "refinement_metric_not_native:" + direction + ":" + key)
                require(value is None or value <= .05, "refinement_threshold:" + direction + ":" + key)
        require(refinement.get("same_physics") is True and refinement.get("passed") is True, "refinement_summary_inconsistent")
        detector = _json(output / "detector_sensitivity.json")
        require(detector["run_id"] == identity["run_id"], "detector_sensitivity_run_id")
        all_indices = np.arange(len(d["times"]))
        require([v["record_floor_m"] for v in detector["floor_variants"]] == [0.0, .00005, .0001, .0002], "missing_detector_floor_variants")
        for variant in detector["floor_variants"]:
            require(variant["records"] == _detected_records(d, all_indices, variant["record_floor_m"]), "detector_floor_records_not_native")
        require([v["stride"] for v in detector["cadence_variants"]] == [1, 2, 4], "missing_detector_cadence_variants")
        for variant in detector["cadence_variants"]:
            selected = all_indices[::variant["stride"]]
            sampled = np.maximum.accumulate(d["eta"][selected][:, [0, -1]], axis=0)
            missed = float(np.max(d["H"][selected] - sampled))
            require(abs(variant["sample_dt_s"] - c["dt"] * variant["stride"]) < 1e-12
                    and abs(variant["max_missed_H_m"] - missed) < 1e-12
                    and variant["records"] == _detected_records(d, selected, .0001), "detector_cadence_not_native")
        stability = _json(output / "record_stability.json")
        require(stability["same_physics"] is True, "record_stability_mixed_physics")
        left_records, gains = [], []
        for role in ("baseline", "time_refined", "space_refined"):
            item = stability["runs"][role]
            require(item["run_id"] in natives, "record_stability_unknown_run")
            nd = load_native(natives[item["run_id"]][0])
            require(nd["identity"]["physical"] == identity["physical"], "record_stability_mixed_physics")
            records = [{**marker, "peak_time_s": float(nd["times"][int(marker["state_id"].rsplit(":", 1)[1])])}
                       for marker in nd["coatings"][-1]["markers"] if marker["side"] == "L" and marker["marker_id"] != "P1"]
            gain = records[1]["height_m"] - records[0]["height_m"] if len(records) >= 2 else None
            require(item["left_records"] == records and item["second_record_gain_m"] == gain, "record_stability_peaks_not_native")
            left_records.append(records)
            gains.append(gain)
        all_second = all(value is not None for value in gains)
        require(stability["second_record_present_all_resolutions"] == all_second, "record_stability_presence_claim")
        require(stability["second_record_gain_range_m"] == ([min(gains), max(gains)] if all_second else None), "record_stability_gain_range")
        if all_second:
            delta = [abs(left_records[0][k]["height_m"] - left_records[1][k]["height_m"])
                     + abs(left_records[1][k]["height_m"] - left_records[2][k]["height_m"]) for k in (0, 1)]
            empirical = stability["empirical_resolution_margin"]
            resolved_gain = gains[-1] > sum(delta) + .0001
            require(np.allclose(empirical["peak_delta_m"], delta, rtol=0, atol=1e-12)
                    and abs(empirical["sum_peak_deltas_m"] - sum(delta)) < 1e-12
                    and empirical["selected_gain_m"] == gains[-1] and empirical["record_floor_m"] == .0001
                    and empirical["gain_exceeds_empirical_margin_plus_floor"] == resolved_gain, "record_gain_margin_not_native")
            require(resolved_gain, "P4_not_resolved_above_empirical_margin")
        else:
            require(not any(marker["marker_id"] == "P4" for marker in markers), "P4_identity_unresolved_across_refinement")
        channel = _channel_reference(c)
        benchmark = _json(output / "boundary_benchmark.json")
        require(channel["independent_channel_relative_error"] < 1e-9
                and benchmark["relative_error"] < 1e-9
                and benchmark["benchmark"] == "steady_vertical_channel_with_Navier_sides_only", "independent_channel_reference_failed")
        # The two physical b variants must be native runs, separate from mesh/time pairs.
        sensitivity = _json(output / "sensitivity.json")
        variants = [v for v in natives.values() if v[1].get("physical", {}).get("slip_length_m") in (.25, 1.0)]
        require({v[1]["physical"]["slip_length_m"] for v in variants} == {.25, 1.0}, "missing_slip_sensitivity_runs")
        require(sensitivity["distinct_physics"] is True and sensitivity["base_slip_length_m"] == .5
                and sensitivity["passed"] is True, "sensitivity_interpretation")
        comparisons = {item["run_id"]: item for item in sensitivity["comparisons"]}
        for path, checked in variants:
            vc = load_native(path)["config"]
            require(checked["run_id"] in comparisons, "sensitivity_report_missing_native_run")
            require(all(vc[k] == c[k] for k in ("nx", "nz", "dt", "x_grading", "z_grading")), "sensitivity_resolution_mismatch")
            other = dict(checked["physical"])
            other["slip_length_m"] = c["slip_length_m"]
            require(other == identity["physical"], "undeclared_sensitivity_physics")
            item = comparisons[checked["run_id"]]
            require(item["base_run_id"] == identity["run_id"] and item["base_slip_length_m"] == .5
                    and item["slip_length_m"] == vc["slip_length_m"] and item["same_mesh_and_dt"] is True,
                    "sensitivity_source_metadata")
            actual = comparison_metrics(primary, path, same_physics=False)
            for key, value in actual.items():
                saved = item["metrics"][key]
                require((value is None and saved is None) or
                        (value is not None and saved is not None and np.isclose(value, saved, rtol=2e-6, atol=1e-9)),
                        "sensitivity_metric_not_native:" + key)
        runup = _json(output / "runup_summary.json")
        require(runup["run_id"] == identity["run_id"] and np.allclose(runup["H_final_m"], d["H"][-1], rtol=0, atol=1e-12)
                and np.allclose(runup["R_final_m"], d["eta"][-1,[0,-1]], rtol=0, atol=1e-12), "runup_summary_not_native")
        expected_peaks = [{**m, "peak_time_s": float(d["times"][int(m["state_id"].rsplit(":", 1)[1])])}
                          for m in markers if m["marker_id"] not in ("P1", "P2")]
        require(runup["records"] == expected_peaks, "runup_records_not_native")
        render = _json(output / "render_evidence.json")
        require(render["model"] == LABEL and render["time_range_s"] == [0, 5] and render["frame_dt_s"] == .01, "render_model_or_time")
        require(render["film_label"] == "Остаточный слой показан условно; толщина не рассчитывалась", "missing_schematic_film_label")
        frame_times = np.arange(501) * .01
        indices = np.rint(frame_times / c["dt"]).astype(int)
        expected_tracks = [(float(t), m["marker_id"], m["side"], m["height_m"])
                           for t in frame_times for m in markers if m["created_at_s"] <= t + 1e-12]
        tracks = [(float(r["time_s"]), r["marker_id"], r["side"], float(r["height_m"])) for r in _csv(output / "marker_tracks.csv")]
        require(tracks == expected_tracks, "marker_tracks_not_native")
        require({v["video"] for v in render["videos"]} == {"final_animation.mp4", "wall_detail.mp4"}, "missing_final_videos")
        decodes = []
        for video in render["videos"]:
            frames = video["frame_map"]
            require(len(frames) == 501, "frame_count_not_5s")
            for j, frame in enumerate(frames):
                if j >= 501:
                    break
                i = indices[j]
                expected_markers = [m for m in markers if m["created_at_s"] <= frame_times[j] + 1e-12]
                require(frame["frame"] == j and abs(frame["time_s"] - frame_times[j]) < 1e-12 and frame["state_id"] == d["ids"][i], "frame_state_time_mismatch")
                require(np.allclose(frame["H_m"], d["H"][i], atol=1e-12, rtol=0)
                        and np.allclose(frame["wet_intervals"], d["intervals"][i], atol=1e-12, rtol=0), "render_coating_not_native")
                require(frame["markers"] == expected_markers, "moving_or_disappearing_marker")
                require(frame["axes"] == frames[0]["axes"], "moving_render_axes")
            probe = json.loads(subprocess.check_output(["ffprobe", "-v", "error", "-count_frames", "-select_streams", "v:0",
                "-show_entries", "stream=width,height,avg_frame_rate,nb_read_frames", "-of", "json", str(output / video["video"])], text=True))["streams"][0]
            subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-i", str(output / video["video"]), "-f", "null", "-"],
                           check=True, capture_output=True, timeout=120)
            require(probe["avg_frame_rate"] == "25/1" and int(probe["nb_read_frames"]) == 501
                    and (probe["width"], probe["height"]) == (1920, 1080), "video_cadence_or_resolution")
            decodes.append({"video": video["video"], "full_decode": True, "probe": probe})
        report.update(run_id=identity["run_id"], physical_model_label=LABEL, computed_refinement=computed,
                      channel_reference=channel, record_gain_m=gains[-1],
                      native_reports={k: v[1] for k,v in natives.items()}, video_decodes=decodes,
                      independent_review={"status": "REQUIRED_SEPARATE_VISUAL_SIGNOFF", "evidence_paths": []})
    except (KeyError, IndexError, ValueError, TypeError, OSError, FloatingPointError, subprocess.SubprocessError) as error:
        errors.append("output_schema_or_read_error:" + str(error))
    report["scientific_checks_passed"] = not errors
    if not errors:
        bound_artifacts = required - {"index.html", "FINAL_REPORT.md"}
        review = _review_signoff(output, identity["run_id"], sha256(primary), bound_artifacts)
        report["independent_review"] = review
        errors.extend(review["failures"])
    report["passed"], report["failures"] = not errors, errors
    return _report_sections(report, output)
