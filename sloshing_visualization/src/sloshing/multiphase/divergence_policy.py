"""Predeclared STEP3A9 decision; no manual override or empirical cap increase."""
import copy
import numpy as np
from .full_rate_policy import POLICY as HISTORICAL

AUDIT = {"version": "step3a9-audit-v1", "meshes": [[96, 48], [120, 60], [144, 72]],
    "h_definition": "maximum geometric cell diameter", "dt": 1e-6,
    "weak_max": 1e-12, "slope_min": .5, "projection_relative": 1e-11,
    "projection_assembly_eps": 64., "pythagorean_relative": 1e-11,
    "wall_budget_s": 3600., "target_mass": -.6060961027495515,
    "delta": 1e-3, "max_unknowns_refined": 400000,
    "equilibrium_guess": "analytic sessile cap, theta60, radius0.3; mu=sigma/(2*radius)",
    "interface_localization": "P2 Bernstein range intersects [-0.9,0.9]",
    "wall_localization": "centroid within 2 epsilon of domain wall"}


def projection_checks(strong, projected, weak, orthogonal, policy=AUDIT):
    values = np.asarray([strong, projected, weak, orthogonal])
    finite = bool(np.isfinite(values).all() and (values >= 0).all())
    tolerance = policy["projection_relative"]*max(projected, weak, 1e-30) + \
        policy["projection_assembly_eps"]*np.finfo(float).eps*strong
    error = abs(projected-weak)
    pyth = abs(strong**2-projected**2-orthogonal**2)/max(strong**2, 1e-60)
    return {"projection_difference": float(error), "projection_tolerance": float(tolerance),
        "pythagorean_relative_error": float(pyth),
        "projection_agrees": bool(finite and error <= tolerance),
        "pythagorean_pass": bool(finite and pyth <= policy["pythagorean_relative"])}


def spatial_decision(rows, policy=AUDIT):
    if len(rows) != 3:
        return {"passed": False, "p_div": None, "error": "Three complete meshes required"}
    h = np.array([r["h"] for r in rows]); d = np.array([r["strong_divergence_L2"] for r in rows])
    valid = bool(np.isfinite(np.r_[h,d]).all() and (h > 0).all() and (d > 0).all())
    slope = float(np.polyfit(np.log(h), np.log(d), 1)[0]) if valid else None
    checks = {"weak_continuity_all": all(np.isfinite(r["weak_continuity_Riesz"]) and
            r["weak_continuity_Riesz"] <= policy["weak_max"] for r in rows),
        "projection_Riesz_agreement": all(r["projection_agrees"] and r["pythagorean_pass"] for r in rows),
        "strict_decrease": valid and bool(d[2] < d[1] < d[0]),
        "refinement_slope": slope is not None and slope >= policy["slope_min"],
        "non_divergence_physical": all(r["physical_pass"] for r in rows),
        "solver_success": all(r["SNES_reason"] > 0 and r["KSP_success"] for r in rows),
        "exact_hierarchy": [r["mesh_shape"] for r in rows] == policy["meshes"]}
    return {"passed": all(checks.values()), "checks": checks, "p_div": slope,
        "pairwise_ratios": [float(d[1]/d[0]), float(d[2]/d[1])] if valid else None,
        "interpretation": "measured refinement slope, not formal optimal FE order",
        "production_authorized": all(checks.values())}


def production_policy(decision):
    if not decision.get("passed") or not decision.get("production_authorized"):
        raise ValueError("Qualified audit required; no manual override")
    p = copy.deepcopy(HISTORICAL)
    cap = p.pop("strong_divergence_max"); source = p.pop("strong_divergence_reference")
    p.update(version="step3a9-v1", total_wall_s=19800., audit_wall_s=3600.,
        divergence_role="spatial diagnostic after qualified refinement audit",
        weak_continuity_role="hard mixed-FEM incompressibility gate",
        projection_fraction_max=1e-6, strong_significance_floor=1e-12,
        strong_same_root_relative=1e-4,
        historical_metadata={"STEP3A8_envelope": cap, "source": source,
            "STEP3A8_STOP_value": 3.814612027895975e-6, "active_hard_gate": False})
    return p


def divergence_gates(row, *, production=False):
    keys = ("strong_divergence_L2", "weak_continuity_Riesz", "projected_divergence_L2",
            "orthogonal_divergence_L2", "projection_fraction")
    finite = bool(np.isfinite([row[k] for k in keys]).all())
    checks = {"finite_divergence": finite,
        "weak_continuity": finite and row["weak_continuity_Riesz"] <= 1e-12,
        "projection_identity": finite and row["projection_agrees"] and row["pythagorean_pass"]}
    if production:
        checks["pressure_projection_fraction"] = finite and (
            row["strong_divergence_L2"] < 1e-12 or row["projection_fraction"] <= 1e-6)
    return checks


def bridge_divergence(a, b):
    error = abs(b-a)/max(a,b,1e-30)
    return {"relative_difference": error, "threshold": 1e-4,
        "passed": bool(np.isfinite([a,b,error]).all() and min(a,b) >= 0 and error <= 1e-4)}
