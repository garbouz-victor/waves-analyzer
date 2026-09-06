"""Predeclared engineering qualification policy, separate from solver tolerances.

The numeric field-comparison limits below operationalize 'small/coherent';
they are not error theorems and must be reported with the measured values.
Changing amplitude cannot improve relative spatial or divergence errors.
"""

POLICY = {
    "slope_conservative": .1, "slope_failure": .3,
    "eta_bulk_relative_limit": .01, "velocity_relative_limit": .02,
    "omega_wall_relative_limit": .1,
    "strong_divergence_particle_caution": .05,
    "significant_fraction": .1,
    "particle_bulk_separation_limit_m": .005,
    "particle_snapshot_separation_limit_m": .00001,
    "particle_rk_step_limit_m": .000001,
    "symmetry_absolute_limit": 1e-12,
}


def fine_short_gate(comparison):
    problems = []
    for key, limit in (("eta_bulk_relative_max", POLICY["eta_bulk_relative_limit"]),
                       ("velocity_relative_max", POLICY["velocity_relative_limit"]),
                       ("omega_wall_left_relative_max", POLICY["omega_wall_relative_limit"]),
                       ("omega_wall_right_relative_max", POLICY["omega_wall_relative_limit"])):
        if comparison[key] > limit:
            problems.append(f"{key}={comparison[key]:.6g} exceeds predeclared {limit:g}; requires analysis")
    for key in ("divergence_l2", "divergence_bulk_l2"):
        if comparison[key+"_fine_max"] >= comparison[key+"_medium_max"]:
            problems.append(key+" does not decrease on fine")
    return {"numerical_gate_passed": not problems, "problems": problems, "policy": POLICY}


def qualification_decision(medium, fine, comparison, particles):
    failures, conditions = [], []
    config = medium["parameters"]
    runs = (medium, fine)
    no_slip = max(r[f"{wall}_{c}_max_max"] for r in runs for wall in ("left", "right", "bottom") for c in ("u", "w"))
    volume = max(r["volume_error_max"] for r in runs)
    contact = max(r["contact_error_max"] for r in runs)
    budget = max(r["energy_balance_relative_max"] for r in runs)
    for value, limit, name in ((no_slip, config["wall_tolerance"], "no-slip"), (volume, config["volume_tolerance"], "volume"),
                               (contact, config["contact_tolerance"], "contact points"), (budget, config["energy_relative_tolerance"], "energy budget")):
        if value > limit:
            failures.append(name+" tolerance violated")
    for r in runs:
        if r["max_step_energy_increase_max"] > config["energy_relative_tolerance"]*r["total_energy_max"]:
            failures.append("Unexplained physical energy growth")
        if r["weak_divergence_l2_max"] > config["weak_divergence_tolerance"]:
            failures.append("Weak divergence tolerance violated")
    slope_run = max(runs, key=lambda r: r["slope_global_max"])
    slope = slope_run["slope_global_max"]
    if slope > POLICY["slope_failure"]:
        failures.append("Local small-slope linearization exceeds 0.3 (model qualification failure, not solver instability)")
    elif slope > POLICY["slope_conservative"]:
        conditions.append("Contact slope exceeds conservative 0.1 policy")
    failures.extend(fine_short_gate(comparison)["problems"])
    if medium["relative_div_above_5_percent_active_fraction"] > POLICY["significant_fraction"]:
        conditions.append("Medium relative strong divergence exceeds 5% for significant active time; particle interpretation requires caution")
    if max(comparison["velocity_symmetry_medium_max"], comparison["velocity_symmetry_fine_max"],
           medium["eta_antisymmetry_max"], fine["eta_antisymmetry_max"],
           medium["omega_wall_symmetry_max"], fine["omega_wall_symmetry_max"]) > POLICY["symmetry_absolute_limit"]:
        failures.append("Unexplained left/right asymmetry")
    if particles["mesh"]["bulk_separation_max_m"] > POLICY["particle_bulk_separation_limit_m"]:
        failures.append("Bulk particle mesh sensitivity exceeds 0.005 m")
    if particles["snapshot"]["bulk_separation_max_m"] > POLICY["particle_snapshot_separation_limit_m"]:
        failures.append("Bulk particle snapshot interpolation error exceeds 10 micrometres")
    if particles["rk_step"]["bulk_separation_max_m"] > POLICY["particle_rk_step_limit_m"]:
        failures.append("Bulk particle ODE-step sensitivity exceeds 1 micrometre")
    if particles["wall_velocity_max_m_per_s"] > config["wall_tolerance"]:
        failures.append("Particle FEM evaluation creates a nonzero wall velocity")
    if any(particles[k]["bulk_invalid_count"] for k in ("mesh", "snapshot", "rk_step")):
        failures.append("A bulk trajectory becomes invalid")
    if any(particles[k]["near_wall_invalid_count"] for k in ("mesh", "snapshot", "rk_step")):
        conditions.append("A near-wall trajectory becomes invalid")
    status = "failed" if failures else "conditional" if conditions else "qualified"
    phrases = {"failed": "NOT QUALIFIED FOR PHYSICAL ANIMATION",
               "conditional": "CONDITIONALLY QUALIFIED: BULK FLOW YES, CONTACT REGION NOT RESOLVED",
               "qualified": "QUALIFIED FOR BULK PHYSICAL ANIMATION"}
    return {"status": status, "decision": phrases[status], "parameters": config, "policy": POLICY,
            "medium_run_complete": medium["complete"], "fine_run_complete": fine["complete"],
            "max_slope_global": slope, "max_slope_time": slope_run["slope_global_time_at_max"],
            "max_slope_mesh": slope_run["parameters"]["mesh"],
            "max_slope_bulk": max(r["slope_bulk_max"] for r in runs),
            "max_slope_intermediate": max(r["slope_intermediate_max"] for r in runs),
            "max_slope_contact": max(r["slope_contact_max"] for r in runs),
            "max_relative_strong_divergence_medium": medium["divergence_relative_max"],
            "max_relative_strong_divergence_fine": fine["divergence_relative_max"],
            "medium_fine_eta_bulk_relative_max": comparison["eta_bulk_relative_max"],
            "medium_fine_velocity_relative_max": comparison["velocity_relative_max"],
            "medium_fine_omega_wall_relative_max": max(comparison["omega_wall_left_relative_max"], comparison["omega_wall_right_relative_max"]),
            "particle_bulk_separation_max_m": particles["mesh"]["bulk_separation_max_m"],
            "particle_near_wall_separation_max_m": particles["mesh"]["near_wall_separation_max_m"],
            "particle_snapshot_sampling_error_max_m": max(particles["snapshot"]["bulk_separation_max_m"], particles["snapshot"]["near_wall_separation_max_m"]),
            "no_slip_max": no_slip, "contact_error_max": contact, "volume_error_max": volume, "energy_balance_relative_max": budget,
            "extrapolated_alpha_for_slope_0_1_deg": slope_run["extrapolated_alpha_for_slope_0_1_deg"],
            "extrapolated_alpha_for_slope_0_05_deg": slope_run["extrapolated_alpha_for_slope_0_05_deg"],
            "extrapolated_angles_have_not_been_run": True, "failures": failures, "conditions": conditions,
            "limitations": ["Pointwise corner omega_max and limiting slope are not spatially resolved",
                            "No moving contact line, wall film, wetting/dewetting or drainage physics",
                            "Taylor-Hood is not pointwise divergence-free",
                            "Particle paths advect the linear velocity field; tiny second-order net drift is not a validated nonlinear transport prediction",
                            "Only two production mesh levels in STEP 1.6; fine is a numerical reference, not exact",
                            "Particle snapshot refinement covers 0...0.5 s; RK invalidity time is stage-resolved",
                            "PDE temporal refinement on 0...5 s is not repeated; dt is inherited from STEP 1.5"]}
