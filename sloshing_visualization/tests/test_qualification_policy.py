from copy import deepcopy

from sloshing.config import SimulationConfig
from dataclasses import asdict
from sloshing.validation.qualification.policy import qualification_decision


def measured_fixture():
    # Explicit synthetic METRICS solely to test decision logic, never report data.
    m = {"parameters": asdict(SimulationConfig(alpha_deg=.02)), "complete": True,
         "volume_error_max": 1e-18, "contact_error_max": 0., "energy_balance_relative_max": 1e-13,
         "max_step_energy_increase_max": 0., "total_energy_max": 4e-7, "weak_divergence_l2_max": 1e-15,
         "slope_global_max": .08, "slope_global_time_at_max": 1., "slope_bulk_max": .001,
         "slope_intermediate_max": .01, "slope_contact_max": .08, "divergence_relative_max": .03,
         "relative_div_above_5_percent_active_fraction": 0., "eta_antisymmetry_max": 1e-16,
         "omega_wall_symmetry_max": 1e-16, "extrapolated_alpha_for_slope_0_1_deg": .025,
         "extrapolated_alpha_for_slope_0_05_deg": .0125}
    m.update({f"{wall}_{c}_max_max": 0. for wall in ("left", "right", "bottom") for c in ("u", "w")})
    f = deepcopy(m)
    f["parameters"]["mesh"] = "fine"
    comparison = {"eta_bulk_relative_max": .001, "velocity_relative_max": .002,
                  "omega_wall_left_relative_max": .01, "omega_wall_right_relative_max": .01,
                  "divergence_l2_medium_max": 1e-4, "divergence_l2_fine_max": 6e-5,
                  "divergence_bulk_l2_medium_max": 1e-6, "divergence_bulk_l2_fine_max": 6e-7,
                  "velocity_symmetry_medium_max": 1e-16, "velocity_symmetry_fine_max": 1e-16}
    particle = {"wall_velocity_max_m_per_s": 0., "rk_step": {"bulk_separation_max_m": 1e-12,
                                                          "bulk_invalid_count": 0, "near_wall_invalid_count": 0},
                "mesh": {"bulk_separation_max_m": 1e-6, "near_wall_separation_max_m": 2e-6,
                         "bulk_invalid_count": 0, "near_wall_invalid_count": 0},
                "snapshot": {"bulk_separation_max_m": 1e-9, "near_wall_separation_max_m": 2e-9,
                             "bulk_invalid_count": 0, "near_wall_invalid_count": 0}}
    return m, f, comparison, particle


def test_slope_policy_cannot_be_overruled_by_passing_invariants():
    args = measured_fixture()
    assert qualification_decision(*args)["status"] == "qualified"
    args[1]["slope_global_max"] = .2
    assert qualification_decision(*args)["status"] == "conditional"
    args[1]["slope_global_max"] = .31
    assert qualification_decision(*args)["status"] == "failed"


def test_divergence_and_particles_can_fail_otherwise_good_dataset():
    args = measured_fixture()
    args[2]["divergence_l2_fine_max"] = 2e-4
    assert qualification_decision(*args)["status"] == "failed"
    args = measured_fixture()
    args[3]["mesh"]["bulk_separation_max_m"] = .006
    assert qualification_decision(*args)["status"] == "failed"
    args = measured_fixture()
    args[3]["rk_step"]["bulk_invalid_count"] = 1
    assert qualification_decision(*args)["status"] == "failed"
