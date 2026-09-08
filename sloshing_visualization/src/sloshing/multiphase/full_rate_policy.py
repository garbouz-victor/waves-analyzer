"""STEP3A8 predeclared, validation-only controls and acceptance constants."""
import numpy as np
from .nonlinear_accuracy import NonlinearControls

CONTROLS={"P0":NonlinearControls(),"P1":NonlinearControls(1e-12,1e-12),
          "P2":NonlinearControls(3e-13,3e-13)}
L0_SHA="21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31"
POLICY={"version":"step3a8-v1","moderate_dt":1e-6,"tiny_dt":4.657657028534679e-12,
    "moderate_levels":["P1","P2"],"controls":{k:v.record() for k,v in CONTROLS.items()},
    "bridge_relative":1e-11,"FD_final_relative":1e-7,
    "moderate_phi_mu_D_relative":1e-10,"moderate_energy_mass_absolute":1e-12,
    "moderate_u_pi_scaled":1e-9,"BC_scaled":1e-12,"weak_continuity":1e-12,
    "strong_divergence_max":3.6880704671415807e-6,
    "strong_divergence_reference":"STEP3A2 stiff bump dt=2.5e-6, same mesh; no relative strong-divergence hard policy exists",
    "phi_relative":1e-3,"F_relative":1e-3,"D_relative":1e-2,
    "kinetic_relative":1e-4,"hydro_relative":1e-4,"power_floor_fraction":1e-8,
    "current_field_floor_fraction":1e-6,"energy_precision_floor":1e-13,
    "final_budget_relative":.05,"initial_excess":5.542888467657825e-8,
    "initial_D":.11400313905309394,"initial_cells":8.328131075834218,
    "restart_field_relative":1e-10,"restart_u_pi_scaled":1e-10,
    "restart_D_relative":1e-10,"restart_energy_absolute":1e-12,
    "preliminary_wall_s":1800.,"full_L0_wall_s":14400.,"total_wall_s":16200.,
    "cost_safety_factor":1.25,"first_guess":"semidiscrete CH at u=0",
    "next_guess":"previous accepted phase_rate, unscaled at dt changes",
    "optional_absolute_tiny":False,"full_CHNS_L1_L2_authorized":False,
    "isolated_rerun_authorized":False,"moving_contact_authorized":False}


def cost_forecast(timings,remaining,used_full,used_total,initialization=0.,policy=POLICY):
    values=np.asarray([r["timing"]["total_s"] for r in timings],dtype=float)
    if not len(values) or not np.isfinite(values).all() or np.min(values)<=0:
        raise ValueError("Measured full accepted-step timings required")
    mean=float(np.mean(values)); p95=float(np.percentile(values,95))
    rate=policy["cost_safety_factor"]*max(mean,p95)
    future=remaining*rate+initialization
    result={"mean_s_per_accepted":mean,"p95_s_per_accepted":p95,
        "conservative_s_per_step":rate,"remaining_steps":remaining,"remaining_wall_s":future,
        "projected_full_L0_wall_s":used_full+future,"projected_total_wall_s":used_total+future}
    result["authorized"]=(result["projected_full_L0_wall_s"]<=policy["full_L0_wall_s"] and
                           result["projected_total_wall_s"]<=policy["total_wall_s"])
    return result


def session_result(error):
    return "failed" if error is not None else "complete"
