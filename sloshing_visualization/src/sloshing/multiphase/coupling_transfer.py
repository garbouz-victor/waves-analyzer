"""Pure measured-coupling gates, separated from isolated temporal uncertainty."""
import numpy as np


def scalar_coupling(full,isolated,initial_D,excess,policy):
    if not np.isfinite([*full.values(),*isolated.values(),initial_D,excess]).all() or min(initial_D,excess)<=0:
        raise ValueError("Finite positive reference scales and scalar inputs required")
    D=isolated["D"]; df=full["D"]; hydro=full["visc"]+full["slip"]
    significant=D>=policy["power_floor_fraction"]*initial_D
    full_significant=df>=policy["power_floor_fraction"]*initial_D
    return {"F_initial_excess_relative":abs(full["F"]-isolated["F"])/excess,
        "CH_power_relative_significant":abs(df-D)/D if significant else None,
        "CH_power_difference_initial_scaled":abs(df-D)/initial_D,
        "Ekin_initial_excess_relative":full["kinetic"]/excess,
        "hydro_over_CH_significant":hydro/df if full_significant else None,
        "mass_full_minus_isolated":full["mass"]-isolated["mass"],
        "CH_significant":significant,"full_CH_significant":full_significant}


def field_coupling(coupling,temporal,direct,current_disturbance,A_phi,mu_gap,current_mu,A_mu,policy):
    if not np.isfinite([coupling,temporal,direct,current_disturbance,A_phi,mu_gap,current_mu,A_mu]).all() or min(A_phi,A_mu)<=0:
        raise ValueError("Valid initial disturbance scales required")
    current_ok=current_disturbance>=policy["current_field_floor_fraction"]*A_phi
    return {"phi_coupling_L2":coupling,"iso_phi_temporal_02_L2":temporal,
        "full0_vs_iso2_L2":direct,"phi_coupling_initial_relative":coupling/A_phi,
        "phi_coupling_current_relative":coupling/current_disturbance if current_ok else None,
        "current_disturbance_significant":current_ok,
        "mu_coupling_L2":mu_gap,"mu_coupling_relative_diagnostic":mu_gap/max(current_mu,policy["current_field_floor_fraction"]*A_mu),
        "combined_phi_proxy":coupling+temporal,"rigorous_continuum_bound":False}


def transfer_gates(scalars,fields,final,policy):
    if not scalars or not fields: raise ValueError("Full scalar and common field observations required")
    def maximum(key):
        values=[r[key] for r in scalars if r[key] is not None]
        return max(values) if values else None
    measured={"phi_initial":max(r["phi_coupling_initial_relative"] for r in fields),
        "F_initial":maximum("F_initial_excess_relative"),"D_significant":maximum("CH_power_relative_significant"),
        "kinetic":maximum("Ekin_initial_excess_relative"),"hydro_power":maximum("hydro_over_CH_significant"),
        "hydro_integral":final["I_hydro"]/final["excess"],
        "phi_coupling":max(r["phi_coupling_L2"] for r in fields),
        "phi_temporal":max(r["iso_phi_temporal_02_L2"] for r in fields),
        "D_coupling":abs(final["I_CH_full"]-final["I_iso0"]),
        "D_temporal":abs(final["I_iso2"]-final["I_iso0"]),
        "F_coupling":abs(final["F_CH_full"]-final["F_iso0"]),
        "F_temporal":abs(final["F_iso2"]-final["F_iso0"]),
        "budget_over_change":abs(final["budget"]/final["DeltaE"]) if final["DeltaE"] else None}
    limits={"phi_initial":policy["phi_relative"],"F_initial":policy["F_relative"],
        "D_significant":policy["D_relative"],"kinetic":policy["kinetic_relative"],
        "hydro_power":policy["hydro_relative"],"hydro_integral":policy["hydro_relative"],
        "phi_coupling":measured["phi_temporal"],"D_coupling":measured["D_temporal"],
        "F_coupling":max(measured["F_temporal"],policy["energy_precision_floor"]),
        "budget_over_change":policy["final_budget_relative"]}
    checks={k:measured[k] is not None and bool(np.isfinite(measured[k]) and measured[k]<=v) for k,v in limits.items()}
    return {"measured":measured,"limits":limits,"gates":checks,"passed":all(checks.values()),
        "combined_D_proxy":measured["D_coupling"]+measured["D_temporal"],
        "max_combined_phi_proxy":max(r["combined_phi_proxy"] for r in fields),
        "rigorous_full_continuum_bound":False}
