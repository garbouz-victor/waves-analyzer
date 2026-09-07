"""Closed-system energy qualification from measured powers, not energy fitting.

Every interval (including the first) is integrated independently of E. NaN/inf
raises a named error: it must never disappear in min/max or pass comparisons.
`qualified` checks a single history; temporal evidence is a separate mandatory
study gate, recorded explicitly as not_measured if absent.
"""
import csv
from dataclasses import asdict
from pathlib import Path
import numpy as np
from .validation_policy import EnergyPolicy

COMPONENTS = ("viscous", "CH", "slip")
ENERGIES = ("E_kin", "E_gravity", "E_interface", "E_wall")


def read_history(path):
    with Path(path).open() as stream:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]


def positive_energy_scale(row):
    return sum(abs(row[key]) for key in ENERGIES)


def validate_energy(history, policy=EnergyPolicy(), temporal_convergence_info=None):
    if len(history)<2:
        raise ValueError("Energy validation requires an initial sample and an accepted interval")
    required=("time", "E_total", "mass_error_relative_to_domain", *ENERGIES,
              *(k+"_dissipation" for k in COMPONENTS))
    for index,row in enumerate(history):
        for key in required:
            if key not in row:
                raise ValueError(f"Missing energy metric {key} at sample {index}")
        for key,value in row.items():
            if (key.startswith("E_") or "energy" in key or "dissipation" in key
                    or key in required) and not np.isfinite(value):
                raise ValueError(f"Nonfinite energy/dissipation metric {key} at sample {index}")
    times=np.array([r["time"] for r in history])
    dt=np.diff(times)
    if np.any(dt<=0):
        raise ValueError("Energy history times must be strictly increasing")
    energy=np.array([r["E_total"] for r in history])
    scales=np.array([positive_energy_scale(r) for r in history])
    scale=max(float(scales[0]),policy.scale_floor)
    legacy_scale=max(abs(float(energy[0])),policy.scale_floor)
    components={}
    total_intervals=np.zeros(len(dt))
    cumulative_consistent=True
    nonnegative=True
    warnings=[]
    for name in COMPONENTS:
        power=np.array([r[name+"_dissipation"] for r in history])
        nonnegative &= bool(np.min(power)>=-policy.negative_power_tolerance)
        intervals=dt*(power[:-1]+power[1:])/2
        cumulative=np.r_[0.,np.cumsum(intervals)]
        total_intervals+=intervals
        saved_key="cumulative_"+name+"_dissipation"
        if any(saved_key in r for r in history):
            if not all(saved_key in r for r in history):
                raise ValueError(f"Incomplete saved {saved_key}")
            saved=np.array([r[saved_key] for r in history])
            cumulative_consistent &= bool(np.allclose(saved,cumulative,rtol=1e-10,atol=1e-12))
        fraction=float(np.max(intervals)/cumulative[-1]) if cumulative[-1]>0 else 0.
        if fraction>policy.interval_fraction_warning:
            warnings.append(f"{name}: >50% of its time integral in one interval; temporally suspicious")
        # Simpson is diagnostic only; never used to replace the primary budget.
        simpson=None
        if len(times)>=3 and np.allclose(dt,dt[0],rtol=1e-10,atol=1e-15):
            from scipy.integrate import simpson as integrate_simpson
            simpson=float(integrate_simpson(power,x=times))
        components[name]={"initial_power_W_per_m":float(power[0]),
            "first_interval_J_per_m":float(intervals[0]),"integral_J_per_m":float(cumulative[-1]),
            "largest_interval_fraction":fraction,"simpson_uniform_J_per_m":simpson,
            "trapezoid_minus_simpson_J_per_m":None if simpson is None else float(cumulative[-1]-simpson)}
    change=energy-energy[0]
    cumulative=np.r_[0.,np.cumsum(total_intervals)]
    defect=change+cumulative
    local=np.diff(energy)+total_intervals
    floor=max(policy.change_floor_absolute,policy.change_floor_fraction*scale)
    active=np.abs(change)>floor
    relative_change=np.abs(defect[active]/change[active])
    relative_to_change=float(relative_change.max()) if relative_change.size else None
    absolute=float(np.max(np.abs(defect)))
    legacy_relative=absolute/legacy_scale
    characteristic_relative=absolute/scale
    change_ok=relative_to_change is None or relative_to_change<=policy.change_relative_tolerance
    energy_components_consistent=bool(np.allclose(energy,
        np.array([sum(r[k] for k in ENERGIES) for r in history]),rtol=1e-10,atol=1e-12))
    # Keep the legacy |E0| metric, but never let cancellation of the arbitrary
    # total-energy reference create a false failure. The positive scale gates.
    budget_ok=(characteristic_relative<=policy.budget_relative_tolerance and change_ok
        and cumulative_consistent and energy_components_consistent and nonnegative)
    mass_max=max(abs(r["mass_error_relative_to_domain"]) for r in history)
    growth=max(0.,float(np.max(np.diff(energy))))
    mass_ok=mass_max<=policy.mass_relative_tolerance
    growth_ok=growth<=policy.growth_absolute_tolerance
    index=int(np.argmax(np.abs(local)))
    return {"finite_metrics":True,"mass_ok":bool(mass_ok),"energy_growth_ok":bool(growth_ok),
        "energy_budget_ok":bool(budget_ok),"qualified":bool(mass_ok and growth_ok and budget_ok),
        "scope":"single-history closure; does not establish temporal convergence",
        "policy":asdict(policy),"nonnegative_dissipation":bool(nonnegative),
        "saved_cumulative_consistent":bool(cumulative_consistent),
        "energy_components_consistent":energy_components_consistent,
        "mass_error_relative_max":float(mass_max),"max_positive_energy_step_J_per_m":growth,
        "absolute_budget_defect":absolute,"final_budget_defect_J_per_m":float(defect[-1]),
        "max_local_budget_defect_J_per_m":float(abs(local[index])),
        "max_local_budget_defect_time_s":float(times[index+1]),
        "energy_change_final_J_per_m":float(change[-1]),
        "legacy_relative_budget_defect":legacy_relative,
        "legacy_metric_gates":False,
        "absolute_budget_limit_from_initial_scale_J_per_m":policy.budget_relative_tolerance*scale,
        "relative_budget_defect_to_initial_scale":characteristic_relative,
        "relative_budget_defect_to_energy_change":relative_to_change,
        "energy_change_floor_J_per_m":floor,"energy_change_applicable_samples":int(active.sum()),
        "E_scale_initial_J_per_m":scale,"E_scale_peak_J_per_m":float(scales.max()),
        "dissipation_components":components,"sampling_warnings":warnings,
        "temporal_convergence_info":temporal_convergence_info or {"status":"not_measured"}}
