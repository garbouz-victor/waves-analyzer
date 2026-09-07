"""Physical-window settling; invariant to numerical step size for the same signal."""
from dataclasses import asdict
import numpy as np
from .validation_policy import SettlingPolicy


def window_series(times, values, window_s):
    times=np.asarray(times,dtype=float)
    values=np.asarray(values,dtype=float)
    if len(times)<2 or not np.isfinite(times).all() or np.any(np.diff(times)<=0) or not np.isfinite(values).all():
        raise ValueError("Settling requires finite values at strictly increasing times")
    start=times[-1]-window_s
    if start<times[0]-1e-12:
        return None
    selected=times>start+1e-12
    return np.r_[start,times[selected]],np.r_[np.interp(start,times,values),values[selected]]


def physical_window_rate(times,values,window_s):
    window=window_series(times,values,window_s)
    if window is None:
        return None
    t,y=window
    t=t-t[0];dt=np.diff(t)
    # Exact integrals of the piecewise-linear signal, not equal sample weights.
    integral_y=np.sum(dt*(y[:-1]+y[1:])/2)
    integral_ty=np.sum(dt*((2*t[:-1]+t[1:])*y[:-1]+(t[:-1]+2*t[1:])*y[1:])/6)
    return float((integral_ty-window_s*integral_y/2)/(window_s**3/12))


def validate_settling(observations,history,initial_energy_scale,policy=SettlingPolicy()):
    times=np.array([r["time"] for r in observations])
    theta=np.array([r["fits"][1]["theta_deg"] for r in observations])
    crossings=[sorted(p["coordinate"] for p in r["crossings"]) for r in observations]
    if any(len(row)!=2 for row in crossings):
        return {"qualified":False,"qualification_status":"failed","reason":"not exactly two crossings",
                "policy":asdict(policy)}
    crossings=np.array(crossings)
    rate=physical_window_rate(times,theta,policy.window_s)
    if rate is None:
        return {"qualified":False,"qualification_status":"not_qualified","reason":"history shorter than physical window",
                "policy":asdict(policy)}
    left=physical_window_rate(times,crossings[:,0],policy.window_s)
    right=physical_window_rate(times,crossings[:,1],policy.window_s)
    ht=np.array([r["time"] for r in history])
    metrics={}
    missing=[]
    for key in ("E_kin","speed_max_dof_sample","relative_mu_variation"):
        if not all(key in r for r in history):
            missing.append(key);metrics[key]=None
            continue
        window=window_series(ht,[r[key] for r in history],policy.window_s)
        metrics[key]=float(np.max(window[1])) if window is not None else None
    kinetic=None if metrics["E_kin"] is None else metrics["E_kin"]/max(initial_energy_scale,1e-30)
    checks={"angle_rate":abs(rate)<=policy.max_angle_rate_deg_per_s,
            "contact_speed":max(abs(left),abs(right))<=policy.max_contact_line_speed_m_per_s,
            "kinetic_fraction":kinetic is not None and kinetic<=policy.max_kinetic_energy_fraction,
            "speed":metrics["speed_max_dof_sample"] is not None and metrics["speed_max_dof_sample"]<=policy.max_speed_m_per_s,
            "chemical_equilibrium":metrics["relative_mu_variation"] is not None and metrics["relative_mu_variation"]<=policy.max_relative_mu_variation}
    qualified=bool(all(checks.values()))
    return {"qualified":qualified,"qualification_status":"passed" if qualified else "failed",
            "policy":asdict(policy),"checks":{k:bool(v) for k,v in checks.items()},
            "window_start_s":float(times[-1]-policy.window_s),"window_end_s":float(times[-1]),
            "angle_rate_deg_per_s":rate,"contact_left_speed_m_per_s":left,
            "contact_right_speed_m_per_s":right,"kinetic_fraction_peak":kinetic,
            "speed_peak_m_per_s":metrics["speed_max_dof_sample"],
            "relative_mu_variation_peak":metrics["relative_mu_variation"],"missing_metrics":missing}
