"""Static apparent-angle and actual moving-crossing measurements from CG roots."""
import json
from pathlib import Path
import numpy as np
from ..initialization import sessile_drop
from ..interface import contour_segments, connected_polylines
from ..contact_line import all_wall_crossings,sessile_apparent_angle
from ..diagnostics import interface_resolution
from .common import run_case


def observe(solver):
    segments=contour_segments(solver)
    components=connected_polylines(segments)
    crossings=all_wall_crossings(components,solver.config.z_min,axis=1)
    fits=[sessile_apparent_angle(segments,solver.config.z_min,solver.config.epsilon,window)
          for window in ((2.,6.),(2.,10.),(3.,8.))]
    return {"time":solver.time,"crossings":crossings,"fits":fits,
            "interface_resolution":interface_resolution(solver)}


def summarize(solver,result,output,initial_angle):
    folder=Path(output)
    observations=json.loads((folder/"observations.json").read_text())
    final=observations[-1]
    c=solver.config
    angle=final["fits"][1]["theta_deg"]
    errors=[abs(f["theta_deg"]-c.theta_equilibrium_deg) for f in final["fits"]]
    before=sorted(p["coordinate"] for p in observations[0]["crossings"])
    after=sorted(p["coordinate"] for p in final["crossings"])
    if len(before)!=2 or len(after)!=2:
        raise RuntimeError("Sessile benchmark needs two crossings; all crossings saved, no silent selection")
    half_width_change=((after[1]-after[0])-(before[1]-before[0]))/2
    expected_direction=np.sign(initial_angle-c.theta_equilibrium_deg)
    movement_correct=abs(initial_angle-c.theta_equilibrium_deg)<1e-6 or half_width_change*expected_direction>0
    theta_history=[r["fits"][1]["theta_deg"] for r in observations]
    last_change=abs(theta_history[-1]-theta_history[-2])
    checks={"apparent_angle_within_3_degrees":max(errors)<3.,
        "fit_window_spread_below_1_degree":np.ptp([f["theta_deg"] for f in final["fits"]])<1.,
        "angle_settling_below_0_05_degree_per_step":last_change<.05,
        "correct_motion_direction":bool(movement_correct),
        "interface_resolved_at_all_samples":all(r["interface_resolution"]["qualified_resolution"] for r in observations),
        "mass_conservation":result["max_mass_error_relative"]<1e-8,
        "no_unexplained_energy_growth":result["max_energy_increase_step"]<1e-9}
    checks={key:bool(value) for key,value in checks.items()}
    summary={"status":"passed" if all(checks.values()) else "failed","checks":checks,
        "initial_theta_deg":initial_angle,"target_theta_deg":c.theta_equilibrium_deg,
        "measured_theta_deg":angle,"fit_windows":final["fits"],
        "contact_half_width_change_m":half_width_change,
        "initial_crossings":before,"final_crossings":after,
        "late_angle_step_change_deg":last_change,"run":result}
    components=connected_polylines(contour_segments(solver))
    if solver.comm.rank==0:
        (folder/"contact.json").write_text(json.dumps(summary,indent=2)+"\n")
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
        times=[o["time"] for o in observations]
        axes[0].plot(times,theta_history);axes[0].axhline(c.theta_equilibrium_deg,color="k",ls="--")
        axes[0].set(xlabel="time [s]",ylabel="apparent angle [degree]")
        axes[1].plot(times,[[p["coordinate"] for p in o["crossings"]] for o in observations])
        axes[1].set(xlabel="time [s]",ylabel="ALL contact crossings x [m]")
        for line in components:
            axes[2].plot(line[:,0],line[:,1])
        axes[2].set(xlabel="x [m]",ylabel="z [m]",aspect="equal")
        axes[2].axhline(c.z_min,color="k")
        for ext in ("pdf","png"):
            fig.savefig(folder/f"contact.{ext}",dpi=140)
        plt.close(fig)
        print(json.dumps({key:value for key,value in summary.items() if key not in ("run","fit_windows")},indent=2),flush=True)
    return summary


def run(config,output,initial_angle=90.):
    radius=config.refinement_circle_radius
    center_z=config.z_min-radius*np.cos(np.deg2rad(initial_angle))
    c=config.changed(refinement_circle_z=float(center_z))
    solver,rows,result=run_case(c,sessile_drop(c.epsilon,radius,initial_angle,c.z_min),output,
                                observer=observe)
    return solver,summarize(solver,result,output,initial_angle)
