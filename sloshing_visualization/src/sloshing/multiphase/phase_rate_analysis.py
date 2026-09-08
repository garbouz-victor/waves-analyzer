"""Read-only measured three-level comparisons; incomplete levels never qualify."""
from pathlib import Path
import numpy as np
from .phase_rate_history import read_json,Journal,load_archive,validate_complete
from .phase_rate_schedule import scalar_convergence
from .energy_validation import validate_energy


def level_metrics(folder, schedule):
    folder=Path(folder); path=folder/"scalar_history.jsonl"
    if not path.exists(): return {"execution_status":"not_run","expected_steps":schedule.nsteps,"accepted_steps":0}
    rows=Journal.read(path); accepted=rows[1:]; last=rows[-1]
    complete=(folder/"COMPLETE.json").exists()
    if complete:
        validate_complete(folder,read_json(folder/"COMPLETE.json")["identity"],schedule)
    times=Journal.read(folder/"timing_history.jsonl") if (folder/"timing_history.jsonl").exists() else []
    sessions=[read_json(p) for p in folder.glob("sessions/*/timing.json")]
    work=[r["work"] for r in accepted]
    I=last["cumulative_CH_dissipation"]; change=last["E_total"]-rows[0]["E_total"]
    worst=min(rows,key=lambda r:r["cells_across_transition_certified_min"])
    maxit=max(accepted,key=lambda r:r["solver"]["snes_iterations"],default=None)
    boundaries=[]
    for a,b in zip(accepted,accepted[1:]):
        if a["dt"]!=b["dt"]:
            boundaries.append({"first_changed_step":b["step"],"time":b["time"],
                "dt_before":a["dt"],"dt_after":b["dt"],"iterations_before":a["solver"]["snes_iterations"],
                "iterations_after":b["solver"]["snes_iterations"],
                "iterations_second_after":rows[b["step"]+1]["solver"]["snes_iterations"] if b["step"]+1<len(rows) else None})
    failed=(folder/"FAILED.json").exists()
    result={"execution_status":"failed" if failed else "complete" if complete else "incomplete",
        "expected_steps":schedule.nsteps,"accepted_steps":last["step"],"final_time":last["time"],
        "rejected_steps":int(failed),"physical_all_accepted":all(all(r["physical_gates"].values()) for r in rows),
        "min_dt":min((r["dt"] for r in accepted),default=None),"max_dt":max((r["dt"] for r in accepted),default=None),
        "wall_s":sum(s["wall_s"] for s in sessions),"accepted_interval_wall_s":sum(r["timing"]["total_s"] for r in times),
        "sec_per_accepted_interval":float(np.mean([r["timing"]["total_s"] for r in times])) if times else None,
        "max_mass_domain":max(abs(r["mass_error_relative_to_domain"]) for r in rows),
        "max_target_mass_domain":max(abs(r["target_mass_error_relative_to_domain"]) for r in rows),
        "final_mass_error":last["mass_error_abs"],"min_certified_cells":worst["cells_across_transition_certified_min"],
        "resolution_min_time":worst["time"],"resolution_min_witness":worst["resolution"]["worst_cell"],
        "F_initial":rows[0]["E_total"],"final_F":last["E_total"],"DeltaF":change,"integrated_D_CH":I,
        "final_budget_defect":last["energy_budget_defect"],"budget_over_DeltaF":abs(last["energy_budget_defect"]/change) if change else None,
        "final_shared_budget_normalizations":{k:v for k,v in last.items() if k.startswith("energy_budget") or k.startswith("energy_change")},
        "first_interval_fraction":accepted[0]["interval_CH_dissipation"]/I if I and accepted else None,
        "largest_interval_fraction":max((r["interval_CH_dissipation"] for r in accepted),default=0)/I if I else None,
        "max_weak_work_defect":max((abs(r["weak_work_defect"]) for r in work),default=None),
        "max_decomposition_roundoff":max((abs(r["decomposition_roundoff"]) for r in work),default=None),
        "max_local_continuum_defect":max((abs(r["continuum_local_defect"]) for r in work),default=None),
        "cumulative_discrete_work":last.get("work_sums"),"boundary_statistics":boundaries,
        "max_Newton":maxit["solver"]["snes_iterations"] if maxit else None,
        "max_Newton_step":maxit["step"] if maxit else None,"max_Newton_time":maxit["time"] if maxit else None,
        "Newton_distribution":{k:float(fn([r["solver"]["snes_iterations"] for r in accepted])) for k,fn in
            (("min",np.min),("mean",np.mean),("p50",np.median),("p95",lambda v:np.percentile(v,95)))} if accepted else None,
        "solver_all_positive":all(r["solver"]["reason"]>0 for r in accepted),
        "KSP_failures":sum(int(any(t["KSP_reason"] is not None and t["KSP_reason"]<0 for t in r["solver"]["trace"])) for r in accepted),
        "first_work":work[:5],"final_step":last.get("solver")}
    # Shared validation is applied to scalar numeric diagnostic rows only.
    numeric=[{k:v for k,v in r.items() if isinstance(v,(float,int))} for r in rows]
    result["shared_energy_validation"]=validate_energy(numeric) if len(numeric)>1 else None
    return result


def compare_complete_levels(root, schedules, M, policy):
    root=Path(root); levels=[level_metrics(root/f"level{i}",s) for i,s in enumerate(schedules)]
    if not all(r["execution_status"]=="complete" for r in levels):
        return {"levels":levels,"overall_temporal_qualification":False,"convergence":"not_run: three complete levels required"}
    identities=[read_json(root/f"level{i}/COMPLETE.json")["identity"] for i in range(3)]
    comparable=[{k:v for k,v in r.items() if k not in ("level","schedule_sha256")} for r in identities]
    if not all(r==comparable[0] for r in comparable): raise ValueError("Mixed-source/physics/stack temporal series forbidden")
    def norm(v): return float(np.sqrt(max(0.,v @ (M@v))))
    histories=[Journal.read(root/f"level{i}/scalar_history.jsonl") for i in range(3)]
    comparisons=[]
    for index in schedules[0].common_parent_indices:
        fields=[]; metas=[]; rows=[]
        for i,s in enumerate(schedules):
            arrays,meta=load_archive(root/f"level{i}/fields/{index*s.factor:06d}")
            fields.append(arrays); metas.append(meta); rows.append(histories[i][index*s.factor])
        if not all(m["time"]==metas[0]["time"] for m in metas): raise ValueError("Common time differs; interpolation forbidden")
        pair=[]
        for a,b in ((0,1),(1,2)):
            e={"phi_L2":norm(fields[b]["phi"]-fields[a]["phi"]),"mu_L2":norm(fields[b]["mu"]-fields[a]["mu"])}
            e.update(phi_relative=e["phi_L2"]/max(norm(fields[b]["phi"]),1e-30),
                mu_relative=e["mu_L2"]/max(norm(fields[b]["mu"]),1e-30))
            e.update({key:abs(rows[b][key]-rows[a][key]) for key in ("E_total","CH_dissipation","cumulative_CH_dissipation")})
            pair.append(e)
        comparisons.append({"parent_step":index,"time":metas[0]["time"],"pairs":pair})
    last=comparisons[-1]; phi=[r["phi_L2"] for r in last["pairs"]]
    convergence=scalar_convergence([r["integrated_D_CH"] for r in levels],[r["DeltaF"] for r in levels],
        [r["final_budget_defect"] for r in levels],phi,[r["final_F"] for r in levels],policy)
    convergence["gates"]["physical_all_levels"]=all(r["physical_all_accepted"] and r["solver_all_positive"] for r in levels)
    convergence["gates"]["shared_finest_energy_gate"]=levels[-1]["shared_energy_validation"]["energy_budget_ok"]
    convergence["passed"]=all(convergence["gates"].values())
    return {"levels":levels,"convergence":convergence,"endpoint_phi_gaps":phi,"common_times":comparisons,
        "max_common_phi_gaps":[max(r["pairs"][i]["phi_L2"] for r in comparisons) for i in (0,1)],
        "overall_temporal_qualification":convergence["passed"]}
