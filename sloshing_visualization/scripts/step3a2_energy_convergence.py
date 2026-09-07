"""Independent comparison of measured BE perturbation histories and FE endpoints."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from dolfinx import fem
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.equilibrium import _integral,_measures
from sloshing.multiphase.provenance import run_provenance


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("runs",nargs=3)
    parser.add_argument("--output",default="validation_results/step3a2/equilibrium_perturbation/convergence.json")
    args=parser.parse_args()
    paths=[Path(p) for p in args.runs]
    records=[json.loads((p/"summary.json").read_text()) for p in paths]
    if any(r.get("status")!="complete" or r.get("kind")!="perturbation" for r in records):
        raise ValueError("Three complete nonzero perturbation runs required")
    if len({r["equilibrium_source_fingerprint"] for r in records})!=1:
        raise ValueError("Different prepared initial states")
    for key in ("psi_coefficient_sha256","initial_phi_coefficient_sha256","delta"):
        if len({r["perturbation"][key] for r in records})!=1:
            raise ValueError("Different discrete perturbations")
    if len({r["t_end"] for r in records})!=1 or any(r["scheme"]!="be" for r in records):
        raise ValueError("Same physical horizon and pure BE required")
    if not np.allclose([r["dt"] for r in records],records[0]["dt"]/np.array([1,2,4]),rtol=1e-12,atol=0):
        raise ValueError("Require dt, dt/2, dt/4 in that order")
    s=CHNSSolver(ModelConfig(**records[-1]["config"]))
    provenance=run_provenance(s)
    Q,_=s.space.sub(2).collapse()
    fields=[]
    for path,record in zip(paths,records):
        if record["provenance"]["mesh_fingerprint"]!=provenance["mesh_fingerprint"]:
            raise ValueError("Different endpoint mesh/partition")
        field=fem.Function(Q)
        with h5py.File(path/"checkpoint/rank_0000.h5","r") as h:
            if h.attrs["status"]!="complete" or h.attrs["time"]!=record["t_end"]:
                raise ValueError("Incomplete or wrong-time endpoint checkpoint")
            if not np.array_equal(h["geometry"][:],s.mesh.geometry.x):
                raise ValueError("Endpoint geometry differs from the measured mesh")
            field.x.array[:]=h["fields/phi/coefficients"][:]
        fields.append(field)
    dx,_=_measures(s)
    def distance(a,b):
        return float(np.sqrt(max(0.,_integral(s,(a-b)**2*dx))))
    to_finest=[distance(f,fields[-1]) for f in fields]
    consecutive=[distance(fields[i],fields[i+1]) for i in (0,1)]
    table=[]
    for path,r,difference in zip(paths,records,to_finest):
        energy=r["energy_validation"]
        delta=energy["energy_change_final_J_per_m"]
        table.append({"run":str(path),"dt":r["dt"],
            "integrated_D_CH":energy["dissipation_components"]["CH"]["integral_J_per_m"],
            "Delta_E":delta,"final_budget_defect":energy["final_budget_defect_J_per_m"],
            "max_local_defect":energy["max_local_budget_defect_J_per_m"],
            "end_phi_difference_to_finest":difference,"end_energy":r["final"]["E_total"],
            "final_defect_over_actual_change":abs(energy["final_budget_defect_J_per_m"])/abs(delta) if delta else None,
            "mass_error":r["metrics"]["mass_error"],"certified_cells":r["minimum_certified_cells"]})
    diss=np.array([r["integrated_D_CH"] for r in table])
    diss_diff=abs(np.diff(diss))
    defects=abs(np.array([r["final_budget_defect"] for r in table]))
    energy_diff=abs(np.diff([r["end_energy"] for r in table]))
    gates={"nonzero_relaxation":all(r["Delta_E"]<0 and r["integrated_D_CH"]>0 for r in table),
        "budget_defect_decreases":bool(np.all(np.diff(defects)<0)),
        "integrated_CH_converges":bool(diss_diff[1]<diss_diff[0] and diss_diff[1]/abs(diss[-1])<=.05),
        "end_phase_converges":consecutive[1]<consecutive[0],
        "end_energy_converges":bool(energy_diff[1]<energy_diff[0]),
        "mass":all(r["mass_error"]<=1e-10 for r in table),
        "resolution":all(r["certified_cells"]>=8 for r in table),
        "no_energy_growth":all(r["energy_validation"]["energy_growth_ok"] for r in records),
        "finest_energy_closure":records[-1]["energy_validation"]["energy_budget_ok"],
        # Small delta-E is NOT allowed to hide this nonzero-transient closure test.
        "finest_final_defect_over_actual_change":table[-1]["final_defect_over_actual_change"] is not None and table[-1]["final_defect_over_actual_change"]<=.05}
    result={"status":"complete","qualification_status":"passed" if all(gates.values()) else "failed",
        "gates":gates,"table":table,"end_phi_successive_L2_differences":consecutive,
        "CH_integral_successive_differences":diss_diff.tolist(),
        "CH_last_relative_difference":float(diss_diff[1]/abs(diss[-1])),
        "budget_observed_orders":np.log2(defects[:-1]/defects[1:]).tolist(),
        "equilibrium_source_fingerprint":records[0]["equilibrium_source_fingerprint"],
        "provenance":provenance,"policy_source":"STEP3A2_DESIGN.md",
        "energy_accounting":"all accepted samples including t=0, independent trapezoidal power integral"}
    with Path(args.output).open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
    for path,record,row in zip(paths,records,table):
        review={"scope":"STEP3A.2 independent nonzero-transient qualification; original summary preserved",
            "original_reported_single_history_status":record["qualification_status"],
            "actual_change_final_ratio":row["final_defect_over_actual_change"],
            "actual_change_limit":.05,"actual_change_floor_applied":False,
            "nonzero_final_closure_pass":row["final_defect_over_actual_change"] is not None and row["final_defect_over_actual_change"]<=.05,
            "temporal_series_qualification":result["qualification_status"],
            "note":"Legacy relative-to-change metric can be N/A below its significance floor; this predeclared nonzero-transient final gate remains mandatory."}
        with (path/"qualification_review.json").open("x") as stream:
            json.dump(review,stream,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2),flush=True)
    if not all(gates.values()):
        raise SystemExit(2)


if __name__=="__main__":
    main()
