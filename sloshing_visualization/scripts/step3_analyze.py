"""Measured STEP 3A.1 aggregation. Historical STEP 3 data are read-only evidence."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/"src"))
os.environ.setdefault("MPLCONFIGDIR",str(root/"validation_results/step3a1/cache/matplotlib"))
from sloshing.multiphase.energy_validation import read_history,validate_energy
from sloshing.multiphase.qualification_summary import repaired_core_summary,read_qualification
from sloshing.multiphase.provenance import source_hash


def aggregate_compatible(base):
    folder=base/"energy_compatible"
    runs=[]
    for path in sorted(folder.glob("*/summary.json")):
        run=json.loads(path.read_text())
        history=path.parent/"history.csv"
        gate=validate_energy(read_history(history)) if history.exists() else None
        passed=(run["status"]=="complete" and gate is not None and gate["qualified"]
                and run.get("interface_resolved_all_times",False))
        runs.append({"path":str(path.parent.relative_to(base)),"execution_status":run["status"],
            "qualification_status":"passed" if passed else "failed",
            "scheme":run["config"]["time_scheme"],"energy_validation":gate,
            "minimum_cells":run.get("minimum_interface_cells_all_times"),
            "provenance":run.get("provenance")})
    if not runs:
        return
    comparison=read_qualification(folder/"temporal_comparison.json")
    schemes={r["scheme"] for r in runs}
    checks={"single_histories_passed":all(r["qualification_status"]=="passed" for r in runs),
            "BE_and_BDF2_measured":{"be","bdf2"}<=schemes,
            "temporal_comparison_passed":comparison["qualification_status"]=="passed"}
    status="passed" if all(checks.values()) else ("failed" if not checks["single_histories_passed"] else "conditional")
    result={"evidence_kind":"compatible_energy_series","qualification_status":status,
        "checks":checks,"runs":runs,"temporal_comparison":comparison,
        "scope":"compatible-angle energy foundation; failed single history blocks later phases"}
    (folder/"qualification.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,default=root/"validation_results/step3a1")
    args=parser.parse_args()
    base=args.results
    historical=root/"validation_results/step3"
    if base.resolve()==historical.resolve() or historical.resolve() in base.resolve().parents:
        parser.error("Historical STEP 3 results must not be overwritten; select the step3a1 output tree")
    base.mkdir(parents=True,exist_ok=True)
    old=read_history(historical/"contact/theta60_resolved/history.csv")
    energy=validate_energy(old)
    ch=energy["dissipation_components"]["CH"]
    audit={"source":"historical STEP 3, different source and time scheme; not a new convergence-series member",
        "energy_validation":energy,"first_CH_interval_J_per_m":ch["first_interval_J_per_m"],
        "first_CH_interval_percent":100*ch["first_interval_J_per_m"]/ch["integral_J_per_m"],
        "actual_energy_decrease_J_per_m":old[0]["E_total"]-old[-1]["E_total"]}
    (base/"historical_energy_audit.json").write_text(json.dumps(audit,indent=2,allow_nan=False)+"\n")
    aggregate_compatible(base)
    summary=repaired_core_summary(base)
    summary["historical_energy_audit"]=audit
    summary["analysis_source_sha256"]=source_hash()
    (base/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    with (base/"summary.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["gate","passed"]);writer.writerows(summary["gates"].items())
    print(summary["status"])
    print(json.dumps(audit,indent=2))


if __name__=="__main__":
    main()
