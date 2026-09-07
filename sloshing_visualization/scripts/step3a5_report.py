#!/usr/bin/env python3
"""Render the measured STOP checkpoint; never run a solver or change policy."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
from sloshing.multiphase.benchmarks.phase_rate_ch import (
    RESULTS, PRIOR, POLICY, HASHES, TINY_DT, read_json, write_json, freeze_policy)
from sloshing.multiphase.provenance import source_hash


def test_counts(path):
    cases = ET.parse(path).findall(".//testcase")
    skipped = sum(c.find("skipped") is not None for c in cases)
    failed = sum(c.find("failure") is not None or c.find("error") is not None for c in cases)
    return {"passed": len(cases)-skipped-failed, "skipped": skipped, "failed_or_error": failed}


def main():
    frozen = freeze_policy()
    moderate = read_json(RESULTS/"moderate_dt/status.json")
    audit = read_json(RESULTS/"algebraic_equivalence/moderate_postmortem/status.json")
    if moderate["qualification_status"] != "failed" or audit["PDE_steps"] != 0:
        raise ValueError("This report renders only the measured moderate STOP checkpoint")
    blocked = ("tiny_step_linear", "tiny_step_nonlinear", "pilot", "tiny_step_full_ch_from_isolated")
    if any((RESULTS/p/"status.json").exists() for p in blocked):
        raise ValueError("Unexpected downstream artifact: manually audit the STOP boundary")
    dif = moderate["differences"]
    timing = read_json(RESULTS/"moderate_dt/timing.json")
    auditing = read_json(RESULTS/"algebraic_equivalence/moderate_postmortem/timing.json")
    old_snes = read_json(PRIOR/"performance/isolated_pilot/status.json")
    tests = {name: test_counts(RESULTS/f"{name}_tests_final.xml") for name in ("pure", "docker")}
    summary = {"iteration": "STEP3A5", "scientific_verdict": "MODEL NOT YET VALIDATED",
        "iteration_verdict": "MODERATE-dt EQUIVALENCE GATE FAILED; PRODUCTION TINY-dt PHASE-RATE PROOF NOT RUN",
        "status": "stopped_at_moderate_equivalence_gate", "moderate": moderate,
        "read_only_postmortem": audit, "required_hashes": HASHES, "policy": frozen,
        "root_cause_qualification": "NOT OBTAINED: mandatory same-dt production rate experiment not authorized",
        "downstream": {stage: {"status": "NOT RUN", "authorized": False,
            "reason": "frozen moderate mu and D_CH equivalence thresholds failed"} for stage in blocked},
        "full_CHNS_rate_implementation": "NOT IMPLEMENTED: phase 10 blocked by phase 4",
        "tests": tests, "pure_deselected_validation_tests": 2,
        "GitHub_Actions": "workflow extended; remote run NOT RUN, no push",
        "newton_monitor": "newtonls/bt unchanged; raw trace full steps (lambda=1); MUMPS INFOG(1)=0",
        "schedule_provenance_note": "dt_schedule in the physical observer provenance is its configured historical schedule, NOT executed here. Actual one-off moderate dt=1e-6 is recorded in each step log. Optimized plan was hash-checked, not executed.",
        "cost": {"moderate_wall_s": timing["wall_s"], "moderate_CPU_s": timing["cpu_s"],
            "absolute_moderate_SNES_s": timing["records"]["nonlinear_SNES"][0]["wall_s"],
            "rate_moderate_SNES_s": timing["records"]["nonlinear_SNES"][1]["wall_s"],
            "postmortem_wall_s": auditing["wall_s"], "pilot_seconds_per_accepted_interval": None},
        "historical_absolute_tiny_source": str(PRIOR/"performance/isolated_pilot/status.json"),
        "historical_absolute_tiny_provenance": read_json(PRIOR/"performance/isolated_pilot/provenance.json"),
        "historical_absolute_tiny_SNES": old_snes}
    write_json(RESULTS/"summary.json", summary)
    with (RESULTS/"summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("scope", "metric", "absolute", "phase_rate", "difference", "threshold", "status"))
        for key, value in dif.items():
            limit = POLICY["moderate_energy_mass_absolute"] if key.endswith("absolute") else POLICY["moderate_phi_mu_D_relative"]
            writer.writerow(("moderate_dt", key, "", "", value, limit, "PASS" if moderate["gates"][key] else "FAIL"))
        for stage in blocked:
            writer.writerow((stage, "execution", "", "", "", "", "NOT RUN: moderate gate failed"))
    performance = RESULTS/"performance"; performance.mkdir(exist_ok=True)
    write_json(performance/"timing_breakdown.json", {"moderate_comparison": timing,
        "read_only_postmortem": auditing, "tiny_production_step": None, "pilot": None,
        "note": "PETSc events are nested/inclusive; initialization includes JIT and weak mu projection"})
    inputs = [RESULTS/"moderate_dt/status.json", RESULTS/"moderate_dt/provenance.json",
        RESULTS/"algebraic_equivalence/moderate_postmortem/status.json", RESULTS/"pure_tests_final.xml",
        RESULTS/"docker_tests_final.xml", PRIOR/"performance/failed_interval_audit/status.json"]
    write_json(RESULTS/"report_provenance.json", {
        "git_commit_at_report_generation": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "current_multiphase_source_sha256": source_hash(),
        "report_generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scientific_run_provenance_is_separate": True,
        "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}})
    print(json.dumps({"verdict": summary["iteration_verdict"], "tests": tests, "cost": summary["cost"]}, indent=2))


if __name__ == "__main__":
    main()
