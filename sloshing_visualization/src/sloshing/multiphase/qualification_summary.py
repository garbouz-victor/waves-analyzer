"""Measured STEP 3A.1 gate dependencies; absence/complete is never scientific PASS."""
import json
from pathlib import Path

CORE_VERDICT="CORE CONTACT MODEL VALIDATED THROUGH REPAIRED ENERGY/RESOLUTION GATES — REMAINING STEP 3A BENCHMARKS REQUIRED"


def read_qualification(path):
    path=Path(path)
    if not path.exists():
        return {"qualification_status":"not_run","source":str(path)}
    result=json.loads(path.read_text())
    result.setdefault("qualification_status","not_qualified")
    return result


def repaired_core_summary(base):
    base=Path(base)
    compatible=read_qualification(base/"energy_compatible/qualification.json")
    laplace=read_qualification(base/"laplace_time/qualification.json")
    startup=read_qualification(base/"startup_energy_qualification.json")
    contact=read_qualification(base/"contact60_resolved/qualification.json")
    parts=contact.get("parts",{})
    def measured_pass(record,kind):
        return record["qualification_status"]=="passed" and record.get("evidence_kind")==kind
    gates={"compatible_energy_baseline_passed":measured_pass(compatible,"compatible_energy_series"),
           "laplace_pressure_and_energy_passed":measured_pass(laplace,"laplace_time_series"),
           "contact_60_geometry_passed":parts.get("geometry")=="passed",
           "contact_60_direction_passed":parts.get("direction")=="passed",
           "contact_60_mass_passed":parts.get("mass")=="passed",
           "contact_60_resolution_passed":parts.get("resolution")=="passed",
           "contact_60_startup_energy_passed":measured_pass(startup,"startup_energy_series"),
           "contact_60_energy_budget_passed":parts.get("energy_budget")=="passed",
           "contact_60_settling_passed":parts.get("settling")=="passed",
           "contact_60_scientific_run_passed":measured_pass(contact,"contact_scientific_run")}
    core_passed=all(gates.values())
    partials={k:parts.get(k,"not_run") for k in ("geometry","direction","mass","resolution","energy_budget","settling")}
    return {"status":CORE_VERDICT if core_passed else "MODEL NOT YET VALIDATED",
            "qualification_status":"passed_repaired_core_only" if core_passed else "failed_or_incomplete",
            "gates":gates,"contact_60":{**partials,"overall":"passed" if core_passed else
                ("failed" if "failed" in partials.values() else "not_qualified")},
            "compatible_energy":compatible,"laplace_time":laplace,"startup_energy":startup,
            "contact_result":contact,
            "remaining_step3a":["theta 90/120 PDE cases","epsilon convergence","slip sensitivity",
                "mobility sensitivity","falling-film PDE benchmark"],
            "tank_bridge_allowed":False,"film_demo_allowed":False}
