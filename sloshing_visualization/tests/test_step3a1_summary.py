import json
from sloshing.multiphase.qualification_summary import repaired_core_summary,CORE_VERDICT


def save(path,payload):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(payload))


def test_complete_execution_or_missing_data_never_passes(tmp_path):
    save(tmp_path/"energy_compatible/qualification.json",{"status":"complete"})
    assert not repaired_core_summary(tmp_path)["gates"]["compatible_energy_baseline_passed"]
    assert repaired_core_summary(tmp_path)["status"]=="MODEL NOT YET VALIDATED"


def test_summary_reads_all_measured_partial_gates(tmp_path):
    for path,kind in (("energy_compatible/qualification.json","compatible_energy_series"),
                      ("laplace_time/qualification.json","laplace_time_series"),
                      ("startup_energy_qualification.json","startup_energy_series")):
        save(tmp_path/path,{"qualification_status":"passed","evidence_kind":kind})
    parts={k:"passed" for k in ("geometry","direction","mass","resolution","energy_budget","settling")}
    save(tmp_path/"contact60_resolved/qualification.json",{"parts":parts,
        "qualification_status":"passed","evidence_kind":"contact_scientific_run"})
    result=repaired_core_summary(tmp_path)
    assert result["status"]==CORE_VERDICT
    assert not result["tank_bridge_allowed"]
    parts["energy_budget"]="failed"
    save(tmp_path/"contact60_resolved/qualification.json",{"parts":parts})
    result=repaired_core_summary(tmp_path)
    assert result["contact_60"]["geometry"]=="passed"
    assert result["contact_60"]["overall"]!="passed"


def test_single_history_pass_does_not_replace_a_temporal_series(tmp_path):
    save(tmp_path/"energy_compatible/qualification.json",{
        "qualification_status":"passed","evidence_kind":"compatible_single_history"})
    assert not repaired_core_summary(tmp_path)["gates"]["compatible_energy_baseline_passed"]
