"""New numerical identity and provenance for the corrective PW1 trajectory."""
import hashlib
import json
from pathlib import Path
import shutil
import time

import yaml

from .asymmetric import MODEL_ID, LABEL, SOURCE_KIND, CONTACT_METHOD, LeftNavierFEM, LeftNavierDiagnostics, save_operator_evidence
from .extension_run import execute_case as shared_execute
from .mission import doctor, read_json, write_json, sha256


def physical_identity(c):
    return {"model_id":MODEL_ID,"a_m":c.a,"d_m":c.d,"g_m_s2":c.g,"nu_m2_s":c.nu,"alpha_deg":c.alpha_deg,
        "left_slip_length_m":c.left_slip_length_m,"sigma_N_m":0.,"bottom":"no_slip",
        "left_wall":"impermeable_Navier_tangential_slip","right_wall":"no_slip",
        "geometry":"linearized_free_surface","layer":"irreversible_lower_connected_zero_volume_one_way"}


def numerical_identity(root):
    names = ["config.py","mesh.py","fem_spaces.py","problems.py","time_integrator.py","diagnostics.py",
        "pinned_wetting/coating.py","pinned_wetting/navier.py","pinned_wetting/extension_run.py","pinned_wetting/asymmetric.py",
        "pinned_wetting/asymmetric_run.py"]
    src = root/"sloshing_visualization/src/sloshing"
    hashes = {name:sha256(src/name) for name in names}
    return hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest(), hashes


def execute_case(root,run_dir,config,identity,mission_state,heartbeat,stop_after=None):
    return shared_execute(root,run_dir,config,identity,mission_state,heartbeat,stop_after,
        fem_factory=LeftNavierFEM,diagnostics_factory=LeftNavierDiagnostics,
        contact_method=CONTACT_METHOD,initialize_native=save_operator_evidence,record_sides=("L",))


def compatible_identity(identity):
    """Commit provenance is deliberately not a numerical compatibility field."""
    return {k:v for k,v in identity.items() if k not in ("execution_HEAD","execution_HEAD_at_creation","execution_HEAD_history")}


def run_case(root,cfg,controls,role,state,heartbeat,stop_after=None):
    source_hash,source_files = numerical_identity(root)
    payload = {"config":json.loads(controls.to_json()),"source_hash":source_hash,"role":role,"model_id":MODEL_ID}
    digest = hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()[:12]
    run_id = f"leftnavier-{role}-{digest}"
    run_dir = root/cfg["runs_root"]/run_id
    physical = physical_identity(controls)
    head = doctor(root)["execution_HEAD"]
    identity = {"model_id":MODEL_ID,"run_id":run_id,"source_hash":source_hash,"source_files":source_files,
        "physical_config_hash":hashlib.sha256(json.dumps(physical,sort_keys=True).encode()).hexdigest(),
        "physical":physical,"physical_model_label":LABEL,"source_kind":SOURCE_KIND,
        "synthetic":False,"role":role,"execution_HEAD":head,"execution_HEAD_at_creation":head}
    if not run_dir.exists():
        run_dir.mkdir(parents=True)
        write_json(run_dir/"identity.json",identity)
        write_json(run_dir/"resolved_case.json",{"model_id":MODEL_ID,"physical_model_label":LABEL,
            "physical":physical,"controls":json.loads(controls.to_json()),"extension_contract":cfg,
            "original_case":yaml.safe_load((root/"mission/pinned_wetting/CONFIG.yaml").read_text()),
            "origins":{"boundary_conditions":"explicit_user_corrective_contract", "slip":"user_selected_left_baseline_0.50_m",
                       "geometry_alpha_nu_g":"explicit_user_corrective_baseline","density":"illustrative_energy_scale_only"},
            "headroom":"Unbounded linear graph eta; z=0 is not a lid.",
            "film_thickness_m":None,"retained_layer_volume_feedback":"neglected_at_leading_order"})
        for name in source_files:
            target = run_dir/"source_snapshot"/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(root/"sloshing_visualization/src/sloshing"/name,target)
    else:
        previous = read_json(run_dir/"identity.json")
        if compatible_identity(previous) != compatible_identity(identity):
            raise ValueError("Existing asymmetric numerical identity differs")
        identity = previous  # immutable native identity retains its creation commit
    provenance_path = run_dir/"provenance.json"
    provenance = read_json(provenance_path) if provenance_path.exists() else {
        "execution_HEAD_at_creation":identity["execution_HEAD"],"execution_HEAD_history":[],"run_id":run_id}
    if not provenance["execution_HEAD_history"] or provenance["execution_HEAD_history"][-1]["HEAD"] != head:
        provenance["execution_HEAD_history"].append({"HEAD":head,"observed_epoch_s":time.time(),"action":"create_or_resume"})
        write_json(provenance_path,provenance)
    if run_id not in state.setdefault("run_ids",[]): state["run_ids"].append(run_id)
    state.update(current_model_id=MODEL_ID,current_run_id=run_id,current_run_scope="LEFT_NAVIER_RIGHT_NOSLIP")
    entry = state.setdefault("corrected_runs",{}).setdefault(role,{})
    entry.update(run_id=run_id,path=str(run_dir),config=json.loads(controls.to_json()))
    heartbeat()
    started = time.monotonic()
    summary = execute_case(root,run_dir,controls,identity,state,heartbeat,stop_after)
    entry.update(last_call_wall_s=time.monotonic()-started,completed=summary["completed"])
    state.update(accepted_step=summary.get("accepted_states",summary.get("accepted_step",0)+1)-1,
                 physical_time_s=summary["time_s"],last_checkpoint=summary["native"])
    heartbeat()
    return summary
