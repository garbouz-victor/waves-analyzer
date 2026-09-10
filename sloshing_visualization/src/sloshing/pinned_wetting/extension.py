"""Bounded continuation of PW1 under the explicit Navier-slip authorization."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import yaml

from .mission import doctor, job, read_json, write_json
from .navier import LABEL, NavierConfig, channel_benchmark
from .extension_run import execute_case, numerical_identity, physical_identity


def run_case(root, cfg, controls, role, state, heartbeat, stop_after=None):
    source_hash, source_files = numerical_identity(root)
    payload = {"config": json.loads(controls.to_json()), "source_hash": source_hash, "role": role}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:12]
    run_id = f"navier-{role}-{digest}"
    run_dir = root / cfg["runs_root"] / run_id
    physical = physical_identity(controls)
    phys_hash = hashlib.sha256(json.dumps(physical, sort_keys=True).encode()).hexdigest()
    identity = {"run_id": run_id, "source_hash": source_hash, "source_files": source_files,
                "physical_config_hash": phys_hash, "physical": physical,
                "physical_model_label": LABEL, "synthetic": False,
                "source_kind": "declared_Navier_linear_FEM", "execution_HEAD": doctor(root)["execution_HEAD"], "role": role}
    if not run_dir.exists():
        run_dir.mkdir(parents=True)
        write_json(run_dir / "identity.json", identity)
        write_json(run_dir / "resolved_case.json", {"model_id": cfg["model_id"], "physical_model_label": LABEL,
            "physical": physical, "controls": json.loads(controls.to_json()), "extension_contract": cfg,
            "original_case": yaml.safe_load((root / "mission/pinned_wetting/CONFIG.yaml").read_text()),
            "origins": {"slip": cfg["slip_origin"], "geometry": "inherited_task", "alpha_and_nu": "package_proposals_accepted_for_mission",
                        "density": "illustrative_energy_scale_only"},
            "headroom": "Unbounded graph eta; z=0 is linearization plane, not a closed upper lid.",
            "film_thickness_m": None, "retained_layer_volume_feedback": "neglected_at_leading_order"})
        for name in source_files:
            dest = run_dir / "source_snapshot" / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / "sloshing_visualization/src/sloshing" / name, dest)
    elif read_json(run_dir / "identity.json") != identity:
        raise ValueError("Existing extension run identity differs")
    if run_id not in state["run_ids"]:
        state["run_ids"].append(run_id)
    state.update(current_run_id=run_id, current_run_scope="DECLARED_EXTENSION", target_run_executed=controls.t_end==5.)
    state.setdefault("extension_runs", {})[role] = {"run_id": run_id, "path": str(run_dir), "config": json.loads(controls.to_json())}
    heartbeat()
    started = time.monotonic()
    summary = execute_case(root, run_dir, controls, identity, state, heartbeat, stop_after=stop_after)
    state["extension_runs"][role]["last_call_wall_s"] = time.monotonic()-started
    state["extension_runs"][role]["completed"] = summary["completed"]
    state.update(accepted_step=summary.get("accepted_states", summary.get("accepted_step", 0)+1)-1,
                 physical_time_s=summary["time_s"], last_checkpoint=summary["native"])
    heartbeat()
    return summary


def extension_main(root, original_cfg, args):
    config_dir = root / "mission/pinned_wetting/extension"
    cfg = yaml.safe_load((config_dir / "CONFIG.yaml").read_text())
    if args.mesh_level != 0 or args.dt_scale != 1.:
        cfg["output_root"] += f"/extra_m{args.mesh_level}_dt{args.dt_scale:g}"
    if args.command == "doctor":
        health=doctor(root)
        if health["status"]=="READY_FOR_REFERENCE":
            health["status"]="READY_FOR_DECLARED_EXTENSION"
        print(json.dumps({**health, "extension": cfg}, indent=2, ensure_ascii=False)); return 0
    if args.command == "status":
        from .mission import lock_busy
        print(json.dumps({**read_json(root / "mission/pinned_wetting/state.json"),
                          "kernel_lock_busy": lock_busy(root / "mission/pinned_wetting/job.lock")}, indent=2)); return 0
    if args.command == "verify":
        from .extension_verify import verify_extension
        # A full independent reconstruction is a heavy local job too. The same
        # lock and measured budget apply; scientific artifacts remain read-only.
        with job(root, " ".join(sys.argv), original_cfg["resources"]) as (state, heartbeat):
            report = verify_extension(root / cfg["output_root"])
            heartbeat()
            print(json.dumps(report, indent=2, ensure_ascii=False))
            state["_job_exit_code"]=0 if report["passed"] else 2
            return state["_job_exit_code"]
    base = original_cfg
    def controls(level, dt_scale, T, slip=cfg["slip_length_m"]):
        return NavierConfig(a=base["geometry"]["half_width_m"], d=base["geometry"]["depth_m"],
             g=base["physics"]["gravity_m_s2"], nu=base["physics"]["kinematic_viscosity_m2_s"],
             alpha_deg=base["initial"]["alpha_deg"], nx=cfg["baseline_mesh"][0]*2**level,
             nz=cfg["baseline_mesh"][1]*2**level, dt=cfg["baseline_dt_s"]*dt_scale,
             t_end=T, snapshot_dt=.01, integrator="sdirk2", slip_length_m=slip)
    with job(root, " ".join(sys.argv), original_cfg["resources"]) as (state, heartbeat):
        if not state["budget"].get("extension_auxiliary_reserved_s"):
            # Preserve the prior reference allowance; charge an additional bounded
            # allowance for this extension's short tests/read-only diagnostics.
            state["budget"]["auxiliary_wall_time_upper_bound_s"] += 600.
            state["budget"]["extension_auxiliary_reserved_s"] = 600.
        state.setdefault("reference_artifact_paths", dict(state.get("artifact_paths", {})))
        state["previous_STOP"] = "mission/pinned_wetting/history/no_slip_STOP_state.json"
        state.update(status="MODEL_DECISION", model_decision=str(config_dir / "model_decision.md"),
                     last_error=None, current_run_scope="DECLARED_EXTENSION")
        state.setdefault("reference_milestones",dict(state.get("milestones",{})))
        state["milestones"]={"M0":"EXTENSION_PREVIEW_COMPLETE","M1":"DECLARED_NAVIER_AND_COATING_IMPLEMENTED",
                             "M2":"TARGET_AND_REFINEMENT_IN_PROGRESS","M3":"FINAL_VERIFICATION_PENDING"}
        state["next_action"] = "Navier pilot -> target5s -> independent dt/mesh -> b sensitivity -> final render/verify"
        heartbeat()
        output = root / cfg["output_root"]
        output.mkdir(parents=True, exist_ok=True)
        if args.command == "package":
            from .extension_release import finalize
            state.update(status="VERIFYING",next_action="verify reviewed final package without restarting PDE")
            heartbeat()
            return finalize(root, output, cfg, state, heartbeat)
        if args.command == "render":
            from .extension_output import render_bundle
            if state.get("selected_extension_native") and (output/"case_summaries.json").exists():
                from .extension_release import write_report, make_manifest, freeze_native_audit
                native=Path(state["selected_extension_native"])
                state.update(status="RENDERING",next_action="renderer-only revision from identical selected native data")
                state["milestones"]["M2"]="TARGET_5S_TIME_SPACE_AND_SLIP_SENSITIVITY_COMPLETE"
                heartbeat()
                freeze_native_audit(root,output)
                render_bundle(root,native,output,final=True,heartbeat=heartbeat)
                summaries=read_json(output/"case_summaries.json")
                write_report(root,output,cfg,state,summaries,read_json(output/"refinement.json"),read_json(output/"sensitivity.json"))
                make_manifest(root,output,cfg,state,summaries)
                state.update(status="VERIFYING",next_action="independent visual/native signoff, then package --extension")
                heartbeat()
                return 0
            native = Path(state["last_checkpoint"])
            render_bundle(root, native, output / "preview", final=False, heartbeat=heartbeat)
            state["status"] = "PREVIEW"
            return 0
        if args.command in ("run","resume") and (output/"manifest.json").exists():
            if read_json(output/"manifest.json").get("status")=="COMPLETE_FOR_DECLARED_MODEL":
                from .extension_release import finalize
                print(json.dumps({"reuse_completed_bundle":str(output),"PDE_restarted":False}),flush=True)
                return finalize(root,output,cfg,state,heartbeat)
        if not (output / "boundary_benchmark.json").exists():
            benchmark = channel_benchmark(NavierConfig(nx=8, nz=12, dt=.005, t_end=.02, snapshot_dt=.01, integrator="sdirk2"))
            if benchmark["relative_error"] > 1e-9:
                raise ValueError("Navier channel benchmark failed")
            write_json(output / "boundary_benchmark.json", benchmark)
        state["status"] = "PREVIEW"
        pilot = run_case(root, cfg, controls(0, 1., cfg["pilot_t_end_s"]), "pilot", state, heartbeat, args.stop_after)
        if not pilot["completed"]:
            return 0
        write_json(output / "pilot_summary.json", pilot)
        print(json.dumps({"pilot": pilot}, ensure_ascii=False, indent=2), flush=True)
        if args.command == "pilot":
            return 0
        from .extension_output import render_bundle
        render_bundle(root, Path(pilot["native"]), output / "preview", final=False, heartbeat=heartbeat)
        if pilot["max_eta_over_a"] > .05 or pilot["max_slope"] > .30:
            raise ValueError("LINEAR_APPLICABILITY: predeclared amplitude/slope gate failed; b is not retuned")
        result=run_to_release(root, original_cfg, cfg, controls, args, state, heartbeat)
        state["_job_exit_code"]=result
        return result


def run_to_release(root, original_cfg, cfg, controls, args, state, heartbeat):
    # Kept separate from the solver for read-only verification and render reuse.
    from .extension_release import complete_release
    return complete_release(root, original_cfg, cfg, controls, args, state, heartbeat, run_case)
