"""Corrective PW1 CLI: bounded pilot qualification, resume, and independent audit.

The actual baseline failed its predeclared linear-geometry gate. These commands
reproduce/refine that scientific blocker, not a substitute five-second release.
No flag bypasses the applicability screen or changes the physical baseline.
"""
import json
import math
from pathlib import Path
import sys

import yaml

from .asymmetric import MODEL_ID, LeftNavierConfig
from .asymmetric_run import numerical_identity, run_case
from .mission import doctor, job, lock_busy, read_json, write_json

BLOCKER = "PW1_PINNED_RIGHT_LINEAR_APPLICABILITY"


def controls(original, cfg, mesh_level=0, dt_scale=1., *, horizon=None):
    if mesh_level not in (0, 1, 2) or not math.isfinite(dt_scale) or dt_scale <= 0:
        raise ValueError("Finite positive dt-scale and mesh-level 0/1/2 required")
    dt = cfg["baseline_dt_s"] * dt_scale
    if dt > cfg["render_sample_dt_s"] or not math.isclose(
            cfg["render_sample_dt_s"] / dt, round(cfg["render_sample_dt_s"] / dt), abs_tol=1e-10):
        raise ValueError("dt must divide the fixed 0.01 s export cadence")
    return LeftNavierConfig(a=original["geometry"]["half_width_m"], d=original["geometry"]["depth_m"],
        g=original["physics"]["gravity_m_s2"], nu=original["physics"]["kinematic_viscosity_m2_s"],
        alpha_deg=original["initial"]["alpha_deg"], left_slip_length_m=cfg["left_slip_length_m"],
        nx=cfg["baseline_mesh"][0]*2**mesh_level, nz=cfg["baseline_mesh"][1]*2**mesh_level,
        dt=dt, t_end=cfg["pilot_t_end_s"] if horizon is None else horizon,
        snapshot_dt=cfg["render_sample_dt_s"], integrator="sdirk2")


def pilot_case(root, original, cfg, args, state, heartbeat):
    """Same entry used by CLI and the post-commit checkpoint regression test."""
    c = controls(original, cfg, args.mesh_level, args.dt_scale)
    requested = json.loads(c.to_json())
    source_hash = numerical_identity(root)[0]
    role = "pilot"
    # A previously computed spatial-confirmation case is also a valid finer
    # pilot. Roles are provenance, not a reason to duplicate the same solve.
    for path in sorted((root / cfg["runs_root"]).glob("*/resolved_case.json")):
        resolved = read_json(path)
        identity = read_json(path.parent / "identity.json")
        if (resolved.get("controls") == requested and identity.get("model_id") == MODEL_ID
                and identity.get("source_hash") == source_hash):
            role = identity["role"]
            # Preserve the saved JSON representation as well as its numerical
            # value (e.g. t_end=1 versus 1.0) in the content-addressed payload.
            c = LeftNavierConfig(**resolved["controls"])
            break
    return run_case(root, cfg, c, role, state, heartbeat, args.stop_after)


def mark_blocked(state, output, witness):
    state.update(status="NEEDS_MODEL_DECISION", current_model_id=MODEL_ID,
        current_run_scope="LEFT_NAVIER_RIGHT_NOSLIP", corrected_complete=False,
        target_run_executed=False, selected_corrected_native=None,
        last_error={"code":BLOCKER, "detail":"Exact wall-adjacent P2 slope exceeds the unchanged 0.30 linear screen.",
                    "evidence":witness,
                    "reproduce":"python sloshing_visualization/scripts/pinned_wetting_mission.py blocker-check --corrected"},
        next_action="A nonlinear free-surface/corner geometry decision is needed; do not retune b or unpin P2.",
        corrected_output_directory=str(output))


def corrected_main(root, original, args):
    config_dir = root / "mission/pinned_wetting/left_navier_right_noslip"
    cfg = yaml.safe_load((config_dir / "CONFIG.yaml").read_text())
    controls(original, cfg, args.mesh_level, args.dt_scale)  # validate even before a job
    output = root / cfg["output_root"]
    extra = args.mesh_level != 0 or args.dt_scale != 1.
    if extra:
        output /= f"extra_m{args.mesh_level}_dt{args.dt_scale:g}"
    if args.command == "doctor":
        health = doctor(root)
        health["status"] = "READY_FOR_CORRECTED_QUALIFICATION" if health["status"] == "READY_FOR_REFERENCE" else health["status"]
        print(json.dumps({**health, "corrected_model":cfg}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "status":
        print(json.dumps({**read_json(root / "mission/pinned_wetting/state.json"),
            "kernel_lock_busy":lock_busy(root / "mission/pinned_wetting/job.lock")}, indent=2))
        return 0
    if args.run_id:
        raise ValueError("Corrected runs are content-addressed; resume with identical mesh/dt flags, not --run-id")
    with job(root, " ".join(sys.argv), original["resources"]) as (state, heartbeat):
        from .asymmetric_verify import surface_applicability, verify_native
        if args.command in ("blocker-check", "verify"):
            summary = read_json(output / "pilot_summary.json")
            if args.command == "blocker-check":
                report = surface_applicability(summary["native"])
                report.update(blocker_code=BLOCKER, target_time_range_s=[0., 5.], target_status="NOT_RUN")
                code = 0 if report["linear_applicability_passed"] else 2
            else:
                from .asymmetric_package import verify_package
                report = verify_package(root, output, heartbeat)
                code = 0 if report["passed"] else 2
            print(json.dumps(report, indent=2, ensure_ascii=False))
            state["_job_exit_code"] = code
            return code
        if args.command == "package":
            from .asymmetric_package import package_blocker
            package_blocker(root, output, state, heartbeat)
            state["_job_exit_code"] = 2
            return 2
        if args.command == "render":
            from .extension_output import render_bundle
            summary = read_json(output / "pilot_summary.json")
            render_bundle(root, Path(summary["native"]), output / "preview", final=False, heartbeat=heartbeat)
            return 0
        output.mkdir(parents=True, exist_ok=True)
        state.update(status="PREVIEW", current_model_id=MODEL_ID, current_run_scope="LEFT_NAVIER_RIGHT_NOSLIP")
        summary = pilot_case(root, original, cfg, args, state, heartbeat)
        if not summary["completed"]:
            state.update(next_action="Resume --corrected with identical mesh/dt flags; accepted checkpoint saved.")
            return 0
        write_json(output / "pilot_summary.json", summary)
        screen = surface_applicability(summary["native"])
        write_json(output / "applicability.json", screen)
        # pilot is a deliberately explicit diagnostic command. run/resume also
        # reuse the already completed fixed-b confirmation and package, if present.
        if not screen["linear_applicability_passed"]:
            mark_blocked(state, output, screen)
            if args.command in ("run", "resume") and (output / "native_audit.json").exists():
                from .asymmetric_package import package_blocker
                package_blocker(root, output, state, heartbeat)
            print(json.dumps({"status":state["status"], "summary":summary, "screen":screen,
                "qualified_5s_result":False}, indent=2, ensure_ascii=False))
            state["_job_exit_code"] = 2
            return 2
        # This branch is for a diagnostic mesh probe, not a release designation.
        # A pilot command never claims target/refinement/video acceptance.
        if args.command == "pilot":
            print(json.dumps({"pilot":summary, "screen":screen, "target_status":"NOT_RUN"}, indent=2))
            return 0
        # This continuation has an independently confirmed baseline blocker.
        # An extra diagnostic mesh that passes by itself does not erase it and
        # must not silently launch an unqualified five-second release.
        baseline = read_json(root / cfg["output_root"] / "native_audit.json")
        failed = [r for r in baseline["native_reports"].values() if not r["linear_applicability_passed"]]
        if not failed:
            raise ValueError("Stored baseline blocker is absent; this diagnostic-only continuation needs a fresh model audit")
        mark_blocked(state, output, failed[0]["first_slope_screen_failure"])
        print(json.dumps({"diagnostic_pilot":summary, "baseline_blocker":BLOCKER,
                          "target_status":"NOT_RUN"}, indent=2))
        state["_job_exit_code"] = 2
        return 2
