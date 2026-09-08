"""External administrative session; frozen numerical classes are never patched."""
import json
import time
from datetime import datetime, timezone
from step3a10_admin import (ROOT, OLD, TRAJECTORY, numerical_guard, prefix_guard,
                            require_authorization, read, sha, lineage)


def main():
    # Fail closed before imports that could create/call a numerical engine.
    auth = require_authorization()
    status_path = ROOT / "continuation/controller_status.json"
    if status_path.exists():
        raise ValueError("One continuation session only; inspect previous result, never automatically retry")
    f = numerical_guard(); manifest = prefix_guard()
    if (auth["STEP3A9_production_policy_sha256"] != f["policy_sha256"] or
            auth["STEP3A9_source_archive_sha256"] != f["source_archive_sha256"] or
            auth["schedule_sha256"] != f["L0_schedule"]["schedule_sha256"]):
        raise ValueError("Authorization does not bind actual frozen numerical identity")

    # All following numerical imports are the unchanged STEP3A9 classes/helpers.
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    from sloshing.multiphase.phase_rate_history import atomic_json
    from sloshing.multiphase.phase_rate_trajectory import CostStop
    from sloshing.multiphase.divergence_trajectory import DivergenceTrajectory
    from sloshing.multiphase.full_phase_rate import CHNSPhaseRateBE
    from sloshing.multiphase.equilibrium import equilibrium_config_fingerprint
    from sloshing.multiphase.linearized_ch import array_hash
    from sloshing.multiphase.performance import petsc_events
    from dolfinx import fem
    from petsc4py import PETSc

    used_total, used_full = d.used(), d.used("full_path")
    calibration_wall = read(ROOT/"cost_audit/forecast.json")["calibration_administrative_wall_s"]
    remaining = min(f["policy"]["total_wall_s"]-used_total-calibration_wall,
                    f["policy"]["full_L0_wall_s"]-used_full)
    if remaining <= 0:
        raise CostStop("Unchanged actual wall cap exhausted before continuation")
    session_root = TRAJECTORY/"sessions"
    session = session_root/f"{len(list(session_root.iterdir())):03d}"
    session.mkdir(exist_ok=False)
    record = {"status":"running", "preliminary":False, "full_path":True,
              "implementation":d.implementation(), **lineage(f),
              "historical_trajectory_identity":manifest["identity"],
              "authorization_sha256":sha(ROOT/"continuation_authorization.json"),
              "external_continuation_controller_sha256":sha(__file__)}
    controller = {"execution_status":"running", "authorization_sha256":record["authorization_sha256"],
                  "start_accepted_step":127, "expected_final_step":1068,
                  "start_checkpoint_sha256":manifest["checkpoint_semantic_sha256"],
                  "numerical_guard":"passed", "source_guard_role":"new external numerical guard, not old HEAD guard",
                  "prefix_integrity":"passed", "full_L0_cap_s":14400., "total_cap_s":19800.,
                  "external_continuation_controller_sha256":sha(__file__), **lineage(f),
                  "start_timestamp_UTC":datetime.now(timezone.utc).isoformat()}
    atomic_json(session/"session_status.json",record); atomic_json(status_path,controller)
    started=time.perf_counter(); runner=None; initialization=None; error=None; PETSc.Log.begin()
    try:
        # Identical operations to frozen input_engine, but no call to its old HEAD guard.
        s,eq,psi,initial,provenance = d.historical.scientific_input()
        if (array_hash(initial)!=f["input_hashes"]["initial_phi_sha256"] or
                array_hash(psi)!=f["input_hashes"]["perturbation_coefficient_sha256"] or
                provenance["equilibrium_fingerprint"]!=f["input_hashes"]["equilibrium_fingerprint"] or
                equilibrium_config_fingerprint(s.config)!=f["physical_config_fingerprint"] or
                provenance["mesh_fingerprint"]!=manifest["identity"]["mesh_fingerprint"]):
            raise ValueError("Historical physical inputs mismatch")
        phi=fem.Function(eq.function_space); phi.x.array[:]=initial; phi.x.scatter_forward(); s.initialize(phi)
        engine=CHNSPhaseRateBE(s)
        if array_hash(s.state.x.array[engine.physical_phi_map])!=f["input_hashes"]["initial_phi_sha256"]:
            raise ValueError("Initial physical phi mismatch")
        runner=DivergenceTrajectory(engine,d.schedules()[0],TRAJECTORY,manifest["identity"],
            f["inherited_policy"]["policy"],f["target_mass"],full_policy=d.POLICY,
            resume=True,source_guard=numerical_guard)
        d.initial_check(runner.journal.rows[0]); initialization=time.perf_counter()-started
        if runner.s.step_number!=127 or runner.schedule.intervals[127].step!=128:
            raise ValueError("Continuation must start at exact interval128")
        numerical_guard(); prefix_guard()
        if time.perf_counter()-started>=remaining:
            raise CostStop("Actual wall budget exhausted in initialization")
        result=runner.run_until(1068,wall_remaining_s=remaining-(time.perf_counter()-started))
        atomic_json(TRAJECTORY/"status.json",result)
        return result
    except Exception as exc:
        error=exc; raise
    finally:
        final=runner.status() if runner is not None else None
        if runner is not None:
            atomic_json(TRAJECTORY/"status.json",final); runner.close()
        wall=time.perf_counter()-started
        atomic_json(session/"session_status.json",{**record,"status":"failed" if error else "complete",
            "error":str(error) if error else None,"wall_s":wall,"initialization_JIT_s":initialization,
            "start_accepted_step":127,"target_steps":1068,"accepted_steps":final["accepted_steps"] if final else 127,
            "PETSc":petsc_events()})
        prefix_guard(after_continuation=True)
        atomic_json(status_path,{**controller,"execution_status":"failed" if error else "complete",
            "error":str(error) if error else None,"wall_s":wall,"initialization_JIT_s":initialization,
            "final_accepted_step":final["accepted_steps"] if final else 127,
            "interval128_executed":bool(final and final["accepted_steps"]>=128 or
                (TRAJECTORY/"FAILED.json").exists() and read(TRAJECTORY/"FAILED.json")["failed_interval"]>=128),
            "COMPLETE_exists":(TRAJECTORY/"COMPLETE.json").exists(),
            "historical_prefix_bytes_preserved":True,"end_timestamp_UTC":datetime.now(timezone.utc).isoformat()})


if __name__=="__main__":
    print(json.dumps(main(),indent=2))
