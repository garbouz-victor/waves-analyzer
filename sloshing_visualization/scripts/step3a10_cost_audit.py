"""Frozen, schedule-aware, read-only STEP3A10 cost experiment."""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile
import time
import numpy as np
from step3a10_admin import ROOT, OLD, TRAJECTORY, MANIFEST, numerical_guard, prefix_guard, read, sha, lineage

OUT = ROOT / "cost_audit"
POLICY = {
    "version": "step3a10-schedule-aware-cost-v1", "full_L0_cap_s": 14400., "total_cap_s": 19800.,
    "observed_last_step": 127, "expected_steps": 1068,
    "core_classes": "is_audit only", "archive_classes": ["is_field", "is_checkpoint", "is_audit"],
    "negative_core_roundoff_s": 1e-9, "negative_session_overhead_stop_s": -1.,
    "n_ge_5": "1.25 * empirical_p95", "n_2_to_4": "1.50 * max", "n_1": "2.00 * sample",
    "percentile_method": "linear", "unseen_core": "fail", "unseen_archive": "5 replays; 1.50 * max",
    "calibration_repetitions": 5, "initialization": "1.25 * max completed primary initialization_JIT_s",
    "session_overhead": "1.50 * max(overhead, 0)", "future_sessions": 1,
    "backtests": [[50, 51, 64], [64, 65, 96], [96, 97, 127], [50, 1, 127]],
    "backtest_gate": "actual <= predicted for ALL windows; no tolerance",
    "calibration_filesystem": "same /work bind, temporary cost_audit subdirectory",
    "full_runner_separate_residual_json_writes": False,
    "full_runner_audit_storage": "full observer audit values are in scalar journal, not isolated residual_audits",
    "calibration_wall_charged_to_total_forecast": True,
    "source_changes_or_factor_retuning_after_measurement": False,
}
ADMIN_FILES = ["scripts/step3a10_admin.py", "scripts/step3a10_cost_audit.py",
               "scripts/step3a10_continue.py", "scripts/step3a10_report.py",
               "STEP3A10_RESUMPTION_DESIGN.md", "tests/test_step3a10_resumption.py"]


class CostAuditFailure(ValueError):
    pass


def event_class(field, checkpoint, audit):
    return f"F{int(bool(field))}_C{int(bool(checkpoint))}_A{int(bool(audit))}"


def events(schedule):
    return [{"step": i.step, "block": i.block, "dt": i.dt, "time": i.end,
             "is_field": i.step in schedule.field_indices,
             "is_checkpoint": i.step in schedule.checkpoint_indices,
             "is_audit": i.step in schedule.audit_indices,
             "block_boundary": i.step in schedule.block_ends,
             "first_post_boundary": i.step-1 in schedule.block_ends,
             "archive_class": event_class(i.step in schedule.field_indices,
                 i.step in schedule.checkpoint_indices, i.step in schedule.audit_indices),
             "core_class": "CORE_AUDIT" if i.step in schedule.audit_indices else "CORE_NORMAL"}
            for i in schedule.intervals]


def counts(schedule_events, split=127):
    result = {}
    for row in schedule_events:
        c = result.setdefault(row["archive_class"], {"total": 0, "observed": 0, "future": 0})
        c["total"] += 1
        c["observed" if row["step"] <= split else "future"] += 1
    return dict(sorted(result.items()))


def estimate(samples, *, calibration=None):
    a = np.asarray(samples, dtype=float)
    if not np.isfinite(a).all() or np.any(a < 0):
        raise CostAuditFailure("Invalid timing sample")
    n = len(a)
    if n >= 5:
        # Linear is the default in both the older host and pinned NumPy APIs.
        value, rule = 1.25*float(np.percentile(a, 95)), "1.25*p95"
    elif n >= 2:
        value, rule = 1.5*float(max(a)), "1.50*max"
    elif n == 1:
        value, rule = 2.*float(a[0]), "2.00*sample"
    elif calibration is not None:
        b = np.asarray(calibration, dtype=float)
        if len(b) != 5 or not np.isfinite(b).all() or np.any(b < 0):
            raise CostAuditFailure("Exactly five finite calibration timings required")
        value, rule = 1.5*float(max(b)), "1.50*max(5 calibration replays)"
    else:
        raise CostAuditFailure("Unseen required class has no qualified calibration")
    return {"observed_n": n, "rule": rule, "upper_s": value,
            "samples_s": a.tolist(), "p95_s": float(np.percentile(a,95)) if n else None}


def decompose(timing_rows, schedule_events):
    result = []
    for r in timing_rows:
        t = r["timing"]
        total, archival = float(t["total_s"]), float(t["archival_s"])
        if not np.isfinite([total, archival]).all() or min(total, archival) < 0:
            raise CostAuditFailure("Nonfinite or negative observed cost")
        core = total-archival
        if core < -POLICY["negative_core_roundoff_s"]:
            raise CostAuditFailure("Materially negative core timing")
        event = schedule_events[r["step"]-1]
        if event["step"] != r["step"]:
            raise CostAuditFailure("Timing/event sequence mismatch")
        result.append({**event, "total_s": total, "archival_s": archival,
                       "core_s": max(0.,core), "core_roundoff_clamped_s": min(0.,core)})
    return result


def predict(training, target, calibration):
    estimates, components = {"core": {}, "archive": {}}, []
    for group, label, sample in (("core", "core_class", "core_s"), ("archive", "archive_class", "archival_s")):
        for cls in sorted({r[label] for r in target}):
            samples = [r[sample] for r in training if r[label] == cls]
            replay = calibration.get(cls) if group == "archive" else None
            est = estimate(samples, calibration=replay)
            count = sum(r[label] == cls for r in target)
            estimates[group][cls] = est
            components.append({"component": cls, "kind": group, "future_count": count,
                               "upper_s": est["upper_s"], "forecast_s": count*est["upper_s"],
                               "observed_n": est["observed_n"], "rule": est["rule"]})
    return {"estimates": estimates, "components": components,
            "predicted_interval_s": sum(c["forecast_s"] for c in components)}


def backtest(training, target, calibration, actual):
    prediction = predict(training,target,calibration)
    upper = prediction["predicted_interval_s"]
    if not np.isfinite(actual) or actual < 0:
        raise CostAuditFailure("Invalid actual backtest cost")
    return {**prediction, "actual_interval_s": actual,
            "forecast_over_actual": upper/actual if actual else None,
            "passed": bool(actual <= upper)}


def session_allowances(sessions, observed):
    previous, rows = 0, []
    for s in sessions:
        end = s["target_steps"]
        if s["status"] != "complete" or not s["full_path"] or end <= previous:
            raise CostAuditFailure("Invalid primary session segmentation")
        init, wall = s["initialization_JIT_s"], s["wall_s"]
        if init is None or not np.isfinite([init, wall]).all() or min(init, wall) < 0:
            raise CostAuditFailure("Missing/invalid initialization wall")
        interval = sum(r["total_s"] for r in observed if previous < r["step"] <= end)
        overhead = wall-init-interval
        if not np.isfinite(overhead) or overhead < -1.:
            raise CostAuditFailure("Session timing accounting inconsistency")
        rows.append({"first_step": previous+1, "last_step": end, "wall_s": wall,
                     "initialization_JIT_s": init, "interval_total_s": interval, "overhead_s": overhead})
        previous = end
    if not rows or previous != max(r["step"] for r in observed):
        raise CostAuditFailure("Incomplete primary session timing coverage")
    return {"sessions": rows, "C_init_s": 1.25*max(r["initialization_JIT_s"] for r in rows),
            "C_session_s": 1.5*max(max(r["overhead_s"],0.) for r in rows),
            "W_prefix_actual_s": sum(r["wall_s"] for r in rows)}


def authorization_gate(complete, total, tests):
    return bool(np.isfinite([complete,total]).all() and 0 <= complete <= POLICY["full_L0_cap_s"] and
                0 <= total <= POLICY["total_cap_s"] and len(tests)==4 and all(r["passed"] for r in tests))


def administrative_bindings():
    return {p: sha(p) for p in ADMIN_FILES}


def freeze():
    from sloshing.multiphase.phase_rate_history import atomic_json
    f = numerical_guard(); prefix_guard()
    if (OUT/"policy.json").exists():
        raise ValueError("Cost policy already frozen; no rewrite")
    record = {**POLICY, "frozen_at_UTC": datetime.now(timezone.utc).isoformat(), **lineage(f),
              "prefix_manifest_sha256": sha(MANIFEST), "numerical_policy_sha256": f["policy_sha256"],
              "numerical_source_archive_sha256": f["source_archive_sha256"],
              "administrative_files_sha256": administrative_bindings()}
    atomic_json(OUT/"policy.json",record)
    # Plain digest file is a generated policy artifact, not a source edit.
    (OUT/"policy_sha256").write_text(sha(OUT/"policy.json")+"\n")
    print("COST POLICY FROZEN",sha(OUT/"policy.json"),flush=True)


def calibration_replay(required, scalar_rows):
    from sloshing.multiphase.phase_rate_history import archive, load_archive, Journal, atomic_json
    from sloshing.multiphase.linearized_ch import json_hash, array_hash
    f = numerical_guard()
    a, cp = load_archive(TRAJECTORY/"checkpoints/000127")
    old_phi = a["older"][a["physical_phi_map"]].copy(); dt = cp["current_dt"]
    fields = {k:a[k].copy() for k in ("u","pi","phi","mu","phase_rate")}
    fields.update(phi_old=old_phi, increment_from_rate=dt*a["phase_rate"],
                  increment_from_physical_states=a["phi"]-old_phi,
                  phase_rate_recovered_from_rounded_phi=(a["phi"]-old_phi)/dt)
    _, fm = load_archive(TRAJECTORY/"fields/000127")
    if {k:array_hash(v) for k,v in fields.items()} != fm["array_sha256"]:
        raise CostAuditFailure("Calibration field payload differs from frozen archive layout")
    strip=lambda m:{k:v for k,v in m.items() if k not in ("archive_version","array_sha256","semantic_sha256")}
    records={}; started=time.perf_counter()
    for cls in sorted(required):
        flag_f, flag_c, flag_a = [bool(int(s[1])) for s in cls.split("_")]
        scalar = scalar_rows[1] if flag_a else scalar_rows[-1]
        scalar = {k:v for k,v in scalar.items() if k not in ("record_sha256","previous_record_sha256")}
        description={"class":cls,"checkpoint_semantic_sha256":cp["semantic_sha256"],
                     "checkpoint_array_hashes":cp["array_sha256"],"field_array_hashes":fm["array_sha256"],
                     "scalar_payload_sha256":json_hash(scalar),"separate_residual_json":False,
                     "archive_function":"frozen phase_rate_history.archive; np.savez_compressed/fsync unchanged",
                     "checkpoint_guard":"new external numerical guard; old guard not modified"}
        times=[]; details=[]
        temp=Path(tempfile.mkdtemp(prefix=".calibration-",dir=OUT))
        try:
            for rep in range(5):
                folder=temp/f"replay-{rep}"; folder.mkdir()
                journal_path=folder/"scalar_history.jsonl"
                shutil.copyfile(TRAJECTORY/"scalar_history.jsonl",journal_path)  # setup outside measured write
                journal=Journal(journal_path)
                tick=time.perf_counter(); journal.append(scalar)
                if flag_f:
                    archive(folder/"fields/000127",fields,strip(fm))
                if flag_c:
                    numerical_guard()
                    meta={**strip(cp),"journal_prefix":journal.prefix()}
                    record=archive(folder/"checkpoints/000127",a,meta)
                    atomic_json(folder/"LATEST.json",{"checkpoint":"checkpoints/000127","checkpoint_sha256":record["semantic_sha256"]})
                elapsed=time.perf_counter()-tick; times.append(elapsed)
                details.append({"repetition":rep+1,"wall_s":elapsed,
                                "files":{str(p.relative_to(folder)):{"bytes":p.stat().st_size,"sha256":sha(p)}
                                         for p in sorted(folder.rglob("*")) if p.is_file()}})
            records[cls]={"times_s":times,"upper_s":1.5*max(times),"payload_description":description,
                          "payload_description_sha256":json_hash(description),"repetitions":details,
                          "filesystem_device":temp.stat().st_dev,"temporary_path":str(temp)}
            atomic_json(OUT/"calibration.json",{"records":records,"wall_s":time.perf_counter()-started,
                        "rule":"exactly 5 per required unseen class; 1.50*max; payload records saved before cleanup"})
        finally:
            # Only this explicit mkdtemp-created directory; never a scientific path.
            # On a failure keep it for diagnosis instead of deleting unrecorded timings.
            if cls in records:
                shutil.rmtree(temp)
    if not required:
        atomic_json(OUT/"calibration.json",{"records":{},"wall_s":0.,"status":"not_needed"})
    return records,time.perf_counter()-started


def run():
    from sloshing.multiphase.phase_rate_history import atomic_json, Journal
    from sloshing.multiphase.benchmarks import divergence_transfer as d
    f=numerical_guard(); prefix_guard(); policy=read(OUT/"policy.json")
    if any(policy[k]!=v for k,v in POLICY.items()) or policy["administrative_files_sha256"]!=administrative_bindings():
        raise CostAuditFailure("Frozen cost policy/source changed")
    if (OUT/"decision.json").exists() or (OUT/"event_schedule.json").exists():
        raise ValueError("Refusing a repeated cost experiment")
    schedule=d.schedules()[0]; ev=events(schedule); table=counts(ev)
    atomic_json(OUT/"event_schedule.json",{"schedule_sha256":schedule.sha256,"events":ev,"class_counts":table})
    print("EXACT EVENT COUNTS",json.dumps(table,sort_keys=True),flush=True)  # before estimates
    observed=decompose(Journal.read(TRAJECTORY/"timing_history.jsonl"),ev)
    future=ev[127:]; required=set()
    for train,target in [(observed,future)]+[(observed[:n],ev[a-1:b]) for n,a,b in POLICY["backtests"]]:
        missing_core={r["core_class"] for r in target}-{r["core_class"] for r in train}
        if missing_core: raise CostAuditFailure("Unseen future core class: "+str(missing_core))
        required |= {r["archive_class"] for r in target}-{r["archive_class"] for r in train}
    replays,calibration_wall=calibration_replay(required,Journal.read(TRAJECTORY/"scalar_history.jsonl"))
    calibration={k:v["times_s"] for k,v in replays.items()}
    tests=[]
    for n,a,b in POLICY["backtests"]:
        test=backtest(observed[:n],ev[a-1:b],calibration,sum(r["total_s"] for r in observed[a-1:b]))
        tests.append({"training_end":n,"validation_start":a,"validation_end":b,**test})
    atomic_json(OUT/"backtests.json",tests)
    prediction=predict(observed,future,calibration)
    all_estimates=predict(observed,ev,calibration)["estimates"]
    primary=[read(p) for p in sorted(TRAJECTORY.glob("sessions/*/session_status.json"))]
    allowances=session_allowances(primary,observed)
    W_total=d.used(); W_full=d.used("full_path")
    if W_full!=allowances["W_prefix_actual_s"]: raise CostAuditFailure("Primary wall mismatch")
    C_future=allowances["C_init_s"]+allowances["C_session_s"]+prediction["predicted_interval_s"]
    forecast={**prediction,**allowances,"future_continuation_s":C_future,
              "projected_complete_L0_s":W_full+C_future,
              "actual_historical_total_scientific_s":W_total,"calibration_administrative_wall_s":calibration_wall,
              "projected_total_scientific_s":W_total+calibration_wall+C_future,
              "future_steps":len(future),"full_L0_cap_s":14400.,"total_cap_s":19800.}
    passed=authorization_gate(forecast["projected_complete_L0_s"],forecast["projected_total_scientific_s"],tests)
    decision={"event_counts":table,"class_estimates":all_estimates,"backtests":tests,
              "forecast":forecast,"old_forecast":read(OLD/"policy/cost_before_full_L0.json"),
              "full_L0_cap_s":14400.,"total_cap_s":19800.,"authorized":passed,
              "policy_sha256":sha(OUT/"policy.json"),"prefix_manifest_sha256":sha(MANIFEST),
              "schedule_sha256":schedule.sha256,"timestamp_UTC":datetime.now(timezone.utc).isoformat(),
              "read_only_cost_audit":True,"numerical_guard_passed":True,**lineage(f)}
    atomic_json(OUT/"observed_classes.json",{"rows":observed,"class_counts":table,"estimates":all_estimates})
    atomic_json(OUT/"forecast.json",forecast); atomic_json(OUT/"decision.json",decision)
    auth={"authorized":passed,"cost_policy_sha256":sha(OUT/"policy.json"),"cost_decision_sha256":sha(OUT/"decision.json"),
          "prefix_manifest_sha256":sha(MANIFEST),"schedule_sha256":schedule.sha256,
          "STEP3A9_production_policy_sha256":f["policy_sha256"],"STEP3A9_source_archive_sha256":f["source_archive_sha256"],
          "external_continuation_controller_sha256":sha("scripts/step3a10_continue.py"),
          "administrative_files_sha256":administrative_bindings(),"timestamp_UTC":datetime.now(timezone.utc).isoformat(),
          **lineage(f)}
    atomic_json(ROOT/"continuation_authorization.json",auth)
    prefix_guard()
    print("COST DECISION",json.dumps({"authorized":passed,"forecast":forecast,
          "backtests":[{k:r[k] for k in ("validation_start","validation_end","predicted_interval_s","actual_interval_s","forecast_over_actual","passed")} for r in tests]},indent=2),flush=True)
    return 0 if passed else 2


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--stage",choices=("freeze","audit"),required=True)
    args=parser.parse_args()
    if args.stage=="freeze": freeze()
    else:
        try: raise SystemExit(run())
        except CostAuditFailure as exc:
            from sloshing.multiphase.phase_rate_history import atomic_json
            atomic_json(OUT/"FAILED.json",{"authorized":False,"error":str(exc),"read_only_cost_audit":True})
            raise
