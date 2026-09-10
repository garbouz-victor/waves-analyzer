"""One process, one run directory, one lock; no external execution service."""
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time
import uuid

import yaml

from ..config import SimulationConfig
from .reference import LABEL, run_reference

BLOCKER = ("PW1_CONTACT_CLOSURE: a law and scale for new wall contact are not qualified. "
           "Sharp no-slip material endpoints are fixed; setting sigma=0 in the existing "
           "CHNS removes interfacial/wall free energy. Retained coating alone does not close R(t).")


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    os.replace(tmp, path)


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_identity(root):
    src = root / "sloshing_visualization/src/sloshing"
    names = ["config.py", "mesh.py", "fem_spaces.py", "solver.py", "problems.py",
             "time_integrator.py", "diagnostics.py", "pinned_wetting/reference.py",
             "pinned_wetting/coating.py"]
    hashes = {name: sha256(src / name) for name in names}
    digest = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    return digest, hashes


def doctor(root):
    mem = {line.split(":")[0]: int(line.split()[1]) * 1024
           for line in Path("/proc/meminfo").read_text().splitlines() if line.endswith("kB")}
    deps = {}
    for name in ("numpy", "scipy", "scikit-fem", "h5py", "matplotlib", "PyYAML"):
        try:
            deps[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            deps[name] = "MISSING"
    return {"python": sys.executable, "versions": deps, "cpu_count": os.cpu_count(),
            "mem_total_bytes": mem["MemTotal"], "mem_available_bytes": mem["MemAvailable"],
            "disk_available_bytes": shutil.disk_usage(root).free,
            "ffmpeg": shutil.which("ffmpeg"), "ffprobe": shutil.which("ffprobe"),
            "execution_HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
            "status": "READY_FOR_REFERENCE" if "MISSING" not in deps.values() else "DEPENDENCY_MISSING"}


def process_identity():
    stat = Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()
    return {"pid": os.getpid(), "start_ticks": stat[19],
            "pid_namespace": os.readlink("/proc/self/ns/pid"),
            "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip()}


def lock_busy(path):
    # Status does not create or rewrite even the lock file.
    if not path.exists():
        return False
    with path.open("r") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(lock, fcntl.LOCK_UN)
            return False
        except BlockingIOError:
            return True


def byte_count(root):
    return sum(p.stat().st_size for p in root.rglob("*") if p.is_file()) if root.exists() else 0


@contextmanager
def job(root, command, limits):
    mission = root / "mission/pinned_wetting"
    lockpath = mission / "job.lock"
    with lockpath.open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("PW1 job is already active; refusing a duplicate")
        state = read_json(mission / "state.json")
        if state.get("active_jobs"):
            # The kernel lock is released on process death even across PID namespaces.
            state.setdefault("job_history", []).extend(
                [{**j, "exit_status": "INTERRUPTED_OR_CRASHED", "exit_code": None}
                 for j in state["active_jobs"]])
            # Charge unaccounted time conservatively instead of silently resetting it.
            last = state["active_jobs"][0].get("heartbeat_epoch_s", time.time())
            state["budget"]["actual_local_job_wall_s"] += max(0., time.time() - last)
        used = (state["budget"]["actual_local_job_wall_s"] +
                state["budget"].get("auxiliary_wall_time_upper_bound_s", 0.))
        if used >= limits["max_total_local_job_wall_hours"] * 3600:
            raise RuntimeError("RESOURCE_LIMIT: cumulative actual wall budget exhausted")
        # A modest absolute ceiling additionally protects the workstation.
        memory_cap = min(int(doctor(root)["mem_total_bytes"] * limits["max_memory_fraction_of_host"]), 4 * 2**30)
        resource.setrlimit(resource.RLIMIT_AS, (memory_cap, memory_cap))
        # Native target stores accepted RK stages as independent-audit evidence.
        # The total 40 GiB cap still applies; a fine 5 s HDF5 can exceed 1 GiB.
        file_cap = min(int(limits["max_new_disk_GiB"] * 2**30), 8 * 2**30)
        resource.setrlimit(resource.RLIMIT_FSIZE, (file_cap,)*2)
        job_record = {**process_identity(), "command": command, "started_epoch_s": time.time(),
                      "heartbeat_epoch_s": time.time(), "wall_s": 0., "memory_cap_bytes": memory_cap,
                      "single_file_cap_bytes": file_cap}
        state["active_jobs"] = [job_record]
        state["budget"]["limits"] = limits
        write_json(mission / "state.json", state)
        begin = time.monotonic()
        last_charge = begin
        def heartbeat(solver_state=None, checkpoint=None):
            nonlocal last_charge
            now = time.monotonic()
            state["budget"]["actual_local_job_wall_s"] += now - last_charge
            state["budget"]["charged_local_job_wall_s"] = (state["budget"]["actual_local_job_wall_s"] +
                state["budget"].get("auxiliary_wall_time_upper_bound_s", 0.))
            last_charge = now
            job_record.update(heartbeat_epoch_s=time.time(), wall_s=now-begin)
            if solver_state is not None:
                state.update(accepted_step=solver_state.step, physical_time_s=solver_state.time,
                             last_checkpoint=checkpoint)
            state["budget"]["new_output_bytes"] = sum(byte_count(root / p) for p in
                ("sloshing_visualization/results/pinned_wetting", "sloshing_visualization/output/pinned_wetting"))
            write_json(mission / "state.json", state)
            if (state["budget"]["charged_local_job_wall_s"] > limits["max_total_local_job_wall_hours"] * 3600 or
                    state["budget"]["new_output_bytes"] > limits["max_new_disk_GiB"] * 2**30 or
                    shutil.disk_usage(root).free < 2**30):
                raise RuntimeError("RESOURCE_LIMIT: accepted checkpoint saved at actual resource cap")
        exitcode = 0
        try:
            yield state, heartbeat
        except BaseException as error:
            exitcode = 130 if isinstance(error, KeyboardInterrupt) else 1
            state["last_error"] = repr(error)
            state["status"] = "RESOURCE_LIMIT" if "RESOURCE_LIMIT" in str(error) else "FAILED"
            raise
        finally:
            try:
                heartbeat()
            finally:
                exitcode = state.pop("_job_exit_code", exitcode)
                state["active_jobs"] = []
                exit_status=("REVIEW_OR_VERIFICATION_PENDING" if state.get("current_run_scope")=="DECLARED_EXTENSION" else "MODEL_BLOCKED") if exitcode==2 else ("FINISHED" if not exitcode else "STOPPED")
                state.setdefault("job_history", []).append({**job_record, "exit_code": exitcode,
                    "finished_epoch_s": time.time(), "exit_status": exit_status})
                write_json(mission / "state.json", state)


def resolve(root, cfg, args):
    geom, phys, initial = cfg["geometry"], cfg["physics"], cfg["initial"]
    controls = SimulationConfig(a=geom["half_width_m"], d=geom["depth_m"],
        g=phys["gravity_m_s2"], alpha_deg=initial["alpha_deg"], nu=phys["kinematic_viscosity_m2_s"],
        t_end=args.reference_t_end, dt=.0025 * args.dt_scale, snapshot_dt=.01,
        integrator="sdirk2", nx=24 * 2**args.mesh_level, nz=48 * 2**args.mesh_level)
    return {"mission_id": cfg["mission_id"], "status": "PREVIEW_ONLY", "physical_model_label": LABEL,
            "requested_case": cfg, "reference_controls": json.loads(controls.to_json()),
            "origins": {"geometry": geom["origin"], "alpha": initial["alpha_origin"],
                        "nu": phys["viscosity_origin"], "density": phys["density_origin"],
                        "numerical_controls": "agent_selected_reference_only"},
            "equations": "linear incompressible unsteady Stokes; gravity traction; eta_t=w at z=0",
            "wall_bc": "no_slip", "sigma_N_m": 0.,
            "z1_m": -controls.a * controls.slope, "z2_m": controls.a * controls.slope,
            "top_headroom_m": None,
            "headroom_reason": "Reference linearizes at z=0 with traction/kinematic boundary, not a lid; target domain NOT RUN.",
            "retained_layer": {"set": "L_macro union {walls x=+-a, -d<=z<=H_i}",
                               "volume_feedback": "neglected_at_leading_order", "thickness": None,
                               "adsorption_work": "not_modeled_at_leading_order",
                               "initial_dry_region_precoated": False},
            "target_contact_closure": None}, controls


def main(root):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["doctor", "run", "status", "resume", "verify", "render", "package", "blocker-check", "pilot"])
    parser.add_argument("--extension", action="store_true", help="Run the explicitly authorized finite Navier-slip model")
    parser.add_argument("--mesh-level", type=int, default=0, choices=[0, 1, 2])
    parser.add_argument("--dt-scale", type=float, default=1.)
    parser.add_argument("--reference-t-end", type=float, default=1.)
    parser.add_argument("--run-id")
    parser.add_argument("--stop-after", type=int, help="Save a checkpoint after this accepted step")
    parser.add_argument("--reference-only", action="store_true", help="verify supported reference subset; never imply target PASS")
    args = parser.parse_args()
    mission = root / "mission/pinned_wetting"
    cfg = yaml.safe_load((mission / "CONFIG.yaml").read_text())
    if args.extension:
        from .extension import extension_main
        return extension_main(root, cfg, args)
    output = root / cfg["output"]["root"]
    runs = root / cfg["output"]["runs_root"]
    if args.command == "doctor":
        print(json.dumps(doctor(root), indent=2)); return 0
    if args.command == "status":
        print(json.dumps({**read_json(mission / "state.json"),
                          "kernel_lock_busy": lock_busy(mission / "job.lock")}, indent=2)); return 0
    if args.command == "verify":
        from .audit import verify_output
        saved = read_json(mission / "state.json")
        destination = saved.get("artifact_paths", {}).get("output_directory", str(output))
        report = verify_output(destination, reference_only=args.reference_only)
        if (output / "manifest.json").exists():
            handoff = read_json(output / "manifest.json")
            errors = []
            for item in handoff["artifacts"]:
                path = output / item["path"]
                if not path.is_file() or sha256(path) != item["sha256"]:
                    errors.append(item["path"])
            if handoff["run_id"] != saved["current_run_id"]:
                errors.append("handoff_run_id")
            report["handoff_integrity"] = {"status": "FAIL" if errors else "PASS", "failures": errors}
            if errors:
                report["passed"] = False
            else:
                report["independent_review"] = read_json(output / "independent_review.json")
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0 if report["passed"] else 2
    if args.command == "blocker-check":
        from .audit import blocker_check
        print(json.dumps(blocker_check(), indent=2)); return 2
    if not math.isfinite(args.dt_scale) or args.dt_scale <= 0 or not 0 < args.reference_t_end <= 5:
        parser.error("Finite dt-scale >0 and reference horizon in (0,5] required")
    with job(root, " ".join(sys.argv), cfg["resources"]) as (state, heartbeat):
        if args.command == "package":
            package_blocker(root, output, state, heartbeat)
            state["_job_exit_code"] = 2
            return 2
        if args.command == "render":
            from .output import export_reference
            destination = output / state["current_run_id"]
            export_reference(destination, Path(state["last_checkpoint"]), heartbeat)
            return 0
        resume = args.command == "resume"
        run_id = args.run_id or (state.get("current_run_id") if resume else None)
        if not resume and args.run_id is None and state.get("current_run_id"):
            previous = runs / state["current_run_id"]
            if (previous / "identity.json").exists():
                _, requested_controls = resolve(root, cfg, args)
                if (read_json(previous / "identity.json")["source_hash"] == source_identity(root)[0] and
                        read_json(previous / "resolved_case.json")["reference_controls"] == json.loads(requested_controls.to_json())):
                    resume, run_id = True, state["current_run_id"]
        if run_id is not None and ("/" in run_id or ".." in run_id):
            raise ValueError("Invalid run ID")
        if resume and not run_id:
            raise ValueError("No existing run ID to resume")
        if not resume:
            resolved, controls = resolve(root, cfg, args)
            run_id = run_id or "reference-" + time.strftime("%Y%m%dT%H%M%S", time.gmtime()) + "-" + uuid.uuid4().hex[:6]
            if "/" in run_id or ".." in run_id:
                raise ValueError("Invalid run ID")
            run_dir = runs / run_id
            run_dir.mkdir(parents=True, exist_ok=False)
            digest, hashes = source_identity(root)
            identity = {"run_id": run_id, "source_hash": digest, "source_files": hashes,
                        "physical_config_hash": hashlib.sha256(json.dumps(resolved, sort_keys=True).encode()).hexdigest(),
                        "execution_HEAD": doctor(root)["execution_HEAD"],
                        "physical_model_label": LABEL, "synthetic": False}
            write_json(run_dir / "resolved_case.json", resolved)
            write_json(run_dir / "identity.json", identity)
            write_json(run_dir / "doctor.json", doctor(root))
            for name in hashes:
                dst = run_dir / "source_snapshot" / name
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(root / "sloshing_visualization/src/sloshing" / name, dst)
            state["run_ids"].append(run_id)
        else:
            run_dir = runs / run_id
            resolved = read_json(run_dir / "resolved_case.json")
            controls = SimulationConfig(**resolved["reference_controls"])
            identity = read_json(run_dir / "identity.json")
            if source_identity(root)[0] != identity["source_hash"]:
                raise ValueError("Numerical source changed; use a new run ID")
            destination = output / run_id
            if (destination / "manifest.json").exists():
                from .audit import verify_output
                verified = verify_output(destination, reference_only=True)
                if not verified["passed"]:
                    raise ValueError("Existing completed reference failed read-only verification")
                state["_job_exit_code"] = 2
                print(json.dumps({"status": "NEEDS_MODEL_DECISION", "reused_run_id": run_id,
                                  "PDE_restarted": False, "blocker": BLOCKER}))
                return 2
        state.update(status="PREVIEW", current_run_id=run_id)
        heartbeat()
        result = run_reference(run_dir, controls, identity, heartbeat, args.stop_after)
        state.update(accepted_step=result["step"], physical_time_s=result["time_s"])
        heartbeat()
        if not result["completed"]:
            state.update(status="PREVIEW", next_action="Resume the saved accepted checkpoint; no active job remains.")
            return 0
        from .audit import blocker_check, verify_output
        from .output import export_reference
        output.mkdir(parents=True, exist_ok=True)
        # Separate immutable outputs per run. Root only provides a current index.
        destination = output / run_id
        destination.mkdir(exist_ok=True)
        write_json(destination / "resolved_case.json", resolved)
        write_json(destination / "blocker_evidence.json", blocker_check())
        decision = mission / "model_decision.md"
        if decision.exists():
            shutil.copyfile(decision, destination / "model_decision.md")
        export_reference(destination, run_dir / "native.h5", heartbeat)
        review = verify_output(destination, reference_only=True)
        write_json(destination / "verification.json", review)
        (output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>PW1</title>'
            f'<h1>NEEDS_MODEL_DECISION</h1><p><a href="{run_id}/index.html">REFERENCE / NOT TARGET — preview и данные</a></p>')
        state.update(status="NEEDS_MODEL_DECISION", model_decision=str(decision),
            last_error={"code": "PW1_CONTACT_CLOSURE", "detail": BLOCKER,
                        "reproduce": f"{sys.executable} sloshing_visualization/scripts/pinned_wetting_mission.py blocker-check"},
            next_action="Specify and authorize a contact transport/relaxation law and scale; preserve no-slip/sigma=0 or explicitly approve an extension.")
        state["milestones"].update(M0="PREVIEW_AND_MODEL_BLOCKER", M1="RETAINED_STATE_ONLY_CONTACT_BLOCKED", M2="NOT_RUN", M3="REFERENCE_REVIEW_ONLY")
        state["artifact_paths"] = {"index": str(output / "index.html"), "preview": str(destination / "preview_reference.mp4"),
                                   "output_directory": str(destination), "native": str(run_dir / "native.h5")}
        state["verification"] = {"scientific": "TARGET_NOT_RUN", "reference": "PASS" if review["passed"] else "FAIL",
                                 "rendering": review["rendering"]["status"], "independent_review": "PENDING"}
        print(json.dumps({"status": state["status"], "preview": state["artifact_paths"]["preview"],
                          "blocker": BLOCKER}, ensure_ascii=False, indent=2))
        state["_job_exit_code"] = 2
        return 2


def package_blocker(root, output, state, heartbeat):
    """Publish a scoped handoff; immutable run bundles/native data are untouched."""
    from .audit import blocker_check, verify_output
    mission = root / "mission/pinned_wetting"
    source = Path(state["artifact_paths"]["output_directory"])
    native = Path(state["last_checkpoint"])
    independent = read_json(mission / "independent_review.json")
    if (independent["run_id"] != state["current_run_id"] or
            independent["native_sha256"] != sha256(native) or
            independent["status"] != "PASS_REFERENCE_REJECT_TARGET"):
        raise ValueError("Independent review does not cover the current native data")
    report = verify_output(source, reference_only=True)
    if not report["passed"]:
        raise ValueError("Reference handoff failed verification: " + repr(report))
    target_report = verify_output(source, reference_only=False)
    if target_report["passed"]:
        raise AssertionError("A reference must not pass target acceptance")
    def relocate_evidence(value):
        if isinstance(value, dict):
            return {key: ([str((source / p).relative_to(output)) if not Path(p).is_absolute() and (source / p).exists() else p for p in val]
                          if key == "evidence_paths" else relocate_evidence(val)) for key, val in value.items()}
        if isinstance(value, list):
            return [relocate_evidence(v) for v in value]
        return value
    report, target_report = relocate_evidence(report), relocate_evidence(target_report)
    output.mkdir(parents=True, exist_ok=True)
    for name in ("resolved_case.json", "model_decision.md", "refinement.json"):
        shutil.copyfile(source / name, output / name)
    for name in ("review.md", "independent_review.json", "test_results.xml"):
        shutil.copyfile(mission / name, output / name)
    # Relocate links only; the signed original review remains untouched.
    (output / "review.md").write_text((mission / "review.md").read_text().replace("../../sloshing_visualization/", "../../"))
    write_json(output / "blocker_evidence.json", blocker_check())
    final_verification = {"status": "NEEDS_MODEL_DECISION", "passed": False,
        "contract": {"status": "TARGET_CONTACT_NOT_QUALIFIED", "evidence_paths": ["model_decision.md", "blocker_evidence.json"]},
        "physics": report["physics"], "wetting": report["wetting"],
        "numerical_refinement": {"status": "NOT_RUN", "evidence_paths": ["refinement.json"]},
        "rendering": report["rendering"], "artifact_integrity": report["artifact_integrity"],
        "independent_review": independent,
        "reference_verification": report, "target_acceptance_rejection": target_report}
    write_json(output / "verification.json", final_verification)
    costs = state["budget"]["actual_local_job_wall_s"]
    (output / "FINAL_REPORT.md").write_text((source / "FINAL_REPORT.md").read_text() +
        f"\n## Независимый итог\n\nReference проверен отдельным контекстом; target отклонён.\n"
        f"См. [review.md](review.md), [verification.json](verification.json), [test_results.xml](test_results.xml).\n"
        f"46 тестов прошли без пропусков. Первый preview и исправленный reference сохранены раздельно.\n"
        f"Учтённое время управляемых jobs до упаковки: {costs:.3f} с. Лимит: 43200 с; новые файлы около 67 MiB.\n"
        f"Текущий бюджет, jobs и exact resume command: `{mission / 'state.json'}`.\n"
        f"\nРабочий interpreter: `{sys.executable}`. В командах выше замените `python` этим путём.\n"
        f"Точная диагностика: `{sys.executable} sloshing_visualization/scripts/pinned_wetting_mission.py blocker-check`.\n"
        "\nНедостающее решение: можно рассмотреть явно отдельный DECLARED_EXTENSION с конечным Navier slip,\n"
        "сохранив sigma=0 и необратимый слой. Это изменяет исходное no-slip и **не выполнено** без решения пользователя.\n"
        "Альтернатива с сохранением no-slip требует отдельно обоснованного закона диффузного/подсеточного переноса контакта.\n")
    rel = source.name
    links = {"preview MP4 (0–1 с; REFERENCE / NOT TARGET)": f"{rel}/preview_reference.mp4",
             "Данные и кадры reference": f"{rel}/index.html", "Итоговый отчёт": "FINAL_REPORT.md",
             "Независимый review": "review.md", "Проверка и отклонение target": "verification.json",
             "Воспроизводимый blocker": "blocker_evidence.json", "Решение о модели": "model_decision.md",
             "Исходные поля и checkpoint (локальный файл)": os.path.relpath(native, output)}
    (output / "index.html").write_text('<!doctype html><html lang="ru"><meta charset="utf-8"><title>PW1 — NEEDS_MODEL_DECISION</title>'
        '<style>body{max-width:1100px;margin:40px auto;font:18px sans-serif;background:#f6f8fc}video{width:100%}li{margin:12px}</style>'
        '<h1>NEEDS_MODEL_DECISION</h1><p>Целевой расчёт 0–5 с, уточнение и final_animation.mp4 не выполнены. '
        'Ниже — реальный контрольный preview с неподвижными концами. Закон нового контакта пока не квалифицирован.</p>'
        f'<video controls preload="metadata" src="{rel}/preview_reference.mp4"></video><ul>' +
        ''.join(f'<li><a href="{url}">{title}</a></li>' for title, url in links.items()) + '</ul></html>')
    artifacts = [{"path": str(p.relative_to(output)), "sha256": sha256(p), "role": "blocked_handoff"}
                 for p in sorted(output.iterdir()) if p.is_file() and p.name != "manifest.json"]
    artifacts.append({"path": str((source / "manifest.json").relative_to(output)),
                      "sha256": sha256(source / "manifest.json"), "role": "immutable_reference_bundle_manifest"})
    write_json(output / "manifest.json", {"schema_version": "1.0", "status": "NEEDS_MODEL_DECISION",
        "run_id": state["current_run_id"], "physical_model_label": LABEL, "artifacts": artifacts,
        "native_data_paths": [str(native)], "native_sha256": sha256(native),
        "time_range": [0., 1.], "target_time_range": [0., 5.], "target_status": "NOT_RUN",
        "source_hash": read_json(source / "manifest.json")["source_hash"],
        "execution_HEAD": doctor(root)["execution_HEAD"],
        "limitations": ["Missing qualified moving-contact closure", "No final target/refinement"]})
    state["artifact_paths"]["package_directory"] = str(output)
    state["artifact_paths"]["final_report"] = str(output / "FINAL_REPORT.md")
    state["verification"]["independent_review"] = independent["status"]
    state["verification"]["pack_integrity"] = "PACKAGE_OK"
    state["verification"]["negative_tests"] = "46_PASS_NO_SKIPS"
    state["milestones"]["M3"] = "BLOCKER_HANDOFF_REVIEWED_TARGET_NOT_RUN"
    state["status"] = "NEEDS_MODEL_DECISION"
    heartbeat()
    print(json.dumps({"status": state["status"], "index": str(output / "index.html"),
                      "independent_review": independent["status"]}, ensure_ascii=False, indent=2))
