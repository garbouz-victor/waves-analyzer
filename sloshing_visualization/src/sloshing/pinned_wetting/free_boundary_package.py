"""Small PW2 release adapter: native evidence first, movies only after qualification."""
import json
import csv
import itertools
import re
from pathlib import Path
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import h5py
import numpy as np

from .mission import doctor, read_json, write_json, sha256
from .free_boundary_run import REL_RUNS, REL_OUTPUT
from .free_boundary_output import load, export, render
from .free_boundary_verify import MODEL_ID, verify_native
from .free_boundary_compare import compare_runs


def tested_source_hashes(root):
    files=[]
    for directory in ("sloshing_visualization/src", "sloshing_visualization/tests",
                      "mission/pinned_wetting/solve_to_animation_v2/tools", "mission/pinned_wetting/solve_to_animation_v2/tests"):
        files.extend((root/directory).rglob("*.py"))
    files.extend(root/"sloshing_visualization"/name for name in ("requirements.txt", "pyproject.toml"))
    return {p.relative_to(root).as_posix():sha256(p) for p in sorted(files)}


def run_complete_tests(root, output, heartbeat):
    """Run within the existing managed job, not a second scheduler."""
    output.mkdir(parents=True,exist_ok=True)
    # Override the project's quick-suite default so validation-marked tests
    # are included in the release evidence rather than silently deselected.
    command=[sys.executable,"-m","pytest","-q","-o","addopts=-ra","tests",
        "../mission/pinned_wetting/solve_to_animation_v2/tests",f"--junitxml={output/'tests.xml'}"]
    before=tested_source_hashes(root)
    with (output/"pytest.log").open("w") as log:
        process=subprocess.Popen(command,cwd=root/"sloshing_visualization",stdout=log,stderr=subprocess.STDOUT)
        try:
            while process.poll() is None:
                heartbeat();time.sleep(2.)
        finally:
            if process.poll() is None:
                process.terminate()
                try:process.wait(timeout=10.)
                except subprocess.TimeoutExpired:
                    process.kill();process.wait()
    record={"command":command,"cwd":str(root/"sloshing_visualization"),"exit_code":process.returncode,
        "full_project_scope":True,"source_hashes":before,"source_unchanged_during_test":before==tested_source_hashes(root),
        "execution_HEAD":doctor(root)["execution_HEAD"],"tests_xml_sha256":sha256(output/"tests.xml")}
    write_json(output/"pytest_evidence.json",record)
    return record


def test_evidence(path,root):
    if not path.exists():
        return {"passed": False, "reason": "complete local pytest XML missing"}
    tree = ET.parse(path).getroot()
    suites = list(tree.iter("testsuite"))
    totals = {key: sum(int(s.get(key, 0)) for s in suites)
              for key in ("tests", "failures", "errors", "skipped")}
    provenance=read_json(path.parent/"pytest_evidence.json") if (path.parent/"pytest_evidence.json").exists() else {}
    current=(provenance.get("source_hashes")==tested_source_hashes(root)
        and provenance.get("source_unchanged_during_test") is True
        and provenance.get("full_project_scope") is True and provenance.get("exit_code")==0
        and provenance.get("tests_xml_sha256")==sha256(path))
    return {**totals, "passed": bool(current and totals["tests"] > 0 and totals["failures"] == totals["errors"] == 0),
            "current_complete_test_provenance":current,
            "path": str(path), "sha256": sha256(path), "remote_CI_claimed": False}


def _guard(root, output, name):
    script = root/"mission/pinned_wetting/solve_to_animation_v2/tools"/name
    target = output/"resolved_case.json" if name == "contract_guard.py" else output
    process = subprocess.run([sys.executable, str(script), str(target)], capture_output=True, text=True)
    result = json.loads(process.stdout)
    return {**result, "exit_code": process.returncode, "checker_sha256": sha256(script)}


def _cached_native(native, heartbeat):
    dest = native.parent/"native_verification.json"
    from . import free_boundary_verify as verifier
    old = read_json(dest) if dest.exists() else {}
    if not (old.get("native_sha256") == sha256(native)
            and old.get("verifier_sha256") == sha256(Path(verifier.__file__))
            and old.get("verifier_dependencies")==verifier.verifier_dependencies()
            and old.get("require_target") is True):
        old = verify_native(native, require_target=True, heartbeat=heartbeat)
        write_json(dest, old)
    return old


def frame_map_matches_native(mapping,times):
    if len(mapping)!=501:return False
    for index,frame in enumerate(mapping):
        step=frame.get("accepted_step",-1)
        if not (frame.get("frame")==index and abs(frame.get("time_s",-1)-index*.01)<1e-12
                and type(step) is int and 0<=step<len(times) and abs(times[step]-index*.01)<1e-10):
            return False
    return True


def verify_full_curve_export(native,exported):
    """Independently reconstruct every CSV curve directly from native P2 DOFs.

    No renderer/exporter sampling helper is used. Merely preserving endpoints
    is insufficient: all ordered points of the thin/overhanging branch count.
    """
    errors=[];maximum=0.;checked=0;seen=[]
    with h5py.File(native,"r") as h,Path(exported).open(newline="") as stream:
        ids=h["topology/interface_nodes"][:]
        times={round(float(group.attrs["time_s"]),12):key for key,group in h["states"].items()}
        for time_text,items in itertools.groupby(csv.DictReader(stream),key=lambda row:row["time_s"]):
            rows=list(items);time_s=float(time_text);key=times.get(round(time_s,12))
            if key is None or (seen and time_s<=seen[-1]):
                errors.append("unordered_or_non_native_export_time");continue
            seen.append(time_s)
            points=h["states"][key]["geometry"][:][:,ids].T
            pieces=[];q=np.arange(32)[:,None]/32.
            for edge in range(0,len(points)-2,2):
                left,mid,right=points[edge:edge+3]
                pieces.append(left+(4*mid-3*left-right)*q+2*(left+right-2*mid)*q*q)
            expected=np.vstack(pieces+[points[-1:]])
            if len(rows)!=len(expected) or [int(row["vertex_id"]) for row in rows]!=list(range(len(expected))):
                errors.append("full_parametric_branch_missing_or_reordered");continue
            actual=np.array([[float(row["x_m"]),float(row["z_m"])] for row in rows])
            defect=float(np.max(abs(actual-expected)))
            if not np.isfinite(defect) or defect>1e-12:
                errors.append("export_geometry_differs_from_full_native_curve")
            if np.isfinite(defect):maximum=max(maximum,defect)
            checked+=1
    return {"passed":not errors and checked>0,"failures":sorted(set(errors)),
        "checked_curves":checked,"maximum_coordinate_error_m":maximum,
        "native_sha256":sha256(native),"export_sha256":sha256(exported),
        "method":"Independent quadratic coefficient reconstruction, all ordered samples; no eta(x) reduction"}


def further_refinement_evidence(root, output, selected, heartbeat):
    """A documented command must have made an actual finer native step."""
    path = output/"further_refinement_command.json"
    if not path.exists():
        return {"passed":False,"reason":"tested further-refinement command missing"}
    record = read_json(path)
    try:
        native = Path(record["native"])
        command = record["command"]
        if (record.get("exit_code") != 0 or not isinstance(command,list)
                or not all(isinstance(part,str) for part in command)
                or "--free-boundary" not in command or "--stop-after" not in command):
            raise ValueError("refinement smoke command did not complete")
        with h5py.File(selected,"r") as first, h5py.File(native,"r") as second:
            one = json.loads(first.attrs["numerical_identity"])
            two = json.loads(second.attrs["numerical_identity"])
            mesh_keys = {"nx","nz","local_right_levels","x_grading","z_grading"}
            same = all(one[key] == two[key] for key in ("model_id","physical_contract","source_hashes","libraries","method"))
            same &= all(value == two["controls"].get(key) for key,value in one["controls"].items()
                        if key not in mesh_keys|{"dt"})
            finer = (second["topology/triangles"].shape[1] > first["topology/triangles"].shape[1]
                     and two["controls"]["dt"] <= one["controls"]["dt"])
            stepped = int(second.attrs["last_accepted_step"]) >= 1
        if not same or not finer or not stepped or record.get("native_sha256") != sha256(native):
            raise ValueError("refinement smoke is stale, incompatible, or not actually finer")
        audit = verify_native(native, require_target=False, heartbeat=heartbeat)
        write_json(output/"further_refinement_native_verification.json",audit)
        return {**record,"passed":bool(audit["passed"]),
                "scope":"CLI smoke only; not an additional 5s convergence result",
                "audit":"further_refinement_native_verification.json"}
    except (KeyError,ValueError,OSError) as error:
        return {"passed":False,"reason":str(error)}


def package(root, args, state, heartbeat, *, rerender=True):
    output = root/REL_OUTPUT
    output.mkdir(parents=True, exist_ok=True)
    runs = state.get("free_boundary_runs", {})
    selected = args.run_id or runs.get("spacefine", {}).get("run_id") or state["current_run_id"]
    native = root/REL_RUNS/selected/"native.h5"
    data = load(native)
    export(native, output, data)
    audits = {selected: _cached_native(native, heartbeat)}
    comparisons = {}
    previous_comparisons=read_json(output/"refinement.json") if (output/"refinement.json").exists() else {}
    missing = []
    if selected != runs.get("spacefine", {}).get("run_id"):
        missing.append("selected_animation_run_is_not_the_spatially_refined_case")
    for kind, roles in (("temporal", ("target", "timefine")), ("spatial", ("timefine", "spacefine"))):
        if any(role not in runs for role in roles):
            missing.append(kind+"_runs_missing")
            continue
        paths = [root/REL_RUNS/runs[role]["run_id"]/"native.h5" for role in roles]
        for path in paths:
            audits.setdefault(path.parent.name, _cached_native(path, heartbeat))
        try:
            from .free_boundary_compare import comparison_dependencies
            previous=previous_comparisons.get(kind,{})
            if (previous.get("first_sha256")==sha256(paths[0]) and previous.get("second_sha256")==sha256(paths[1])
                    and previous.get("comparison_dependencies")==comparison_dependencies()
                    and previous.get("time_end_s")==5. and previous.get("kind")==kind
                    and previous.get("common_time_policy")=="union_of_all_accepted_times"):
                comparisons[kind]=previous
            else:
                comparisons[kind] = compare_runs(*paths, kind=kind, heartbeat=heartbeat)
        except (ValueError, OSError) as error:
            comparisons[kind] = {"passed": False, "error": str(error)}
    write_json(output/"refinement.json", comparisons)
    write_json(output/"native_verification.json", audits)
    write_json(output/"discretization_diagnostics.json",{
        "note":"Diagnostic values, not extra physical-error gates; complete h/dt geometry and raw-energy gates remain mandatory.",
        "native_maxima":{run_id:audit.get("maxima",{}) for run_id,audit in audits.items()}})
    guards = {name: _guard(root, output, name) for name in ("contract_guard.py", "check_exports.py")}
    guards["full_native_curve_export"]=verify_full_curve_export(native,output/"interface_samples.csv")
    write_json(output/"export_verification.json", guards)
    tests = test_evidence(output/"tests.xml",root)
    refinement_command = further_refinement_evidence(root,output,native,heartbeat)
    native_passed = bool(audits) and all(a.get("passed") and a.get("actual_5s_trajectory") for a in audits.values())
    refinement_passed = all(comparisons.get(k, {}).get("passed") for k in ("temporal", "spatial"))
    numerical_candidate = native_passed and refinement_passed and not missing
    rendered = read_json(output/"render_evidence.json") if (output/"render_evidence.json").exists() else {}
    from . import free_boundary_output as renderer
    renderer_hash=sha256(Path(renderer.__file__))
    if numerical_candidate and rerender and not (rendered.get("native_sha256") == sha256(native)
            and rendered.get("renderer_sha256")==renderer_hash and rendered.get("final_candidate")):
        rendered = render(native, output, heartbeat, final=True)
    review = read_json(output/"independent_review.json") if (output/"independent_review.json").exists() else {}
    review_current = (review.get("native_sha256") == sha256(native)
                      and (output/"render_evidence.json").exists()
                      and review.get("render_evidence_sha256") == sha256(output/"render_evidence.json")
                      and review.get("refinement_sha256") == sha256(output/"refinement.json")
                      and review.get("native_verification_sha256") == sha256(output/"native_verification.json")
                      and bool(review.get("reviewer")))
    video_passed = bool(rendered.get("final_candidate") and rendered.get("native_sha256") == sha256(native)
                        and rendered.get("renderer_sha256")==renderer_hash and len(rendered.get("videos", [])) == 2)
    if video_passed:
        for item in rendered["videos"]:
            video = output/item["video"]
            mapping=item.get("frame_map",[])
            map_valid=frame_map_matches_native(mapping,data["times"])
            decode = subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-i", str(video), "-f", "null", "-"],
                                    capture_output=True, timeout=120)
            video_passed &= (decode.returncode == 0 and item.get("full_decode") is True
                             and map_valid and video.exists() and item.get("sha256")==sha256(video))
    result = {
        "model_id": MODEL_ID, "run_id": selected, "native_sha256": sha256(native),
        "contract_passed": guards["contract_guard.py"].get("passed", False),
        "full_native_curve_export_passed":guards["full_native_curve_export"]["passed"],
        "independent_native_equations_passed": native_passed,
        "boundary_conditions_passed": native_passed,
        "geometry_resolution_passed": bool(numerical_candidate and review_current and review.get("geometry_resolution_passed")),
        "mass_energy_passed": native_passed,
        "spatial_refinement_passed": bool(comparisons.get("spatial", {}).get("passed")),
        "temporal_refinement_passed": bool(comparisons.get("temporal", {}).get("passed")),
        "retained_coating_passed": native_passed and guards["check_exports.py"].get("passed", False),
        "video_decode_passed": bool(video_passed),
        "independent_visual_review_passed": bool(review_current and review.get("independent_visual_review_passed")),
        "local_tests_passed": tests["passed"], "tests": tests, "missing": missing,
        "further_refinement_command_passed":refinement_command["passed"],
        "further_refinement_command":refinement_command,
        "evidence": {"native": "native_verification.json", "refinement": "refinement.json",
                     "exports": "export_verification.json", "review": "independent_review.json", "video": "render_evidence.json"},
        "verification_code_sha256": sha256(Path(__file__)), "continuous_model_error_bound_claimed": False}
    result["passed"] = all(value is True for key, value in result.items() if key.endswith("_passed"))
    result["status"] = "COMPLETE_FOR_DECLARED_FREE_BOUNDARY_MODEL" if result["passed"] else "INCOMPLETE_NOT_FINAL"
    write_json(output/"verification.json", result)
    worktree = {"HEAD":doctor(root)["execution_HEAD"],
        "branch":subprocess.check_output(["git","branch","--show-current"],cwd=root,text=True).strip(),
        "porcelain_status":subprocess.check_output(["git","status","--porcelain=v1"],cwd=root,text=True).splitlines(),
        "commit_performed_by_package":False,"push_performed":False}
    write_json(output/"worktree_provenance.json",worktree)
    write_json(output/"runtime.json", {"budget_snapshot": state["budget"], "environment": doctor(root),
        "worktree":worktree,
        "note": "Current package wall interval is closed in mission/state.json after this snapshot.", "new_paid_compute": False})
    status = result["status"]
    p2_error=float(max(abs(data["surface"][:,1,-1]-data["surface"][0,1,-1])))
    controls=read_json(output/"resolved_case.json")["numerics"]
    runup = read_json(output/"runup_summary.json")
    runup["qualified_final"] = result["passed"]
    write_json(output/"runup_summary.json", runup)
    diagnostic_names=("momentum_residual","weak_divergence","accepted_endpoint_weak_divergence",
        "stage_strong_divergence_Linf_s_inv","stage_strong_divergence_L2","maximum_step_dilation",
        "volume_relative","continuous_energy_relative","split_energy_relative",
        "rezone_L2_velocity_projection_relative","rezone_surface_velocity_change_m_s",
        "surface_traction_L2_per_density","surface_augmented_traction_L2_per_density",
        "surface_grad_div_traction_L2_per_density","left_traction_L2_per_density")
    diagnostic_summary={run_id:{name:audit.get("maxima",{}).get(name) for name in diagnostic_names}
                        for run_id,audit in audits.items()}
    report = [f"# PW2 — {status}", "", "LEFT WALL: Navier-slip, b=0.50 m. RIGHT WALL: no-slip. BOTTOM: no-slip. sigma=0.",
        "Full nonlinear incompressible Navier–Stokes on the moving liquid domain; not the old fixed-domain linear result.",
        ("P2 remains the actual current right contact point at the same height for the entire saved trajectory."
         if p2_error<=1e-12 else "WARNING: saved trajectory fails the required fixed-P2 invariant."),
        f"max_t |eta_right(t)-z2| = {p2_error:.17g} m.", "",
        f"Native interval: 0–{data['times'][-1]:.12g} s. Selected run: {selected}.",
        f"Global mesh controls: {controls['nx']}×{controls['nz']}, local RIGHT refinement level {controls['local_right_levels']}; dt cap {controls['dt']} s. Actual topology is native evidence.",
        f"Numerical grad-div gamma={controls.get('grad_div_gamma_m2_s',0.):g} m²/s, fixed in each dt/h comparison; not physical viscosity. Its separately saved signed work is numerical loss, not physical heat. The raw physical continuous-energy gate is not compensated by this work.",
        f"HEAD at packaging: {doctor(root)['execution_HEAD']}; source snapshots and numerical identity are stored with native data.",
        f"Worktree branch: {worktree['branch']}; uncommitted paths present: {bool(worktree['porcelain_status'])}. Exact status: worktree_provenance.json. No automatic commit or push.",
        "", "Retained coating is a separate irreversible one-way subgrid state. Exact thickness, volume, inertia, separate dynamics and feedback are not calculated. Navier slip does not itself cause this memory.",
        "Any resolved thin liquid region remains in the computational domain, full interface and mass/energy balances.",
        "Signed geometric/remapping/grad-div work is numerical, not physical heat. Two-level comparisons do not establish formal convergence order or a physical-error bound.", "",
        "With grad-div enabled the finite-h natural weak stress includes gamma div(v) I. Physical free-surface traction, its numerical normal addition, and augmented traction are diagnosed separately; augmented weak balance is not a pointwise physical traction certificate. This consistent numerical term vanishes for exactly incompressible fields and is not a new physical boundary law.", "",
        "Incompressibility and natural traction are enforced weakly in the finite-element space. A small mixed-system residual does not imply zero pointwise divergence or zero pointwise traction at finite resolution. The diagnostics below include the full wall-adjacent region, accepted-end divergence and actual remap velocity changes; none is silently equated with a rigorous physical-error estimate.", "",
        "## Verification", "", "```json", json.dumps(result, indent=2), "```", "",
        "## Discretization and transfer diagnostics", "", "```json", json.dumps(diagnostic_summary,indent=2), "```", "",
        "## Runup and historical markers", "", "```json", (output/"runup_summary.json").read_text().strip(), "```", "",
        "## Complete boundary refinement", "", "```json", json.dumps(comparisons, indent=2), "```", "",
        "## Tested further accuracy command", "", "```json", json.dumps(refinement_command,indent=2), "```", "",
        "Native fields and all accepted coating states: native_manifest.json. Review: review.md. All old reference/FAIL bundles remain separate."]
    (output/"FINAL_REPORT.md").write_text("\n".join(report)+"\n")
    if (output/"index.html").exists():
        page = (output/"index.html").read_text()
        page=re.sub(r'<p id="qualification_status">.*?</p>',f'<p id="qualification_status">{status}</p>',page)
        if 'href="FINAL_REPORT.md"' not in page:
            page=page.replace('</body>','').replace('</html>','')+'<h2>Проверка и данные</h2><ul>'+''.join(
                f'<li><a href="{name}">{name}</a></li>' for name in
                ('FINAL_REPORT.md','verification.json','refinement.json','review.md','native_manifest.json','contacts.csv','marker_catalog.csv'))+'</ul></html>'
        (output/"index.html").write_text(page)
    write_json(output/"manifest.json", {"model_id": MODEL_ID, "status": status, "run_id": selected,
        "native_runs": [{"path": str(root/REL_RUNS/run_id/"native.h5"), "sha256": audit.get("native_sha256")}
                        for run_id, audit in audits.items()],
        "artifacts": [{"path": p.relative_to(output).as_posix(), "sha256": sha256(p)} for p in sorted(output.iterdir())
                      if p.is_file() and p.name != "manifest.json"], "push_performed": False})
    state.update(status=status, free_boundary_verification=str(output/"verification.json"),
                 current_model_id=MODEL_ID, free_boundary_final_candidate=str(output),
                 current_run_id=selected, current_run_scope="PW2_FREE_BOUNDARY_FIXED_CONDITIONS",
                 target_run_executed=bool(audits[selected].get("actual_5s_trajectory")),
                 verification=result,
                 next_action="Mission complete for declared free-boundary model" if result["passed"] else "Resolve failed gates; no final claim")
    state["artifact_paths"]={"output_directory":str(output),"native":str(native),
        "index":str(output/"index.html"),"final_report":str(output/"FINAL_REPORT.md"),
        "final_animation":str(output/"final_animation.mp4") if video_passed else None,
        "wall_detail":str(output/"wall_detail.mp4") if video_passed else None}
    if result["passed"]:
        state["milestones"]={"M0_contract_method":"PASS","M1_actual_preview":"PASS",
            "M2_nonlinear_pilot":"PASS","M3_target_dt_h":"PASS","M4_final_review_package":"PASS"}
    state["_job_exit_code"] = 0 if result["passed"] else 2
    return result
