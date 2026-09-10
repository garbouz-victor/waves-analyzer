"""Small, honest blocker handoff for the independently rejected corrected case."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

import h5py
import numpy as np

from .asymmetric import MODEL_ID
from .asymmetric_verify import verify_native
from .mission import doctor, read_json, write_json, sha256


def test_counts(path):
    suites = ET.parse(path).getroot()
    if suites.tag == "testsuite":
        suites = [suites]
    return {key:sum(int(s.get(key, 0)) for s in suites) for key in ("tests", "failures", "errors", "skipped")}


def energy_comparisons(audit):
    """Common-time series already checked against native fields by the auditor."""
    from .reference import load_native
    result = {}
    fields = ("kinetic_energy", "gravitational_potential_energy", "cumulative_bulk_dissipation",
              "cumulative_wall_dissipation", "rk_energy_correction")
    for name, first, second in (("temporal", "pilot", "pilot_time"), ("spatial", "pilot_time", "pilot_space")):
        ra, rb = (audit["native_reports"][role] for role in (first, second))
        if ra["physical"] != rb["physical"] or not all(r["numerical_checks_passed"] for r in (ra, rb)):
            raise ValueError("Energy comparison requires the same independently audited corrected physics")
        a, b = (load_native(r["native"]) for r in (ra, rb))
        common, ia, ib = np.intersect1d(np.round(a["times"], 12), np.round(b["times"], 12), return_indices=True)
        scale = a["diagnostics"][0]["gravitational_potential_energy"]
        norms = {key:float(np.max(abs(np.array([r[key] for r in a["diagnostics"]])[ia] -
                                      np.array([r[key] for r in b["diagnostics"]])[ib])))/scale for key in fields}
        result[name] = {"run_ids":[ra["run_id"], rb["run_id"]], "time_range_s":[float(common[0]), float(common[-1])],
            "common_times":len(common), "max_absolute_differences_normalized_by_initial_potential":norms}
    return {"scope":"UNQUALIFIED_PILOT_ONLY", "normalization":"initial potential energy per density",
        "wall_term":"LEFT_Navier_only", "signed_RK_is_not_heat":True,
        "source":"diagnostic series independently reconstructed from accepted native fields/stages", "comparisons":result}


def verify_package(root, output, heartbeat=lambda:None, *, check_manifest=True):
    """Always reconstruct native asymmetric operators; never trust passed=true."""
    output = Path(output)
    review = read_json(output / "independent_review.json")
    audit = read_json(output / "native_audit.json")
    integrity = []
    for name, expected in review["artifact_hashes"].items():
        path = output / name
        if not path.is_file() or sha256(path) != expected:
            integrity.append("review_binding:" + name)
    contract = root / "mission/pinned_wetting/left_navier_right_noslip"
    for name, expected in review["contract_hashes"].items():
        if sha256(contract / name) != expected:
            integrity.append("contract_binding:" + name)
    if review["model_id"] != MODEL_ID or audit["model_id"] != MODEL_ID:
        integrity.append("wrong_model_id")
    reports = {}
    right_error = right_H_error = 0.
    for role, original in audit["native_reports"].items():
        native = Path(original["native"])
        if sha256(native) != original["native_sha256"]:
            integrity.append("native_binding:" + role)
        reports[role] = verify_native(native, heartbeat=heartbeat)
        with h5py.File(native) as h:
            z2 = float(h["accepted/0/eta"][-1])
            for group in h["accepted"].values():
                right_error = max(right_error, abs(float(group["eta"][-1])-z2))
                right_H_error = max(right_H_error, abs(float(group["H"][1])-z2))
        heartbeat()
    videos = []
    evidence = read_json(output / "preview/render_evidence.json")
    from ..visualization.render import probe_video
    for video in evidence["videos"]:
        path = output / "preview" / video["video"]
        subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-threads", "1", "-i", str(path),
                        "-f", "null", "-"], check=True, capture_output=True, timeout=60)
        probe = probe_video(path)
        frames = video["frame_map"]
        if (int(probe["nb_read_frames"]) != len(frames) or probe["avg_frame_rate"] != "25/1"
                or (probe["width"], probe["height"]) != (1920,1080)):
            integrity.append("video_frames:" + video["video"])
        expected_times = np.arange(101)*.01
        if len(frames) != 101 or not np.allclose([f["time_s"] for f in frames], expected_times, rtol=0, atol=1e-12):
            integrity.append("video_time_map:" + video["video"])
        with h5py.File(audit["native_reports"]["pilot"]["native"]) as h:
            for frame in frames:
                step = frame["state_id"].rsplit(":", 1)[1]
                group = h["accepted"][step]
                if (group.attrs["state_id"] != frame["state_id"] or
                        abs(group.attrs["time_s"]-frame["time_s"]) > 1e-12 or
                        not np.array_equal(group["H"][:], frame["H_m"])):
                    integrity.append("native_frame_map:" + video["video"])
                    break
        videos.append({"file":"preview/"+video["video"], "full_decode_passed":True, "probe":probe})
        heartbeat()
    if check_manifest:
        manifest = read_json(output / "manifest.json")
        if manifest["model_id"] != MODEL_ID or manifest["status"] != "NEEDS_MODEL_DECISION":
            integrity.append("manifest_identity_or_status")
        for artifact in manifest["artifacts"]:
            path = output / artifact["path"]
            if not path.is_file() or sha256(path) != artifact["sha256"]:
                integrity.append("manifest:" + artifact["path"])
    numerical = all(r["numerical_checks_passed"] for r in reports.values())
    confirmed = numerical and all(r["failures"] == ["limit:max_slope"] for r in reports.values())
    tests = test_counts(output / "tests.xml")
    return {"status":"NEEDS_MODEL_DECISION", "passed":False, "model_id":MODEL_ID,
        "target_5s":"NOT_RUN", "qualified_final_animation":"NOT_CREATED",
        "blocker":review["blocker"], "scientific_blocker_confirmed":confirmed,
        "discrete_equations_and_BCs_passed":numerical,
        "linear_applicability_passed":False, "native_reports":reports,
        "max_t_abs_eta_right_minus_z2_m":right_error, "max_t_abs_H_right_minus_z2_m":right_H_error,
        "pilot_refinement":audit["refinement"], "videos":videos,
        "artifact_integrity_passed":not integrity, "integrity_failures":integrity,
        "local_tests":tests, "local_tests_passed":tests["failures"]==tests["errors"]==0,
        "verifier_source_sha256":sha256(Path(__file__).with_name("asymmetric_verify.py")),
        "verifier_source_snapshot":"asymmetric_verifier_source.py",
        "independent_review_status":review["status"]}


def package_blocker(root, output, state, heartbeat):
    """Keep the signed pilot evidence intact and clearly reject final acceptance."""
    from .asymmetric_mission import mark_blocked
    output = Path(output)
    audit = read_json(output / "native_audit.json")
    verification = verify_package(root, output, heartbeat, check_manifest=False)
    if not (verification["scientific_blocker_confirmed"] and verification["artifact_integrity_passed"]
            and verification["local_tests_passed"]):
        raise ValueError("Blocker handoff itself failed: " + repr(verification))
    preview = output / "preview"
    shutil.copyfile(Path(__file__).with_name("asymmetric_verify.py"), output / "asymmetric_verifier_source.py")
    for name in ("resolved_case.json", "surface_history.csv", "wetting_history.csv", "marker_catalog.csv",
                 "marker_tracks.csv", "energy_history.csv", "key_frames.png", "model_decision.md"):
        shutil.copyfile(preview / name, output / name)
    runup = read_json(preview / "runup_summary.json")
    runup.update(qualification="UNQUALIFIED_LINEAR_PILOT_NOT_PHYSICAL_PREDICTION", time_range_s=[0., 1.],
        H_at_5s_m=None, P4_at_5s="UNKNOWN_NOT_RUN",
        uncertainty="See refinement.json. Spatial whole-surface engineering gate FAIL; no converged physical record claim.")
    write_json(output / "runup_summary.json", runup)
    write_json(output / "refinement.json", {"scope":"PILOT_0_TO_1S_ONLY", "target_5s":"NOT_RUN",
        "same_left_b_m":.5, "engineering_limit":.05, "comparisons":audit["refinement"],
        "energy_observables":energy_comparisons(audit), "passed":False, "formal_order_claimed":False})
    write_json(output / "sensitivity.json", {"status":"NOT_RUN_DUE_TO_LINEAR_APPLICABILITY_BLOCKER",
        "requested_left_b_m":[.25, .5, 1.], "right_wall":"no_slip", "is_mesh_convergence":False,
        "b_was_not_retuned":True})
    write_json(output / "verification.json", verification)
    mark_blocked(state, output, verification["blocker"])
    state.update(model_decision=str(root / "mission/pinned_wetting/left_navier_right_noslip/model_decision.md"),
        corrected_milestones={"implementation":"PASS", "pilot_preview":"SAVED_NOT_FINAL",
            "independent_native_audit":"PASS_DISCRETE_EQUATIONS_REJECT_TARGET",
            "linear_applicability":"FAIL", "pilot_spatial_refinement":"FAIL",
            "target_5s":"NOT_RUN", "final_package":"BLOCKER_EVIDENCE_ONLY"})
    state["milestones"] = dict(state["corrected_milestones"])
    for role,r in audit["native_reports"].items():
        native=Path(r["native"])
        state.setdefault("corrected_runs",{}).setdefault(role,{}).update(
            run_id=r["run_id"], path=str(native.parent), completed=True,
            config=read_json(native.parent/"resolved_case.json")["controls"])
    if (output / "interrupted_cli_attempt.json").exists():
        state["corrected_interrupted_attempts"] = [read_json(output / "interrupted_cli_attempt.json")]
    state["artifact_paths"] = {"output_directory":str(output), "index":str(output / "index.html"),
        "preview":str(preview / "preview_extension.mp4"), "wall_detail_preview":str(preview / "wall_detail_preview.mp4"),
        "final_report":str(output / "FINAL_REPORT.md"), "native":audit["native_reports"]["pilot"]["native"]}
    state["verification"] = {"target":False, "discrete_equations_and_BCs":True,
        "linear_applicability":False, "pilot_spatial_refinement":False, "blocker_evidence":True}
    heartbeat()
    write_json(output / "runtime.json", {"budget_snapshot":state["budget"], "environment":doctor(root),
        "authoritative_job_close_ledger":"mission/pinned_wetting/state.json", "forecast_was_not_a_blocker":True,
        "new_external_expenses":0, "push_performed":False})
    write_report(root, output, verification)
    links = {"Отчёт и единственный blocker":"FINAL_REPORT.md", "Native/operator verify":"verification.json",
        "Независимый review":"review.md", "Исходные данные и кадры":"preview/index.html",
        "Уточнение (только pilot)":"refinement.json", "Тесты чистого окружения":"tests.xml"}
    links.update({"Native " + role:os.path.relpath(r["native"], output) for role,r in audit["native_reports"].items()})
    (output / "index.html").write_text('<!doctype html><html lang="ru"><meta charset="utf-8">'
        '<title>PW1 corrected — NEEDS_MODEL_DECISION</title><style>body{max-width:1100px;margin:40px auto;'
        'font:18px sans-serif}video{width:100%}li{margin:12px}</style><h1>NEEDS_MODEL_DECISION</h1>'
        '<p>PW1_LEFT_NAVIER_RIGHT_NOSLIP_V1. Слева Навье b=0,50 м, справа и на дне no-slip, σ=0.</p>'
        '<p>Только реальный линейный pilot 0–1 с. Наклон 9,445 / 15,701 превышает предел 0,30; '
        'пространственное различие поверхности 18,24% &gt; 5%. Это НЕ конечный результат 0–5 с.</p>'
        '<p>P2 — фактический неподвижный endpoint. Остаточное покрытие задано отдельно и подсеточно; '
        'толщина, объём, инерция, отдельная динамика и обратное влияние не вычислялись.</p>'
        '<video controls preload="metadata" src="preview/preview_extension.mp4"></video>'
        '<video controls preload="metadata" src="preview/wall_detail_preview.mp4"></video><ul>' +
        ''.join(f'<li><a href="{url}">{label}</a></li>' for label,url in links.items()) + '</ul></html>')
    artifacts = [{"path":str(p.relative_to(output)), "sha256":sha256(p)} for p in sorted(output.rglob("*"))
                 if p.is_file() and p.name != "manifest.json" and not any(part.startswith("extra_") for part in p.relative_to(output).parts)]
    write_json(output / "manifest.json", {"model_id":MODEL_ID, "status":"NEEDS_MODEL_DECISION",
        "physical":audit["native_reports"]["pilot"]["physical"], "time_range_s":[0., 1.], "target_time_range_s":[0., 5.],
        "target_status":"NOT_RUN", "artifacts":artifacts,
        "native_data":[{"role":role, "path":r["native"], "sha256":r["native_sha256"]} for role,r in audit["native_reports"].items()],
        "execution_HEAD":doctor(root)["execution_HEAD"], "numerical_identity_excludes_execution_HEAD":True})
    heartbeat()


def write_report(root, output, report):
    audit = read_json(output / "native_audit.json")
    runs = audit["native_reports"]
    rows = []
    for role, r in runs.items():
        c = read_json(Path(r["native"]).parent / "resolved_case.json")["controls"]
        m = r["maxima"]
        rows.append(f'| {role} | {c["nx"]}×{c["nz"]} | {c["dt"]} | {m["max_slope"]:.9f} | '
                    f'{m["convective_indicator"]:.9f} | {m["kinematic_indicator"]:.9f} |')
    git_state = subprocess.check_output(["git", "status", "--short"], cwd=root, text=True)
    counts = report["local_tests"]
    passed_tests = counts["tests"]-counts["skipped"]-counts["failures"]-counts["errors"]
    text = f'''# PW1 corrected — NEEDS_MODEL_DECISION

LEFT WALL: Navier-slip, b_left=0.50 m. RIGHT WALL: full no-slip.
BOTTOM: full no-slip. sigma=0. Model ID: PW1_LEFT_NAVIER_RIGHT_NOSLIP_V1.
P2 remains the actual current right contact point at the same height for the entire
computed trajectory (0–1 s); no qualified 0–5 s trajectory was produced.

## Один воспроизводимый blocker

`PW1_PINNED_RIGHT_LINEAR_APPLICABILITY`: у закреплённого правого контакта точный
наклон P2-поверхности превышает заранее принятый предел 0.30. На основной сетке
первое превышение при t=.090 s (|eta_x|=.307613890456), на мелкой при t=.0775 s
(|eta_x|=.314260270008). Максимумы ниже. Уменьшение dt нарушение не устраняет;
пространственное различие полной поверхности 18.241059% также превышает 5%.
Это не доказательство отсутствия слабого continuum-решения; расчётная линейная
траектория существует, но требуемая физическая квалификация не проходит.
Необходима новая постановка нелинейной геометрии поверхности/пристеночной зоны.
Она не реализована без нового физического решения пользователя. b, alpha, P2 и пороги не изменялись.

## Реально выполнено

Три новых расчёта 0–1 s; уточнение времени и сетки при неизменном b_left=.50 m.
Ранний preview: 24×48, dt=.005 s. Обе MP4 полностью декодированы: 101 кадр,
1920×1080,25fps; отображаемое время0–1 s. Это preview, не final_animation.mp4.

| Расчёт | Сетка | dt, s | max exact slope | convective indicator | kinematic indicator |
|---|---|---|---|---|---|
{chr(10).join(rows)}

max exact |eta|/a={runs['pilot_space']['maxima']['max_eta_over_a']:.12g} (<.05).
Все surface edges, в том числе у стенок, включены. Omitted-term indicators не являются
строгими оценками физической ошибки. Причина остановки — hard slope screen, не просто ненулевые indicators.

`max_t |eta_right(t)-z2| = {report['max_t_abs_eta_right_minus_z2_m']:.1f} m` для всех 1003 принятых состояний трёх runs.
`max_t |H_R-z2| = {report['max_t_abs_H_right_minus_z2_m']:.1f} m`.
Максимумы обеих компонент правой velocity trace равны 0; FEM trace rows=[1,0].
Непроницаемость слева, no-slip дна, momentum, weak divergence, mass и
bulk/LEFT-wall/signed-RK energy audit прошли независимо. RK correction — не физическое тепло.
См. реальные residuals и интегралы в native_audit.json и verification.json.

P1=-34.920769492mm, P2=+34.920769492mm, обе заданы приt=0 и неизменны.
В НЕквалифицированном линейном pilot P3=24.025594873mm приt=.785 s,
подтверждение/создание при.790 s; fine P3=24.030200794mm при.785 s,
создание при.7875 s. Эти пики происходят ПОСЛЕ провала screen и не являются
достоверным физическим прогнозом. P4 в 0–1 s не возник; наличие P4 за 0–5 s неизвестно.
H_L(1)=24.025594873mm (fine24.030200794mm), H_R(1)=34.920769492mm.
H_L(5), H_R(5) не вычислены; условие модели предписывает H_R=z2, но данных за 5 s нет.

Все реальные измеренные различия: refinement.json и native_audit.json/refinement.
Временной тест: одна сетка,dt/2. Пространственный: одинаковыйdt=.0025s.
Energy-series comparisons отдельно сохранены в refinement.json/energy_observables
с нормировкой на начальную потенциальную энергию; signed RK не суммируется с физическим теплом.
Formal convergence order по двум уровням не заявлен. Близость левого пика не
компенсирует провал нормы всей поверхности. Target refinement 0–5 s и sensitivity
b_left=.25/1m не запускались после подтверждённого hard blocker; старые symmetric runs не использованы.

## Покрытие и ограничения

R_L — actual FEM endpoint, не extrapolation. H_L обновлён max после каждого
accepted step; lower-connected coating хранится в native/checkpoints и доступен
без renderer. P3 возникает из вычисленного левого рекорда, правая P2 закреплена
физически. Navier и необратимая память — два разных закона.
DECLARED_EXTENSION: слой задан подсеточно, его толщина, объём, инерция,
отдельная динамика и feedback to bulk не рассчитывались. Молекулярное прикрепление не доказано.

## Engineering / воспроизведение

PyYAML>=6,<7 добавлен в requirements.txt и обычные project dependencies.
Чистый venv установлен из requirements; `pip check` и mission doctor выполнены.
Полный локальный pytest: {passed_tests} PASS, {counts['skipped']} SKIP,
{counts['failures']} failures, {counts['errors']} errors; подробности tests.xml и clean_environment.json.
Remote GitHub CI не запускался, push не выполнялся. Дополнительные обязательные
asymmetric negative tests:16 PASS в независимой reviewer_tests.xml.
Регрессия создаёт два реальных git commit во временном fixture repo: после изменения
только документации CLI resume сохраняет численную identity и записывает новый HEAD в provenance.
Первый smoke дополнительно выявил разницу JSON 1/1.0 в digest. Лишний запуск был остановлен;
частичные данные сохранены и исключены из scientific evidence (interrupted_cli_attempt.json).
Исправлен выбор прежнего exact control payload; добавлен отдельный regression.
Отдельно проверена переносимость сетки NumPy1.21.5→2.0.2: различие координат
не более 2.22e−16 m при точном совпадении topology. Guard допускает только
8 machine eps × масштаб координатной оси; constrained DOFs, P2, physical screens
и residual tolerances не ослаблены. Исходный подписанный аудит сохранён;
текущий verifier повторно восстанавливает операторы в чистом окружении.

Команды из корня (python с установленными requirements):

```sh
python sloshing_visualization/scripts/pinned_wetting_mission.py blocker-check --corrected
python sloshing_visualization/scripts/pinned_wetting_mission.py verify --corrected
python sloshing_visualization/scripts/pinned_wetting_mission.py pilot --corrected --mesh-level 1 --dt-scale 0.5
```

Ожидаемый exit2 означает подтверждённый научный blocker, не crash. Последняя команда
проверена CLI regression и на сохранённом mission run; она переиспользует совместимый
48×96 pilot и не создаёт дубликат PDE. Для ещё одной сетки mesh-level2 потребует
нового реального расчёта под теми же ресурсными caps. Никакая команда не подменяет final за 5 s.

HEAD: {doctor(root)['execution_HEAD']}. Изменения corrective iteration оставлены в worktree,
без нового commit/push. Старый no-slip STOP и symmetric Navier/Navier native bundle
не изменены; symmetric COMPLETE сохранён в history и классифицирован как superseded.

Точные пути данных — manifest.json/native_data; видео — preview/preview_extension.mp4
и preview/wall_detail_preview.mp4; HTML — index.html. Состояние и окончательный
учёт jobs/ресурсов — mission/pinned_wetting/state.json. НЕ COMPLETE.

Worktree при упаковке:

```text
{git_state.rstrip()}
```
'''
    (output / "FINAL_REPORT.md").write_text(text)
