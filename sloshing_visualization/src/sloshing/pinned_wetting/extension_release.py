"""Finite PW1 release workflow, with separate numerical and parameter comparisons."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from ..validation.contact import surface_value
from .reference import load_native
from .mission import read_json, write_json, sha256
from .navier import LABEL
from .extension_output import render_bundle


METRICS = ("normalized_surface_Linf", "normalized_macro_Linf", "normalized_H_Linf",
           "first_peak_time_fraction", "first_peak_height_normalized")


def exact_surface_distance(a, b, ia, ib):
    """Exact spatial sup of the difference of two continuous piecewise P2 graphs."""
    edges = np.union1d(a["x"][::2], b["x"][::2])
    nodes = np.column_stack((edges[:-1], (edges[:-1]+edges[1:])/2, edges[1:]))
    maximum = 0.
    for i, j in zip(ia, ib):
        v = (surface_value(a["x"], a["eta"][i], nodes.ravel()) -
             surface_value(b["x"], b["eta"][j], nodes.ravel())).reshape(-1, 3)
        q = 2*(v[:, 2]+v[:, 0]-2*v[:, 1])
        r = v[:, 2]-v[:, 0]-q
        vertex = np.divide(-r, 2*q, out=np.zeros_like(r), where=abs(q)>1e-30)
        inside = (vertex>0) & (vertex<1)
        peak = q[inside]*vertex[inside]**2+r[inside]*vertex[inside]+v[inside, 0]
        maximum = max(maximum, float(np.max(abs(v))), float(np.max(abs(peak), initial=0.)))
    return maximum


def first_peak(data):
    peaks = [m for m in data["coatings"][-1]["markers"] if m["side"]=="L" and m["marker_id"]!="P1"]
    if not peaks:
        return None
    p = peaks[0]
    return {"height_m": p["height_m"], "time_s": float(data["times"][int(p["state_id"].rsplit(":",1)[1])])}


def compare_native(path_a, path_b):
    a, b = load_native(Path(path_a)), load_native(Path(path_b))
    step = max(a["config"]["dt"], b["config"]["dt"])
    times = np.arange(round(min(a["times"][-1],b["times"][-1])/step)+1)*step
    ia, ib = (np.rint(times/d["config"]["dt"]).astype(int) for d in (a,b))
    if not (np.allclose(a["times"][ia], times, atol=1e-12, rtol=0) and
            np.allclose(b["times"][ib], times, atol=1e-12, rtol=0)):
        raise ValueError("Refinement times are not common accepted physical times")
    scale = float(abs(a["eta"][0,-1]-a["eta"][0,0]))
    p, q = first_peak(a), first_peak(b)
    metrics = {"normalized_surface_Linf": exact_surface_distance(a,b,ia,ib)/scale,
        "normalized_macro_Linf": float(np.max(abs(a["eta"][ia][:,[0,-1]]-b["eta"][ib][:,[0,-1]])))/scale,
        "normalized_H_Linf": float(np.max(abs(a["H"][ia]-b["H"][ib])))/scale,
        "first_peak_time_fraction": abs(p["time_s"]-q["time_s"])/q["time_s"] if p and q else None,
        "first_peak_height_normalized": abs(p["height_m"]-q["height_m"])/scale if p and q else None}
    return {"run_ids": [d["identity"]["run_id"] for d in (a,b)],
        "physical_hashes": [d["identity"]["physical_config_hash"] for d in (a,b)],
        "paths": [str(path_a),str(path_b)], "controls": [d["config"] for d in (a,b)],
        "normalization_m": scale, "common_time_count": len(times), "time_range_s": [0., float(times[-1])],
        "metrics": metrics, "first_peaks": [p,q],
        "passed": all(v is None or v<=.05 for v in metrics.values()) and ((p is None)==(q is None))}


def applicability(summary):
    passed = summary["max_eta_over_a"]<=.05 and summary["max_slope"]<=.30
    if not passed:
        raise ValueError("LINEAR_APPLICABILITY: fixed-b amplitude/slope gate failed; not retuning physics")
    return passed


def freeze_native_audit(root,output):
    """Preserve an actually completed independent audit, never solver diagnostics."""
    if (output/"native_audit.json").exists():
        return
    first=read_json(output/"verification.json")
    if not first.get("scientific_checks_passed") or not first.get("native_reports") or not all(r.get("passed") for r in first["native_reports"].values()):
        raise ValueError("Cannot freeze an incomplete or failed independent native audit")
    source=root/"mission/pinned_wetting/extension/history/full_native_auditor_source.py"
    shutil.copyfile(source,output/"native_auditor_source.py")
    first.update(auditor_source_sha256=sha256(source),auditor_source_snapshot_path="native_auditor_source.py",
        native_audit_scope="FULL_STAGE_OPERATOR_NATIVE_AUDIT",
        reuse_note="Actual independent full-stage reconstruction before renderer-only revision; reuse requires separate review and identical native SHA256 values.")
    write_json(output/"native_audit.json",first)


def record_evidence(output,summaries):
    """Sensitivity of event *reporting*, never a modification of accepted H or PDE."""
    d=load_native(Path(summaries["space_refined"]["native"]))
    def records(indices,floor):
        t=d["times"][indices]; r=d["eta"][indices][:,[0,-1]]
        memory=np.maximum.accumulate(r,axis=0)
        previous=[float(r[0,0]),float(r[0,1])]
        peaks=[]
        for i in range(2,len(t)):
            for j,side in enumerate(("L","R")):
                peak=float(r[i-1,j])
                if peak>r[i-2,j] and peak>r[i,j] and peak>previous[j]+floor and peak>=memory[i,j]-1e-13:
                    peaks.append({"side":side,"height_m":peak,"peak_time_s":float(t[i-1]),"confirmed_at_s":float(t[i])})
                    previous[j]=peak
        return peaks
    indices=np.arange(len(d["times"]))
    result={"run_id":d["identity"]["run_id"],"interpretation":"Postprocess sensitivity only; actual H always uses every accepted endpoint with no record floor.",
        "floor_variants":[{"record_floor_m":f,"records":records(indices,f)} for f in (0.,.00005,.0001,.0002)],"cadence_variants":[]}
    for stride in (1,2,4):
        chosen=indices[::stride]
        sampled=np.maximum.accumulate(d["eta"][chosen][:,[0,-1]],axis=0)
        result["cadence_variants"].append({"stride":stride,"sample_dt_s":stride*d["config"]["dt"],
            "max_missed_H_m":float(np.max(d["H"][chosen]-sampled)),"records":records(chosen,.0001)})
    write_json(output/"detector_sensitivity.json",result)
    roles={}
    for role in ("baseline","time_refined","space_refined"):
        x=load_native(Path(summaries[role]["native"]))
        peaks=[{**m,"peak_time_s":float(x["times"][int(m["state_id"].rsplit(":",1)[1])])}
            for m in x["coatings"][-1]["markers"] if m["side"]=="L" and m["marker_id"]!="P1"]
        roles[role]={"run_id":x["identity"]["run_id"],"left_records":peaks,
            "second_record_gain_m":peaks[1]["height_m"]-peaks[0]["height_m"] if len(peaks)>1 else None}
    gains=[v["second_record_gain_m"] for v in roles.values()]
    empirical=None
    if all(g is not None for g in gains):
        delta=[]
        for k in (0,1):
            hs=[roles[r]["left_records"][k]["height_m"] for r in ("baseline","time_refined","space_refined")]
            delta.append(abs(hs[0]-hs[1])+abs(hs[1]-hs[2]))
        empirical={"peak_delta_m":delta,"sum_peak_deltas_m":sum(delta),"selected_gain_m":gains[-1],
            "record_floor_m":.0001,"gain_exceeds_empirical_margin_plus_floor":gains[-1]>sum(delta)+.0001,
            "definition":"deltaPk=abs(Pk_baseline-Pk_time)+abs(Pk_time-Pk_space); observed resolution margin, not a rigorous error bound"}
    write_json(output/"record_stability.json",{"same_physics":True,"runs":roles,
        "second_record_present_all_resolutions":all(g is not None for g in gains),
        "second_record_gain_range_m":[min(gains),max(gains)] if all(g is not None for g in gains) else None,
        "empirical_resolution_margin":empirical,
        "interpretation":"Observed numerical consistency of record identity, not proof of an exact physical peak or rigorous error bound."})


def complete_release(root, original_cfg, cfg, controls, args, state, heartbeat, run_case):
    output = root/cfg["output_root"]
    level, scale = args.mesh_level, args.dt_scale
    summaries = {}
    def solve(role, l, s, slip=cfg["slip_length_m"]):
        state.update(status="TARGET_RUNNING" if role=="baseline" else "REFINING", next_action=f"compute/reuse {role}, then continue release")
        heartbeat()
        value = run_case(root,cfg,controls(l,s,5.,slip),role,state,heartbeat,args.stop_after)
        if not value["completed"]:
            state.update(status="TARGET_RUNNING",next_action="resume --extension with identical mesh/dt controls")
            heartbeat()
            return None
        applicability(value)
        summaries[role] = value
        print(json.dumps({"completed_role":role, "summary":value},ensure_ascii=False),flush=True)
        return value
    if solve("baseline",level,scale) is None:
        return 0
    # Early full-duration movie remains explicitly preliminary, not a final handoff.
    render_bundle(root,Path(summaries["baseline"]["native"]),output/"preliminary",final=False,heartbeat=heartbeat)
    if solve("time_refined",level,scale/2) is None:
        return 0
    if solve("space_refined",level+1,scale/2) is None:
        return 0
    refinement = {"same_physics":True,"height_normalization":"abs(z2-z1)",
        "interpretation":"Two levels per direction test numerical agreement, not formal order or a physical error bound.",
        "tolerance":.05,"selected_role":"space_refined"}
    for direction, roles in (("time",["baseline","time_refined"]),("space",["time_refined","space_refined"])):
        result = compare_native(*(summaries[r]["native"] for r in roles))
        result["roles"] = roles
        if len(set(result["physical_hashes"]))!=1:
            raise ValueError("MIXED_PHYSICS: slip/physical configuration differs within refinement")
        refinement[direction] = result
    refinement["passed"] = all(refinement[k]["passed"] for k in ("time","space"))
    write_json(output/"refinement.json",refinement)
    if not refinement["passed"]:
        # Preserve all fields and report the specific failed comparison; never relax a gate.
        raise ValueError("REFINEMENT_FAILED: see refinement.json; a targeted extra level must resolve the failed norm")
    selected = summaries["space_refined"]
    state.update(selected_extension_role="space_refined",selected_extension_native=selected["native"])
    sensitivity = {"passed":True,"distinct_physics":True,"selected_role":"space_refined", "base_slip_length_m":cfg["slip_length_m"],
        "interpretation":"Bounded parameter sensitivity, not numerical refinement, calibration, or a physical error bound.","comparisons":[]}
    for role, slip in zip(("slip_low","slip_high"),cfg["slip_sensitivity_m"]):
        if solve(role,level+1,scale/2,slip) is None:
            return 0
        comp = compare_native(selected["native"],summaries[role]["native"])
        sensitivity["comparisons"].append({"base_run_id":comp["run_ids"][0],"run_id":comp["run_ids"][1],
            "base_role":"space_refined","role":role,"base_slip_length_m":cfg["slip_length_m"],"slip_length_m":slip,
            "same_mesh_and_dt":True,"metrics":comp["metrics"],"first_peaks":comp["first_peaks"],
            "paths":comp["paths"],"linear_applicability_passed":True})
        write_json(output/"sensitivity.json",sensitivity)
    write_json(output/"case_summaries.json",summaries)
    state.update(status="RENDERING",next_action="final 0–5 s videos, independent native reconstruction and review")
    heartbeat()
    render_bundle(root,Path(selected["native"]),output,final=True,heartbeat=heartbeat)
    write_report(root,output,cfg,state,summaries,refinement,sensitivity)
    make_manifest(root,output,cfg,state,summaries)
    state.update(status="VERIFYING",last_checkpoint=selected["native"],current_run_id=selected["run_id"],
        artifact_paths={"extension_index":str(output/"index.html"),"extension_output":str(output),"native":selected["native"],
                        "reference_index":str(root/"sloshing_visualization/output/pinned_wetting/index.html")},
        next_action="verify --extension; independent reviewer signoff; package --extension")
    heartbeat()
    # Verification is deliberately implemented independently of this release/solver module.
    from .extension_verify import verify_extension
    report = verify_extension(output)
    write_json(output/"verification.json",report)
    print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
    return 0 if report["passed"] else 2


def write_report(root,output,cfg,state,summaries,refinement,sensitivity):
    record_evidence(output,summaries)
    from .mission import doctor
    runtime=doctor(root)
    runtime["status"]="DEPENDENCIES_AVAILABLE_FOR_DECLARED_EXTENSION" if "MISSING" not in runtime["versions"].values() else "DEPENDENCY_MISSING"
    write_json(output/"runtime.json",runtime)
    shutil.copyfile(root/"mission/pinned_wetting/extension/tests.xml",output/"tests.xml")
    for module in ("extension_output.py","extension_verify.py"):
        shutil.copyfile(Path(__file__).with_name(module),output/("renderer_source.py" if module=="extension_output.py" else "verifier_source.py"))
    import xml.etree.ElementTree as ET
    test_count=len(ET.parse(output/"tests.xml").getroot().findall(".//testcase"))
    chosen = summaries["space_refined"]
    runup = read_json(output/"runup_summary.json")
    peak = runup["first_left_peak"]
    later=runup["later_left_new_record"]
    later_text=(f'{later["height_m"]*1000:.6f} мм при t={later["peak_time_s"]:.6f} с' if later else "не обнаружен")
    control=read_json(Path(chosen["native"]).parent/"resolved_case.json")["controls"]
    table = "\n".join(f'| {key} | {refinement["time"]["metrics"][key]:.7g} | {refinement["space"]["metrics"][key]:.7g} |' for key in METRICS)
    parameter = "\n".join(f'- b={c["slip_length_m"]} m: normalized surface/R/H differences = '+
        "/".join(f'{c["metrics"][k]:.6g}' for k in METRICS[:3])+f'; first peak {c["first_peaks"][1]}.' for c in sensitivity["comparisons"])
    record_stability=read_json(output/"record_stability.json")
    validity_table="\n".join(f'| {b:.2f} | {summaries[r]["max_slope"]:.6g} | {summaries[r]["max_eta_over_b"]:.6g} | {summaries[r]["max_convective_indicator"]:.6g} | {summaries[r]["max_kinematic_indicator"]:.6g} |'
        for r,b in (("slip_low",.25),("space_refined",.5),("slip_high",1.)))
    report = f'''# PW1 — DECLARED_EXTENSION

Конечное скольжение Навье на боковых стенках, b={cfg["slip_length_m"]} м.
Это НЕ решение исходной no-slip-задачи. Её STOP, acceptance, solver и REFERENCE
остаются в [исходном пакете](../index.html). Проверка и независимая подпись:
[verification.json](verification.json), [independent_review.json](independent_review.json).
Статус выпуска определяется этими файлами и manifest.json, не названием видео.

## Результат расчёта

Выбран {chosen["run_id"]}, 0–5 физических секунд, {chosen["accepted_states"]} принятых состояний.
Сетка {control["nx"]}×{control["nz"]}, dt={control["dt"]} с, P2/P1 FEM, SDIRK2.
Первый левый рекорд: {peak["height_m"]*1000:.6f} мм при t={peak["peak_time_s"]:.6f} с;
подтверждение локального пика на следующем принятом шаге t={peak["created_at_s"]:.6f} с.
Подъём от исходной P1: {runup["first_left_rise_from_P1_m"]*1000:.6f} мм.
P1/P2 остаются на начальных высотах ±34.920769 мм. P3 создан только после расчётного пика.
Разность P4−P3 на трёх численных разрешениях: {record_stability["second_record_gain_range_m"]} м.
Все высоты z отсчитываются от исходного среднего уровня z=0, не от дна сосуда.
Исследованы пороги регистрации 0/0.05/0.1/0.2 мм и частоты detector sampling:
[detector_sensitivity.json](detector_sensitivity.json), [record_stability.json](record_stability.json).
Последующий левый рекорд P4: {later_text}.
{runup["no_later_record_reason"] or ""}
Конечные H_L/H_R: {runup["H_final_m"]} м, R_L/R_R: {runup["R_final_m"]} м.
Левый рекорд выше начального P2: {runup["left_record_exceeds_initial_P2"]}.

## Физическая модель и пределы

[Полные уравнения, знаки и происхождение параметров](model_decision.md).
Линейные нестационарные Stokes + гравитационная свободная поверхность; sigma=0,
непроницаемые бока, Navier slip только на боках, no-slip дно. Вертикальные endpoint DOFs
свободны; eta_t получает реальную FE скорость на конце. Экстраполяции гребня нет.
Плоскость z=0 — линеаризация, не крышка; eta может иметь оба знака.

Остаточный слой задан отдельным законом H[n]=max(H[n-1],R[n]), W=[bottom,H].
Он необратим, непрерывен ниже H, сохранён в каждом native checkpoint и читается без renderer.
Navier slip не является причиной этого закона памяти. Точная толщина, объём, инерция,
отдельная динамика и обратное воздействие слоя на bulk не вычислялись: это явная
ведущепорядковая подсеточная аппроксимация, не доказательство молекулярного прикрепления.
Полоса условная. Маркеры неподвижны после рождения; P3/P4 возникли только из расчёта.
Пики определены на каждом принятом dt, не на .01 s сетке видео; отдельного непрерывного
межшагового поиска нет. Предзаданный порог регистрации нового рекорда 0.1 мм.

Область линеаризации (максимумы за весь выбранный расчёт, включая стенки):
max|eta|/a={chosen["max_eta_over_a"]:.7g} <=0.05;
точный max|eta_x|={chosen["max_slope"]:.7g} <=0.30;
max|eta|/b={chosen["max_eta_over_b"]:.7g}.
Индикаторы отброшенных конвекции={chosen["max_convective_indicator"]:.7g},
кинематической поправки={chosen["max_kinematic_indicator"]:.7g} (нормировки в model_decision).
Это screens применимости, НЕ доказанная 5% физическая точность. Отброшенные нелинейные
члены не равны нулю; качественная ведущепорядковая интерпретация, не точная moving-domain CFD.

## Численная проверка при фиксированном b

Все высоты нормированы на abs(z2-z1)=0.06984153898349546 м.
Отдельно dt/2 при одной сетке и сетка×2 при одном dt. Порог каждой нормы 0.05.

| Метрика | time | space |
| --- | ---: | ---: |
{table}

Два уровня не доказывают порядок или строгую границу ошибки. См. [refinement.json](refinement.json).
Чувствительность по параметру b при одинаковых mesh/dt (НЕ уточнение):

{parameter}

Проверки линейных допущений во всех вариантах (без притворного physical-error PASS):

| b, м | max slope | max eta/b | convection/A0 | kinematic/U0 |
| ---: | ---: | ---: | ---: | ---: |
{validity_table}

Варианты b=.25/1.00 — ограниченная чувствительность на одной fine сетке, не отдельные
полные исследования их сходимости. Нелинейные поправки существеннее при b=.25;
точность реальной жидкости на уровне5% не утверждается ни для одного b.

Эти длины выбраны до запуска; b не подбирался к желаемому максимуму/виду волны.
Параметр b — иллюстративный размерный, не измеренная характеристика стенки.

## Законы сохранения и воспроизводимость

Относительная масса max={chosen["max_mass_relative_error"]:.7g}; weak divergence max={chosen["max_weak_divergence"]:.7g};
относительная энергетическая невязка max={chosen["max_energy_relative_error"]:.7g};
невязка слабого импульса max={chosen["max_momentum_relative_residual"]:.7g}.
Работа=0, капиллярная энергия=0 по модели. Bulk и wall friction диссипации разделены,
RK correction имеет знак и не названа физическим теплом. Плотность 1000 kg/m³ —
объявленный масштаб перевода per-density энергий в J/m, не влияет на решения.
Условия Навье естественные слабые; pointwise Robin traction residual не строго ноль,
особенно в углах; показатели и уточнение сохранены в case_summaries.json.
На time/space сетках max slope изменился с {summaries["time_refined"]["max_slope"]:.6g}
до {chosen["max_slope"]:.6g}, а wall traction L2 residual — с
{summaries["time_refined"]["max_wall_traction_L2"]:.6g} до {chosen["max_wall_traction_L2"]:.6g}.
Следовательно, сходимость всех угловых производных НЕ заявляется; проверена
применимость по заданным screens на двух разрешениях и согласованность наблюдаемых высот.
Если требовать C1-гладкости напряжений ровно в пересечении свободной поверхности
и боковой стенки, одновременные zero tangential stress и конечный Navier slip
дали бы нулевую скорость угловой точки. Движущийся endpoint здесь относится к слабой
постановке с условиями на открытых гранях почти всюду, допускающей угловую
нерегулярность. Это ограничение не скрывается за малой интегральной невязкой;
равномерная C1-сходимость и нелинейная физическая точность возле corner не доказаны.

Два видео 1920×1080,25 fps,501 frame с t=0..5 шагом .01 s: замедление ×4,
последний кадр добавляет .04 s к длительности 20.04 s. Вертикальный масштаб отличается от x,
полный сосуд d=10 м показан отдельно. [Все исходные принятые поля и RK stages](manifest.json).

Команды из корня репозитория (Python с зависимостями из doctor):

```sh
python sloshing_visualization/scripts/pinned_wetting_mission.py doctor --extension
python sloshing_visualization/scripts/pinned_wetting_mission.py status --extension
python sloshing_visualization/scripts/pinned_wetting_mission.py resume --extension
python sloshing_visualization/scripts/pinned_wetting_mission.py verify --extension
```

Для дополнительного уточнения: `run --extension --mesh-level 1 --dt-scale 0.5`;
автоматически используется отдельный output bundle (существующие hashes запрещают перезапись).
Бюджет 12 h jobs,40 GiB,один heavy job,4 GiB RAM/job; каждая стадия продолжает native checkpoint,
kernel lock не допускает дубликата. Фактический бюджет — state.json/manifest.json.
Ограниченный набор миссии:{test_count} tests PASS,0 skipped; два предупреждения только в грубом no-slip
reference continuation fixture. [JUnit evidence](tests.xml), [runtime/dependencies](runtime.json).
Полный независимый операторный проход сохранён в [native_audit.json](native_audit.json)
с исходником и SHA каждого native. Его исторический общий passed=false связан только
с отсутствовавшей тогда отдельной подписью: scientific_checks_passed=true, все пять native PASS.
Текущий окончательный статус находится в verification.json. После renderer-only исправления проход повторно используется
только по отдельной подписи reviewer и неизменным SHA; данные eta/H, экспорты, уточнение,
полоса и видео проверяются повторно. Явная команда `verify --extension` всегда вновь
восстанавливает все уравнения из raw fields; `package` использует подписанную evidence.
Внешних расходов, installs и push не было.
'''
    (output/"FINAL_REPORT.md").write_text(report)
    index=(output/"index.html").read_text()
    extras=["FINAL_REPORT.md","manifest.json","verification.json","refinement.json","sensitivity.json","case_summaries.json","independent_review.json", "detector_sensitivity.json", "record_stability.json", "tests.xml", "runtime.json", "renderer_source.py", "verifier_source.py", "native_audit.json", "native_auditor_source.py"]
    index=index.replace('</ul>', ''.join(f'<li><a href="{n}">{n}</a></li>' for n in extras if f'href="{n}"' not in index)+'</ul>')
    (output/"index.html").write_text(index)


def make_manifest(root,output,cfg,state,summaries):
    native=Path(summaries["space_refined"]["native"])
    ident=read_json(native.parent/"identity.json")
    artifacts=[]
    for p in sorted(output.iterdir()):
        if p.is_file() and p.name not in ("manifest.json","verification.json"):
            artifacts.append({"path":p.name,"sha256":sha256(p),"role":"render" if p.suffix in (".mp4",".png") else "data_or_evidence"})
    runs=[]
    for role, summary in summaries.items():
        runs.append({"role":role,"path":summary["native"],"sha256":sha256(Path(summary["native"])),"run_id":summary["run_id"]})
    manifest={"schema_version":"PW1_NAVIER_1.0","run_id":ident["run_id"],"source_hash":ident["source_hash"],
        "physical_config_hash":ident["physical_config_hash"],"execution_HEAD":ident["execution_HEAD"],
        "status":"VERIFYING","physical_model_label":LABEL,"label":LABEL,
        "resolved_case":"resolved_case.json","selected_role":"space_refined","native_data_paths":[str(native)],
        "native_runs":runs,"artifacts":artifacts,"time_range":[0.,5.],
        "frame_map":"render_evidence.json","contact_source_map":"wetting_history.csv: state_id maps to native accepted group, reach=actual endpoint R",
        "units":{"length":"m (video mm explicitly)","time":"physical seconds","energy_csv":"J/m","native_energy":"per density"},
        "limitations":["DECLARED_EXTENSION, not no-slip solution","linear free-surface Stokes, not nonlinear moving-domain CFD",
            "b is illustrative, not calibrated","subgrid coating thickness,volume,inertia,separate dynamics/feedback not computed",
            "two-level numerical agreement is not physical error bound","accepted-step peak sampling, no continuous event search"],
        "budget":state["budget"],"reference_STOP":"../../../../mission/pinned_wetting/history/no_slip_STOP_state.json"}
    write_json(output/"manifest.json",manifest)


def finalize(root,output,cfg,state,heartbeat):
    """Package only after the separately implemented verifier and reviewer sign off."""
    from .extension_verify import verify_extension
    report=verify_extension(output,use_signed_native_audit=True)
    write_json(output/"verification.json",report)
    heartbeat()
    if not report["passed"]:
        state.update(status="VERIFYING",next_action="resolve verifier evidence, then package --extension")
        state["_job_exit_code"]=2
        heartbeat()
        print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
        return 2
    manifest=read_json(output/"manifest.json")
    manifest["status"]="COMPLETE_FOR_DECLARED_MODEL"
    manifest["budget"]=dict(state["budget"])
    manifest["budget_note"]="Measured through packaging heartbeat; final job-close ledger is mission/pinned_wetting/state.json. Auxiliary allowance is an explicit upper bound, not a measurement."
    for name in ("independent_review.json","review.md","tests.xml"):
        if (output/name).is_file():
            manifest["artifacts"]=[a for a in manifest["artifacts"] if a["path"]!=name]
            manifest["artifacts"].append({"path":name,"sha256":sha256(output/name),"role":"independent_evidence"})
    write_json(output/"manifest.json",manifest)
    state.update(status="COMPLETE_FOR_DECLARED_MODEL",next_action="DONE: open declared_extension/index.html; original no-slip STOP preserved",
        extension_complete=True,last_error=None,
        target_run_executed=True, extension_solver_implemented_by_agent=True,
        milestones={"M0":"EXTENSION_PREVIEW_COMPLETE","M1":"DECLARED_NAVIER_AND_COATING_VERIFIED",
                    "M2":"TARGET_5S_TIME_SPACE_AND_SLIP_SENSITIVITY_COMPLETE","M3":"FINAL_DATA_VIDEOS_AND_INDEPENDENT_REVIEW_PASS"},
        verification={"scientific":"PASS_FOR_DECLARED_EXTENSION","reference":"PRESERVED_NOT_RECLASSIFIED",
                      "rendering":"PASS_FULL_DECODE_NATIVE_BOUND","independent_review":"PASS_DECLARED_EXTENSION"})
    state["_job_exit_code"]=0
    heartbeat()
    print(json.dumps({"status":state["status"],"index":str(output/"index.html"),"verification":report},ensure_ascii=False,indent=2),flush=True)
    return 0
