"""PW1 reference adapter: native states -> CSV and visibly labelled MP4.

Reuses the project's exact P2 surface evaluator and STEP2 ffprobe helper.
No solving, wall extrapolation, or contact memory construction happens here.
"""
import csv
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Rectangle
import numpy as np

from ..validation.contact import surface_value
from ..visualization.render import probe_video
from .reference import LABEL, load_native

FILM_LABEL = "Остаточный слой показан условно; толщина не рассчитывалась"


def csv_write(path, columns, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(columns)
        writer.writerows(rows)


def export_reference(output, native, heartbeat=lambda: None):
    from .mission import sha256, write_json
    output, native = Path(output), Path(native)
    output.mkdir(parents=True, exist_ok=True)
    if (output / "manifest.json").exists():
        old = json.loads((output / "manifest.json").read_text())
        if old["native_sha256"] == sha256(native) and (output / "preview_reference.mp4").exists():
            return  # Never overwrite an existing completed bundle.
        raise ValueError("Existing preview bundle differs; select a new run directory")
    data = load_native(native)
    c = data["config"]
    x, times, eta, H = data["x"], data["times"], data["eta"], data["H"]
    run_id = data["identity"]["run_id"]
    rho = json.loads((native.parent / "resolved_case.json").read_text())["requested_case"]["physics"]["density_kg_m3"]
    frame_times = np.arange(round(times[-1] / .01) + 1) * .01
    frame_ids = np.rint(frame_times / c["dt"]).astype(int)
    if np.max(abs(times[frame_ids] - frame_times)) > 1e-10:
        raise ValueError("Reference frames must match accepted states exactly")
    csv_write(output / "surface_history.csv", ["time_s", "point_index", "x_m", "z_m", "component_id", "state_id"],
              ((t, j, xx, zz, "main", data["ids"][i]) for i, t in enumerate(times)
               for j, (xx, zz) in enumerate(zip(x, eta[i]))))
    csv_write(output / "wetting_history.csv",
              ["time_s", "R_L_m", "R_R_m", "reach_L_m", "reach_R_m", "H_L_m", "H_R_m", "state_id"],
              ((t, *eta[i, [0, -1]], *eta[i, [0, -1]], *H[i], data["ids"][i]) for i, t in enumerate(times)))
    markers = data["coatings"][-1]["markers"]
    csv_write(output / "marker_catalog.csv", ["marker_id", "side", "height_m", "created_at_s", "state_id"],
              ([m[k] for k in ("marker_id", "side", "height_m", "created_at_s", "state_id")] for m in markers))
    csv_write(output / "marker_tracks.csv", ["time_s", "marker_id", "side", "height_m"],
              ((float(t), m["marker_id"], m["side"], m["height_m"]) for t in frame_times
               for m in markers if m["created_at_s"] <= t + 1e-12))
    energy_cols = ["time_s", "kinetic_J_per_m", "potential_J_per_m", "capillary_J_per_m", "wall_J_per_m",
                   "physical_dissipation_integral_J_per_m", "external_work_J_per_m",
                   "numerical_energy_remainder_J_per_m", "budget_defect_J_per_m", "state_id"]
    csv_write(output / "energy_history.csv", energy_cols,
              ((r["time"], rho*r["kinetic_energy"], rho*r["gravitational_potential_energy"], 0., 0.,
                rho*r["cumulative_viscous_dissipation"], rho*r["cumulative_work"], rho*r["rk_energy_correction"],
                rho*r["energy_balance_residual"], data["ids"][i]) for i, r in enumerate(data["diagnostics"])))
    write_json(output / "runup_summary.json", {"status": "TARGET_NOT_RUN", "physical_model_label": LABEL,
        "first_target_left_peak": None, "subsequent_target_record": None,
        "reason": "Moving-contact closure not qualified; reference endpoints are fixed by construction.",
        "reference_only": {"H_final_m": H[-1].tolist(), "peak_time_s": None,
                           "marker_ids": [m["marker_id"] for m in markers]}})
    write_json(output / "refinement.json", {"status": "NOT_RUN", "same_physics": None,
        "time_refinement": None, "mesh_refinement": None, "contact_regularizer_sensitivity": None,
        "reason": "No qualified target contact law; reference refinement cannot validate advancing contact.",
        "reference_resolution_controls": "run --mesh-level 1 --dt-scale 0.5 (new reference run only)"})
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor="#f6f8fc")
    ax = fig.add_axes([.09, .32, .68, .50])
    inset = fig.add_axes([.83, .30, .09, .51])
    history = fig.add_axes([.09, .12, .68, .12])
    extent = max(.06, float(np.max(abs(eta))) * 1.5)
    ax.set(xlim=(-c["a"]-0.15, c["a"]+0.15), ylim=(-extent*1000, extent*1000),
           xlabel="x, м", ylabel="z, мм", title="Расчётная поверхность; разные масштабы x и z")
    ax.grid(alpha=.2)
    for side in (-c["a"], c["a"]):
        ax.axvline(side, color="#24354b", lw=3)
    dense_x = np.unique(np.concatenate([np.linspace(l, r, 9) for l, r in zip(x[::2][:-1], x[::2][1:])]))
    line, = ax.plot([], [], color="#1268ac", lw=3, zorder=4)
    ax.plot(x, eta[0]*1000, ls="--", color="#8391a5", lw=1.2, label="Начальная прямая")
    strips = []
    for i, side in enumerate((-c["a"], c["a"])):
        patch = Rectangle((side if i == 0 else side-.018, -extent*1000), .018, 0,
                          color="#15b6c8", alpha=.85, zorder=5)
        ax.add_patch(patch)
        strips.append(patch)
    for m in markers:
        xx = -c["a"] if m["side"] == "L" else c["a"]
        ax.scatter([xx], [m["height_m"]*1000], s=90, c="#7d2b65", zorder=10)
        ax.annotate(f'{m["marker_id"]}: {m["height_m"]*1000:+.3f} мм',
                    (xx, m["height_m"]*1000), xytext=(20 if xx < 0 else -20, 18), textcoords="offset points",
                    ha="left" if xx < 0 else "right", color="#7d2b65", fontsize=13)
        ax.axhline(m["height_m"]*1000, lw=.8, color="#7d2b65", ls=":")
    inset.set(xlim=(-1.2*c["a"], 1.2*c["a"]), ylim=(-c["d"], .5), title="Сосуд, м")
    inset.set_aspect("equal")
    inset.fill_between([-c["a"], c["a"]], -c["d"], 0, color="#bddff3")
    inset.plot([-c["a"], -c["a"], c["a"], c["a"]], [.3, -c["d"], -c["d"], .3], c="#24354b")
    inset.axhspan(-extent, extent, facecolor="#ed9a43", alpha=.6)
    inset.set_xticks([-c["a"], c["a"]])
    history.plot(times, H[:, 0]*1000, c="#168d9d", label="H слева")
    history.plot(times, H[:, 1]*1000, c="#7d2b65", label="H справа")
    history.set(xlim=(0, times[-1]), ylim=(-extent*1000, extent*1000), xlabel="Физическое время, с", ylabel="H, мм")
    history.legend(loc="center right", ncol=2)
    cursor = history.axvline(0, color="black", lw=1)
    fig.text(.09, .93, LABEL, fontsize=28, weight="bold", color="#b1372b")
    fig.text(.09, .885, "Контрольный линейный FEM: концы закреплены. Новое касание здесь не вычисляется.", fontsize=17)
    timestamp = fig.text(.78, .89, "", fontsize=19, weight="bold")
    fig.text(.09, .035, FILM_LABEL + ". Объём и обратное влияние слоя пренебрежимо малы в ведущем порядке.", fontsize=13)
    fig.text(.80, .21, "no-slip; σ = 0\nν = 0.01 м²/с\nα = 2°\nГлубина: 10 м", fontsize=14)
    fill = [None]
    frame_map = []
    representative = set(np.rint(np.linspace(0, len(frame_ids)-1, 5)).astype(int))
    preview_paths = []
    writer = FFMpegWriter(fps=25, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "20", "-threads", "2"])
    video = output / "preview_reference.mp4"
    with writer.saving(fig, str(video), 100):
        for frame_number, (t, i) in enumerate(zip(frame_times, frame_ids)):
            yy = surface_value(x, eta[i], dense_x)*1000
            line.set_data(dense_x, yy)
            if fill[0] is not None:
                fill[0].remove()
            fill[0] = ax.fill_between(dense_x, -extent*1000, yy, color="#b3d8f2", zorder=1)
            for side in (0, 1):
                strips[side].set_height((H[i, side]+extent)*1000)
            cursor.set_xdata([t, t])
            timestamp.set_text(f"t = {t:.2f} с")
            writer.grab_frame()
            frame_map.append({"frame": frame_number, "time_s": float(t), "state_id": data["ids"][i],
                "H_m": H[i].tolist(), "wet_intervals": data["intervals"][i].tolist(),
                "markers": [{"marker_id": m["marker_id"], "side": m["side"], "height_m": m["height_m"]} for m in markers],
                "xlim": list(ax.get_xlim()), "ylim_mm": list(ax.get_ylim())})
            if frame_number in representative:
                path = output / f"frame_{frame_number:04d}.png"
                fig.savefig(path, dpi=100)
                preview_paths.append(path)
            if frame_number % 10 == 0:
                heartbeat()
    plt.close(fig)
    from PIL import Image
    thumbs = []
    for path in preview_paths:
        with Image.open(path) as original:
            thumbs.append(original.resize((960, 540)))
    montage = Image.new("RGB", (1920, 540*3), "white")
    for i, thumb in enumerate(thumbs):
        montage.paste(thumb, ((i % 2)*960, (i // 2)*540))
    montage.save(output / "key_frames.png")
    probe = probe_video(video)
    subprocess.run(["ffmpeg", "-v", "error", "-xerror", "-i", str(video), "-f", "null", "-"], check=True,
                   capture_output=True, timeout=120)
    write_json(output / "render_evidence.json", {"status": "DECODED_REFERENCE", "probe": probe,
        "full_decode": True, "frame_count": len(frame_ids), "frame_sample_dt_s": .01,
        "physical_time_range_s": [float(times[0]), float(times[-1])], "film_label": FILM_LABEL,
        "label": LABEL, "fixed_axes": True, "frame_map": frame_map})
    report = f"""# PW1 — NEEDS_MODEL_DECISION

Получен реальный контрольный FEM preview {times[0]:g}–{times[-1]:g} с, **REFERENCE / NOT TARGET**.
Это не final_animation.mp4 и не квалифицированный расчёт нового смачивания.

Посчитано: линейное вязкое течение, P2/P1, SDIRK2; a={c['a']} м, d={c['d']} м,
g={c['g']} м/с², ν={c['nu']} м²/с, α={c['alpha_deg']}°, σ=0, no-slip,
сетка {c['nx']}×{c['nz']}, dt={c['dt']} с. Начальная поверхность прямая, скорость нулевая.
Каждый из {len(times)} accepted states содержит velocity/eta/H и связные интервалы W.
Слой входит в физическое состояние как L_macro ∪ W; его объём, инерция,
обратное влияние и adsorption work пренебрегаются в ведущем порядке. Толщина неизвестна.

Контрольные концы: P1={eta[0,0]:.15g} м, P2={eta[0,-1]:.15g} м, неизменны по уравнениям.
P3/P4 не созданы. Это свойство fixed-endpoint reference, а не вывод о максимумах целевой модели.
Новый целевой максимум/время: NOT RUN. Целевой горизонт 5 с и уточнение: NOT RUN.
Сильная точечная дивергенция не тождественно нулевая в mixed FEM; проверяется слабая норма.
Наибольший наклон в reference: {max(r['max_slope'] for r in data['diagnostics']):.5g}.
При больших пристеночных наклонах линейная геометрия не квалифицируется локально.

Один blocker: **PW1_CONTACT_CLOSURE**. Требуется закон переноса/релаксации нового контакта
и его масштаб, согласованные с no-slip, σ=0, массой и энергией. Память уже смоченных
точек этот закон не заменяет. Доказательства двух проверенных путей — model_decision.md
и blocker_evidence.json. Диффузное смачивание при no-slip не объявляется невозможным.

Воспроизведение: `python sloshing_visualization/scripts/pinned_wetting_mission.py blocker-check`.
Продолжение существующего reference: `python sloshing_visualization/scripts/pinned_wetting_mission.py resume`.
Новый reference с более мелкой сеткой/временем: `python sloshing_visualization/scripts/pinned_wetting_mission.py run --mesh-level 1 --dt-scale 0.5`.
Последняя команда уточняет только reference и не снимает физический blocker.

Native данные: `{native}`. Run ID: `{run_id}`.
Новые внешние расходы и push отсутствуют. Исторические результаты не изменялись.
"""
    (output / "FINAL_REPORT.md").write_text(report)
    links = ["surface_history.csv", "wetting_history.csv", "marker_catalog.csv", "marker_tracks.csv", "energy_history.csv",
             "resolved_case.json", "model_decision.md", "blocker_evidence.json", "runup_summary.json", "refinement.json",
             "key_frames.png", "FINAL_REPORT.md", "verification.json", "manifest.json"]
    (output / "index.html").write_text('<!doctype html><html lang="ru"><meta charset="utf-8">'
        '<title>PW1: контрольное preview</title><style>body{max-width:1100px;margin:40px auto;font:18px sans-serif;background:#f6f8fc}video{width:100%}li{margin:8px}</style>'
        '<h1>NEEDS_MODEL_DECISION</h1><p>REFERENCE / NOT TARGET. Новое касание не квалифицировано; финальный результат отсутствует.</p>'
        '<video controls preload="metadata" src="preview_reference.mp4"></video><p>' + FILM_LABEL + '</p><ul>' +
        ''.join(f'<li><a href="{name}">{name}</a></li>' for name in links) + '</ul></html>')
    artifacts = [{"path": str(p.relative_to(output)), "sha256": sha256(p), "role": "reference_evidence"}
                 for p in sorted(output.iterdir()) if p.is_file() and p.name not in ("manifest.json", "verification.json")]
    write_json(output / "manifest.json", {"schema_version": "1.0", **data["identity"],
        "status": "PREVIEW_ONLY", "contact_method": data["contact_method"], "source_kind": "real_reference_solver",
        "artifacts": artifacts, "native_data_paths": [str(native)], "native_sha256": sha256(native),
        "time_range": [float(times[0]), float(times[-1])], "resolved_case": "resolved_case.json",
        "frame_map": "render_evidence.json#frame_map", "contact_source_map": "native.h5:/accepted/<step>/eta endpoints",
        "units": {"length": "m", "time": "s", "energy": "J/m span, converted from per density with declared rho"},
        "limitations": ["Fixed endpoints; no advancing wall-contact prediction", "Linear geometry; local corner slopes may exceed small-slope validity",
                        "Zero-volume one-way retained coating", "Target horizon/refinement/final release NOT RUN"]})
