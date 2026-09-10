"""The existing PW1 native/CSV/video path adapted to moving Navier endpoints."""
import json
from pathlib import Path
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image

from ..validation.contact import surface_value
from ..visualization.render import probe_video
from .output import csv_write, FILM_LABEL
from .reference import load_native
from .navier import LABEL
from .mission import read_json, write_json, sha256


def render_bundle(root, native, output, *, final=False, heartbeat=lambda: None):
    native, output = Path(native), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    signature = {"native_sha256": sha256(native), "renderer_sha256": sha256(Path(__file__)), "final": final}
    if (output / "render_manifest.json").exists():
        old = read_json(output / "render_manifest.json")["signature"]
        if old == signature:
            return
        if old["native_sha256"] != signature["native_sha256"] or old["final"] != final:
            raise ValueError("Different native data/profile requires another output directory")
        # Only a renderer revision on identical immutable data may replace a
        # candidate. Preserve every previous root artifact before regeneration.
        archive=output/"render_history"/old["renderer_sha256"][:12]
        archive.mkdir(parents=True,exist_ok=True)
        for previous in output.iterdir():
            if previous.is_file() and not (archive/previous.name).exists():
                shutil.copy2(previous,archive/previous.name)
        write_json(archive/"ARCHIVE_NOTICE.json",{"status":"SUPERSEDED_RENDER_VERSION_NOT_CURRENT_RELEASE",
            "original_output":str(output),"signature":old,"native_data_changed":False,
            "reason":"Separate nearby historical marker captions without moving world coordinates."})
    d = load_native(native)
    c, times, x, eta, H = d["config"], d["times"], d["x"], d["eta"], d["H"]
    resolved = read_json(native.parent / "resolved_case.json")
    rho = resolved["original_case"]["physics"]["density_kg_m3"]
    frame_times = np.arange(round(times[-1]/.01)+1)*.01
    indices = np.rint(frame_times/c["dt"]).astype(int)
    if not np.allclose(times[indices], frame_times, rtol=0, atol=1e-12):
        raise ValueError("Render cadence is not on accepted solver states")
    for name in ("resolved_case.json",):
        shutil.copyfile(native.parent / name, output / name)
    shutil.copyfile(root / "mission/pinned_wetting/extension/model_decision.md", output / "model_decision.md")
    shutil.copyfile(root / "mission/pinned_wetting/extension/ACCEPTANCE_ADDENDUM.yaml", output / "ACCEPTANCE_ADDENDUM.yaml")
    header = ["time_s", "point_index", "x_m", "z_m", "component_id", "state_id"]
    def surface_rows(selected):
        return ((times[i], j, xx, zz, "main", d["ids"][i]) for i in selected for j, (xx, zz) in enumerate(zip(x, eta[i])))
    csv_write(output / "surface_history.csv", header, surface_rows(range(len(times))))
    slice_ids = [int(round(t/c["dt"])) for t in np.arange(0., times[-1]+1e-12, .5)]
    csv_write(output / "surface_snapshots.csv", header, surface_rows(slice_ids))
    csv_write(output / "wetting_history.csv",
        ["time_s", "R_L_m", "R_R_m", "reach_L_m", "reach_R_m", "H_L_m", "H_R_m", "state_id"],
        ((t, *eta[i,[0,-1]], *eta[i,[0,-1]], *H[i], d["ids"][i]) for i, t in enumerate(times)))
    markers = d["coatings"][-1]["markers"]
    csv_write(output / "marker_catalog.csv", ["marker_id", "side", "height_m", "created_at_s", "state_id"],
        ([m[k] for k in ("marker_id", "side", "height_m", "created_at_s", "state_id")] for m in markers))
    csv_write(output / "marker_tracks.csv", ["time_s", "marker_id", "side", "height_m"],
        ((float(t), m["marker_id"], m["side"], m["height_m"]) for t in frame_times for m in markers if m["created_at_s"] <= t+1e-12))
    energy_keys = ["time_s", "kinetic_J_per_m", "potential_J_per_m", "capillary_J_per_m", "wall_J_per_m",
        "physical_dissipation_integral_J_per_m", "external_work_J_per_m", "numerical_energy_remainder_J_per_m",
        "budget_defect_J_per_m", "state_id", "bulk_viscous_dissipation_integral_J_per_m", "wall_friction_dissipation_integral_J_per_m"]
    csv_write(output / "energy_history.csv", energy_keys,
        ((r["time"], rho*r["kinetic_energy"], rho*r["gravitational_potential_energy"], 0., 0.,
          rho*r["cumulative_viscous_dissipation"], rho*r["cumulative_work"], rho*r["rk_energy_correction"],
          rho*r["energy_balance_residual"], d["ids"][i], rho*r["cumulative_bulk_dissipation"], rho*r["cumulative_wall_dissipation"])
         for i, r in enumerate(d["diagnostics"])))
    peaks = [{**m, "peak_time_s": float(times[int(m["state_id"].rsplit(":",1)[1])])} for m in markers if m["marker_id"] not in ("P1","P2")]
    left_peaks = [p for p in peaks if p["side"] == "L"]
    write_json(output / "runup_summary.json", {"model": LABEL, "run_id": d["identity"]["run_id"],
        "first_left_peak": left_peaks[0] if left_peaks else None,
        "later_left_new_record": left_peaks[1] if len(left_peaks)>1 else None,
        "no_later_record_reason": "No later resolved local peak exceeded the retained record plus the predeclared 0.1 mm detector floor." if len(left_peaks)<2 else None,
        "H_final_m": H[-1].tolist(), "R_final_m": eta[-1,[0,-1]].tolist(),
        "first_left_rise_from_P1_m": left_peaks[0]["height_m"]-eta[0,0] if left_peaks else None,
        "left_record_exceeds_initial_P2": bool(H[-1,0] > eta[0,-1]), "records": peaks,
        "uncertainty": "See independent mesh/time comparisons in refinement.json; two levels are not an error bound."})
    names = ["final_animation.mp4", "wall_detail.mp4"] if final else ["preview_extension.mp4"]
    render_records = []
    for name in names:
        detail = name == "wall_detail.mp4"
        record = render_video(output/name, d, indices, frame_times, detail=detail, heartbeat=heartbeat)
        render_records.append(record)
    expected_frames = len(indices)
    for video in names:
        probe = probe_video(output/video)
        subprocess.run(["ffmpeg","-v","error","-xerror","-i",str(output/video),"-f","null","-"],
                       check=True, capture_output=True, timeout=120)
        if (probe["avg_frame_rate"] != "25/1" or int(probe["nb_read_frames"]) != expected_frames or
                (probe["width"], probe["height"]) != (1920,1080)):
            raise ValueError("Video cadence/size mismatch")
        record = next(r for r in render_records if r["video"] == video)
        record["probe"], record["full_decode"] = probe, True
    write_json(output / "render_evidence.json", {"model": LABEL, "film_label": FILM_LABEL,
        "time_range_s": [float(times[0]), float(times[-1])], "frame_dt_s": .01, "videos": render_records})
    frames = sorted(output.glob("main_frame_*.png"))
    montage = Image.new("RGB", (1920, 540*((len(frames)+1)//2)), "white")
    for n, p in enumerate(frames):
        with Image.open(p) as im:
            montage.paste(im.resize((960,540)), ((n%2)*960,(n//2)*540))
    montage.save(output / "key_frames.png")
    write_json(output / "render_manifest.json", {"signature": signature, "run_id": d["identity"]["run_id"],
        "native": str(native), "videos": names, "frame_count": expected_frames,
        "status": "RENDERED_NOT_YET_FINAL_VERIFIED" if final else "PREVIEW_ONLY"})
    links = ["surface_history.csv", "surface_snapshots.csv", "wetting_history.csv", "marker_catalog.csv", "marker_tracks.csv",
             "energy_history.csv", "runup_summary.json", "key_frames.png", "model_decision.md", "resolved_case.json"]
    (output/"index.html").write_text('<!doctype html><html lang="ru"><meta charset="utf-8"><title>PW1 — Navier</title>'
        '<style>body{max-width:1150px;margin:40px auto;font:18px sans-serif;background:#f7f9fc}video{width:100%}li{margin:10px}</style>'
        f'<h1>{LABEL}</h1><p>Основная жидкость: линейная модель, b={c["slip_length_m"]} м, σ=0. '
        'Остаточное покрытие — отдельный необратимый подсеточный закон; толщина и отдельная динамика не вычислялись.</p>'
        ''.join(f'<video controls preload="metadata" src="{n}"></video>' for n in names) + '<ul>' +
        ''.join(f'<li><a href="{n}">{n}</a></li>' for n in links) + '</ul></html>')


def render_video(path, d, indices, times, *, detail, heartbeat):
    c, x, eta, H = d["config"], d["x"], d["eta"], d["H"]
    markers = d["coatings"][-1]["markers"]
    fig = plt.figure(figsize=(19.2,10.8), dpi=100, facecolor="#f7f9fc")
    axes = ([fig.add_axes([.08,.38,.25,.43]),fig.add_axes([.41,.38,.25,.43]),
             fig.add_axes([.735,.38,.12,.43])] if detail else [fig.add_axes([.08,.38,.73,.43])])
    history = fig.add_axes([.08,.14,.73,.16])
    inset = fig.add_axes([.91,.38,.055,.43] if detail else [.87,.35,.07,.45])
    limits = (-.048, .048)
    if np.max(abs(eta)) > .045:
        raise ValueError("Fixed display headroom exceeded; preserve fields and choose a new renderer range")
    dense = np.unique(np.concatenate([np.linspace(l,r,9) for l,r in zip(x[::2][:-1],x[::2][1:])]))
    lines, fills, strips, artists = [], [], [], []
    for axis_index, ax in enumerate(axes):
        ylimits=np.array(limits)*1000
        if detail:
            xl = (.72,1.06) if axis_index==1 else (-1.06,-.72)
            title = "Правая стенка" if axis_index==1 else "Левая стенка"
            if axis_index==2:
                peaks=[m["height_m"]*1000 for m in markers if m["side"]=="L" and m["marker_id"]!="P1"]
                xl=(-1.02,-.94)
                if peaks:
                    ylimits=np.array([min(peaks)-1.,max(peaks)+1.])
                title="Рекорды слева\nувеличение по z"
        else:
            xl, title = (-1.12,1.12), "Вычисленная поверхность; разные масштабы x и z"
        ax.set(xlim=xl, ylim=ylimits, xlabel="x, м", ylabel="z, мм", title=title)
        if detail and axis_index==2:
            ax.title.set_fontsize(11)
            ax.tick_params(labelsize=10)
        ax.grid(alpha=.2)
        for initial_height in eta[0,[0,-1]]:
            ax.axhline(initial_height*1000,c="#84235e",ls="--",lw=.8,alpha=.35,zorder=2)
        for side in (-c["a"],c["a"]):
            ax.axvline(side,c="#223653",lw=3)
        line, = ax.plot([],[],c="#1265ab",lw=3,zorder=5)
        lines.append(line); fills.append(None)
        ax.plot(x,eta[0]*1000,"--",c="#8395ad",lw=1)
        side_strips = []
        for j,side in enumerate((-c["a"],c["a"])):
            strip = Rectangle((side if j==0 else side-.012, limits[0]*1000), .012, 0, color="#15bac3",zorder=6)
            ax.add_patch(strip); side_strips.append(strip)
        strips.append(side_strips)
        ma = []
        for m in markers:
            xx = -c["a"] if m["side"]=="L" else c["a"]
            shape="D" if m["marker_id"]=="P4" else "o"
            point, = ax.plot([xx],[m["height_m"]*1000],marker=shape,linestyle="none",c="#84235e",
                            mfc="white" if m["marker_id"]=="P3" else "#84235e",ms=8,zorder=9)
            neighbors=sorted([q for q in markers if q["side"]==m["side"] and abs(q["height_m"]-m["height_m"])<.004],
                             key=lambda q:q["height_m"])
            dy=12
            if len(neighbors)>1 and not (detail and axis_index==2):
                rank=next(k for k,q in enumerate(neighbors) if q["marker_id"]==m["marker_id"])
                dy=24*(2*rank-len(neighbors)+1)
            text = ax.annotate(f'{m["marker_id"]}: {m["height_m"]*1000:+.2f} мм', (xx,m["height_m"]*1000),
                               xytext=(12 if xx<0 else -12,dy),textcoords="offset points",ha="left" if xx<0 else "right",
                               fontsize=10 if detail else 12,c="#84235e",zorder=10,
                               arrowprops={"arrowstyle":"-","lw":.8,"color":"#84235e"} if len(neighbors)>1 else None)
            ma.append((m,point,text))
        artists.append(ma)
    history.plot(d["times"],eta[:,0]*1000,"--",c="#176aab",alpha=.6,label="R слева — bulk")
    history.plot(d["times"],eta[:,-1]*1000,"--",c="#d78329",alpha=.6,label="R справа — bulk")
    history.plot(d["times"],H[:,0]*1000,c="#147e92",lw=2,label="H слева — верх покрытия")
    history.plot(d["times"],H[:,1]*1000,c="#8a3264",lw=2,label="H справа — верх покрытия")
    history.set(xlim=(0,d["times"][-1]),ylim=np.array(limits)*1000,xlabel="Физическое время, с",ylabel="Высота, мм")
    history.grid(alpha=.2); history.legend(loc="lower right",ncol=2,fontsize=10)
    cursor = history.axvline(0,c="black",lw=1)
    inset.set(xlim=(-1.2,1.2),ylim=(-c["d"],.2),title="Сосуд, м")
    inset.set_aspect("equal"); inset.fill_between([-1,1],-c["d"],0,color="#bcdbef")
    inset.plot([-1,-1,1,1],[.1,-c["d"],-c["d"],.1],c="#243b55")
    inset.axhspan(*limits,color="#dd9543",alpha=.7); inset.set_xticks([-1,1])
    fig.text(.08,.94,LABEL,fontsize=22,weight="bold",c="#184265")
    fig.text(.08,.89,f'Линейные гравитационные колебания • b = {c["slip_length_m"]:.2f} м • σ = 0 • ν = {c["nu"]} м²/с',fontsize=16)
    if path.name not in ("final_animation.mp4", "wall_detail.mp4"):
        fig.text(.08,.845,"PREVIEW / NOT FINAL — уточнение и итоговая проверка ещё не завершены",fontsize=13,c="#a14b20")
    else:
        fig.text(.08,.845,"Пунктир: начальная прямая и уровни P1/P2. В увеличениях масштабы x (м) и z (мм) различны.",fontsize=12,c="#526277")
    timer = fig.text(.82,.89,"",fontsize=20,weight="bold")
    values = fig.text(.08,.325,"",fontsize=14)
    fig.text(.08,.065,FILM_LABEL,fontsize=14,c="#28596b")
    fig.text(.08,.032,"Покрытие задано подсеточно: объём, инерция и отдельная динамика не вычислялись. Навье и память покрытия — разные законы.",fontsize=12)
    fig.text(.86,.23,"Боковые стенки:\nнепроницаемы,\nNavier slip.\nДно: no-slip.\nЗамедление ×4.",fontsize=12)
    representative = set(np.rint(np.linspace(0,len(indices)-1,6)).astype(int))
    for marker in markers[2:]:
        representative.add(min(len(indices)-1,int(np.ceil(marker["created_at_s"]/.01-1e-10))))
    frame_map=[]
    writer=FFMpegWriter(fps=25,codec="libx264",extra_args=["-pix_fmt","yuv420p","-crf","20","-threads","2"])
    with writer.saving(fig,str(path),100):
        for frame,(t,i) in enumerate(zip(times,indices)):
            y=surface_value(x,eta[i],dense)*1000
            for a,ax in enumerate(axes):
                lines[a].set_data(dense,y)
                if fills[a] is not None: fills[a].remove()
                fills[a]=ax.fill_between(dense,limits[0]*1000,y,color="#b6d9f1",zorder=1)
                for j in (0,1): strips[a][j].set_height((H[i,j]-limits[0])*1000)
                for m,point,text in artists[a]:
                    point.set_visible(m["created_at_s"]<=t+1e-12); text.set_visible(m["created_at_s"]<=t+1e-12)
            timer.set_text(f"t = {t:.2f} с"); cursor.set_xdata([t,t])
            values.set_text(f'Слева: R = {eta[i,0]*1000:+.2f} мм; H = {H[i,0]*1000:+.2f} мм     Справа: R = {eta[i,-1]*1000:+.2f} мм; H = {H[i,1]*1000:+.2f} мм')
            writer.grab_frame()
            if frame in representative:
                prefix="detail" if detail else "main"
                fig.savefig(path.parent/f"{prefix}_frame_{frame:04d}.png",dpi=100)
            frame_map.append({"frame":frame,"time_s":float(t),"state_id":d["ids"][i],"H_m":H[i].tolist(),
                "wet_intervals":d["intervals"][i].tolist(),"markers":[m for m in markers if m["created_at_s"]<=t+1e-12],
                "axes":[{"xlim":list(ax.get_xlim()),"ylim_mm":list(ax.get_ylim())} for ax in axes]})
            if frame%10==0: heartbeat()
    plt.close(fig)
    return {"video":path.name,"frame_map":frame_map,"fixed_world_markers":True,"fixed_axes":True}
