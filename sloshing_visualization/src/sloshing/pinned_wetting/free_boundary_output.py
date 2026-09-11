"""Native parametric-interface export and the existing matplotlib/FFmpeg path."""
import json
from pathlib import Path
import shutil
import subprocess
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Polygon
from PIL import Image
from .mission import read_json, write_json, sha256
from .output import csv_write
from ..visualization.render import probe_video


def curve(points, count=16):
    t = np.linspace(0.,1.,count+1)
    phi = np.array([(1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)])
    pieces = [points[:,j:j+3]@phi for j in range(0,points.shape[1]-2,2)]
    return np.hstack([p[:,:-1] for p in pieces]+[pieces[-1][:,-1:]])


def curve_bounds(surfaces):
    """Exact coordinate extrema of every native quadratic edge, all times."""
    left=surfaces[:,:,0:-2:2];mid=surfaces[:,:,1::2];right=surfaces[:,:,2::2]
    aa=2*(left+right-2*mid);bb=-3*left+4*mid-right
    with np.errstate(divide="ignore",invalid="ignore"):
        location=-bb/(2*aa)
        value=aa*location**2+bb*location+left
    valid=(location>0)&(location<1)&np.isfinite(value)
    lower=np.minimum(left,right);upper=np.maximum(left,right)
    lower=np.minimum(lower,np.where(valid,value,np.inf))
    upper=np.maximum(upper,np.where(valid,value,-np.inf))
    return lower.min(axis=(0,2)),upper.max(axis=(0,2))


def tessellation_bound(surfaces,count=32):
    deviation=surfaces[:,:,1::2]-.5*(surfaces[:,:,0:-2:2]+surfaces[:,:,2::2])
    return float(np.max(np.sqrt(np.sum(deviation**2,axis=1))))/count**2


def load(native):
    with h5py.File(native,"r") as h:
        ids = h["topology/interface_nodes"][:]
        keys = sorted(h["states"])
        times=[]; surface=[]; coatings=[]; rows=[]; trace=[]
        components = h["topology/component_dofs"][:]
        X0 = h["topology/initial_geometry"][:]
        right = np.flatnonzero(X0[0]==1.); left=np.flatnonzero(X0[0]==-1.)
        bottom = np.flatnonzero(X0[1]==-10.)
        for key in keys:
            s=h["states"][key]; times.append(float(s.attrs["time_s"]))
            surface.append(s["geometry"][:][:,ids])
            coatings.append(json.loads(s.attrs["coating"]))
            rows.append(json.loads(s.attrs["diagnostics"]))
            v=s["velocity"][:]
            trace.append([np.max(abs(v[components[0,right]])),np.max(abs(v[components[1,right]])),
                np.max(np.sqrt(np.sum(v[components[:,bottom]]**2,axis=0))),np.max(abs(v[components[0,left]]))])
        return {"times":np.array(times),"surface":np.array(surface),"coatings":coatings,
                "rows":rows,"trace":np.array(trace),"run_id":str(h.attrs["run_id"]),
                "model_id":str(h.attrs["model_id"]),"native":str(native)}


def export(native, output, d):
    output.mkdir(parents=True,exist_ok=True)
    for name in ("resolved_case.json","model_notes.md"):
        shutil.copyfile(Path(native).parent/name,output/name)
    times, surfaces, coatings = d["times"],d["surface"],d["coatings"]
    csv_write(output/"contacts.csv",["step","time_s","R_left_m","R_right_m","H_left_m","H_right_m",
        "right_u_max_m_s","right_w_max_m_s","bottom_speed_max_m_s","left_normal_max_m_s"],
        ((i,t,*surfaces[i,1,[0,-1]],*coatings[i]["H"],*d["trace"][i]) for i,t in enumerate(times)))
    csv_write(output/"retained_intervals.csv",["step","side","z_bottom_m","z_top_m"],
        ((i,side,-10.,coat["H"][j]) for i,coat in enumerate(coatings) for j,side in enumerate(("L","R"))))
    # Full quadratic parametric curve; no projection to a single-valued eta(x).
    selected = np.unique(np.r_[0,np.arange(0,len(times),max(1,len(times)//100)),len(times)-1,
        [np.argmin(abs(times-t)) for t in np.arange(0.,times[-1]+1e-12,.5)]]).astype(int)
    csv_write(output/"interface_samples.csv",["time_s","vertex_id","x_m","z_m"],
        ((times[i],j,*point) for i in selected for j,point in enumerate(curve(surfaces[i],32).T)))
    markers = coatings[-1]["markers"]
    csv_write(output/"marker_catalog.csv",["marker_id","side","x_m","z_m","created_at_s","source_step"],
        ((m["marker_id"],m["side"],-1. if m["side"]=="L" else 1.,m["height_m"],m["created_at_s"],
          int(m["state_id"].rsplit(":",1)[1])) for m in markers))
    csv_write(output/"marker_tracks.csv",["step","time_s","marker_id","side","x_m","z_m"],
        ((i,t,m["marker_id"],m["side"],-1. if m["side"]=="L" else 1.,m["height_m"])
         for i,t in enumerate(times) for m in coatings[i]["markers"]))
    energy_keys=["energy","kinetic","potential","volume_m2","cumulative_bulk_loss","cumulative_left_wall_loss",
        "cumulative_geometry_kinetic_work","cumulative_geometry_potential_work","continuous_energy_error",
        "bulk_loss","left_wall_loss","raw_geometry_kinetic_work","raw_geometry_potential_work",
        "geometric_divergence_work","hydrostatic_split_work","geometry_kinetic_work","geometry_potential_work",
        "rezone_work","cumulative_rezone_work","numerical_grad_div_work",
        "cumulative_numerical_grad_div_work","split_energy_error"]
    csv_write(output/"energy_history.csv",["step","time_s",*energy_keys],
        ((i,t,*[d["rows"][i].get(k,0.) for k in energy_keys]) for i,t in enumerate(times)))
    write_json(output/"runup_summary.json",{"run_id":d["run_id"],"model_id":d["model_id"],
        "time_s":float(times[-1]),"R_final_m":surfaces[-1,1,[0,-1]].tolist(),"H_final_m":coatings[-1]["H"],
        "markers":[{**m,"peak_time_s":float(times[int(m["state_id"].rsplit(":",1)[1])])} for m in markers],
        "P4_exists":any(m["marker_id"]=="P4" for m in markers),"qualified_final":False})
    with h5py.File(native,"r") as source:
        pressure_variable=str(source.attrs.get("pressure_variable","gauge_p"))
    write_json(output/"native_manifest.json",{"model_id":d["model_id"],"run_id":d["run_id"],
        "path":str(native),"sha256":sha256(native),"format":"HDF5 full P2 material geometry/velocity at every accepted step; pressure_midpoint for each physical step, startup pressure at state0",
        "time_range_s":[0.,float(times[-1])],"accepted_states":len(times),"coating":"each accepted state's coating attribute",
        "film_thickness":None,"retained_layer_dynamics":"zero-volume one-way; separate feedback not computed",
        "pressure_variable":pressure_variable,"physical_pressure_per_density":"p=pi-g*z" if pressure_variable=="pi=p+gz" else "stored gauge_p"})
    write_json(output/"export_geometry.json",{"parametric_P2_curve_retained":True,"subdivisions_per_quadratic_edge":32,
        "exact_chord_error_bound_m":tessellation_bound(surfaces),
        "coordinate_extrema_m":[q.tolist() for q in curve_bounds(surfaces)],
        "energy_units":"per unit density and per out-of-plane depth; multiply by rho=1000 for J/m",
        "signed_numerical_work_is_physical_heat":False})


def render(native, output, heartbeat=lambda:None, final=False):
    d=load(native); output=Path(output)
    if final and abs(d["times"][-1]-5.)>1e-12: raise ValueError("Cannot label partial native as final")
    if np.max(abs(d["surface"][:,1,-1]-d["surface"][0,1,-1]))>1e-12:
        raise ValueError("Native P2 moved; renderer refuses repair")
    export(native,output,d)
    numerical_controls=read_json(output/"resolved_case.json")["numerics"]
    gamma=float(numerical_controls.get("grad_div_gamma_m2_s",0.))
    frame_times=np.arange(int(np.floor(d["times"][-1]/.01+1e-10))+1)*.01
    indices=np.array([int(np.argmin(abs(d["times"]-t))) for t in frame_times])
    if np.max(abs(d["times"][indices]-frame_times))>1e-10:
        raise ValueError("Video samples must correspond to actual accepted times")
    H=np.array([c["H"] for c in d["coatings"]]); S=d["surface"]
    lower,upper=curve_bounds(S)
    low=min(float(lower[1])-.008,-.043);high=max(float(upper[1])+.008,.043)
    reports=[]
    for detail in (False,True):
        name=("wall_detail.mp4" if detail else "final_animation.mp4") if final else ("wall_detail_preview.mp4" if detail else "preview.mp4")
        fig=plt.figure(figsize=(19.2,10.8),dpi=100,facecolor="#f7f9fc")
        positions=([(.075,.40,.39,.40),(.535,.40,.39,.40)] if detail else [(.075,.40,.73,.40)])
        axes=[fig.add_axes(pos) for pos in positions]
        panels=[]
        for i,ax in enumerate(axes):
            ax.set(xlim=((-1.02,-.75) if i==0 else (.75,1.02)) if detail else (-1.07,1.07),
                ylim=(low*1000,high*1000),xlabel="x, м",ylabel="z, мм",
                title=("LEFT — Navier" if i==0 else "RIGHT — no-slip, P2 fixed") if detail else "Полная вычисленная граница; разные масштабы x/z")
            ax.grid(alpha=.2)
            ax.axvline(-1,color="#263f54",lw=3);ax.axvline(1,color="#263f54",lw=3)
            ax.plot(S[0,0],S[0,1]*1000,"--",color="#8595a4",lw=1)
            line,=ax.plot([],[],color="#086aae",lw=2.5,zorder=5)
            area=Polygon(np.zeros((3,2)),facecolor="#b8dbee",edgecolor="none",zorder=1);ax.add_patch(area)
            retained,=ax.plot([],[],color="#17bbc2",lw=7,zorder=7,solid_capstyle="butt")
            dotH,=ax.plot([],[],marker="_",ms=17,color="#00949f",zorder=8)
            marks=[]
            for m in d["coatings"][-1]["markers"]:
                x=-1 if m["side"]=="L" else 1; z=m["height_m"]*1000
                point,=ax.plot([x],[z],"o",color="#8b2862",ms=6,zorder=9)
                text=ax.annotate(m["marker_id"],(x,z),xytext=(10 if x<0 else -10,8),
                    textcoords="offset points",ha="left" if x<0 else "right",color="#8b2862",fontsize=12)
                marks.append((m,point,text))
            panels.append((line,area,retained,dotH,marks))
        hist=fig.add_axes([.075,.16,.73,.15]); hist.grid(alpha=.2)
        hist.plot(d["times"],S[:,1,0]*1000,label="R_L — текущий контакт",color="#086aae")
        hist.plot(d["times"],H[:,0]*1000,label="H_L — история",color="#16a6af")
        hist.plot(d["times"],H[:,1]*1000,label="R_R = H_R = P2",color="#8b2862")
        hist.set(xlim=(0,max(.01,d["times"][-1])),ylim=(low*1000,high*1000),xlabel="Физическое время, с",ylabel="Высота, мм")
        hist.legend(fontsize=11,ncol=3);cursor=hist.axvline(0,color="black")
        vessel=fig.add_axes([.86,.39,.065,.4]) if not detail else None
        if vessel:
            vessel.set(xlim=(-1.2,1.2),ylim=(-10.1,.15),title="Сосуд, м");vessel.set_aspect("equal")
            vessel.plot([-1,-1,1,1],[.1,-10,-10,.1],color="#263f54")
            whole=Polygon(np.zeros((3,2)),facecolor="#b8dbee",edgecolor="#086aae",lw=.7);vessel.add_patch(whole)
        fig.text(.075,.945,"PW2 — FULL NONLINEAR FREE-BOUNDARY NAVIER–STOKES",fontsize=21,weight="bold",color="#21455f")
        fig.text(.075,.89,"DECLARED_EXTENSION • LEFT Navier b=0.50 м • RIGHT/BOTTOM no-slip • σ=0",fontsize=17)
        fig.text(.075,.845,"PREVIEW — НЕ ФИНАЛЬНЫЙ РЕЗУЛЬТАТ" if not final else "Полная материальная граница; P2 — неподвижный текущий правый контакт",fontsize=14,color="#ad5427" if not final else "#21455f")
        clock=fig.text(.82,.32,"",fontsize=18);values=fig.text(.075,.34,"",fontsize=13)
        fig.text(.075,.075,"Остаточное покрытие задано подсеточно; ширина полосы условна, толщина не рассчитывалась.",fontsize=14,color="#197882")
        fig.text(.075,.04,"Объём, инерция, отдельная динамика и feedback покрытия не вычислялись. Разрешённая жидкость учитывается целиком. Замедление ×4.",fontsize=11)
        if gamma:
            fig.text(.075,.012,f"Численная grad-div стабилизация γ={gamma:g} м²/с; её работа учитывается отдельно от физического тепла.",fontsize=10,color="#596d7c")
        frame_map=[]; key=set(np.rint(np.linspace(0,len(indices)-1,min(6,len(indices)))).astype(int))
        required_times=[0.,1.,2.5,5.]
        first_peak=next((m for m in d["coatings"][-1]["markers"] if m["marker_id"]=="P3"),None)
        if first_peak is not None:
            required_times.append(float(d["times"][int(first_peak["state_id"].rsplit(":",1)[1])]))
        right_mask=S[0,0]>.75
        deformation=np.max(np.linalg.norm(S[:,:,right_mask]-S[0][:,right_mask][None,:,:],axis=1),axis=1)
        required_times.append(float(d["times"][int(deformation.argmax())]))
        for required_time in required_times:
            if required_time<=frame_times[-1]+1e-12:
                key.add(int(np.argmin(abs(frame_times-required_time))))
        writer=FFMpegWriter(fps=25,codec="libx264",extra_args=["-pix_fmt","yuv420p","-crf","20","-threads","2"])
        with writer.saving(fig,str(output/name),100):
            for frame,(t,i) in enumerate(zip(frame_times,indices)):
                q=curve(S[i],32); poly=np.vstack(([-1,low],[1,low],q.T[::-1]))
                for line,area,retained,dotH,marks in panels:
                    line.set_data(q[0],q[1]*1000);area.set_xy(poly*np.array([1,1000]))
                    retained.set_data([-1,-1],[S[i,1,0]*1000,H[i,0]*1000])
                    dotH.set_data([-1],[H[i,0]*1000])
                    for m,point,text in marks:
                        visible=m["created_at_s"]<=t+1e-12;point.set_visible(visible);text.set_visible(visible)
                if vessel: whole.set_xy(np.vstack(([-1,-10],[1,-10],q.T[::-1])))
                cursor.set_xdata([t,t]);clock.set_text(f"t={t:.2f} с")
                values.set_text(f"R_L={S[i,1,0]*1000:+.3f} мм  H_L={H[i,0]*1000:+.3f} мм    R_R=H_R=P2={H[i,1]*1000:+.3f} мм")
                if frame in key: fig.savefig(output/f"{'wall' if detail else 'main'}_frame_{frame:04d}.png",dpi=100)
                writer.grab_frame();frame_map.append({"frame":frame,"time_s":float(t),"accepted_step":int(i)})
                if frame%25==0: heartbeat()
        plt.close(fig)
        subprocess.run(["ffmpeg","-v","error","-xerror","-i",str(output/name),"-f","null","-"],check=True,capture_output=True,timeout=120)
        probe=probe_video(output/name)
        if int(probe["nb_read_frames"])!=len(indices) or (probe["width"],probe["height"])!=(1920,1080) or probe["avg_frame_rate"]!="25/1":
            raise ValueError("Video frame count, dimensions or cadence mismatch")
        decoded_checks=[]
        for frame in sorted(key):
            decoded=output/f"{'wall' if detail else 'main'}_decoded_{frame:04d}.png"
            subprocess.run(["ffmpeg","-v","error","-y","-i",str(output/name),"-vf",f"select=eq(n\\,{frame})",
                "-frames:v","1",str(decoded)],check=True,capture_output=True,timeout=120)
            with Image.open(decoded) as actual,Image.open(output/f"{'wall' if detail else 'main'}_frame_{frame:04d}.png") as expected:
                difference=np.asarray(actual.convert("RGB"),dtype=float)-np.asarray(expected.convert("RGB"),dtype=float)
                rms=float(np.sqrt(np.mean(difference**2))/255.)
            if rms>.03:raise ValueError("Decoded key frame does not match the native-rendered frame")
            decoded_checks.append({"frame":int(frame),"decoded_file":decoded.name,"normalized_RGB_RMS":rms})
        reports.append({"video":name,"sha256":sha256(output/name),"full_decode":True,"probe":probe,"frame_map":frame_map,"decoded_key_frames":decoded_checks})
    frames=[output/f"main_frame_{frame:04d}.png" for frame in sorted(key)]
    montage=Image.new("RGB",(1920,540*((len(frames)+1)//2)),"white")
    for i,p in enumerate(frames):
        with Image.open(p) as im:montage.paste(im.resize((960,540)),((i%2)*960,(i//2)*540))
    montage.save(output/"key_frames.png")
    evidence={"native_sha256":sha256(native),"renderer_sha256":sha256(Path(__file__)),"final_candidate":final,
        "videos":reports,"time_range_s":[0.,float(frame_times[-1])],"model_id":d["model_id"],
        "key_frame_indices":sorted(int(q) for q in key),
        "right_deformation_keyframe_policy":"max displacement of native surface nodes initially x>0.75; nodal diagnostic, not exact continuum maximum",
        "P2_tessellation_error_bound_m":tessellation_bound(S)}
    write_json(output/"render_evidence.json",evidence)
    (output/"index.html").write_text('<!doctype html><html lang="ru"><meta charset="utf-8"><title>PW2</title><style>body{max-width:1200px;margin:30px auto;font:18px sans-serif;background:#f7f9fc}video,img{width:100%}li{margin:10px}</style>'
        '<h1>PW2 — full nonlinear moving-domain NS</h1><p id="qualification_status">'+('Кандидат, ожидающий итоговой проверки.' if final else 'PREVIEW ONLY — НЕ финальный результат.')+'</p>'
        '<p>LEFT Navier b=0.50 м; RIGHT и BOTTOM no-slip; σ=0. P2 — неподвижный текущий контакт.</p>'
        '<p>Остаточный слой задан подсеточно: толщина, объём, инерция и отдельная динамика/feedback не вычислялись.</p>'+
        (f'<p>Численная grad-div стабилизация γ={gamma:g} м²/с; работа отдельно от физического тепла, не изменение ν.</p>' if gamma else '')+
        ''.join(f'<video controls src="{r["video"]}"></video>' for r in reports)+'<img src="key_frames.png"><ul>'+
        ''.join(f'<li><a href="{p.name}">{p.name}</a></li>' for p in sorted(output.iterdir()) if p.suffix in (".csv",".json",".md"))+'</ul></html>')
    return evidence


def output_main(root,args,state,heartbeat):
    from .free_boundary_run import REL_RUNS,REL_OUTPUT
    run_id=args.run_id or state["current_run_id"]
    native=root/REL_RUNS/run_id/"native.h5"
    if args.command=="render":
        out=root/REL_OUTPUT/"preview"/run_id
        result=render(native,out,heartbeat,final=False)
        state.update(status="IN_PROGRESS",free_boundary_preview=str(out/"index.html"))
        print(json.dumps({"preview":str(out),"time_range_s":result["time_range_s"]}));return 0
    if args.command=="verify":
        if not args.run_id and (root/REL_OUTPUT/"verification.json").exists():
            from .free_boundary_package import package
            result=package(root,args,state,heartbeat,rerender=False)
            print(json.dumps(result,indent=2));return 0 if result["passed"] else 2
        from .free_boundary_verify import verify_native
        result=verify_native(native,heartbeat=heartbeat)
        write_json(native.parent/"native_verification.json",result)
        code=0 if result["passed"] else 2
        state["_job_exit_code"]=code
        print(json.dumps(result,indent=2));return code
    if args.command=="package":
        from .free_boundary_package import package
        result=package(root,args,state,heartbeat)
        print(json.dumps(result,indent=2));return 0 if result["passed"] else 2
    raise ValueError("Unsupported PW2 output command")
