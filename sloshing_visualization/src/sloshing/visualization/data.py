"""Read-only FEM -> bounded display cache. No solver calls, moving mesh or FD curl."""

import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np

from ..validation.qualification.dataset import dataset_fingerprint, validate_cache
from ..validation.qualification.fem_evaluation import FEMGeometry, evaluate_map
from ..validation.qualification.signals import modal_events, horizontal_quadrature
from ..validation.contact import surface_value
from ..validation.particles import particle_seeds, FEMVelocityHistory, advect
from .scope import load_scope, compatible_cache, event_frames, LIMITATIONS
from .tracers import linear_lagrangian_displacement, choose_magnification
from .wall_profiles import profile_points,wall_normal_profiles,trace_speed_max,arrow_geometry

SCHEMA = "step2-display-v2-no-slip"
HISTORY_FIELDS = {key:key for key in ("kinetic_energy", "total_energy", "modal_eta", "slope_bulk", "slope_global")}
HISTORY_FIELDS["potential_energy"] = "gravitational_potential_energy"


def surface_extrema(x, eta):
    """Exact P2 range, including stationary points inside each surface edge."""
    a = 2*eta[::2][:-1]-4*eta[1::2]+2*eta[2::2]
    b = -3*eta[::2][:-1]+4*eta[1::2]-eta[2::2]
    nonzero = a != 0
    s = np.zeros_like(a)
    s[nonzero] = -b[nonzero]/(2*a[nonzero])
    hit = nonzero & (s > 0) & (s < 1)
    values = np.r_[eta, (a*s*s+b*s+eta[::2][:-1])[hit]]
    return float(values.min()), float(values.max())


def grid_map(geometry, x, z, degree=2):
    xx, zz = np.meshgrid(x, z)
    return geometry.interpolation(np.array([xx.ravel(), zz.ravel()]), degree)


def prepare_data(dataset, medium, output, scope=None):
    """Exclusive cache creation; COMPLETE matching caches reused, never overwritten."""
    started = time.perf_counter()
    dataset, medium, output = Path(dataset), Path(medium), Path(output)
    scope = scope or load_scope()
    for path in (dataset, medium):
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}. Restore explicitly with qualify_animation_dataset.py; renderer never runs PDE.")
    params = validate_cache(dataset)
    mc = validate_cache(medium)
    if params["mesh"] != "fine" or mc["mesh"] != "medium":
        raise ValueError("STEP 2 expects fine rendering source and medium comparison source")
    for k in ("a", "d", "g", "alpha_deg", "nu", "dt", "snapshot_dt", "t_end", "integrator"):
        if params[k] != mc[k]:
            raise ValueError("Physical source mismatch: "+k)
    if (params["a"], params["d"], params["t_end"]) != (1., 10., 5.):
        raise ValueError("This STEP 2 scope is defined for a=1, d=10, t_end=5")
    sources = {"fine": dataset_fingerprint(dataset), "medium": dataset_fingerprint(medium)}
    expected = {"schema": SCHEMA, "sources": sources, "scope": scope}
    output.mkdir(parents=True, exist_ok=True)
    cache = output/"display_data.h5"
    if cache.exists():
        with h5py.File(cache, "r") as h:
            if h.attrs.get("status") != "complete":
                raise RuntimeError("Incomplete display cache; use a new output directory")
            compatible_cache(json.loads(h.attrs["identity"]), expected)
            metadata = json.loads(h.attrs["metadata"])
            for key in ("u", "w", "omega", "eta", "linear_paths", "advected_paths", "wall_left_speed_max", "wall_right_speed_max", "profile_left_speed", "profile_right_speed"):
                if h[key].shape[0] != metadata["frame_count"]:
                    raise RuntimeError("Truncated display cache: "+key)
        (output/"animation_metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
        print("REUSE display cache", cache, flush=True)
        return cache, metadata
    with h5py.File(dataset, "r") as source, h5py.File(cache, "x") as out:
        out.attrs.update(status="running", identity=json.dumps(expected, sort_keys=True))
        try:
            geo = FEMGeometry.from_hdf5(source)
            times = source["times"][:]
            x = np.linspace(-1, 1, scope["display_nx"])
            z = np.linspace(scope["near_surface_z_min_m"], 0, scope["display_nz"])
            xf = np.linspace(-1, 1, scope["full_depth_nx"])
            zf = np.linspace(-10, 0, scope["full_depth_nz"])
            # Resolve every surface edge for drawing; never hide contact features
            # behind a uniform 161-point surface curve.
            sx = source["mesh/surface_x"][:]
            ex = np.unique(np.concatenate([np.linspace(l, r, 9) for l, r in zip(sx[::2][:-1], sx[::2][1:])]))
            seeds, groups = particle_seeds()
            seeds = seeds[groups == "bulk"]  # 28, including seven deep tracers
            mappings = {"near": grid_map(geo, x, z), "omega": grid_map(geo, x, z, 1),
                        "full": grid_map(geo, xf, zf), "full_omega": grid_map(geo, xf, zf, 1),
                        "seeds": geo.interpolation(seeds.T)}
            qx, qz = np.meshgrid(np.linspace(-scope['arrow_seed_abs_x_max_m'], scope['arrow_seed_abs_x_max_m'], 25), np.linspace(-1.45, -.15, 18))
            mappings["arrows"] = geo.interpolation(np.array([qx.ravel(), qz.ravel()]))
            eps=np.unique(np.r_[np.linspace(0,scope['profile_epsilon_max_m'],101),1e-5,1e-4,5e-4,1e-3,2e-3,5e-3,1e-2,2e-2,5e-2,1e-1])
            profile_depths=np.array(scope['boundary_layer_depths_m'])
            for side in ('left','right'):
                points,profile_shape=profile_points(geo,side,profile_depths,eps)
                mappings['profile_'+side]=geo.interpolation(points)
            wall_z=np.sort(np.r_[geo.z,(geo.z[:-1]+geo.z[1:])/2])
            wallmaps={side:geo.interpolation(np.array([np.full_like(wall_z,edge),wall_z])) for side,edge in (('left',-1.),('right',1.))}
            zoomx,zoomz=np.linspace(-1,-.8,161),np.linspace(-1.2,0,121)
            mappings['wall_zoom']=grid_map(geo,zoomx,zoomz)
            zx,zz=np.meshgrid(np.linspace(-.94,-.82,5),np.linspace(-1.15,-.16,16))
            mappings['zoom_arrows']=geo.interpolation(np.array([zx.ravel(),zz.ravel()]))
            # Exact line quadrature; depth uses a fixed log range later.
            depths = np.unique(np.r_[np.linspace(-10, -1.5, 25), np.linspace(-1.5, 0, 41)])
            depth_maps = []
            for depth in depths:
                points, weights = horizontal_quadrature(geo, depth)
                depth_maps.append((geo.interpolation(points), weights))
            wallpoints = np.array([(side, depth) for side in (-1., 1.) for depth in np.linspace(-10, 0, 101)]).T
            wallmap = geo.interpolation(wallpoints)
            for key, values in (("times", times), ("x", x), ("z", z), ("full_x", xf), ("full_z", zf),
                                ("eta_x", ex), ("seeds", seeds), ("arrow_x", qx), ("arrow_z", qz), ("depths", depths),
                                ('profile_epsilon',eps),('profile_depths',profile_depths),('wall_zoom_x',zoomx),('wall_zoom_z',zoomz),('zoom_arrow_x',zx),('zoom_arrow_z',zz)):
                out.create_dataset(key, data=values)
            shapes = {"u": (len(z), len(x)), "w": (len(z), len(x)), "omega": (len(z), len(x)),
                      "full_u": (len(zf), len(xf)), "full_w": (len(zf), len(xf)), "full_omega": (len(zf), len(xf)),
                      "eta": (len(ex),), "arrow_u": qx.shape, "arrow_w": qx.shape, "seed_v": seeds.shape,
                      "depth_q": (len(depths),), "ranges": (7,), 'wall_zoom_u':(len(zoomz),len(zoomx)), 'wall_zoom_w':(len(zoomz),len(zoomx)),
                      'zoom_arrow_u':zx.shape,'zoom_arrow_w':zx.shape}
            for side in ('left','right'):
                out.create_dataset('wall_'+side+'_speed_max',shape=(len(times),),dtype='f8')
                for key in ('u','w','speed'):shapes['profile_'+side+'_'+key]=profile_shape
            for key, shape in shapes.items():
                out.create_dataset(key, shape=(len(times),)+shape, dtype="f8" if key in ("eta", "seed_v", "ranges") or key.startswith('profile_') else "f4",
                                   chunks=(1,)+shape, compression="gzip", compression_opts=1, shuffle=True, fletcher32=True)
            for key, source_key in HISTORY_FIELDS.items():
                out.create_dataset(key, data=source["diagnostics/"+source_key][:])
            speed_max, wall_max = 0., 0.
            percentile_samples = []  # <= 20 MB; not all source coefficients in RAM
            band = np.abs(x) <= scope["transition_abs_x_max_m"]
            bulk = np.abs(x) <= scope["bulk_abs_x_max_m"]
            for i, t in enumerate(times):
                u, w, om, eta = (source["fields/"+key][i] for key in ("u", "w", "omega", "eta"))
                ud, wd = (evaluate_map(c, mappings["near"]).reshape(shapes["u"]) for c in (u, w))
                od = evaluate_map(om, mappings["omega"]).reshape(shapes["omega"])
                out["u"][i], out["w"][i], out["omega"][i] = ud, wd, od
                for key, coeff, mapping in (("full_u",u,"full"),("full_w",w,"full"),("full_omega",om,"full_omega"),
                                            ("arrow_u",u,"arrows"),("arrow_w",w,"arrows")):
                    out[key][i] = evaluate_map(coeff, mappings[mapping]).reshape(shapes[key])
                out["eta"][i] = surface_value(sx, eta, ex)
                out["seed_v"][i] = np.column_stack([evaluate_map(c, mappings["seeds"]) for c in (u, w)])
                out["depth_q"][i] = [np.sqrt(np.dot(weights, evaluate_map(u, mp)**2+evaluate_map(w, mp)**2)) for mp, weights in depth_maps]
                speed = np.hypot(ud, wd)
                speed_max = max(speed_max, float(speed[:, band].max()))
                wall_max = max(wall_max, float(np.hypot(evaluate_map(u, wallmap), evaluate_map(w, wallmap)).max()))
                for side in ('left','right'):
                    out['wall_'+side+'_speed_max'][i]=trace_speed_max(evaluate_map(u,wallmaps[side]),evaluate_map(w,wallmaps[side]))
                    for key,value in wall_normal_profiles(geo,u,w,side,profile_depths,eps,mappings['profile_'+side]).items():out['profile_'+side+'_'+key][i]=value
                for key,coeff,mp in (('wall_zoom_u',u,'wall_zoom'),('wall_zoom_w',w,'wall_zoom'),('zoom_arrow_u',u,'zoom_arrows'),('zoom_arrow_w',w,'zoom_arrows')):
                    out[key][i]=evaluate_map(coeff,mappings[mp]).reshape(shapes[key])
                percentile_samples.append(np.asarray(abs(od[::2, band][:, ::2]), dtype=np.float32).ravel())
                lo, hi = surface_extrema(sx, eta)
                out["ranges"][i] = [lo, hi, speed.min(), speed.max(), speed[:, bulk].max(), abs(om).max(), abs(od[:, band]).max()]
                if i % 100 == 0:
                    out.flush()
                    print(f"DISPLAY {i+1}/{len(times)} t={t:.3f}; elapsed={time.perf_counter()-started:.1f}s", flush=True)
            if wall_max > params["wall_tolerance"]:
                raise RuntimeError("Rendered FEM interpolation violates exact-wall tolerance")
            xi = linear_lagrangian_displacement(times, out["seed_v"][:])
            linear = seeds[None]+xi
            # Actual-domain validity is separate from display magnification.
            valid = np.ones(linear.shape[:2], dtype=bool)
            for i in range(len(times)):
                p = linear[i]
                valid[i] = (abs(p[:, 0]) <= 1) & (p[:, 1] >= -10) & (p[:, 1] <= 0)
                ids = np.flatnonzero(valid[i])
                valid[i, ids] &= p[ids, 1] <= surface_value(sx, source["fields/eta"][i], p[ids, 0])
            valid = np.logical_and.accumulate(valid, axis=0)
            if not valid.all():
                raise RuntimeError("Linear tracer becomes invalid; inspect seeds before rendering")
            out.create_dataset("linear_paths", data=linear)
            out.create_dataset("linear_valid", data=valid)
            history = FEMVelocityHistory.open(dataset)
            try:
                paths = advect(history, seeds)
            finally:
                history.close()
            if np.isfinite(paths["first_invalid_time"]).any():
                raise RuntimeError("Advected diagnostic tracer invalid; inspect before rendering")
            out.create_dataset("advected_paths", data=paths["paths"])
            mag, p95 = choose_magnification(xi, scope["particle_target_p95_display_m"])
            # Same factor for BOTH tracer models and BOTH coordinates. Ensure
            # amplified glyphs remain inside reference domain; no clipping.
            combined = np.concatenate((xi, paths["paths"]-seeds[None]), axis=0)
            for axis, low, high in ((0, -1., 1.), (1, -10., 0.)):
                delta = combined[:, :, axis]
                allowed = np.where(delta > 0, (high-seeds[:, axis])/np.maximum(delta, 1e-300),
                                   (low-seeds[:, axis])/np.minimum(delta, -1e-300))
                mag = min(mag, .85*float(allowed.min()))
            mag = float(np.floor(mag)) if mag >= 1 else float(mag)
            modal = modal_events(times, out["modal_eta"][:], params["a"], params["d"], params["g"])
            omega_limit = float(np.percentile(np.concatenate(percentile_samples), scope["omega_percentile"]))
            eta_range = [float(out["ranges"][:, 0].min()), float(out["ranges"][:, 1].max())]
            repo = Path(__file__).resolve().parents[4]
            sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
            metadata = {"schema": SCHEMA, "source_commit_sha": sha, "sources": sources, "parameters": params,
                        "scope": scope, "interpretation": scope["interpretation"], "limitations": LIMITATIONS,
                        "frame_count": len(times), "fps": scope["fps"], "physical_time_range_s": [float(times[0]), float(times[-1])],
                        "playback_slowdown": 1/(scope["fps"]*params["snapshot_dt"]), "video_duration_s": len(times)/scope["fps"],
                        "velocity_color_range_m_per_s": [0., speed_max],
                        "velocity_scale_method": "Maximum on fixed 161x121 near-surface display grid, |x|<=0.98, ALL snapshots; a display statistic, not a norm",
                        "vorticity_color_range_per_s": [-omega_limit, omega_limit],
                        "vorticity_scale_method": "99.5 percentile |DG1 omega| on every second near-surface display row/eligible column, |x|<=0.98, ALL snapshots; corners excluded only from scale",
                        "color_saturation": "Out-of-range values remain unchanged; fixed colorbar extend indicators denote saturation, not clipping of fields",
                        "particle_displacement_magnification": mag, "particle_p95_displacement_m": p95,
                        "particle_displacement_max_m": float(np.linalg.norm(xi, axis=-1).max()),
                        "tracer_model_separation_max_m": float(np.linalg.norm(linear-paths["paths"], axis=-1).max()),
                        "surface_display_magnification": 1., "surface_axis_units": "mm; conversion only, no amplitude magnification",
                        "true_eta_range_m": eta_range, "surface_ylim_mm": [-1.12*max(abs(np.array(eta_range)))*1000, 1.12*max(abs(np.array(eta_range)))*1000],
                        "fixed_arrow_seconds": .12/speed_max, "arrow_key_speed_m_per_s": .0005,
                        "wall_interpolation_speed_max_m_per_s": wall_max, "modal": modal, "key_frames": event_frames(times, modal),
                        "depth_q_max": float(out["depth_q"][:].max()), "depth_log_display_min": 1e-10,
                        "cache_quantization": "float32 display velocity/omega/depth; float64 eta, tracers and scientific histories",
                        "dependency_versions": {n: importlib.metadata.version(n) for n in ("numpy", "scipy", "matplotlib", "h5py", "plotly")},
                        "precompute_runtime_seconds": time.perf_counter()-started}
            arrow=arrow_geometry(qx,qz,out['arrow_u'][:],out['arrow_w'][:],metadata['fixed_arrow_seconds'],scope)
            profile_speed_max=max(float(out['profile_'+s+'_speed'][:].max()) for s in ('left','right'))
            component_max={k:max(float(abs(out['profile_'+s+'_'+k][:]).max()) for s in ('left','right')) for k in ('u','w')}
            wall_residual={s:float(out['wall_'+s+'_speed_max'][:].max()) for s in ('left','right')}
            if max(wall_residual.values())>params['wall_tolerance']:raise RuntimeError('Exact P2 wall-trace maximum exceeds tolerance')
            symmetry=float(abs(out['profile_left_speed'][:]-out['profile_right_speed'][:]).max())
            metadata.update(arrow)
            metadata['no_slip_visualization']={'wall_velocity_source':'actual P2 FEM; exact polynomial trace norm maximum',
                'wall_speed_max':wall_residual,'arrow_seed_x_range':[float(qx.min()),float(qx.max())],
                'arrow_wall_margin_m':arrow['wall_visual_margin'],'quiver_pivot':'tail',
                'boundary_layer_profile_depths_m':profile_depths.tolist(),'profile_epsilon_range_m':[0.,float(eps[-1])],
                'profile_speed_ylim_m_per_s':[0.,1.06*profile_speed_max],
                'profile_component_ylim_m_per_s':{k:[-1.06*v,1.06*v] for k,v in component_max.items()},
                'profile_left_right_speed_max_difference_m_per_s':symmetry,'profile_relative_symmetry_error':symmetry/max(profile_speed_max,1e-300),
                'zoom_fixed_arrow_seconds':20.,'corner_patch_z_range_m':[scope['contact_corner_z_min_m'],scope['contact_corner_z_max_m']]}
            out.attrs.update(metadata=json.dumps(metadata, allow_nan=False), status="complete")
        except BaseException as exc:
            out.attrs.update(status="failed", failure=repr(exc))
            raise
    (output/"animation_metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n")
    return cache, metadata
