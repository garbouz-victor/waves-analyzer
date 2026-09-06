"""Medium/fine norms evaluated on exact common-element quadrature."""

import json
from pathlib import Path
import time

import h5py
import numpy as np

from ..artifacts import write_tables
from ..contact import surface_regions
from .dataset import dataset_fingerprint, validate_cache
from .fem_evaluation import FEMGeometry, evaluate_map
from .overlay import common_quadrature


def compare_datasets(medium, fine, output, t_end=None):
    validate_cache(medium)
    validate_cache(fine)
    output = Path(output)
    start = time.perf_counter()
    with h5py.File(medium) as m, h5py.File(fine) as f:
        cm, cf = json.loads(m.attrs["config_json"]), json.loads(f.attrs["config_json"])
        for k in cm:
            if k not in ("mesh", "nx", "nz", "t_end") and cm[k] != cf[k]:
                raise ValueError(f"Comparison parameters differ: {k}")
        end = min(m["times"][-1], f["times"][-1]) if t_end is None else t_end
        name = "comparison_short" if end <= 1 else "comparison_full"
        sources = {"medium": dataset_fingerprint(medium), "fine": dataset_fingerprint(fine),
                   "t_end": float(end), "method": "exact-overlay-P4-v1"}
        existing = output/(name+".json")
        if existing.exists() and (output/"raw"/(name+".json")).exists():
            cached = json.loads(existing.read_text())["summary"]
            if cached.get("sources") == sources:
                print("REUSE comparison", name, flush=True)
                return cached
        fm, ff = FEMGeometry.from_hdf5(m), FEMGeometry.from_hdf5(f)
        points, weights, triangles = common_quadrature(fm, ff)
        print(f"OVERLAY {triangles} integration triangles, {len(weights)} degree-4 points", flush=True)
        mm, mf = fm.interpolation(points), ff.interpolation(points)
        om, of = fm.interpolation(points, 1), ff.interpolation(points, 1)
        masks = {"full": np.ones(len(weights), dtype=bool), "wall_left": points[0] <= -.9,
                 "wall_right": points[0] >= .9, "surface_layer": points[1] >= -1.}
        regional_weights = {key: weights*mask for key, mask in masks.items()}
        xm, xf = m["mesh/surface_x"][:], f["mesh/surface_x"][:]
        mirror = []
        for h in (m, f):
            xz = h["mesh/velocity_dof_coordinates"][:]
            order = np.lexsort((xz[1], xz[0]))
            reflected = np.lexsort((xz[1], -xz[0]))
            index = np.empty_like(order)
            index[order] = reflected
            np.testing.assert_allclose(xz[:, index], np.vstack((-xz[0], xz[1])), atol=1e-12)
            mirror.append(index)
        rows = []
        for i, t in enumerate(m["times"]):
            if t > end+1e-12:
                break
            j = int(np.searchsorted(f["times"][:], t-1e-12))
            if not np.isclose(f["times"][j], t, atol=1e-12, rtol=0):
                raise ValueError("Datasets do not have common snapshots")
            se = surface_regions(xm, m["fields/eta"][i], xf, f["fields/eta"][j])
            sn = surface_regions(xf, np.zeros_like(xf), xf, f["fields/eta"][j])
            row = {"time": float(t), **{f"eta_{key}_difference": value for key, value in se.items()},
                   **{f"eta_{key}_reference_norm": value for key, value in sn.items()}}
            dv2, vf2 = np.zeros(len(weights)), np.zeros(len(weights))
            symmetry = [0., 0.]
            for key in ("u", "w"):
                coefficients = [m[f"fields/{key}"][i], f[f"fields/{key}"][j]]
                for mesh_index, coeff in enumerate(coefficients):
                    symmetry[mesh_index] = max(symmetry[mesh_index], float(np.max(abs(coeff-(1 if key == "u" else -1)*coeff[mirror[mesh_index]]))))
                mv = evaluate_map(coefficients[0], mm)
                fv = evaluate_map(coefficients[1], mf)
                dv2 += (mv-fv)**2
                vf2 += fv**2
            row.update(velocity_difference=float(np.sqrt(weights@dv2)), velocity_reference_norm=float(np.sqrt(weights@vf2)))
            row.update(velocity_symmetry_medium=symmetry[0], velocity_symmetry_fine=symmetry[1])
            wo_m, wo_f = evaluate_map(m["fields/omega"][i], om), evaluate_map(f["fields/omega"][j], of)
            for key, rw in regional_weights.items():
                row[f"omega_{key}_difference"] = float(np.sqrt(rw@(wo_m-wo_f)**2))
                row[f"omega_{key}_reference_norm"] = float(np.sqrt(rw@wo_f**2))
            for key in ("total_energy", "divergence_l2", "divergence_relative", "divergence_bulk_l2", "divergence_bulk_relative", "max_slope"):
                row[key+"_medium"] = float(m[f"diagnostics/{key}"][i])
                row[key+"_fine"] = float(f[f"diagnostics/{key}"][j])
                row[key+"_difference"] = abs(row[key+"_medium"]-row[key+"_fine"])
            for component in ("u", "w"):
                for k, value in enumerate(abs(m[f"probes/{component}"][i]-f[f"probes/{component}"][j])):
                    row[f"probe_{component}_{k}_difference"] = float(value)
            rows.append(row)
            if i % 100 == 0:
                print(f"COMPARE t={t:.3f}, delta_v={row['velocity_difference']:.4g}, delta_omega_wall={row['omega_wall_left_difference']:.4g}", flush=True)
        times = np.array([r["time"] for r in rows])
        summary = {"medium": str(medium), "fine": str(fine), "t_end": float(end),
                   "sources": sources,
                   "snapshots": len(rows), "common_triangles": triangles, "quadrature_points": len(weights),
                   "quadrature_order": 4, "quadrature_is_exact_for_piecewise_fields": True,
                   "integrated_domain_area": float(weights.sum()), "runtime_seconds": time.perf_counter()-start}
        for key in rows[0]:
            if not key.endswith("_difference"):
                continue
            base = key[:-len("_difference")]
            curve = np.array([r[key] for r in rows])
            summary[base+"_max"] = float(curve.max())
            summary[base+"_final"] = float(curve[-1])
            summary[base+"_time_at_max"] = float(times[np.argmax(curve)])
            reference = base+"_reference_norm"
            if reference in rows[0]:
                norm = max(r[reference] for r in rows)
                summary[base+"_reference_peak_norm"] = norm
                summary[base+"_relative_max"] = float(curve.max()/max(norm, 1e-30))
                summary[base+"_relative_final"] = float(curve[-1]/max(norm, 1e-30))
        for key in ("divergence_l2", "divergence_relative", "divergence_bulk_l2", "divergence_bulk_relative"):
            for mesh in ("medium", "fine"):
                summary[f"{key}_{mesh}_max"] = max(r[f"{key}_{mesh}"] for r in rows)
            significant = np.array([r["velocity_reference_norm"] for r in rows]) > .01*max(r["velocity_reference_norm"] for r in rows)
            summary[key+"_fine_below_medium_active_fraction"] = float(np.mean(np.array([r[key+"_fine"] < r[key+"_medium"] for r in rows])[significant])) if np.any(significant) else 1.
        for mesh in ("medium", "fine"):
            summary["velocity_symmetry_"+mesh+"_max"] = max(r["velocity_symmetry_"+mesh] for r in rows)
        selected = [0, .1, .5, 1, 2, 3, 4, 5, summary["velocity_time_at_max"]]
        for h in (m, f):
            tt = h["times"][:]
            valid = tt <= end+1e-12
            for key in ("kinetic_energy", "max_slope", "omega_wall_left_L2"):
                selected.append(float(tt[valid][np.argmax(h[f"diagnostics/{key}"][:][valid])]))
        indices = np.unique([int(np.argmin(abs(times-t))) for t in selected if t <= end+1e-12])
        name = "comparison_short" if end <= 1 else "comparison_full"
        write_tables(output, name, [rows[i] for i in indices], summary)
        # Full raw curves are reproducible and intentionally ignored in Git.
        write_tables(output/"raw", name, rows, summary)
        return summary
