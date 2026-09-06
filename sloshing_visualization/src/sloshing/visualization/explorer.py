"""Small, offline Plotly explorer. Display samples are never scientific norms."""

import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np
from plotly.offline import get_plotlyjs

from .panels import display_cell_edges
from .scope import phase_label


def display_numbers(values):
    """Six significant digits for HTML only; input/source arrays stay unchanged."""
    values = np.asarray(values)
    return np.array([float(f"{v:.6g}") for v in values.ravel()]).reshape(values.shape).tolist()


def explorer_payload(h, meta):
    indices = np.unique(np.r_[np.arange(0, len(h["times"]), 10), len(h["times"])-1])
    data = {"meta": meta, "indices": indices.tolist(), "frames": [], "grids": {}}
    for name, prefix in (("near", ""), ("full", "full_")):
        x, z = h[prefix+"x"][:][::2], h[prefix+"z"][:][::2]
        data["grids"][name] = {"x": display_cell_edges(x).tolist(), "z": display_cell_edges(z).tolist()}
    for key in ("times", "eta_x", "seeds", "kinetic_energy", "potential_energy", "total_energy"):
        data[key] = h[key][:].tolist()
    data["initial_eta"] = h["eta"][0].tolist()
    data["arrow_x"] = h["arrow_x"][:][::2, ::2].ravel().tolist()
    data["arrow_z"] = h["arrow_z"][:][::2, ::2].ravel().tolist()
    for i in indices:
        frame = {"t": float(h["times"][i]), "phase": phase_label(float(h["times"][i]), meta["modal"]),
                 "eta": h["eta"][i].tolist(), "paths": h["linear_paths"][i].tolist()}
        for name, prefix in (("near", ""), ("full", "full_")):
            frame[name] = {key: display_numbers(h[prefix+key][i][::2, ::2]) for key in ("u", "w", "omega")}
        for key in ("arrow_u", "arrow_w"):
            frame[key] = display_numbers(h[key][i][::2, ::2].ravel())
        data["frames"].append(frame)
    return data


def explorer(cache, meta, output):
    started = time.perf_counter()
    with h5py.File(cache, "r") as h:
        payload = explorer_payload(h, meta)
    template = Path(__file__).with_name("explorer_template.html").read_text()
    html = template.replace("__PLOTLY_JS__", get_plotlyjs()).replace(
        "__STEP2_DATA__", json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/"))
    path = Path(output)/"explorer.html"
    path.write_text(html)
    summary = {"status": "complete", "sources": meta["sources"], "frames": len(payload["frames"]),
               "source_snapshot_stride": 10, "physical_sampling_s": .05,
               "near_grid": [81, 61], "full_grid": [41, 81], "bytes": path.stat().st_size,
               "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
               "offline": True, "field_precision": "6 significant digits, display only",
               "normalization": "Same fixed physical limits and arrow/displacement factors as MP4",
               "scope": "Exploration only, not scientific error norms; no new PDE states",
               "runtime_seconds": time.perf_counter()-started}
    (Path(output)/"explorer_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(f"EXPLORER {summary['frames']} snapshots, {summary['bytes']/1e6:.1f} MB, {summary['runtime_seconds']:.1f}s", flush=True)
