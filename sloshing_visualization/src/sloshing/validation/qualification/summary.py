"""Build a measured decision; passing unit tests is not a qualification gate."""

import json
from pathlib import Path

import h5py

from ..artifacts import write_tables
from .comparison import compare_datasets
from .dataset import dataset_fingerprint, dataset_summary
from .plots import comparison_plot, history_plots
from .policy import qualification_decision
from .metrics import pinned_discrete_equilibrium
from .signals import modal_events


def build_summary(output, plots=True):
    output = Path(output)
    m, f = output/"runs/medium-full.h5", output/"runs/fine-full.h5"
    medium, fine = dataset_summary(m), dataset_summary(f)
    comparison = compare_datasets(m, f, output)
    particle_path = output/"particle_summary.json"
    if not particle_path.exists():
        raise RuntimeError("Particle qualification must be completed before the final decision")
    particles = json.loads(particle_path.read_text())
    sources = {"medium": dataset_fingerprint(m), "fine": dataset_fingerprint(f)}
    if particles.get("sources") != sources:
        raise RuntimeError("Particle qualification does not match current completed datasets; rerun particles phase")
    if particles.get("sampling_source") != dataset_fingerprint(output/"runs/sampling-dense.h5"):
        raise RuntimeError("Particle sampling comparison does not match the current dense-snapshot dataset")
    result = qualification_decision(medium, fine, comparison, particles)
    result["sources"] = sources
    result["discrete_static_surface_diagnostic"] = {}
    for mesh, path in (("medium", m), ("fine", f)):
        with h5py.File(path) as h:
            _, info = pinned_discrete_equilibrium(h["mesh/surface_x"][:], h["fields/eta"][0][[0, -1]], medium["parameters"]["g"])
        result["discrete_static_surface_diagnostic"][mesh] = info
    if plots:
        modal = history_plots(m, output)
        comparison_plot(output)
    else:
        modal = json.loads((output/"modal_summary.json").read_text())
    result["modal_diagnostic"] = modal
    with h5py.File(f) as h:
        c = fine["parameters"]
        fine_modal = modal_events(h["times"][:], h["diagnostics/modal_eta"][:], c["a"], c["d"], c["g"])
    result["modal_diagnostic_fine"] = fine_modal
    with (output/"modal_fine_summary.json").open("w") as handle:
        json.dump(fine_modal, handle, indent=2, allow_nan=False)
    for name, data in (("summary", result), ("medium-full_gate", medium), ("fine-full_gate", fine)):
        with (output/(name+".json")).open("w") as handle:
            json.dump(data, handle, indent=2, allow_nan=False)
    keys = ("slope_global_max", "slope_bulk_max", "slope_intermediate_max", "slope_contact_max",
            "divergence_l2_max", "divergence_relative_max", "divergence_bulk_l2_max", "divergence_bulk_relative_max",
            "volume_error_max", "contact_error_max", "energy_balance_relative_max", "max_step_energy_increase_max",
            "omega_L2_max", "omega_wall_left_L2_max", "omega_wall_right_L2_max", "omega_max_max",
            "eta_antisymmetry_max", "first_slope_above_0.1", "first_slope_above_0.3")
    write_tables(output, "dataset_metrics", [{"mesh": r["parameters"]["mesh"], **{k:r[k] for k in keys}} for r in (medium, fine)])
    print(json.dumps(result, indent=2), flush=True)
    return result
