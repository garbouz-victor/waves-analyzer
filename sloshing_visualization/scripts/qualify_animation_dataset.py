#!/usr/bin/env python3
"""Sequential measured gates; expensive simulations are never pytest defaults."""

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from sloshing.config import SimulationConfig
from sloshing.validation.qualification.dataset import dataset_summary, run_dataset, validate_cache
from sloshing.validation.qualification.comparison import compare_datasets
from sloshing.validation.qualification.policy import fine_short_gate


def phase_run(phase, output):
    mesh, duration = phase.split("-")
    filename = "animation_candidate_fine.json" if mesh == "fine" else "animation_candidate.json"
    base = SimulationConfig.from_json(Path(__file__).resolve().parents[1]/"configs"/filename)
    config = replace(base, mesh=mesh, t_end=1. if duration == "short" else 5.)
    path = output/"runs"/(phase+".h5")
    order = ("medium-short", "medium-full", "fine-short", "fine-full")
    index = order.index(phase)
    if index:
        previous = output/"runs"/(order[index-1]+".h5")
        validate_cache(previous)
    if phase == "fine-full":
        comparison = compare_datasets(output/"runs/medium-full.h5", output/"runs/fine-short.h5", output, t_end=1.)
        gate = fine_short_gate(comparison)
        if not gate["numerical_gate_passed"]:
            raise RuntimeError("Fine-short gate requires analysis before fine-full: "+"; ".join(gate["problems"]))
    prefix = output/"runs"/(mesh+"-short.h5") if duration == "full" else None
    run_dataset(config, path, prefix)
    summary = dataset_summary(path)
    output.mkdir(parents=True, exist_ok=True)
    with (output/(phase+"_gate.json")).open("w") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps(summary, indent=2), flush=True)
    if phase == "fine-short":
        comparison = compare_datasets(output/"runs/medium-full.h5", path, output, t_end=1.)
        gate = fine_short_gate(comparison)
        with (output/"fine_short_numerical_gate.json").open("w") as handle:
            json.dump(gate, handle, indent=2, allow_nan=False)
        from sloshing.validation.qualification.plots import comparison_plot
        comparison_plot(output, "comparison_short")
        if not gate["numerical_gate_passed"]:
            raise RuntimeError("Fine-short comparison needs analysis: "+"; ".join(gate["problems"]))


def compare_full_and_gate(output):
    comparison = compare_datasets(output/"runs/medium-full.h5", output/"runs/fine-full.h5", output)
    gate = fine_short_gate(comparison)
    with (output/"full_numerical_gate.json").open("w") as handle:
        json.dump(gate, handle, indent=2, allow_nan=False)
    if not gate["numerical_gate_passed"]:
        raise RuntimeError("Long-window field comparison requires analysis before particles: "+"; ".join(gate["problems"]))
    return comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("medium-short", "medium-full", "fine-short", "fine-full", "compare", "particles", "summary", "all"), required=True)
    parser.add_argument("--output", type=Path, default=Path("validation_results/animation_qualification"))
    args = parser.parse_args()
    if args.phase == "all":
        for phase in ("medium-short", "medium-full", "fine-short", "fine-full"):
            phase_run(phase, args.output)
    elif args.phase in ("medium-short", "medium-full", "fine-short", "fine-full"):
        phase_run(args.phase, args.output)
        return
    if args.phase in ("compare", "all"):
        compare_full_and_gate(args.output)
    if args.phase in ("particles", "all"):
        if args.phase == "particles":
            compare_full_and_gate(args.output)
        from sloshing.validation.qualification.particle_study import particle_study
        particle_study(args.output)
    if args.phase in ("summary", "all"):
        from sloshing.validation.qualification.summary import build_summary
        build_summary(args.output)


if __name__ == "__main__":
    main()
