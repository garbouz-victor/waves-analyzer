"""Read-only historical P2 re-audit: no PDE advance, no historical writes."""
import hashlib
import argparse
import json
from pathlib import Path
import numpy as np
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.diagnostics import Diagnostics, interface_resolution, resolution_cell_data
from sloshing.multiphase.storage import load_checkpoint
from sloshing.multiphase.resolution import polynomial_bounds, normal_widths
from sloshing.multiphase.free_energy import transition_width
from sloshing.multiphase.provenance import run_provenance


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="validation_results/step3a2/p2_resolution/historical_reaudit.json")
    args=parser.parse_args()
    source=Path("validation_results/step3/contact/theta60_resolved")
    output_file=Path(args.output)
    output=output_file.parent
    cfg=ModelConfig(**json.loads((source/"summary.json").read_text())["config"])
    solver=CHNSSolver(cfg)
    if solver.comm.size!=1:
        raise ValueError("Historical rank-local checkpoint requires one MPI rank")
    checkpoint=source/"checkpoint/rank_0000.h5"
    before=hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    load_checkpoint(source/"checkpoint",solver,Diagnostics(solver))
    new=interface_resolution(solver)
    tri,v,g,ranges,indices,cert,directional=resolution_cell_data(solver)
    lo,hi=polynomial_bounds(v,cfg.phase_degree)
    bernstein=(lo<=.9)&(hi>=-.9)
    actual=(ranges["min"]<=.9)&(ranges["max"]>=-.9)
    false=np.flatnonzero(bernstein & ~actual)
    old_width=normal_widths(tri[bernstein],g[bernstein])
    previous=json.loads(Path("validation_results/step3a1/mesh_resolution_tests/historical_final.json").read_text())
    after=hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    assert before==after
    result={"scope":"read-only historical t=0.5 checkpoint; no PDE advance",
        "source":str(source),"checkpoint_sha256_before":before,"checkpoint_sha256_after":after,
        "provenance":run_provenance(solver),"old_centroid_indicator":previous["old_centroid_indicator"],
        "step3a1_bernstein_sampled":previous["new_Bernstein_indicator"],
        "remeasured_bernstein_active_count":int(bernstein.sum()),
        "remeasured_bernstein_sampled_min":float(transition_width(cfg.epsilon)/old_width.max()),
        "bernstein_false_positive_count":len(false),"bernstein_false_positive_cell_ids":false.tolist(),
        "exact_p2_certified_normal":new,
        "false_positive_witnesses":[{"local_cell":int(i),"bernstein_min":float(lo[i]),
            "bernstein_max":float(hi[i]),"actual_min":float(ranges["min"][i]),
            "actual_max":float(ranges["max"][i]),"range_status":str(ranges["status"][i])} for i in false]}
    output.mkdir(parents=True,exist_ok=True)
    with output_file.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
    np.savez_compressed(output/(output_file.stem+"_cell_certificates.npz"),triangles=tri,
        nodal_values=v,gradients=g,active_cell_ids=indices,h_cert=cert["width"],h_directional=directional,
        width_fallback_reason=cert["fallback_reason"],**ranges)
    print(json.dumps({"Bernstein_false_positives":len(false),"new":new},indent=2),flush=True)


if __name__=="__main__":
    main()
