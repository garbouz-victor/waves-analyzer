"""Re-measure historical final checkpoint geometry; absolutely no PDE advance."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from mpi4py import MPI
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.solver import CHNSSolver
from sloshing.multiphase.storage import load_checkpoint
from sloshing.multiphase.diagnostics import Diagnostics,interface_resolution
from sloshing.multiphase.provenance import source_hash


def main():
    source=Path("validation_results/step3/contact/theta60_resolved")
    output=Path("validation_results/step3a1/mesh_resolution_tests")
    run=json.loads((source/"summary.json").read_text())
    if run["status"]!="complete":
        raise ValueError("Only a complete historical execution is eligible for this read-only audit")
    solver=CHNSSolver(ModelConfig(**run["config"]))
    load_checkpoint(source/"checkpoint",solver,Diagnostics(solver))
    measured=interface_resolution(solver)
    record={"scope":"historical final checkpoint only; no timesteps, not an all-time new run",
            "source":str(source),"old_centroid_indicator":run["final_interface_resolution"],
            "new_Bernstein_indicator":measured,"analysis_source_sha256":source_hash()}
    if MPI.COMM_WORLD.rank==0:
        output.mkdir(parents=True,exist_ok=True)
        with (output/"historical_final.json").open("x") as stream:
            json.dump(record,stream,indent=2,allow_nan=False)
        print(json.dumps(record,indent=2))


if __name__=="__main__":
    main()
