"""Check reflection invariance of a candidate triangulation before its solve."""
import json
from pathlib import Path
import numpy as np
from mpi4py import MPI
from sloshing.multiphase.config import ModelConfig
from sloshing.multiphase.mesh import make_mesh

c=ModelConfig.from_json("configs/step3a2/equilibrium60.json").changed(nx=24,nz=12,mesh_diagonal="crossed")
mesh,_,_,cost=make_mesh(c)
tri=mesh.geometry.x[mesh.geometry.dofmap,:2]
def key(points):
    return tuple(sorted(tuple(p) for p in np.rint(points*1e12).astype(np.int64)))
keys={key(t) for t in tri}
reflected=tri.copy();reflected[:,:,0]*=-1
matched=sum(key(t) in keys for t in reflected)
result={"config":c.as_dict(),"cost":cost,"reflection_matched_cells":matched,"total_cells":len(tri),
        "reflection_invariant":matched==len(tri)}
output=Path("validation_results/step3a2/conditioning/crossed_mesh_symmetry.json")
with output.open("x") as stream:
    json.dump(result,stream,indent=2)
print(json.dumps(result,indent=2),flush=True)
