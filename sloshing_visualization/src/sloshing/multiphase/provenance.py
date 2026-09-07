"""Run-time provenance, not a retrospective claim about historical results."""
import hashlib
import os
from pathlib import Path


def source_hash():
    base=Path(__file__).resolve().parent
    digest=hashlib.sha256()
    for path in sorted(base.rglob("*.py")):
        digest.update(str(path.relative_to(base)).encode());digest.update(path.read_bytes())
    return digest.hexdigest()


def run_provenance(solver):
    mesh=solver.mesh
    digest=hashlib.sha256(mesh.geometry.x.tobytes()+mesh.geometry.dofmap.tobytes()).hexdigest()
    rank_hashes=solver.comm.allgather(digest)
    return {"git_commit":os.environ.get("STEP3_GIT_COMMIT","unavailable"),
            "multiphase_source_sha256":source_hash(),
            "config_fingerprint":solver.config.fingerprint(),
            "docker_image_digest":os.environ.get("STEP3_DOCKER_IMAGE_DIGEST","unavailable"),
            "rank_count":solver.comm.size,"dt_schedule":solver.schedule.descriptor,
            "mesh_fingerprint":hashlib.sha256("".join(rank_hashes).encode()).hexdigest(),
            "mesh_fingerprint_scope":"rank-ordered geometry/connectivity; rank/partition sensitive"}
