"""Small validation helpers: controls, retained arrays and candidate publication.

No physical config mutation, residual scaling, field correction or timestep policy.
"""
from dataclasses import dataclass, asdict
import copy
import json
from pathlib import Path
import numpy as np
from .linearized_ch import array_hash, json_hash


@dataclass(frozen=True)
class NonlinearControls:
    atol: float = 1e-10
    rtol: float = 1e-9
    stol: float = 0.
    max_it: int = 30

    def __post_init__(self):
        if not np.isfinite([self.atol, self.rtol, self.stol]).all() or min(self.atol, self.rtol) <= 0:
            raise ValueError("Finite positive nonlinear tolerances required")
        if self.stol != 0 or self.max_it != 30:
            raise ValueError("Validation keeps stol=0, max_it=30")

    def record(self):
        return {**asdict(self), "snes_type": "newtonls", "line_search_type": "bt",
            "ksp_type": "preonly", "pc_type": "lu", "factor_solver_type": "mumps"}


def actual_controls(snes):
    rtol, atol, stol, max_it = snes.getTolerances()
    ksp = snes.getKSP(); pc = ksp.getPC()
    return {"atol": atol, "rtol": rtol, "stol": stol, "max_it": max_it,
        "snes_type": snes.getType(), "line_search_type": snes.getLineSearch().getType(),
        "ksp_type": ksp.getType(), "pc_type": pc.getType(), "factor_solver_type": pc.getFactorSolverType()}


def apply_controls(snes, controls, prefix):
    """Fresh context, explicit setters; no persistent PETSc tolerance options."""
    if not prefix or not prefix.endswith("_"):
        raise ValueError("Unique explicit PETSc prefix ending in '_' required")
    snes.setOptionsPrefix(prefix)
    snes.setTolerances(**asdict(controls))
    if actual_controls(snes) != controls.record():
        raise ValueError("Unexpected solver settings; refusing silent control changes")


def finite_arrays(*arrays):
    return all(np.isfinite(np.asarray(a)).all() for a in arrays)


def synchronized_solution(problem, state):
    """DOLFINx 0.10: x and Function do NOT have to share storage."""
    from dolfinx.fem.petsc import assign
    if problem.solver.getConvergedReason() <= 0:
        raise RuntimeError("Only converged PETSc iterates may be synchronized as solutions")
    assign(problem.x, state)
    state.x.scatter_forward()
    if not np.array_equal(problem.x.array, state.x.array):
        raise RuntimeError("Serial PETSc/Function solution synchronization mismatch")


def check_candidate_material(solver, candidate):
    """Same quadrature and SAME existing checker, on a private observer."""
    import basix
    import ufl
    from dolfinx import fem
    if not finite_arrays(candidate.x.array):
        raise RuntimeError("Nonfinite physical candidate; accepted state unchanged")
    observer = copy.copy(solver)
    observer.state = candidate
    points, _ = basix.make_quadrature(solver.mesh.basix_cell(), solver.config.quadrature_degree)
    observer.phase_quadrature = fem.Expression(ufl.split(candidate)[2], points)
    observer._check_material()


def publish_candidate(solver, candidate, dt, end_time, phase):
    """No admissibility decisions here; caller checks before this transaction."""
    names = ("state", "old", "older")
    previous = {key: getattr(solver, key).x.array.copy() for key in names}
    attrs = ("time", "step_number", "current_dt", "current_phase")
    old_attrs = {key: getattr(solver, key) for key in attrs}
    try:
        solver.older.x.array[:] = previous["state"]
        solver.state.x.array[:] = candidate.x.array
        solver.old.x.array[:] = candidate.x.array
        for key in names:
            getattr(solver, key).x.scatter_forward()
        solver.time, solver.step_number = end_time, old_attrs["step_number"]+1
        solver.current_dt, solver.current_phase = dt, phase
    except Exception:
        for key in names:
            getattr(solver, key).x.array[:] = previous[key]
            getattr(solver, key).x.scatter_forward()
        for key, value in old_attrs.items():
            setattr(solver, key, value)
        raise


def rate_arrays(phi_old, phase_rate, mu, dt, **extra):
    """Independent retained arrays; recovered rate is explicitly diagnostic."""
    if not np.isfinite(dt) or dt <= 0 or not finite_arrays(phi_old, phase_rate, mu):
        raise ValueError("Invalid retained rate snapshot")
    old, rate, chemical = (np.array(a, copy=True) for a in (phi_old, phase_rate, mu))
    new = old+dt*rate
    arrays = {"phi_old": old, "phase_rate": rate, "mu_new": chemical,
        "phi_new": new, "increment_from_rate": dt*rate,
        "increment_from_physical_states": new-old,
        "phase_rate_recovered_from_rounded_phi": (new-old)/dt, "dt": np.array(dt), **extra}
    if not finite_arrays(*arrays.values()):
        raise ValueError("Nonfinite snapshot")
    return {key: np.array(value, copy=True) for key, value in arrays.items()}


def write_snapshot(folder, arrays, metadata):
    folder = Path(folder)
    if (folder/"state.npz").exists() or (folder/"snapshot.json").exists():
        raise ValueError("Refusing to overwrite a snapshot")
    arrays = {key: np.array(value, copy=True) for key, value in arrays.items()}
    record = {**metadata, "snapshot_type": "BE-validation-unknowns", "snapshot_version": 1,
        "array_sha256": {key: array_hash(value) for key, value in arrays.items()}}
    record["snapshot_sha256"] = json_hash(record)
    np.savez_compressed(folder/"state.npz", **arrays)
    (folder/"snapshot.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")
    return record


def read_snapshot(folder):
    folder = Path(folder)
    record = json.loads((folder/"snapshot.json").read_text())
    if record["snapshot_sha256"] != json_hash({k: v for k, v in record.items() if k != "snapshot_sha256"}):
        raise ValueError("Snapshot metadata fingerprint mismatch")
    with np.load(folder/"state.npz") as data:
        arrays = {key: data[key].copy() for key in data.files}
    if {k: array_hash(v) for k, v in arrays.items()} != record["array_sha256"]:
        raise ValueError("Snapshot array fingerprint mismatch")
    return arrays, record
