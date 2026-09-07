"""Measurement only: wall/CPU categories and peak RSS, never solver policy."""
from contextlib import contextmanager
import resource
import time


CATEGORIES = ("assembly", "factorization_setup", "nonlinear_SNES", "linear_solve",
              "diagnostics", "interface_resolution", "contact_geometry",
              "filesystem_write", "BE_work", "initialization")


class Performance:
    def __init__(self):
        self.started = time.perf_counter()
        self.cpu_started = time.process_time()
        self.records = {key: [] for key in CATEGORIES}

    @contextmanager
    def measure(self, category, **metadata):
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            yield
        finally:
            self.records.setdefault(category, []).append({**metadata,
                "wall_s": time.perf_counter()-wall, "cpu_s": time.process_time()-cpu})

    def snapshot(self):
        categories = {}
        for name, rows in self.records.items():
            categories[name] = ({"count": len(rows),
                "wall_s": sum(r["wall_s"] for r in rows),
                "cpu_s": sum(r["cpu_s"] for r in rows)} if rows else None)
        return {"wall_s": time.perf_counter()-self.started,
            "cpu_s": time.process_time()-self.cpu_started,
            "peak_RSS_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            "peak_RSS_scope": "Linux process lifetime high-water mark, not stage increment",
            "categories": categories, "records": self.records,
            "timing_semantics": "application category timers; nested PETSc events reported separately"}


def petsc_events():
    """Optional inclusive PETSc event data; never sum it with application timers."""
    try:
        from petsc4py import PETSc
        names = ("MatLUFactorSym", "MatLUFactorNum", "MatSolve", "KSPSolve",
                 "PCSetUp", "SNESFunctionEval", "SNESJacobianEval", "SNESSolve")
        result = {name: PETSc.Log.Event(name).getPerfInfo() for name in names}
        return {"inclusive_nested_events": result, "available": True}
    except (ImportError, AttributeError, RuntimeError) as exc:
        return {"available": False, "reason": str(exc)}
