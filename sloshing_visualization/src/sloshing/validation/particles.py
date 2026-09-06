"""Deterministic trajectory qualification from saved FEM velocities, no renderer."""

from collections import OrderedDict

import h5py
import numpy as np

from .contact import surface_value
from .qualification.dataset import validate_cache
from .qualification.fem_evaluation import FEMGeometry, evaluate_map


class FEMVelocityHistory:
    def __init__(self, geometry, times, u, w, eta, surface_x):
        self.geometry, self.times = geometry, np.asarray(times)
        if len(self.times) < 2 or np.any(np.diff(self.times) <= 0):
            raise ValueError("Velocity time samples must be strictly increasing")
        self.u, self.w, self.eta, self.surface_x = u, w, eta, surface_x
        self.cache = OrderedDict()
        self.file = None

    @classmethod
    def open(cls, path):
        validate_cache(path)
        h = h5py.File(path, "r")
        try:
            result = cls(FEMGeometry.from_hdf5(h), h["times"][:], h["fields/u"], h["fields/w"],
                         h["fields/eta"], h["mesh/surface_x"][:])
            result.file = h
            return result
        except BaseException:
            h.close()
            raise

    def close(self):
        if self.file is not None:
            self.file.close()

    def bracket(self, t):
        if t < self.times[0]-1e-12 or t > self.times[-1]+1e-12:
            raise ValueError("Particle time outside saved solution")
        i = min(max(np.searchsorted(self.times, t, side="right")-1, 0), len(self.times)-2)
        return i, (t-self.times[i])/(self.times[i+1]-self.times[i])

    def snapshot(self, i):
        if i not in self.cache:
            self.cache[i] = (self.u[i], self.w[i], self.eta[i])
            if len(self.cache) > 3:
                self.cache.popitem(last=False)
        self.cache.move_to_end(i)
        return self.cache[i]

    def velocity(self, t, points):
        i, theta = self.bracket(t)
        mapping = self.geometry.interpolation(points)
        a, b = self.snapshot(i), self.snapshot(i+1)
        return np.stack([(1-theta)*evaluate_map(a[c], mapping)+theta*evaluate_map(b[c], mapping) for c in (0, 1)])

    def surface(self, t, x):
        i, theta = self.bracket(t)
        eta = (1-theta)*self.snapshot(i)[2]+theta*self.snapshot(i+1)[2]
        return surface_value(self.surface_x, eta, x)

    def valid(self, t, points):
        x, z = points
        g = self.geometry
        inside = np.isfinite(x) & np.isfinite(z) & (x >= g.x[0]) & (x <= g.x[-1]) & (z >= g.z[0]) & (z <= g.z[-1])
        if np.any(inside):
            ids = np.flatnonzero(inside)
            inside[ids] &= z[ids] <= self.surface(t, x[ids])
        return inside


def particle_seeds():
    bulk = [(x, z) for z in (-.1, -.3, -1., -3.) for x in (-.75, -.5, -.25, 0., .25, .5, .75)]
    wall = [(x, z) for z in (-.05, -.2, -1.) for x in (-.95, -.9, .9, .95)]
    return np.array(bulk+wall), np.array(["bulk"]*len(bulk)+["near_wall"]*len(wall))


def advect(history, seeds, times=None, max_step=.0025):
    """RK4 with linear snapshot interpolation; no reflection or silent clipping.

    Invalid particles stop permanently. First-invalid time is the first failing
    RK stage/end check, with time resolution at most max_step (not a root solve).
    Steps also stop at every velocity snapshot kink.
    """
    times = history.times if times is None else np.asarray(times)
    if times[0] != history.times[0] or np.any(np.diff(times) <= 0) or max_step <= 0:
        raise ValueError("Invalid particle integration time grid")
    points = np.asarray(seeds, dtype=float).T.copy()
    active = history.valid(times[0], points)
    invalid = np.full(points.shape[1], np.nan)
    invalid[~active] = times[0]
    points[:, ~active] = np.nan
    paths = [points.T.copy()]
    current = float(times[0])
    for target in times[1:]:
        while current < target-1e-13:
            next_saved = history.times[min(np.searchsorted(history.times, current+1e-12), len(history.times)-1)]
            end = min(float(target), current+max_step, float(next_saved))
            step = end-current
            if step <= 0:
                raise RuntimeError("Particle step failed to advance")
            ids = np.flatnonzero(active)
            if len(ids):
                y = points[:, ids]
                def stage(t, trial):
                    good = history.valid(t, trial) & active[ids]
                    bad = ids[~good & active[ids]]
                    invalid[bad] = t
                    active[bad] = False
                    value = np.full_like(trial, np.nan)
                    if np.any(good):
                        value[:, good] = history.velocity(t, trial[:, good])
                    return value
                k1 = stage(current, y)
                k2 = stage(current+step/2, y+step*k1/2)
                k3 = stage(current+step/2, y+step*k2/2)
                k4 = stage(end, y+step*k3)
                final = y+step*(k1+2*k2+2*k3+k4)/6
                good = history.valid(end, final) & active[ids]
                bad = ids[~good & active[ids]]
                invalid[bad] = end
                active[bad] = False
                points[:, ids[good]] = final[:, good]
                points[:, ~active] = np.nan
            current = end
        paths.append(points.T.copy())
    return {"times": np.asarray(times), "paths": np.asarray(paths), "first_invalid_time": invalid}


def compare_paths(a, b, seeds, groups):
    np.testing.assert_allclose(a["times"], b["times"], atol=1e-12, rtol=0)
    distance = np.linalg.norm(a["paths"]-b["paths"], axis=2)
    rows = []
    for i, (seed, group) in enumerate(zip(seeds, groups)):
        finite = np.isfinite(distance[:, i])
        value = distance[finite, i]
        row = {"particle": i, "group": str(group), "x0": float(seed[0]), "z0": float(seed[1]),
               "separation_max_m": float(value.max()) if len(value) else None,
               "separation_final_m": float(distance[-1, i]) if np.isfinite(distance[-1, i]) else None,
               "time_at_max": float(a["times"][finite][np.argmax(value)]) if len(value) else None,
               "invalid_a": float(a["first_invalid_time"][i]) if np.isfinite(a["first_invalid_time"][i]) else None,
               "invalid_b": float(b["first_invalid_time"][i]) if np.isfinite(b["first_invalid_time"][i]) else None}
        rows.append(row)
    return rows
