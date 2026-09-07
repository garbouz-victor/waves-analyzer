"""Mass-inner-product Arnoldi with reorthogonalization and explicit defects.

This analysis module accepts an operator with apply/project/M0/norm methods.
Only the projected Hessenberg exponential is dense, never the full CH matrix.
"""
import numpy as np
import time
from scipy.linalg import eig, expm
from .linearized_ch import array_hash, json_hash

ALGORITHM_VERSION = "M-arnoldi-two-pass-mgs-full-Ritz-residual-v2"


class MassArnoldi:
    def __init__(self, operator, seed, max_dimension=1024):
        self.operator = operator
        self.seed = operator.project(np.asarray(seed, dtype=float))
        self.beta = operator.norm(self.seed)
        if not np.isfinite(self.beta) or self.beta == 0:
            raise ValueError("Nonzero finite mass-zero Arnoldi seed required")
        capacity = min(max_dimension, operator.size-1)
        self.V = np.zeros((operator.size, capacity+1), order="F")
        self.MV = np.zeros_like(self.V, order="F")
        self.hessenberg = np.zeros((capacity+1, capacity))
        self.V[:, 0] = self.seed/self.beta
        self.MV[:, 0] = operator.M0 @ self.V[:, 0]
        self.dimension = 0
        self.breakdown = False

    def extend(self, dimension, monitor=None):
        target = min(dimension, self.hessenberg.shape[1])
        for j in range(self.dimension, target):
            if self.breakdown:
                break
            w = self.operator.apply(self.V[:, j])
            initial_norm = self.operator.norm(w)
            for _ in range(2):
                w = self.operator.project(w)
                for i in range(j+1):
                    hij = float(self.MV[:, i] @ w)
                    self.hessenberg[i, j] += hij
                    w -= hij*self.V[:, i]
            w = self.operator.project(w)
            beta = self.operator.norm(w)
            self.hessenberg[j+1, j] = beta
            self.dimension = j+1
            if beta <= 64*np.finfo(float).eps*max(initial_norm, 1.):
                self.breakdown = True
            else:
                self.V[:, j+1] = w/beta
                self.MV[:, j+1] = self.operator.M0 @ self.V[:, j+1]
            if monitor:
                monitor(self.dimension)
        return self.dimension

    def projected(self, dimension=None):
        m = self.dimension if dimension is None else dimension
        if m > self.dimension or m < 1:
            raise ValueError("Unbuilt Arnoldi dimension")
        return self.hessenberg[:m, :m]

    def diagnostics(self, dimension=None):
        m = self.dimension if dimension is None else dimension
        V, MV = self.V[:, :m], self.MV[:, :m]
        orth = np.linalg.norm(V.T @ MV-np.eye(m), ord=2)
        # Check the actual full-operator relation on deterministic columns,
        # including the terminal residual, not merely the stored recurrence.
        columns = sorted(set((0, m//2, m-1)))
        residuals = []
        for j in columns:
            actual = self.operator.apply(V[:, j])
            expected = self.V[:, :m+1] @ self.hessenberg[:m+1, j]
            residuals.append(self.operator.norm(actual-expected)/max(self.operator.norm(actual), 1e-300))
        return {"dimension": m, "orthogonality_defect": float(orth),
                "max_checked_Arnoldi_relative_residual": float(max(residuals)),
                "checked_columns": columns, "terminal_residual_norm": float(self.hessenberg[m, m-1]),
                "maximum_basis_mass": float(np.max(abs(self.operator.mass @ V))),
                "happy_breakdown": bool(self.breakdown),
                "basis_sha256": array_hash(V), "Hessenberg_sha256": array_hash(self.projected(m))}

    def ritz(self, dimension=None, previous=None, dissipative_significance=False):
        m = self.dimension if dimension is None else dimension
        values, vectors = eig(self.projected(m))
        weights = np.linalg.solve(vectors, np.eye(m)[:, 0])
        previous_values = np.array([complex(r["real"], r["imag"]) for r in previous]) if previous else None
        power0 = self.operator.power(self.seed) if dissipative_significance else None
        rows = []
        for index in np.argsort(values.real)[::-1]:
            value, z = values[index], vectors[:, index]
            recurrence_residual = abs(self.hessenberg[m, m-1]*z[-1])/(max(abs(value), 1.)*np.linalg.norm(z))
            vector = self.V[:, :m] @ z
            # The recurrence estimate can fall far below floating-point error.
            # Scientific Ritz qualification uses the ACTUAL full-operator
            # defect, including accumulated Arnoldi and rounding errors.
            residual = self.operator.norm(self.operator.apply(vector)-value*vector)/(
                max(abs(value), 1.)*self.operator.norm(vector))
            change = (float(min(abs(previous_values-value))/max(abs(value), 1.)) if previous is not None else None)
            significance = float(abs(weights[index])*np.linalg.norm(z))
            power_fraction = None
            if dissipative_significance:
                power_fraction = float(self.beta**2*abs(weights[index])**2*self.operator.power(vector)/max(power0, 1e-300))
            significant = significance >= 1e-8 or (power_fraction is not None and power_fraction >= 1e-8)
            converged = residual <= 1e-6 and change is not None and change <= 1e-6
            rows.append({"mode": int(index), "real": float(value.real), "imag": float(value.imag),
                "tau_s": float(-1/value.real) if value.real < 0 else None,
                "relative_Ritz_residual": float(residual), "last_level_relative_change": change,
                "Arnoldi_recurrence_residual_estimate": float(recurrence_residual),
                "seed_L2_significance": significance, "initial_power_self_fraction": power_fraction,
                "significant": bool(significant), "converged": bool(converged),
                "old_z_finest": float(abs(value.real)*2.5e-6)})
        return rows

    def reference(self, times, dimension=None):
        m = self.dimension if dimension is None else dimension
        h = self.projected(m)
        reduced = np.array([self.beta*expm(float(t)*h)[:, 0] for t in times])
        states = reduced @ self.V[:, :m].T
        energies, powers, defects = [], [], []
        for q, c in zip(states, reduced):
            rate = self.V[:, :m] @ (h @ c)
            energy = self.operator.energy(q)
            power = self.operator.power(q)
            defect = float(q @ (self.operator.H @ rate)+power)
            energies.append(energy); powers.append(power)
            defects.append(abs(defect)/max(abs(power), 1e-300))
        return {"times": np.asarray(times), "states": states, "reduced": reduced,
                "energy": np.array(energies), "power": np.array(powers),
                "projected_energy_identity_relative_defect": np.array(defects)}

    def fingerprint(self, dimension=None):
        return json_hash({"algorithm_version": ALGORITHM_VERSION, "operator": self.operator.fingerprint,
                          "seed": array_hash(self.seed), **self.diagnostics(dimension)})


def compare_references(operator, coarse, fine):
    if not np.array_equal(coarse["times"], fine["times"]):
        raise ValueError("Krylov levels require identical reference times")
    initial = operator.norm(fine["states"][0])
    fields = [operator.norm(a-b)/max(operator.norm(b), 1e-14*initial)
              for a, b in zip(coarse["states"], fine["states"])]
    result = {"max_relative_L2_difference": float(max(fields)), "relative_L2_by_time": fields}
    for key in ("energy", "power"):
        relative = abs(coarse[key]-fine[key])/np.maximum(abs(fine[key]), 1e-14*abs(fine[key][0]))
        result["max_relative_"+key+"_difference"] = float(max(relative))
    return result


class RestartedExponential:
    """Independent SLEPc exponential, certified externally in the FEM M norm.

    MFN 3.24 uses Euclidean normalization internally; setting only BV's matrix
    does NOT implement an M-Arnoldi exponential (a nonidentity-M tiny regression
    catches this). Leave its internal norm Euclidean. The primary Ritz analyzer
    above remains explicit M-Arnoldi. Both compute the same mathematical exp(L).
    """
    def __init__(self, operator, dimension=64, tolerance=1e-11, max_seconds=3600.):
        from petsc4py import PETSc
        from slepc4py import SLEPc
        self.operator = operator
        self.applications = 0
        self.deadline = time.perf_counter()+max_seconds
        owner = self
        class Shell:
            def mult(self, matrix, x, y):
                if time.perf_counter() > owner.deadline:
                    raise RuntimeError("Restarted exponential resource-time guard exceeded")
                owner.applications += 1
                y.array[:] = operator.apply(x.array_r)
        self.A = PETSc.Mat().createPython((operator.size, operator.size), context=Shell(), comm=PETSc.COMM_SELF)
        self.A.setUp()
        self.mfn = SLEPc.MFN().create(PETSc.COMM_SELF)
        self.mfn.setOperator(self.A)
        self.mfn.setType(SLEPc.MFN.Type.KRYLOV)
        self.mfn.setDimensions(min(dimension, operator.size-1))
        self.mfn.setTolerances(tolerance, max_it=10000)
        self.mfn.getFN().setType(SLEPc.FN.Type.EXP)
        bv = self.mfn.getBV()
        bv.setOrthogonalization(SLEPc.BV.OrthogType.MGS, SLEPc.BV.OrthogRefineType.ALWAYS)
        self.x, self.y = self.A.createVecs()

    def action(self, seed, duration):
        if duration == 0:
            return self.operator.project(seed), {"iterations": 0, "operator_applications": 0, "reason": 1}
        before = self.applications
        self.x.array[:] = self.operator.project(seed)
        self.mfn.getFN().setScale(duration, 1.)
        self.mfn.solve(self.x, self.y)
        reason = int(self.mfn.getConvergedReason())
        if reason <= 0:
            raise RuntimeError(f"Restarted exponential did not converge: {reason}")
        return self.operator.project(self.y.array.copy()), {"iterations": self.mfn.getIterationNumber(),
            "operator_applications": self.applications-before, "reason": reason}

    def reference(self, seed, times, monitor=None):
        states, history = [], []
        q = self.operator.project(seed)
        previous = 0.
        for t in times:
            q, log = self.action(q, float(t-previous))
            states.append(q.copy()); history.append({"time": float(t), **log})
            previous = float(t)
            if monitor:
                monitor(history[-1])
        energies = [self.operator.energy(q) for q in states]
        powers = [self.operator.power(q) for q in states]
        return {"times": np.asarray(times), "states": np.asarray(states),
                "energy": np.asarray(energies), "power": np.asarray(powers)}, history

    def close(self):
        self.mfn.destroy(); self.x.destroy(); self.y.destroy(); self.A.destroy()


class RationalMassKrylov(MassArnoldi):
    """Fixed positive shifted-resolvent basis; reference is exp(t V^T M L V).

    Resolvent solves generate the space ONLY. They are not used as a BE-time
    reference. Original L, exact quadratic energies and full Ritz defects are
    retained. This avoids a single huge polynomial space at t*rho >> 1.
    """
    def __init__(self, operator, seed, shifts, max_dimension=128):
        super().__init__(operator, seed, max_dimension)
        self.shifts = tuple(float(x) for x in shifts)
        if not self.shifts or min(self.shifts) <= 0:
            raise ValueError("Positive resolvent times required")
        self.LV = np.zeros((operator.size, self.hessenberg.shape[1]), order="F")

    def extend(self, dimension, monitor=None):
        target = min(dimension, self.hessenberg.shape[1])
        for j in range(self.dimension, target):
            if self.breakdown:
                break
            self.LV[:, j] = self.operator.apply(self.V[:, j])
            shift = self.shifts[j % len(self.shifts)]
            w = self.operator.be_step(self.V[:, j], shift)
            initial_norm = self.operator.norm(w)
            for _ in range(2):
                w = self.operator.project(w)
                for i in range(j+1):
                    hij = float(self.MV[:, i] @ w)
                    self.hessenberg[i, j] += hij
                    w -= hij*self.V[:, i]
            w = self.operator.project(w)
            beta = self.operator.norm(w)
            self.hessenberg[j+1, j] = beta
            self.dimension = j+1
            if beta <= 64*np.finfo(float).eps*max(initial_norm, 1.):
                self.breakdown = True
                break
            self.V[:, j+1] = w/beta
            self.MV[:, j+1] = self.operator.M0 @ self.V[:, j+1]
            # Inspect actual sparse LU storage, not dense N^2 estimates.
            cached_bytes = sum((lu.L.nnz+lu.U.nnz)*12 for lu in self.operator._be_cache.values() if hasattr(lu, "L"))
            if cached_bytes > 2*1024**3 or (self.operator.resolvent_backend == "petsc_mumps" and len(self.operator._be_cache) > 2):
                self.operator.clear_be_cache()
            if monitor:
                monitor(self.dimension)
        return self.dimension

    def projected(self, dimension=None):
        m = self.dimension if dimension is None else dimension
        if not 1 <= m <= self.dimension:
            raise ValueError("Unbuilt rational dimension")
        return self.MV[:, :m].T @ self.LV[:, :m]

    def diagnostics(self, dimension=None):
        m = self.dimension if dimension is None else dimension
        V, MV = self.V[:, :m], self.MV[:, :m]
        return {"dimension": m, "method": "fixed-pole rational M-Krylov; projected ORIGINAL L",
                "orthogonality_defect": float(np.linalg.norm(V.T @ MV-np.eye(m), ord=2)),
                "maximum_basis_mass": float(max(abs(self.operator.mass @ V))),
                "happy_breakdown": bool(self.breakdown), "resolvent_times": list(self.shifts),
                "basis_sha256": array_hash(V), "projected_L_sha256": array_hash(self.projected(m))}

    def ritz(self, dimension=None, previous=None, dissipative_significance=True):
        m = self.dimension if dimension is None else dimension
        values, vectors = eig(self.projected(m))
        weights = np.linalg.solve(vectors, np.eye(m)[:, 0])
        previous_values = np.array([complex(r["real"], r["imag"]) for r in previous]) if previous else None
        rows = []
        for j in np.argsort(values.real)[::-1]:
            value, z = values[j], vectors[:, j]
            vector = self.V[:, :m] @ z
            defect = self.LV[:, :m] @ z-value*vector
            residual = self.operator.norm(defect)/(max(abs(value), 1.)*self.operator.norm(vector))
            change = float(min(abs(previous_values-value))/max(abs(value), 1.)) if previous is not None else None
            significance = float(abs(weights[j])*self.operator.norm(vector))
            power_fraction = float(self.beta**2*abs(weights[j])**2*self.operator.power(vector)/self.operator.power(self.seed))
            rows.append({"mode": int(j), "real": float(value.real), "imag": float(value.imag),
                "tau_s": float(-1/value.real) if value.real < 0 else None,
                "relative_Ritz_residual": residual, "last_level_relative_change": change,
                "seed_L2_significance": significance, "initial_power_self_fraction": power_fraction,
                "significant": bool(significance >= 1e-8 or power_fraction >= 1e-8),
                "converged": bool(residual <= 1e-6 and change is not None and change <= 1e-6),
                "old_z_finest": float(abs(value.real)*2.5e-6)})
        return rows
