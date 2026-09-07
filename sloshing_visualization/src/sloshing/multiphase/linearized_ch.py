"""Analysis-only sparse linearized CH; never replaces production CHNS.

L = -mobility M0^{-1} K M0^{-1} H. H is the exact discrete energy Hessian,
including the wall second variation. All dense arrays below are vectors or
small diagnostic data; the production operator/inverse is never densified.
"""
import hashlib
import json
from contextlib import nullcontext
import numpy as np
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import LinearOperator, splu

ALGORITHM_VERSION = "mass-zero-sparse-linearized-ch-v1"


def array_hash(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


class LinearizedCH:
    def __init__(self, M0, K, H, mobility=1., mass_vector=None, provenance=None):
        self.M0, self.K, self.H = (csr_matrix(a, dtype=float) for a in (M0, K, H))
        self.size = self.M0.shape[0]
        if any(a.shape != (self.size, self.size) for a in (self.M0, self.K, self.H)):
            raise ValueError("Incompatible CH matrix shapes")
        if not all(np.isfinite(a.data).all() for a in (self.M0, self.K, self.H)):
            raise ValueError("Nonfinite CH matrix")
        if not np.isfinite(mobility) or mobility <= 0:
            raise ValueError("Positive mobility required")
        self.mobility = float(mobility)
        self.one = np.ones(self.size)
        self.mass = self.M0 @ self.one if mass_vector is None else np.asarray(mass_vector).copy()
        self.area = float(self.mass @ self.one)
        if self.area <= 0 or self.mass.shape != (self.size,):
            raise ValueError("Invalid discrete mass vector")
        self.mass_lu = splu(self.M0.tocsc())
        self._be_cache = {}
        self.resolvent_backend = "scipy_superlu"
        self.performance = None  # instrumentation has no effect on operator identity
        identity = {"algorithm_version": ALGORITHM_VERSION, "mobility": self.mobility,
                    "mass_sha256": array_hash(self.mass), "provenance": provenance or {}}
        for name in ("M0", "K", "H"):
            a = getattr(self, name)
            a.sort_indices()
            identity[name] = {"shape": list(a.shape), "indptr": array_hash(a.indptr),
                              "indices": array_hash(a.indices), "values": array_hash(a.data)}
        self.identity = identity
        self.fingerprint = json_hash(identity)
        self.operator = LinearOperator((self.size, self.size), matvec=self.apply, dtype=float)

    def mass_solve(self, rhs):
        if np.iscomplexobj(rhs):
            return self.mass_lu.solve(rhs.real) + 1j*self.mass_lu.solve(rhs.imag)
        return self.mass_lu.solve(np.asarray(rhs))

    def project(self, q):
        q = np.asarray(q)
        if q.ndim == 1:
            return q-self.one*(self.mass @ q)/self.area
        return q-self.one[:, None]*(self.mass @ q)[None, :]/self.area

    def inner(self, q, r):
        return np.vdot(q, self.M0 @ r)

    def norm(self, q):
        return float(np.sqrt(max(0., self.inner(q, q).real)))

    def dual_norm(self, f):
        return float(np.sqrt(max(0., np.vdot(f, self.mass_solve(f)).real)))

    def raw(self, q):
        return -self.mobility*self.mass_solve(self.K @ self.mass_solve(self.H @ q))

    def apply(self, q):
        """Projection is ONLY an analysis operation, not a solution correction."""
        return self.project(self.raw(self.project(q)))

    def chemical(self, q):
        return self.mass_solve(self.H @ q)

    def energy(self, q):
        return float(.5*np.vdot(q, self.H @ q).real)

    def power(self, q):
        r = self.chemical(q)
        return float(self.mobility*np.vdot(r, self.K @ r).real)

    def be_step(self, q, dt):
        """Coupled sparse linear BE, no dense inverse or endpoint heat surrogate."""
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Positive finite dt required")
        if dt not in self._be_cache:
            # Algebraic auxiliary scaling r_tilde=sqrt(dt)*r. Eliminating it
            # gives EXACTLY M(q_new-q_old)+dt*mobility*K*M^-1*H*q_new=0.
            # Balanced block magnitudes avoid pathological SuperLU pivot fill
            # for the extremely small resolvent times in the CH spectrum.
            with self.performance.measure("assembly", dt=float(dt)) if self.performance else nullcontext():
                root_dt = np.sqrt(dt)
                a = bmat([[self.M0, root_dt*self.mobility*self.K], [-root_dt*self.H, self.M0]], format="csc")
            # Only a few distinct block steps are normally used. Call clear to
            # avoid retaining LU factors from rejected planner candidates.
            with self.performance.measure("factorization_setup", dt=float(dt)) if self.performance else nullcontext():
                self._be_cache[dt] = (PETScSparseFactor(a) if self.resolvent_backend == "petsc_mumps" else splu(a))
        with self.performance.measure("linear_solve", dt=float(dt)) if self.performance else nullcontext():
            rhs = np.concatenate((self.M0 @ q, np.zeros(self.size)))
            result = self._be_cache[dt].solve(rhs)
        return self.project(result[:self.size])

    def clear_be_cache(self):
        for factor in self._be_cache.values():
            if isinstance(factor, PETScSparseFactor):
                factor.close()
        self._be_cache.clear()

    def matrix_checks(self, seed=3301, count=20):
        rng = np.random.default_rng(seed)
        conservation, projection, bilinear, identity = [], [], {k: [] for k in ("M0", "K", "H")}, []
        for _ in range(count):
            q, r = self.project(rng.normal(size=(2, self.size)).T).T
            rate = self.raw(q)
            scale = max(self.norm(rate), 1e-300)
            conservation.append(abs(self.mass @ rate)/(np.sqrt(self.area)*scale))
            projection.append(self.norm(self.apply(q)-rate)/scale)
            for name in bilinear:
                a = getattr(self, name)
                denominator = np.linalg.norm(q)*np.linalg.norm(a @ r)+np.linalg.norm(r)*np.linalg.norm(a @ q)
                bilinear[name].append(abs(q @ (a @ r)-r @ (a @ q))/max(denominator, 1e-300))
            power = self.power(q)
            identity.append(abs(q @ (self.H @ rate)+power)/max(abs(power), 1e-300))
        skew = {name: float(np.linalg.norm((getattr(self, name)-getattr(self, name).T).data)/
                           max(np.linalg.norm(getattr(self, name).data), 1e-300)) for name in bilinear}
        return {"random_seed": seed, "vectors": count,
                "max_relative_mass_defect": float(max(conservation)),
                "max_projection_relative_defect": float(max(projection)),
                "frobenius_symmetry_defects": skew,
                "bilinear_symmetry_defects": {k: float(max(v)) for k, v in bilinear.items()},
                "max_relative_energy_identity_defect": float(max(identity))}


def assemble_linearized_ch(solver, phi_eq, provenance=None):
    """Use UFL's second derivative of the existing scalar energy, degree 12."""
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from petsc4py import PETSc
    from .equilibrium import energy_form, _measures
    if solver.comm.size != 1:
        raise ValueError("STEP 3A.3 sparse analyzer requires one MPI rank")
    if solver.config.rho_liquid != solver.config.rho_gas or solver.config.g or solver.config.a_x:
        raise ValueError("CH spectral planner requires matched density and no body force")
    dx, _ = _measures(solver)
    test, trial = ufl.TestFunction(phi_eq.function_space), ufl.TrialFunction(phi_eq.function_space)
    first = ufl.derivative(energy_form(solver, phi_eq), phi_eq, test)
    forms = (trial*test*dx, ufl.inner(ufl.grad(trial), ufl.grad(test))*dx,
             ufl.derivative(first, phi_eq, trial))
    matrices = []
    for form in forms:
        a = fp.assemble_matrix(fem.form(form)); a.assemble()
        indptr, indices, values = a.getValuesCSR()
        matrices.append(csr_matrix((values.copy(), indices.copy(), indptr.copy()), shape=a.getSize()))
        a.destroy()
    v = fp.assemble_vector(fem.form(test*dx))
    v.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    mass = v.array.copy(); v.destroy()
    return LinearizedCH(*matrices, mobility=solver.config.mobility, mass_vector=mass, provenance=provenance)


class NonlinearCHRate:
    """Independent chemical weak assembly and rate; NO time advance."""
    def __init__(self, solver, phase_space, linear):
        import ufl
        from dolfinx import fem
        from .equilibrium import chemical_stationary_form
        self.phi = fem.Function(phase_space)
        self.form = fem.form(chemical_stationary_form(solver, self.phi, 0., ufl.TestFunction(phase_space)))
        self.linear = linear

    def residual(self, coefficients):
        from dolfinx.fem import petsc as fp
        from petsc4py import PETSc
        self.phi.x.array[:] = coefficients; self.phi.x.scatter_forward()
        v = fp.assemble_vector(self.form)
        v.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        result = v.array.copy(); v.destroy()
        return result

    def rate(self, coefficients):
        mu = self.linear.mass_solve(self.residual(coefficients))
        return -self.linear.mobility*self.linear.mass_solve(self.linear.K @ mu)

    def finite_difference_checks(self, phi_eq, psi, amplitudes=(1e-3, 5e-4, 2.5e-4, 1.25e-4)):
        linear = self.linear
        hq = linear.H @ psi
        lq = linear.raw(psi)
        rows = []
        for delta in amplitudes:
            rp = self.residual(phi_eq+delta*psi)
            rm = self.residual(phi_eq-delta*psi)
            plus = -linear.mobility*linear.mass_solve(linear.K @ linear.mass_solve(rp))
            minus = -linear.mobility*linear.mass_solve(linear.K @ linear.mass_solve(rm))
            central = (plus-minus)/(2*delta)
            rows.append({"delta": delta,
                "hessian_relative_M_dual_error": linear.dual_norm((rp-rm)/(2*delta)-hq)/linear.dual_norm(hq),
                "central_rate_relative_L2_error": linear.norm(central-lq)/linear.norm(lq),
                "one_sided_rate_relative_L2_error": linear.norm(plus/delta-lq)/linear.norm(lq),
                "nonlinear_mass_rate": float(linear.mass @ plus),
                "relative_nonlinear_mass_defect": float(abs(linear.mass @ plus)/(np.sqrt(linear.area)*linear.norm(plus)))})
        return rows


class PETScSparseFactor:
    """Serial sparse MUMPS factor, using the already pinned PETSc stack."""
    def __init__(self, matrix):
        from petsc4py import PETSc
        a = matrix.tocsr()
        self.A = PETSc.Mat().createAIJ(size=a.shape, csr=(a.indptr, a.indices, a.data), comm=PETSc.COMM_SELF)
        self.A.assemble()
        self.ksp = PETSc.KSP().create(PETSc.COMM_SELF)
        self.ksp.setType("preonly")
        self.ksp.getPC().setType("lu")
        self.ksp.getPC().setFactorSolverType("mumps")
        self.ksp.setOperators(self.A)
        self.ksp.setErrorIfNotConverged(True)
        self.ksp.setUp()
        self.x, self.b = self.A.createVecs()

    def solve(self, rhs):
        self.b.array[:] = rhs
        self.ksp.solve(self.b, self.x)
        if self.ksp.getConvergedReason() <= 0:
            raise RuntimeError("Sparse resolvent factor solve failed")
        return self.x.array.copy()

    def close(self):
        self.x.destroy(); self.b.destroy(); self.ksp.destroy(); self.A.destroy()


def shifted_modes(operator, target_rate, count=8):
    """Sparse shift-invert analysis; full-L residuals certify returned modes.

    For sigma<0 and a=sqrt(-1/sigma), the two-field matrix is
    [[M,-a*mobility*K],[-a*H,M]], with RHS [-M*b/sigma,0].
    Eliminating a*r gives (L-sigma I)q=b. No dense production L.
    """
    from scipy.sparse.linalg import eigs
    sigma = float(target_rate)
    if sigma >= 0 or operator.size <= count+2:
        raise ValueError("Negative target and sufficient phase dimension required")
    a = np.sqrt(-1/sigma)
    matrix = bmat([[operator.M0, -a*operator.mobility*operator.K],
                   [-a*operator.H, operator.M0]], format="csc")
    factor = PETScSparseFactor(matrix)
    def inverse(b):
        rhs = np.concatenate((-operator.M0 @ operator.project(b)/sigma, np.zeros(operator.size)))
        return operator.project(factor.solve(rhs)[:operator.size])
    shell = LinearOperator((operator.size, operator.size), matvec=inverse, dtype=float)
    rng = np.random.default_rng(3304)
    seed = operator.project(rng.normal(size=operator.size))
    try:
        inverse_values, vectors = eigs(shell, k=count, which="LM", v0=seed,
            ncv=min(operator.size, max(32, 2*count+1)), tol=1e-12, maxiter=1000)
    finally:
        factor.close()
    values = sigma+1/inverse_values
    rows, modes = [], []
    for j, value in enumerate(values):
        v = operator.project(vectors[:, j])
        v /= operator.norm(v)
        residual = operator.norm(operator.raw(v)-value*v)/max(abs(value), 1.)
        rows.append({"real": float(value.real), "imag": float(value.imag),
            "relative_full_operator_residual": residual,
            "relative_mass": float(abs(operator.mass @ v)/np.sqrt(operator.area)),
            "tau": float(-1/value.real) if value.real < 0 else None})
        modes.append(v)
    return rows, np.column_stack(modes)
