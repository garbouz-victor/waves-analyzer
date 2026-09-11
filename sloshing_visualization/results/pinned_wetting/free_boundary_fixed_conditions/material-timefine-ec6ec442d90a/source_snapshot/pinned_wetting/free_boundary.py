"""Full material NS on an evolving isoparametric P2 mesh, PW2 only.

No fixed-domain surface extrapolation and no Eulerian convection omission:
ALL geometry nodes move at the P2 material velocity. See MODEL_NOTES.md.
"""
from dataclasses import dataclass, replace, asdict
import numpy as np
from scipy.sparse import bmat, csr_matrix, coo_matrix
from scipy.sparse.linalg import spsolve, splu, gmres, LinearOperator
from skfem import (MeshTri2, Basis, ElementVector, ElementTriP2, ElementTriP1,
                   BilinearForm, LinearForm, asm)
from skfem.helpers import dot, sym_grad, ddot, div
from ..mesh import generate_mesh
from ..config import SimulationConfig

MODEL_ID = "PW2_FULL_NS_LEFT_NAVIER_RIGHT_NOSLIP_MATERIAL_P2_V1"


@dataclass(frozen=True)
class Controls:
    nx: int = 12
    nz: int = 24
    dt: float = .0025
    t_end: float = .12
    x_grading: float = 2.
    z_grading: float = 4.
    local_right_levels: int = 0
    intorder: int = 8
    nonlinear_tolerance: float = 2e-11
    nonlinear_iterations: int = 15
    minimum_dt: float = 1e-6
    max_halvings: int = 8
    rezone_interval_s: float = 0.
    hydrostatic_split: bool = False
    skew_divergence: bool = False
    rezone_strategy: str = "harmonic_p2"
    grad_div_gamma_m2_s: float = 0.

    def __post_init__(self):
        numbers=(self.dt,self.t_end,self.minimum_dt,self.rezone_interval_s,self.nonlinear_tolerance,self.grad_div_gamma_m2_s)
        if not all(np.isfinite(x) for x in numbers) or min(self.dt,self.t_end,self.minimum_dt,self.nonlinear_tolerance)<=0:
            raise ValueError("PW2 requires finite positive time/solver controls")
        if self.t_end>5 or self.dt>.01 or self.rezone_interval_s<0 or min(self.nx,self.nz)<2:
            raise ValueError("PW2 controls exceed declared horizon/cadence or have invalid mesh")
        if self.rezone_strategy not in ("harmonic_p2","straight_interior","quality_optimized"):
            raise ValueError("Unknown PW2 rezoning strategy")
        if self.grad_div_gamma_m2_s < 0:
            raise ValueError("PW2 grad-div coefficient must be nonnegative, in m2/s")


@BilinearForm
def mass(u, v, w): return dot(u, v)

@BilinearForm
def strain(u, v, w): return 2*.01*ddot(sym_grad(u), sym_grad(v))

@BilinearForm
def divergence(u, q, w): return q*div(u)

@BilinearForm
def left_friction(u, v, w): return .01/.50*u[1]*v[1]

@LinearForm
def gravity(v, w): return -9.81*v[1]


@LinearForm
def surface_gravity(v, w): return -9.81*w.x[1]*dot(v,w.n)


@BilinearForm
def divergence_correction(u,v,w): return .5*w.div_vm*dot(u,v)


def _scatter(local, rows, columns, shape):
    """Assemble element matrices without vector-component zero kernels."""
    r = np.broadcast_to(rows.T[:, :, None], local.shape).ravel()
    c = np.broadcast_to(columns.T[:, None, :], local.shape).ravel()
    result = coo_matrix((local.ravel(), (r, c)), shape=shape).tocsr()
    result.eliminate_zeros()
    return result


def volume_operators(scalar, pressure, components, div_vm=None):
    """Same curved-cell quadrature/forms, expanded into scalar P2 blocks.

    This changes assembly cost, not the trial spaces, quadrature, or PDE.
    Independent tests compare every full matrix to separate skfem forms.
    """
    phi = np.array([np.asarray(b[0]) for b in scalar.basis])
    grad = np.array([b[0].grad for b in scalar.basis])
    psi = np.array([np.asarray(b[0]) for b in pressure.basis])
    dx = scalar.dx
    def product(left, right, weight=dx):
        return np.einsum("ieq,jeq,eq->eij", left, right, weight, optimize=True)
    mass_local = product(phi, phi)
    xx, zz = product(grad[:, 0], grad[:, 0]), product(grad[:, 1], grad[:, 1])
    zx = product(grad[:, 1], grad[:, 0])
    nodes = scalar.element_dofs
    dofs = [component[nodes] for component in components]
    shape = (2*scalar.N, 2*scalar.N)
    M = sum((_scatter(mass_local, d, d, shape) for d in dofs), csr_matrix(shape))
    K = (_scatter(.01*(2*xx+zz), dofs[0], dofs[0], shape)
         + _scatter(.01*(xx+2*zz), dofs[1], dofs[1], shape)
         + _scatter(.01*zx, dofs[0], dofs[1], shape)
         + _scatter(.01*zx.transpose(0, 2, 1), dofs[1], dofs[0], shape))
    B = sum((_scatter(product(psi, grad[:, axis]), pressure.element_dofs, dofs[axis],
                      (pressure.N, 2*scalar.N)) for axis in (0, 1)),
            csr_matrix((pressure.N, 2*scalar.N)))
    force = np.zeros(2*scalar.N)
    np.add.at(force, dofs[1].T.ravel(), (-9.81*np.einsum("ieq,eq->ei", phi, dx)).ravel())
    Kg = csr_matrix(shape)
    if div_vm is not None:
        local = product(phi, phi, .5*div_vm*dx)
        Kg = sum((_scatter(local, d, d, shape) for d in dofs), csr_matrix(shape))
    return M, K, B, force, Kg


def grad_div_operator(scalar, components, gamma):
    """Consistent numerical gamma*(div u, div v), NOT physical viscosity.

    This remains a separate full-space operator and separately measured work.
    It does not replace the pressure-space incompressibility constraint.
    """
    shape = (2*scalar.N, 2*scalar.N)
    if gamma == 0.:
        return csr_matrix(shape)
    grad = np.array([b[0].grad for b in scalar.basis])
    dofs = [component[scalar.element_dofs] for component in components]
    result = csr_matrix(shape)
    for test_axis in (0, 1):
        for trial_axis in (0, 1):
            local = gamma*np.einsum("ieq,jeq,eq->eij", grad[:,test_axis],
                grad[:,trial_axis], scalar.dx, optimize=True)
            result += _scatter(local,dofs[test_axis],dofs[trial_axis],shape)
    return result


class PicardLinearSolver:
    """Fresh per material step; lagged LU is only a GMRES preconditioner.

    Every solve uses the CURRENT matrix and RHS. Failed or insufficiently
    accurate Krylov solves refactor the CURRENT matrix with direct SuperLU.
    No tolerance of the outer iteration or accepted PDE checks is relaxed.
    """
    def __init__(self):
        self.factor = None
        self.statistics = {"factorizations": 0, "krylov_iterations": 0,
                           "direct_fallbacks": 0, "maximum_true_block_residual": 0.}

    def solve(self, matrix, rhs, velocity_size):
        def residual(solution):
            if not np.all(np.isfinite(solution)):
                return np.inf
            r = matrix@solution-rhs
            return max(float(np.max(abs(r[:velocity_size]))),
                       float(np.max(abs(r[velocity_size:]))))
        solution = None
        if self.factor is not None:
            def counted(_):
                self.statistics["krylov_iterations"] += 1
            try:
                solution, info = gmres(matrix, rhs,
                    M=LinearOperator(matrix.shape, matvec=self.factor.solve, dtype=matrix.dtype),
                    rtol=2e-13, atol=1e-15, restart=30, maxiter=4,
                    callback=counted, callback_type="pr_norm")
                if info != 0 or residual(solution) > 1e-13:
                    solution = None
            except (RuntimeError, ValueError):
                solution = None
            if solution is None:
                self.statistics["direct_fallbacks"] += 1
        if solution is None:
            self.factor = splu(matrix, permc_spec="COLAMD")
            self.statistics["factorizations"] += 1
            solution = self.factor.solve(rhs)
        error = residual(solution)
        if not np.isfinite(error) or error > 1e-13:
            raise ValueError(f"current mixed linear system residual={error}")
        self.statistics["maximum_true_block_residual"] = max(
            self.statistics["maximum_true_block_residual"], error)
        return solution


def boundary_operators(mesh,components,intorder):
    """Parametric-edge quadrature: no physical→reference Newton inversion.

    n ds=(-z',x') ds_parameter on the complete LEFT→RIGHT free curve.
    Thus this uses the actual curved surface normal, not z=0 normals.
    """
    t,weights=np.polynomial.legendre.leggauss(max(4,(intorder+2)//2))
    t=(t+1)/2;weights=weights/2
    phi=np.array([(1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)])
    derivative=np.array([4*t-3,4-8*t,4*t-1])
    scalar=Basis(mesh,ElementTriP2())
    facets=mesh.boundaries["left"]
    nodes=np.vstack((mesh.facets[0,facets],scalar.dofs.facet_dofs[0,facets],mesh.facets[1,facets])).T
    tangent=np.einsum("ien,eq->inq",mesh.p[:,nodes.T],derivative)
    measure=np.sqrt(np.sum(tangent**2,axis=0))*weights
    local=.02*np.einsum("iq,jq,eq->eij",phi,phi,measure)
    dofs=components[1][nodes]
    row=np.broadcast_to(dofs[:,:,None],local.shape).ravel()
    col=np.broadcast_to(dofs[:,None,:],local.shape).ravel()
    n=2*scalar.N
    wall=coo_matrix((local.ravel(),(row,col)),shape=(n,n)).tocsr()
    ids=interface_nodes(mesh)
    surface_nodes=np.array([ids[i:i+3] for i in range(0,len(ids)-2,2)])
    coordinates=np.einsum("ien,eq->inq",mesh.p[:,surface_nodes.T],phi)
    tangent=np.einsum("ien,eq->inq",mesh.p[:,surface_nodes.T],derivative)
    normal_measure=np.array([-tangent[1],tangent[0]])
    local_force=-9.81*np.einsum("eq,ieq,jq,q->iej",coordinates[1],normal_measure,phi,weights)
    force=np.zeros(n)
    for axis in (0,1):np.add.at(force,components[axis][surface_nodes].ravel(),local_force[axis].ravel())
    return wall,force


def initial_mesh(c, alpha_deg=2.):
    cfg = SimulationConfig(nx=c.nx, nz=c.nz, x_grading=c.x_grading,
                           z_grading=c.z_grading)
    m = generate_mesh(cfg)
    # Local INITIAL spatial refinement is a numerical choice, not a profile.
    for level in range(c.local_right_levels):
        centers = m.p[:, m.t].mean(axis=1)
        marked = np.flatnonzero((centers[0] > 1-.30/(2**level)) &
                               (centers[1] > -.40/(2**level)))
        m = m.refined(marked)
    m = m.with_boundaries({"left": lambda x: x[0] == -1.,
        "right": lambda x: x[0] == 1., "bottom": lambda x: x[1] == -10.,
        "surface": lambda x: x[1] == 0.})
    p = m.p.copy()
    p[1] += (p[1]+10.)/10.*p[0]*np.tan(np.deg2rad(alpha_deg))
    # Straight initial triangles: edge nodes are exact midpoints after tilt.
    return replace(MeshTri2.from_mesh(replace(m, doflocs=p)), _boundaries=m.boundaries)


def surface_has_resolved_crossing(points, subdivisions=64):
    """Early rejection of resolved full-curve crossings, not a certificate.

    Sample every P2 edge parametrically; never sort by x.  The independent
    verifier additionally checks controlled curve geometry and exact contacts.
    Tangencies or sub-tessellation features are not certified by this guard.
    """
    t=np.arange(subdivisions+1,dtype=float)/subdivisions
    phi=np.array([(1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)])
    curves=np.array([(points[:,edge:edge+3]@phi).T for edge in range(0,points.shape[1]-2,2)])
    lower,upper=curves.min(axis=1),curves.max(axis=1)
    def cross(one,two):return one[...,0]*two[...,1]-one[...,1]*two[...,0]
    for first in range(len(curves)):
        for second in range(first+1,len(curves)):
            if np.any(upper[first]<lower[second]) or np.any(upper[second]<lower[first]):continue
            a=curves[first,:-1,None,:];b=curves[second,None,:-1,:]
            r=np.diff(curves[first],axis=0)[:,None,:];s=np.diff(curves[second],axis=0)[None,:,:]
            determinant=cross(r,s)
            with np.errstate(divide="ignore",invalid="ignore"):
                one=cross(b-a,s)/determinant;two=cross(b-a,r)/determinant
            hit=(determinant!=0)&(one>=0)&(one<=1)&(two>=0)&(two<=1)
            if second==first+1:hit[-1,0]=False  # sole shared topological node
            if np.any(hit):return True
    return False


def quadratic_minima(mesh):
    """Exact minimum of signed P2 Jacobian determinant on each triangle."""
    pts = np.array([[0, 1, 0, .5, 0, .5], [0, 0, 1, 0, .5, .5]])
    mapping = mesh.mapping()
    # detDF raises a generic exception at exact zeros, before we can classify
    # the invalid candidate and invoke the finite rejection policy.
    val = (mapping.J(0,0,pts)*mapping.J(1,1,pts)-
           mapping.J(0,1,pts)*mapping.J(1,0,pts))
    # scikit-fem triangles can have either orientation; initial sign is fixed.
    sign = np.sign(val[:, 0])
    val = val*sign[:, None]
    f, f10, f01, fh0, f0h, fhh = val.T
    aa = 2*(f10+f-2*fh0); cc = 2*(f01+f-2*f0h)
    dd = f10-f-aa; ee = f01-f-cc
    bb = 4*(fhh-f-.25*(aa+cc)-.5*(dd+ee))
    mins = np.minimum.reduce((f, f10, f01))
    def take(x, y, valid):
        nonlocal mins
        v = aa*x*x+bb*x*y+cc*y*y+dd*x+ee*y+f
        mins = np.minimum(mins, np.where(valid, v, np.inf))
    with np.errstate(divide="ignore", invalid="ignore"):
        x = -dd/(2*aa); take(x, np.zeros_like(x), (x>0)&(x<1))
        y = -ee/(2*cc); take(np.zeros_like(y), y, (y>0)&(y<1))
        q = aa-bb+cc; r = bb-2*cc+dd-ee
        x = -r/(2*q); take(x, 1-x, (x>0)&(x<1))
        determinant = 4*aa*cc-bb*bb
        x = (bb*ee-2*cc*dd)/determinant
        y = (bb*dd-2*aa*ee)/determinant
        take(x, y, (x>0)&(y>0)&(x+y<1))
    return mins, sign


class Geometry:
    """Quadrature-only geometry work; no redundant momentum assembly."""
    def __init__(self,mesh,c):
        self.mesh, self.controls = mesh, c
        self.scalar = Basis(mesh, ElementTriP2(), intorder=c.intorder)
        self.vbasis = Basis(mesh, ElementVector(ElementTriP2()), intorder=c.intorder)
        self.components = self.vbasis.split_indices()

    def nodal_vector(self,v):
        return np.vstack([v[indices] for indices in self.components])

    def measure(self,v):
        coords=np.asarray(self.scalar.global_coordinates())
        value=np.asarray(self.vbasis.interpolate(v))
        kinetic=.5*float(np.sum(np.sum(value**2,axis=0)*self.scalar.dx))
        potential=9.81*(float(np.sum(coords[1]*self.scalar.dx))+100.)
        return {"volume_m2":float(np.sum(self.scalar.dx)),"kinetic":kinetic,
                "potential":potential,"energy":kinetic+potential}


class Operators(Geometry):
    def __init__(self, mesh, c, advecting_velocity=None):
        super().__init__(mesh,c)
        self.pbasis = Basis(mesh, ElementTriP1(), intorder=c.intorder)
        self.wall = {side: self.vbasis.get_dofs(side).all()
                     for side in ("left", "right", "bottom")}
        left_normal = np.intersect1d(self.wall["left"], self.components[0])
        self.fixed = np.unique(np.r_[left_normal, self.wall["right"], self.wall["bottom"]])
        self.free = np.setdiff1d(np.arange(self.vbasis.N), self.fixed)
        div_vm = None
        if c.skew_divergence and advecting_velocity is not None:
            field = self.vbasis.interpolate(advecting_velocity)
            div_vm = field.grad[0,0]+field.grad[1,1]
        self.M, self.Kb, self.B, self.f_body, self.Kg = volume_operators(
            self.scalar, self.pbasis, self.components, div_vm)
        self.Kd = grad_div_operator(self.scalar,self.components,c.grad_div_gamma_m2_s)
        self.Kl,free_force = boundary_operators(mesh,self.components,c.intorder)
        self.f = free_force if c.hydrostatic_split else self.f_body

    def solve(self, v0, dt, linear_solver=None):
        ids = self.free
        M, K, B = self.M[ids][:, ids], (self.Kb+self.Kl+self.Kg+self.Kd)[ids][:, ids], self.B[:, ids]
        A = bmat([[M/dt+K*.5, -B.T], [B*.5, None]], format="csc")
        rhs = np.r_[(M/dt-K*.5)@v0[ids]+self.f[ids], -.5*B@v0[ids]]
        sol = spsolve(A, rhs) if linear_solver is None else linear_solver.solve(A, rhs, len(ids))
        if not np.all(np.isfinite(sol)): raise ValueError("nonfinite mixed solve")
        v = np.zeros_like(v0); v[ids] = sol[:len(ids)]
        return v, sol[len(ids):]

    def startup(self):
        ids = self.free; B = self.B[:, ids]
        A = bmat([[self.M[ids][:, ids], -B.T], [B, None]], format="csc")
        sol = spsolve(A, np.r_[self.f[ids], np.zeros(B.shape[0])])
        a = np.zeros(self.vbasis.N); a[ids] = sol[:len(ids)]
        return a, sol[len(ids):]

    def nodal_vector(self, v):
        return np.vstack([v[indices] for indices in self.components])

    def measure(self, v):
        coords = np.asarray(self.scalar.global_coordinates())
        volume = float(np.sum(self.scalar.dx))
        kinetic = .5*float(v@(self.M@v))
        # Domain z moment minus flat -100, double precision cancellation <1e-9 E0.
        potential = 9.81*(float(np.sum(coords[1]*self.scalar.dx))+100.)
        return {"volume_m2": volume, "kinetic": kinetic, "potential": potential,
                "energy": kinetic+potential}


def interface_nodes(mesh):
    """Ordered original topology, including every quadratic surface midpoint."""
    scalar = Basis(mesh, ElementTriP2())
    facets = mesh.boundaries["surface"]
    # Ordering is topological, not current x: steep/overhanging branches survive.
    edges = mesh.facets[:, facets]
    adjacency = {}
    for f, (a, b) in zip(facets, edges.T):
        mid = scalar.dofs.facet_dofs[0, f]
        adjacency.setdefault(int(a), []).append((int(b), int(mid)))
        adjacency.setdefault(int(b), []).append((int(a), int(mid)))
    endpoints = [a for a, nbrs in adjacency.items() if len(nbrs)==1]
    start = min(endpoints, key=lambda a: mesh.p[0, a])
    result, previous, current = [start], None, start
    while True:
        neighbors = [(a,m) for a,m in adjacency[current] if a != previous]
        if not neighbors: break
        nxt, mid = neighbors[0]; result += [mid, nxt]
        previous, current = current, nxt
    return np.array(result, dtype=int)


def sample_interface(mesh, ids, subdivisions=8):
    t = np.linspace(0.,1.,subdivisions+1)
    shape = np.vstack(((1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)))
    parts = [mesh.p[:, ids[k:k+3]]@shape for k in range(0,len(ids)-2,2)]
    return np.hstack([q[:,:-1] for q in parts]+[parts[-1][:,-1:]])


def step(mesh0, v0, dt, c, orientation):
    X0 = mesh0.p
    initial = Geometry(mesh0, c)
    vm_guess = v0.copy()
    v1 = v0.copy()
    linear_solver = PicardLinearSolver()
    for iteration in range(c.nonlinear_iterations):
        Xm = X0+.5*dt*initial.nodal_vector(vm_guess)
        mid = replace(mesh0, doflocs=Xm)
        mins, signs = quadratic_minima(mid)
        if np.min(mins) <= 0 or not np.array_equal(signs, orientation):
            raise ValueError("midpoint mesh Jacobian invalid")
        op = Operators(mid,c,vm_guess)
        v1, p = op.solve(v0,dt,linear_solver)
        vm = (v0+v1)*.5
        error = float(np.max(abs(vm-vm_guess)))
        if error < c.nonlinear_tolerance: break
        vm_guess = vm
    else:
        raise ValueError("nonlinear material iteration did not converge: "+str(error))
    X1 = X0+dt*op.nodal_vector(vm)
    mesh1 = replace(mesh0, doflocs=X1)
    mins, signs = quadratic_minima(mesh1)
    if np.min(mins)<=0 or not np.array_equal(signs,orientation):
        raise ValueError("accepted-candidate mesh Jacobian invalid")
    # Reconstruct on the actual converged midpoint, not the last lagged geometry.
    op = Operators(replace(mesh0, doflocs=(X0+X1)*.5),c,vm)
    residual = op.M@((v1-v0)/dt)+(op.Kb+op.Kl+op.Kg+op.Kd)@vm-op.B.T@p-op.f
    momentum = float(np.linalg.norm(residual[op.free],np.inf))
    weakdiv = float(np.linalg.norm(op.B@vm,np.inf))
    if momentum > 1e-8 or weakdiv > 1e-10:
        raise ValueError(f"equation residuals momentum={momentum}, div={weakdiv}")
    end = Geometry(mesh1,c)
    row = end.measure(v1)
    # Explicit geometry quadrature work, computed independently of E1-E0.
    vmq = np.asarray(op.vbasis.interpolate(vm))
    dvq = np.asarray(op.vbasis.interpolate(v1-v0))
    d0, d1, dm = initial.scalar.dx, end.scalar.dx, op.scalar.dx
    da, dj = (d0+d1)*.5, d1-d0
    z0 = np.asarray(initial.scalar.global_coordinates())[1]
    z1 = np.asarray(end.scalar.global_coordinates())[1]
    Wk = np.sum((da-dm)*np.sum(vmq*dvq,axis=0)+dj*(.5*np.sum(vmq**2,axis=0)+.125*np.sum(dvq**2,axis=0)))
    Wp = 9.81*np.sum((da-dm)*(z1-z0)+dj*(z1+z0)*.5)
    geometric_divergence_work = dt*float(vm@(op.Kg@vm))
    hydrostatic_split_work = dt*float(vm@(op.f-op.f_body))
    row.update(dt_s=dt, nonlinear_iterations=iteration+1, nonlinear_error=error,
        linear_solver=linear_solver.statistics,
        momentum_residual=momentum, weak_divergence=weakdiv,
        bulk_loss=dt*float(vm@(op.Kb@vm)), left_wall_loss=dt*float(vm@(op.Kl@vm)),
        numerical_grad_div_work=-dt*float(vm@(op.Kd@vm)),
        raw_geometry_kinetic_work=float(Wk),raw_geometry_potential_work=float(Wp),
        geometric_divergence_work=geometric_divergence_work,
        hydrostatic_split_work=hydrostatic_split_work,
        geometry_kinetic_work=float(Wk)-geometric_divergence_work,
        geometry_potential_work=float(Wp)+hydrostatic_split_work,
        min_jacobian=float(np.min(mins)), min_jacobian_cell=int(np.argmin(mins)))
    if abs(row["volume_m2"]-20.)/20.>1e-6:
        raise ValueError("volume acceptance failed: "+str(row["volume_m2"]))
    return mesh1, v1, p, row
