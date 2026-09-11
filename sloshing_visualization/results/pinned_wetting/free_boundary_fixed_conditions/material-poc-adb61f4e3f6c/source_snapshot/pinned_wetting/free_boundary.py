"""Full material NS on an evolving isoparametric P2 mesh, PW2 only.

No fixed-domain surface extrapolation and no Eulerian convection omission:
ALL geometry nodes move at the P2 material velocity. See MODEL_NOTES.md.
"""
from dataclasses import dataclass, replace, asdict
import numpy as np
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import spsolve
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


def quadratic_minima(mesh):
    """Exact minimum of signed P2 Jacobian determinant on each triangle."""
    pts = np.array([[0, 1, 0, .5, 0, .5], [0, 0, 1, 0, .5, .5]])
    val = mesh.mapping().detDF(pts)
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


class Operators:
    def __init__(self, mesh, c):
        self.mesh, self.controls = mesh, c
        self.scalar = Basis(mesh, ElementTriP2(), intorder=c.intorder)
        self.vbasis = Basis(mesh, ElementVector(ElementTriP2()), intorder=c.intorder)
        self.pbasis = Basis(mesh, ElementTriP1(), intorder=c.intorder)
        self.components = self.vbasis.split_indices()
        self.wall = {side: self.vbasis.get_dofs(side).all()
                     for side in ("left", "right", "bottom")}
        left_normal = np.intersect1d(self.wall["left"], self.components[0])
        self.fixed = np.unique(np.r_[left_normal, self.wall["right"], self.wall["bottom"]])
        self.free = np.setdiff1d(np.arange(self.vbasis.N), self.fixed)
        self.M = asm(mass, self.vbasis).tocsr()
        self.Kb = asm(strain, self.vbasis).tocsr()
        self.Kl = asm(left_friction, self.vbasis.boundary("left", intorder=c.intorder)).tocsr()
        self.B = asm(divergence, self.vbasis, self.pbasis).tocsr()
        self.f = asm(gravity, self.vbasis)

    def solve(self, v0, dt):
        ids = self.free
        M, K, B = self.M[ids][:, ids], (self.Kb+self.Kl)[ids][:, ids], self.B[:, ids]
        A = bmat([[M/dt+K*.5, -B.T], [B*.5, None]], format="csc")
        rhs = np.r_[(M/dt-K*.5)@v0[ids]+self.f[ids], -.5*B@v0[ids]]
        sol = spsolve(A, rhs)
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
    initial = Operators(mesh0, c)
    vm_guess = v0.copy()
    v1 = v0.copy()
    for iteration in range(c.nonlinear_iterations):
        Xm = X0+.5*dt*initial.nodal_vector(vm_guess)
        mid = replace(mesh0, doflocs=Xm)
        mins, signs = quadratic_minima(mid)
        if np.min(mins) <= 0 or not np.array_equal(signs, orientation):
            raise ValueError("midpoint mesh Jacobian invalid")
        op = Operators(mid,c)
        v1, p = op.solve(v0,dt)
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
    op = Operators(replace(mesh0, doflocs=(X0+X1)*.5),c)
    residual = op.M@((v1-v0)/dt)+(op.Kb+op.Kl)@vm-op.B.T@p-op.f
    momentum = float(np.linalg.norm(residual[op.free],np.inf))
    weakdiv = float(np.linalg.norm(op.B@vm,np.inf))
    if momentum > 1e-8 or weakdiv > 1e-10:
        raise ValueError(f"equation residuals momentum={momentum}, div={weakdiv}")
    end = Operators(mesh1,c)
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
    row.update(dt_s=dt, nonlinear_iterations=iteration+1, nonlinear_error=error,
        momentum_residual=momentum, weak_divergence=weakdiv,
        bulk_loss=dt*float(vm@(op.Kb@vm)), left_wall_loss=dt*float(vm@(op.Kl@vm)),
        geometry_kinetic_work=float(Wk), geometry_potential_work=float(Wp),
        min_jacobian=float(np.min(mins)), min_jacobian_cell=int(np.argmin(mins)))
    if abs(row["volume_m2"]-20.)/20.>1e-6:
        raise ValueError("volume acceptance failed: "+str(row["volume_m2"]))
    return mesh1, v1, p, row
