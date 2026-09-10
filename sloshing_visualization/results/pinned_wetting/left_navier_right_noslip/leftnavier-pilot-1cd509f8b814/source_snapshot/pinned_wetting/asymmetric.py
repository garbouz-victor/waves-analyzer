"""PW1 correction: LEFT Navier, RIGHT and BOTTOM full no-slip.

The symmetric solver and its immutable native bundles remain historical.
Only the SDIRK2 integrator/energy-state machinery is shared with that solver.
"""
from dataclasses import dataclass

import numpy as np
from skfem import BilinearForm, asm

from ..config import SimulationConfig
from ..diagnostics import Diagnostics, quadratic_edge_max
from ..fem_spaces import FEMSystem

MODEL_ID = "PW1_LEFT_NAVIER_RIGHT_NOSLIP_V1"
LABEL = "DECLARED_EXTENSION — LEFT NAVIER / RIGHT NO-SLIP"
CONTACT_METHOD = "actual_FE_trace_left_free_right_fixed"
SOURCE_KIND = "declared_left_Navier_right_noslip_linear_FEM"


@dataclass(frozen=True)
class LeftNavierConfig(SimulationConfig):
    left_slip_length_m: float = .50

    def __post_init__(self):
        super().__post_init__()
        if not np.isfinite(self.left_slip_length_m) or self.left_slip_length_m <= 0:
            raise ValueError("Finite positive LEFT slip length required")
        if self.integrator != "sdirk2":
            raise ValueError("PW1 corrected model uses SDIRK2")


@BilinearForm
def left_tangential_mass(u, v, _):
    return u[1] * v[1]


class LeftNavierFEM(FEMSystem):
    def __init__(self, config):
        super().__init__(config)
        left_normal = np.intersect1d(self.wall_dofs["left"], self.component_dofs[0])
        self.fixed = np.unique(np.r_[left_normal, self.wall_dofs["right"], self.wall_dofs["bottom"]])
        self.free = np.setdiff1d(np.arange(self.velocity.N), self.fixed)
        self.left_basis = self.velocity.boundary("left", intorder=6)
        self.K_bulk_full = self.K_full.copy()
        self.K_wall_full = (config.nu/config.left_slip_length_m * asm(left_tangential_mass, self.left_basis)).tocsr()
        self.K_full = self.K_bulk_full + self.K_wall_full
        self.M = self.M_full[self.free][:, self.free].tocsr()
        self.K_bulk = self.K_bulk_full[self.free][:, self.free].tocsr()
        self.K_wall = self.K_wall_full[self.free][:, self.free].tocsr()
        self.K = self.K_bulk + self.K_wall
        self.B = self.B_full[:, self.free].tocsr()
        self.R = self.R_full[:, self.free].tocsr()
        self.C = (self.R.T @ self.S).tocsr()
        self.H = (self.C @ self.R).tocsr()
        self.surface_basis = self.velocity.boundary("surface", intorder=6)
        if [self.R.getrow(i).nnz for i in (0, self.R.shape[0]-1)] != [1, 0]:
            raise ValueError("LEFT endpoint must be free; RIGHT endpoint must be Dirichlet fixed")


class LeftNavierDiagnostics(Diagnostics):
    def measure(self, state):
        row = super().measure(state)
        f, c = self.fem, self.fem.config
        full = f.expand_velocity(state.velocity)
        eta_max = quadratic_edge_max(state.eta[:-2:2], state.eta[1::2], state.eta[2::2])
        wall_loss = float(state.velocity @ (f.K_wall @ state.velocity))
        row.update(max_eta_over_a=eta_max/c.a, max_eta_over_b=eta_max/c.left_slip_length_m,
            right_endpoint_error=abs(float(state.eta[-1]-self.contacts[1])),
            wall_dissipation=wall_loss, left_wall_dissipation=wall_loss,
            bulk_dissipation=float(state.velocity @ (f.K_bulk @ state.velocity)),
            cumulative_wall_dissipation=state.wall_dissipated_energy,
            cumulative_left_wall_dissipation=state.wall_dissipated_energy,
            cumulative_bulk_dissipation=state.bulk_dissipated_energy)
        field = f.velocity.interpolate(full)
        adv = np.einsum("ij...,j...->i...", field.grad, np.asarray(field))
        row["convective_indicator"] = float(np.max(np.sqrt(np.sum(adv**2, axis=0)))) / (c.g*max(abs(c.slope),1e-12))
        trace = f.surface_basis.interpolate(full)
        eta_scalar = np.zeros(f.scalar.N)
        eta_scalar[f.surface_dofs] = state.eta
        eta_field = f.scalar.boundary("surface", intorder=6).interpolate(eta_scalar)
        omitted = np.asarray(eta_field)*trace.grad[1,1] - np.asarray(trace)[0]*eta_field.grad[0]
        row["kinematic_indicator"] = float(np.max(abs(omitted))) / (np.sqrt(c.g*c.a)*max(abs(c.slope),1e-12))
        value = f.left_basis.interpolate(full)
        # LEFT outward normal is -ex: traction + nu/b*w must vanish.
        residual = c.nu*(-value.grad[1,0]-value.grad[0,1]+np.asarray(value)[1]/c.left_slip_length_m)
        row["wall_traction_L2"] = row["left_wall_traction_L2"] = float(np.sqrt(np.sum(residual**2*f.left_basis.dx)))
        wall_power = c.nu/c.left_slip_length_m*float(np.sum(np.asarray(value)[1]**2*f.left_basis.dx))
        row["wall_power_quadrature_error"] = abs(wall_power-wall_loss)
        return row

    def validate(self, row):
        if not all(np.isfinite(v) for v in row.values()):
            raise ValueError("Nonfinite asymmetric fields/diagnostics")
        limits = {"weak_divergence_l2":1e-9, "volume_error":1e-10,
            "energy_balance_relative":1e-8, "left_u_max":1e-12,
            "right_u_max":1e-12, "right_w_max":1e-12,
            "bottom_u_max":1e-12, "bottom_w_max":1e-12,
            "right_endpoint_error":1e-12, "wall_power_quadrature_error":1e-11}
        bad = {k:row[k] for k,tol in limits.items() if abs(row[k])>tol}
        if bad or min(row["wall_dissipation"],row["bulk_dissipation"]) < -1e-15:
            raise ValueError("Asymmetric boundary/energy constraint violation: "+repr(bad))


def save_operator_evidence(h, f):
    """Full, not reduced, wall operator: an erroneous RIGHT term cannot hide."""
    op = h.create_group("operators")
    wall = op.create_group("K_wall_full")
    matrix = f.K_wall_full.copy().tocsr()
    matrix.sum_duplicates(); matrix.sort_indices()
    for name in ("data", "indices", "indptr"):
        wall.create_dataset(name, data=getattr(matrix,name))
    wall.attrs["shape"] = matrix.shape
    op.create_dataset("fixed_velocity_dofs",data=f.fixed)
    op.attrs.update(endpoint_trace_nnz=[1,0],wall_support="left_only")
    h.attrs["model_id"] = MODEL_ID


def channel_benchmark(config):
    """Independent scalar FE channel; does not instantiate the production FEM.

No vessel bottom/top Dirichlet conditions in this separate steady channel test.
"""
    from scipy.sparse.linalg import spsolve
    from skfem import Basis, ElementTriP2, MeshTri, LinearForm
    from skfem.helpers import dot, grad
    c = config
    mesh = MeshTri.init_tensor(np.linspace(-c.a,c.a,9),np.linspace(-1,0,4)).with_boundaries(
        {"left":lambda x:np.isclose(x[0],-c.a), "right":lambda x:np.isclose(x[0],c.a)})
    basis = Basis(mesh,ElementTriP2(),intorder=6)
    @BilinearForm
    def laplace(u,v,_): return c.nu*dot(grad(u),grad(v))
    @BilinearForm
    def robin(u,v,_): return c.nu/c.left_slip_length_m*u*v
    @LinearForm
    def load(v,_): return .001*v
    matrix = asm(laplace,basis)+asm(robin,basis.boundary("left"))
    right = basis.get_dofs("right").all()
    free = np.setdiff1d(np.arange(basis.N),right)
    w = np.zeros(basis.N)
    w[free] = spsolve(matrix[free][:,free],asm(load,basis)[free])
    x, a, b, nu, G = basis.doflocs[0],c.a,c.left_slip_length_m,c.nu,.001
    exact = G/(2*nu)*(a*a-x*x)+G*a*b/(nu*(b+2*a))*(a-x)
    return {"benchmark":"steady_channel_LEFT_Navier_RIGHT_no_slip", "left_slip_length_m":b,
        "relative_error":float(np.max(abs(w-exact))/np.max(abs(exact))),
        "max_absolute_error_m_s":float(np.max(abs(w-exact))),
        "right_velocity_max":float(np.max(abs(w[right]))),
        "analytic_left_robin_residual":abs(2*G*a*a/(nu*(b+2*a))-(2*G*a*a*b/(nu*(b+2*a)))/b),
        "exact_formula":"G/(2*nu)*(a*a-x*x)+G*a*b/(nu*(b+2*a))*(a-x)","G_m_s2":G}
