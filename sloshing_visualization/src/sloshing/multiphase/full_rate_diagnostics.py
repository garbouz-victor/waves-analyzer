"""Independent coupled physical/constraint observers; same quadrature."""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu


def assemble_vector(form):
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from petsc4py import PETSc
    vector=fp.assemble_vector(fem.form(form))
    vector.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
    a=vector.array.copy(); vector.destroy(); return a


class FieldMetric:
    def __init__(self,solver,space,*,riesz=False):
        import ufl
        from dolfinx import fem
        from dolfinx.fem import petsc as fp
        from .equilibrium import _measures
        dx,_=_measures(solver)
        a=fp.assemble_matrix(fem.form(ufl.inner(ufl.TrialFunction(space),ufl.TestFunction(space))*dx))
        a.assemble(); ptr,col,val=a.getValuesCSR()
        self.M=csr_matrix((val.copy(),col.copy(),ptr.copy()),shape=a.getSize()); a.destroy()
        self.lu=splu(self.M.tocsc()) if riesz else None

    def norm(self,q): return float(np.sqrt(max(0.,q @ (self.M@q))))

    def riesz_norm(self,r): return float(np.sqrt(max(0.,r @ self.lu.solve(r))))


def constrained_maps(solver):
    maps=[np.asarray(solver.space.sub(i).collapse()[1]) for i in range(4)]
    constrained=np.unique(np.concatenate([b.dof_indices()[0] for b in solver.bcs]))
    if len(np.intersect1d(constrained,np.r_[maps[2],maps[3]])):
        raise ValueError("BC erroneously applied to phase/rate or chemical block")
    velocity=np.intersect1d(constrained,maps[0]); gauge=np.intersect1d(constrained,maps[1])
    if len(gauge)!=1: raise ValueError("Exactly the original single pressure gauge is required")
    return velocity,gauge


class FullObserver:
    def __init__(self,solver,phase_metric):
        import ufl
        from .equilibrium import _measures
        self.s=solver; self.phase_metric=phase_metric
        self.spaces=[solver.space.sub(i).collapse()[0] for i in range(4)]
        self.maps=[np.asarray(solver.space.sub(i).collapse()[1]) for i in range(4)]
        self.velocity_metric=FieldMetric(solver,self.spaces[0])
        self.pressure_metric=FieldMetric(solver,self.spaces[1],riesz=True)
        self.velocity_bcs,self.gauge=constrained_maps(solver)
        c=solver.config
        self.U_ref=float(np.sqrt(c.sigma/(c.rho_liquid*c.epsilon))); self.P_ref=float(c.sigma/c.epsilon)
        u,pi,phi,mu=ufl.split(solver.state); dx,_=_measures(solver)
        self.continuity=ufl.TestFunction(self.spaces[1])*ufl.div(u)*dx
        self.advection=-phi*ufl.inner(u,ufl.grad(ufl.TestFunction(self.spaces[2])))*dx

    def measure(self,row,policy,*,expanded=False):
        s=self.s; a=s.state.x.array
        velocity=a[self.maps[0]]; pressure=a[self.maps[1]]
        continuity=assemble_vector(self.continuity)
        vbc=float(np.max(abs(a[self.velocity_bcs]),initial=0.)); gauge=float(np.max(abs(a[self.gauge]),initial=0.))
        result={"F_CH":row["E_interface"]+row["E_wall"],
            "phi_L2":self.phase_metric.norm(a[self.maps[2]]),"mu_L2":self.phase_metric.norm(a[self.maps[3]]),
            "velocity_L2":self.velocity_metric.norm(velocity),"pressure_L2":self.pressure_metric.norm(pressure),
            "pressure_min_dof":float(pressure.min()),"pressure_max_dof":float(pressure.max()),
            "pressure_gauge_value":float(a[self.gauge[0]]),"velocity_BC_max":vbc,
            "velocity_BC_scaled":vbc/self.U_ref,"pressure_gauge_scaled":gauge/self.P_ref,
            "weak_continuity_algebraic":float(np.linalg.norm(continuity)),
            "weak_continuity_Riesz":self.pressure_metric.riesz_norm(continuity)}
        checks={"velocity_BC":result["velocity_BC_scaled"]<=policy["BC_scaled"],
            "pressure_gauge":result["pressure_gauge_scaled"]<=policy["BC_scaled"],
            "weak_continuity":result["weak_continuity_Riesz"]<=policy["weak_continuity"],
            "strong_divergence_range":row["strong_divergence_L2"]<=policy["strong_divergence_max"],
            "finite_full":bool(np.isfinite(list(result.values())).all())}
        result["full_physical_checks"]={k:bool(v) for k,v in checks.items()}
        if expanded:
            adv=assemble_vector(self.advection)
            diff=s.config.mobility*(self.phase_metric.K @ a[self.maps[3]])
            an=self.phase_metric.riesz_norm(adv); dn=self.phase_metric.riesz_norm(diff)
            significant=row["CH_dissipation"]>=policy["power_floor_fraction"]*policy["initial_D"]
            result["phase_coupling_audit"]={"advection_Riesz":an,"diffusion_Riesz":dn,
                "advective_over_diffusive":an/dn if significant and dn>0 else None,
                "diffusion_significant":bool(significant),"diagnostic_only":True}
        return result
