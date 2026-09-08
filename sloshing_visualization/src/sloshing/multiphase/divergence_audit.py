"""Independent Q_h projection of divergence, not a renamed observer Riesz norm."""
import numpy as np
from .divergence_policy import projection_checks


class DivergenceAudit:
    def __init__(self, solver):
        import ufl
        from dolfinx import fem
        from dolfinx.fem import petsc as fp
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import splu
        from .equilibrium import _measures
        if solver.comm.size != 1: raise ValueError("Serial audit only")
        self.s = solver
        self.Q = solver.space.sub(1).collapse()[0]
        self.projected = fem.Function(self.Q)
        dx, _ = _measures(solver)
        q = ufl.TestFunction(self.Q); p = ufl.TrialFunction(self.Q)
        # Deliberately independent assembly/storage/factorization from FullObserver.
        matrix = fp.assemble_matrix(fem.form(p*q*dx)); matrix.assemble()
        ptr, col, val = matrix.getValuesCSR()
        self.M = csr_matrix((val.copy(), col.copy(), ptr.copy()), shape=matrix.getSize())
        matrix.destroy(); self.lu = splu(self.M.tocsc())
        u = ufl.split(solver.state)[0]; d = ufl.div(u)
        self.load = fem.form(q*d*dx)
        self.forms = {"strong_divergence_L2": fem.form(d*d*dx),
            "projected_divergence_L2": fem.form(self.projected**2*dx),
            "orthogonal_divergence_L2": fem.form((d-self.projected)**2*dx),
            "grad_u_L2": fem.form(ufl.inner(ufl.grad(u),ufl.grad(u))*dx),
            "symgrad_u_L2": fem.form(ufl.inner(ufl.sym(ufl.grad(u)),ufl.sym(ufl.grad(u)))*dx)}

    def measure(self, existing_weak_riesz):
        from dolfinx import fem
        from dolfinx.fem import petsc as fp
        from petsc4py import PETSc
        vector = fp.assemble_vector(self.load)
        vector.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
        b = vector.array.copy(); vector.destroy()
        self.projected.x.array[:] = self.lu.solve(b); self.projected.x.scatter_forward()
        result = {}
        for name, form in self.forms.items():
            squared = float(fem.assemble_scalar(form))
            result[name] = float(np.sqrt(squared)) if squared >= 0 else float("nan")
        strong = result["strong_divergence_L2"]; projected = result["projected_divergence_L2"]
        result.update(weak_continuity_algebraic_independent=float(np.linalg.norm(b)),
            weak_continuity_Riesz=float(existing_weak_riesz),
            projection_linear_residual=float(np.linalg.norm(self.M@self.projected.x.array-b)),
            projection_fraction=projected/max(strong,1e-30))
        result.update(projection_checks(strong,projected,existing_weak_riesz,result["orthogonal_divergence_L2"]))
        return result

    def cells(self):
        import ufl
        from dolfinx import fem
        from .equilibrium import _measures
        from .full_rate_diagnostics import assemble_vector
        s = self.s; dx, _ = _measures(s)
        D = fem.functionspace(s.mesh, ("DG",0)); d = ufl.div(ufl.split(s.state)[0])
        raw = assemble_vector(ufl.TestFunction(D)*d*d*dx)
        n = s.mesh.topology.index_map(2).size_local
        contributions = raw[np.array([D.dofmap.cell_dofs(i)[0] for i in range(n)])]
        tri = s.mesh.geometry.x[s.mesh.geometry.dofmap[:n], :2]
        centers = tri.mean(axis=1)
        phi = s.state.sub(2).collapse()
        points = np.array([[0.,0.],[1.,0.],[0.,1.],[.5,0.],[.5,.5],[0.,.5]])
        values = fem.Expression(phi,points).eval(s.mesh,np.arange(n,dtype=np.int32)).reshape(n,6)
        edge = 2*values[:,3:]-.5*values[:,[0,1,2]]-.5*values[:,[1,2,0]]
        bern = np.column_stack((values[:,:3],edge))
        interface = (bern.min(axis=1)<=.9)&(bern.max(axis=1)>=-.9)
        c = s.config
        distance = np.column_stack((centers[:,0]-c.x_min,c.x_max-centers[:,0],
                                    centers[:,1]-c.z_min,c.z_max-centers[:,1]))
        wall = distance.min(axis=1)<=2*c.epsilon
        contact = interface & (centers[:,1]-c.z_min<=2*c.epsilon)
        total = float(contributions.sum()); denom = max(total,1e-300)
        ranked = np.sort(contributions)[::-1]
        shares = {str(p):float(ranked[:int(np.ceil(n*p/100))].sum()/denom) for p in (1,5,10)}
        regions = {name:{"cells":int(mask.sum()),"squared_norm_fraction":float(contributions[mask].sum()/denom)}
            for name,mask in {"interface_active":interface,"wall_near":wall,
                "contact_near":contact,"bulk":~(interface|wall)}.items()}
        stats = {"cells":n,"sum_divergence_squared":total,"max_cell_contribution":float(ranked[0]),
            "percentiles":{str(p):float(np.percentile(contributions,p)) for p in (50,90,95,99)},
            "top_percent_squared_norm_fractions":shares,
            "top_percent_norm_fractions":{k:float(np.sqrt(v)) for k,v in shares.items()},
            "localization":regions,"localization_note":"interface/wall/contact regions overlap; bulk excludes interface and wall"}
        return {"triangles":tri,"contributions":contributions,"phi_vertices":values[:,:3],
                "interface_active":interface,"wall_near":wall,"contact_near":contact}, stats


def kernel_demonstration(solver):
    """Tiny-only SVD projects a deterministic velocity into the actual ker(B)."""
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from scipy.linalg import null_space
    from .equilibrium import _measures
    from .full_rate_diagnostics import FullObserver
    from .phase_rate import PhaseMetric
    V, vm = solver.space.sub(0).collapse(); Q = solver.space.sub(1).collapse()[0]
    if V.dofmap.index_map.size_global*V.dofmap.index_map_bs > 2000:
        raise ValueError("Kernel demonstration is tiny-mesh only")
    dx,_ = _measures(solver)
    mat = fp.assemble_matrix(fem.form(ufl.TestFunction(Q)*ufl.div(ufl.TrialFunction(V))*dx)); mat.assemble()
    ptr,col,val = mat.getValuesCSR()
    from scipy.sparse import csr_matrix
    B = csr_matrix((val.copy(),col.copy(),ptr.copy()),shape=mat.getSize()).toarray(); mat.destroy()
    constrained = np.unique(np.concatenate([bc.dof_indices()[0] for bc in solver.bcs]))
    free = np.flatnonzero(~np.isin(np.asarray(vm),constrained))
    N = null_space(B[:,free],rcond=1e-12)
    seed = np.random.default_rng(309).normal(size=len(free)); v = np.zeros(B.shape[1])
    v[free] = N@(N.T@seed); v /= np.linalg.norm(v)
    saved = solver.state.x.array.copy()
    try:
        solver.state.x.array[np.asarray(vm)] = v; solver.state.x.scatter_forward()
        phase = PhaseMetric(solver,solver.space.sub(2).collapse()[0])
        obs = FullObserver(solver,phase)
        from .full_rate_diagnostics import assemble_vector
        weak = obs.pressure_metric.riesz_norm(assemble_vector(obs.continuity))
        result = DivergenceAudit(solver).measure(weak)
        result.update(kernel_dimension=N.shape[1],B_u_algebraic=float(np.linalg.norm(B@v)),
                      velocity_coefficient_norm=float(np.linalg.norm(v)))
        result["passed"] = bool(weak < 1e-12 and result["strong_divergence_L2"] > 1e-3 and
                                result["projection_agrees"] and result["pythagorean_pass"])
        return result
    finally:
        solver.state.x.array[:] = saved; solver.state.x.scatter_forward()
