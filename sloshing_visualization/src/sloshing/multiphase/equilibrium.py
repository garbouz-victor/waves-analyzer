"""Mass-constrained discrete free-energy equilibrium; one true global scalar.

The serial augmented SNES has N phase coefficients plus ONE algebraic mu DOF.
Its last row is the global mass constraint. No P0 multiplier, field clipping,
post-solve mass shift, change of wall physics, or modified functional is used.
The stationary chemical form is written independently of transient assembly;
container tests compare its residual term by term against that weak form.
"""
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import numpy as np

ALGORITHM_VERSION = "constrained-p2-global-scalar-snes-v2"
CHECKPOINT_SCHEMA = "step3a2-prepared-variational-equilibrium-v1"


@dataclass(frozen=True)
class EquilibriumPolicy:
    mass_domain_tolerance: float = 1e-10
    normalized_riesz_tolerance: float = 1e-8
    snes_atol: float = 1e-13
    snes_rtol: float = 1e-12
    snes_max_it: int = 80
    mu_spread_tolerance: float = 1e-12
    initial_ch_tolerance: float = 1e-16
    wall_improvement_ratio: float = .5
    first_variation_tolerance: float = 1e-8
    line_search: str = "nleqerr"


@dataclass
class PreparedEquilibrium:
    phi_coefficients: np.ndarray
    metadata: dict


class EquilibriumFailure(RuntimeError):
    """Diagnostic candidate is explicitly incomplete and cannot initialize CHNS."""
    def __init__(self, message, coefficients, diagnostics):
        super().__init__(message)
        self.coefficients=coefficients
        self.diagnostics=diagnostics


def safeguarded_scalar_root(function, bracket, tolerance=1e-12, max_iterations=60):
    """Bounded alternative for a fixed-mu inner solver; never expands a bracket.

    Function evaluations must themselves reject branch departures. This utility
    is tested independently; the default equilibrium uses augmented SNES.
    """
    a,b=map(float,bracket)
    if not np.isfinite([a,b]).all() or not a<b:
        raise ValueError("Invalid scalar bracket")
    history=[]
    def evaluate(x):
        value=float(function(x))
        if not np.isfinite(value):
            raise RuntimeError("Nonfinite mass residual; branch rejected")
        history.append({"mu_star":x,"mass_residual":value})
        return value
    fa,fb=evaluate(a),evaluate(b)
    if abs(fa)<=tolerance:
        return a,history
    if abs(fb)<=tolerance:
        return b,history
    if np.signbit(fa)==np.signbit(fb):
        raise RuntimeError("Mass root not bracketed; no branch expansion permitted")
    for _ in range(max_iterations):
        x=(a*fb-b*fa)/(fb-fa)
        if not a+.1*(b-a)<x<b-.1*(b-a):
            x=(a+b)/2
        fx=evaluate(x)
        if abs(fx)<=tolerance:
            return x,history
        if np.signbit(fx)==np.signbit(fa):
            a,fa=x,fx
        else:
            b,fb=x,fx
    raise RuntimeError("Bounded mass root did not converge")


def equilibrium_config_fingerprint(config):
    """Same physics/mesh/FE/quadrature, explicitly independent of time schedule.

    Changing dt/scheme is the purpose of preservation/refinement, so these are
    separately recorded as the transient full config fingerprint.
    """
    values=config.as_dict()
    for key in ("name","dt","t_end","time_scheme","startup_dt","startup_t_end",
                "startup_stages","snes_atol","snes_rtol","snes_max_it","max_unknowns"):
        values.pop(key)
    return hashlib.sha256(json.dumps(values,sort_keys=True).encode()).hexdigest()


def algorithm_hash():
    base=Path(__file__).parent
    digest=hashlib.sha256()
    for name in ("equilibrium.py","solver.py","weak_form.py","free_energy.py","mesh.py",
                 "material.py","boundary.py","p2_certification.py"):
        digest.update(name.encode());digest.update((base/name).read_bytes())
    return digest.hexdigest()


def _measures(solver):
    import ufl
    metadata={"quadrature_degree":solver.config.quadrature_degree}
    return (ufl.Measure("dx",domain=solver.mesh,metadata=metadata),
            ufl.Measure("ds",domain=solver.mesh,subdomain_data=solver.tags,metadata=metadata))


def _integral(solver, expression):
    from dolfinx import fem
    from mpi4py import MPI
    return solver.comm.allreduce(float(fem.assemble_scalar(fem.form(expression))),op=MPI.SUM)


def energy_form(solver, phi):
    import ufl
    from .free_energy import bulk_energy, wall_energy, lambda_from_sigma
    from .mesh import WALL_IDS
    c=solver.config
    dx,ds=_measures(solver)
    result=(bulk_energy(phi,c.sigma,c.epsilon)+lambda_from_sigma(c.sigma)*c.epsilon/2*
            ufl.inner(ufl.grad(phi),ufl.grad(phi)))*dx
    for wall in c.wetting_walls:
        result+=wall_energy(phi,c.sigma,c.theta_equilibrium_deg)*ds(WALL_IDS[wall])
    return result


def chemical_stationary_form(solver, phi, mu_star, test):
    import ufl
    from .free_energy import bulk_derivative, wall_derivative, lambda_from_sigma
    from .mesh import WALL_IDS
    c=solver.config
    dx,ds=_measures(solver)
    result=((bulk_derivative(phi,c.sigma,c.epsilon)-mu_star)*test+
            lambda_from_sigma(c.sigma)*c.epsilon*ufl.dot(ufl.grad(phi),ufl.grad(test)))*dx
    for wall in c.wetting_walls:
        result+=wall_derivative(phi,c.sigma,c.theta_equilibrium_deg)*test*ds(WALL_IDS[wall])
    return result


def stationarity(solver, phi, mu_star):
    """L2 Riesz representative, not a raw coefficient-vector norm."""
    import ufl
    from dolfinx.fem.petsc import LinearProblem
    Q=phi.function_space
    trial,test=ufl.TrialFunction(Q),ufl.TestFunction(Q)
    dx,_=_measures(solver)
    problem=LinearProblem(trial*test*dx,chemical_stationary_form(solver,phi,mu_star,test),
        petsc_options_prefix="step3a2_riesz_",petsc_options={"ksp_type":"preonly","pc_type":"lu",
        "pc_factor_mat_solver_type":"mumps","ksp_error_if_not_converged":True})
    r=problem.solve()
    norm=np.sqrt(max(0.,_integral(solver,r*r*dx)))
    area=_integral(solver,1*dx)
    scale=solver.config.sigma/solver.config.length_scale_m
    return {"L2":float(norm),"normalized":float(norm/(scale*np.sqrt(area))),
            "mu_scale_Pa":scale,"domain_area":area,
            "definition":"||L2 Riesz(EL defect)|| / ((sigma/L_reference)*sqrt(area))"}


def solve_equilibrium(solver, target_mass=None, mu_initial=0., policy=None, history_callback=None):
    """Solve from the current phase, with its mass as default reproducible target.

    Serial-only augmented assembly is explicit: the saved partition fingerprint
    prevents reuse in a different MPI partition. No DOLFINx Real-space support
    is assumed. PETSc SNES line search and MUMPS factor the actual bordered
    Jacobian of the chemical and global mass equations.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.fem import petsc as fp
    from petsc4py import PETSc
    from scipy.sparse import bmat, csr_matrix
    from .provenance import run_provenance
    if solver.comm.size!=1:
        raise ValueError("Prepared augmented SNES currently requires one MPI rank")
    c=solver.config
    if c.rho_liquid!=c.rho_gas or c.g!=0 or c.a_x!=0:
        raise ValueError("STEP 3A.2 equilibrium supports matched density and zero body force only")
    policy=policy or EquilibriumPolicy()
    phi=solver.state.sub(2).collapse()
    n=len(phi.x.array)
    Q=phi.function_space
    test,trial=ufl.TestFunction(Q),ufl.TrialFunction(Q)
    dx,_=_measures(solver)
    area=_integral(solver,1*dx)
    target_mass=_integral(solver,phi*dx) if target_mass is None else float(target_mass)
    mu=fem.Constant(solver.mesh,PETSc.ScalarType(mu_initial))
    form=chemical_stationary_form(solver,phi,mu,test)
    residual_form=fem.form(form)
    jac_form=fem.form(ufl.derivative(form,phi,trial))
    mass_vector=fp.assemble_vector(fem.form(test*dx))
    mass_vector.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
    m=mass_vector.array.copy()
    native=fp.assemble_matrix(jac_form);native.assemble()
    def bordered():
        ip,ix,data=native.getValuesCSR()
        a=csr_matrix((data,ix,ip),shape=(n,n))
        return bmat([[a,csr_matrix(-m[:,None])],[csr_matrix(m[None,:]),csr_matrix((1,1))]],format="csr")
    block=bordered()
    jac=PETSc.Mat().createAIJ(size=(n+1,n+1),csr=(block.indptr,block.indices,block.data),comm=solver.comm)
    jac.assemble()
    x=PETSc.Vec().createSeq(n+1,comm=solver.comm)
    x.array[:n]=phi.x.array;x.array[n]=mu_initial
    f=x.duplicate()
    def sync(vector):
        phi.x.array[:]=vector.array_r[:n]
        phi.x.scatter_forward()
        mu.value=vector.array_r[n]
    def residual_callback(snes,vector,out):
        sync(vector)
        r=fp.assemble_vector(residual_form)
        r.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
        out.array[:n]=r.array
        out.array[n]=float(np.dot(m,phi.x.array))-target_mass
        r.destroy()
    def jacobian_callback(snes,vector,A,P):
        sync(vector)
        native.zeroEntries();fp.assemble_matrix(native,jac_form);native.assemble()
        block=bordered()
        A.zeroEntries();A.setValuesCSR(block.indptr,block.indices,block.data);A.assemble()
    history=[]
    energy=fem.form(energy_form(solver,phi))
    def monitor(snes,iteration,norm):
        sync(snes.getSolution())
        mass=float(np.dot(m,phi.x.array))
        row={"iteration":iteration,"mu_star":float(mu.value),"phase_mass":mass,
             "mass_residual":mass-target_mass,"snes_residual":float(norm),
             "energy":float(fem.assemble_scalar(energy)),
             "newton_update_norm":float(snes.getSolutionUpdate().norm()) if iteration else 0.}
        history.append(row)
        if history_callback:
            history_callback(row)
    snes=PETSc.SNES().create(solver.comm)
    snes.setOptionsPrefix("step3a2_equilibrium_")
    options=PETSc.Options()
    for key,value in {"snes_type":"newtonls","snes_linesearch_type":policy.line_search,"snes_stol":0.,
        "ksp_type":"preonly","pc_type":"lu","pc_factor_mat_solver_type":"mumps"}.items():
        options["step3a2_equilibrium_"+key]=value
    snes.setFunction(residual_callback,f);snes.setJacobian(jacobian_callback,jac)
    snes.setFromOptions()
    snes.setTolerances(atol=policy.snes_atol,rtol=policy.snes_rtol,stol=0.,max_it=policy.snes_max_it)
    snes.setMonitor(monitor)
    try:
        snes.solve(None,x)
        sync(x)
        mass=_integral(solver,phi*dx)
        riesz=stationarity(solver,phi,mu)
        if (snes.getConvergedReason()<=0 or abs(mass-target_mass)/area>policy.mass_domain_tolerance or
                riesz["normalized"]>policy.normalized_riesz_tolerance):
            # Inspect the phase Hessian near zero to distinguish poorly resolved
            # stationary modes from a chemical-form/constraint assembly error.
            from scipy.sparse.linalg import eigsh
            ip,ix,data=native.getValuesCSR()
            phase_hessian=csr_matrix((data,ix,ip),shape=(n,n))
            eigenvalues,eigenvectors=eigsh(phase_hessian,k=3,sigma=0.)
            mode=fem.Function(Q)
            mode.x.array[:]=eigenvectors[:,np.argmin(abs(eigenvalues))]
            mode_norm=np.sqrt(_integral(solver,mode*mode*dx))
            translation_norm=np.sqrt(_integral(solver,phi.dx(0)**2*dx))
            overlap=_integral(solver,mode*phi.dx(0)*dx)/max(mode_norm*translation_norm,1e-300)
            details={"status":"incomplete","diagnostic_only":True,"target_mass":target_mass,
                "phase_mass":mass,"mass_residual":mass-target_mass,"mu_star_candidate":float(mu.value),
                "stationarity":riesz,"snes_reason":snes.getConvergedReason(),"snes_history":history,
                "phase_hessian_near_zero_eigenvalues":eigenvalues.tolist(),
                "nearest_mode_L2_overlap_with_translation":float(overlap),"provenance":run_provenance(solver)}
            raise EquilibriumFailure("PREPARED EQUILIBRIUM NOT OBTAINED: direct stationary solve failed",
                                     phi.x.array.copy(),details)
        metadata={"schema":CHECKPOINT_SCHEMA,"status":"complete","algorithm_version":ALGORITHM_VERSION,
            "algorithm_sha256":algorithm_hash(),"target_mass":target_mass,"phase_mass":mass,
            "mass_residual":mass-target_mass,"mass_residual_relative_to_domain":abs(mass-target_mass)/area,
            "free_energy":_integral(solver,energy_form(solver,phi)),"mu_star":float(mu.value),
            "stationarity":riesz,"policy":asdict(policy),"snes_history":history,
            "snes_reason":snes.getConvergedReason(),"snes_iterations":snes.getIterationNumber(),
            "outer_root_history":[],"outer_root_note":"not applicable: one coupled global-scalar augmented SNES",
            "quadrature_degree":c.quadrature_degree,"config":c.as_dict(),
            "equilibrium_config_fingerprint":equilibrium_config_fingerprint(c),
            "provenance":run_provenance(solver)}
        coefficients=phi.x.array.copy()
        metadata["coefficient_sha256"]=hashlib.sha256(coefficients.tobytes()).hexdigest()
        metadata["fingerprint"]=hashlib.sha256(json.dumps(metadata,sort_keys=True).encode()).hexdigest()
        return PreparedEquilibrium(coefficients,metadata)
    finally:
        snes.destroy();x.destroy();f.destroy();jac.destroy();native.destroy();mass_vector.destroy()


def validate_prepared(solver, prepared):
    from .provenance import run_provenance
    meta=prepared.metadata
    if meta.get("schema")!=CHECKPOINT_SCHEMA or meta.get("status")!="complete":
        raise ValueError("Incomplete prepared equilibrium cache")
    if meta["algorithm_version"]!=ALGORITHM_VERSION or meta["algorithm_sha256"]!=algorithm_hash():
        raise ValueError("Prepared equilibrium source algorithm mismatch")
    if meta["equilibrium_config_fingerprint"]!=equilibrium_config_fingerprint(solver.config):
        raise ValueError("Prepared equilibrium physical/FE/quadrature/mesh config mismatch")
    if meta["provenance"]["mesh_fingerprint"]!=run_provenance(solver)["mesh_fingerprint"]:
        raise ValueError("Prepared equilibrium mesh/partition fingerprint mismatch")
    if hashlib.sha256(prepared.phi_coefficients.tobytes()).hexdigest()!=meta["coefficient_sha256"]:
        raise ValueError("Prepared equilibrium coefficient integrity failure")
    data={k:v for k,v in meta.items() if k!="fingerprint"}
    if hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()!=meta["fingerprint"]:
        raise ValueError("Prepared equilibrium metadata integrity failure")
    if (meta["mass_residual_relative_to_domain"]>meta["policy"]["mass_domain_tolerance"] or
            meta["stationarity"]["normalized"]>meta["policy"]["normalized_riesz_tolerance"]):
        raise ValueError("Uncertified stationary equilibrium cache")


def save_prepared(folder, prepared):
    import h5py
    folder=Path(folder)
    folder.mkdir(parents=True,exist_ok=True)
    with h5py.File(folder/"prepared.h5","x") as h:
        h.attrs["schema"]=CHECKPOINT_SCHEMA
        h.attrs["metadata_json"]=json.dumps(prepared.metadata,sort_keys=True)
        h.create_dataset("phi_coefficients",data=prepared.phi_coefficients)
    with (folder/"prepared.json").open("x") as stream:
        json.dump(prepared.metadata,stream,indent=2,allow_nan=False)


def load_prepared(folder, solver, target_mass):
    import h5py
    folder=Path(folder)
    manifest=json.loads((folder/"prepared.json").read_text())
    with h5py.File(folder/"prepared.h5","r") as h:
        if h.attrs["schema"]!=CHECKPOINT_SCHEMA or json.loads(h.attrs["metadata_json"])!=manifest:
            raise ValueError("Incomplete prepared equilibrium checkpoint transaction")
        prepared=PreparedEquilibrium(h["phi_coefficients"][:],manifest)
    validate_prepared(solver,prepared)
    if target_mass!=manifest["target_mass"]:
        raise ValueError("Prepared equilibrium target mass mismatch")
    return prepared


def zero_mass_perturbations(solver, phi):
    """Deterministic FE directions; all mean removal precedes any solution run."""
    from dolfinx import fem
    dx,_=_measures(solver)
    area=_integral(solver,1*dx)
    c=solver.config
    def bump(x,xc,zc,rx,rz):
        r2=((x[0]-xc)/rx)**2+((x[1]-zc)/rz)**2
        inside=r2<1
        result=np.zeros_like(r2)
        result[inside]=np.exp(1-1/(1-r2[inside]))
        return result
    directions={}
    # Bulk bump centers at z=.32,.43 m for the matched-density benchmark.
    for name,centers in {"smooth_bulk":((-.12,.32),(.12,.43)),
                         "near_interface":((-.10,.13),(.10,.13)),
                         "near_wall":((-.18,.03),(.18,.03))}.items():
        one,two=fem.Function(phi.function_space),fem.Function(phi.function_space)
        height=c.z_max-c.z_min;width=c.x_max-c.x_min
        rx,rz=.075*width,(.10 if name!="near_wall" else .065)*height
        for field,(xc,zc) in zip((one,two),centers):
            field.interpolate(lambda x,xc=xc,zc=zc:bump(x,xc/1.2*width, c.z_min+zc/.6*height,rx,rz))
        m1,m2=_integral(solver,one*dx),_integral(solver,two*dx)
        if abs(m2)<1e-16 or abs(m1)<1e-16:
            raise ValueError("Mesh too coarse to represent deterministic perturbation bumps")
        one.x.array[:]-=(m1/m2)*two.x.array
        one.x.array[:]/=np.max(np.abs(one.x.array))
        mean=_integral(solver,one*dx)/area
        one.x.array[:]-=mean
        one.x.scatter_forward()
        directions[name]=(one,{"subtracted_discrete_mean":mean,"mass":_integral(solver,one*dx),
            "construction":"two C-infinity compact bumps with opposite discrete integrals, normalized then zero mean",
            "centers_reference_domain":centers,"radii_m":[rx,rz]})
    return directions


def first_variations(solver, phi, step=1e-3):
    """Independent five-point ENERGY differences, exact for a quartic in delta.

    This never accesses the Newton residual or its Jacobian. The symmetric
    second difference is only a directional local-minimum diagnostic.
    """
    from dolfinx import fem
    dx,_=_measures(solver)
    base=_integral(solver,energy_form(solver,phi))
    area=_integral(solver,1*dx)
    result=[]
    for name,(psi,metadata) in zero_mass_perturbations(solver,phi).items():
        energies={}
        trial=fem.Function(phi.function_space)
        for k in (-2,-1,1,2):
            trial.x.array[:]=phi.x.array+k*step*psi.x.array
            trial.x.scatter_forward()
            energies[k]=_integral(solver,energy_form(solver,trial))
        derivative=(energies[-2]-8*energies[-1]+8*energies[1]-energies[2])/(12*step)
        norm=np.sqrt(_integral(solver,psi*psi*dx))
        normalization=solver.config.sigma/solver.config.length_scale_m*np.sqrt(area)*norm
        result.append({"name":name,**metadata,"energy_difference_step":step,"F0":base,
            "energies":energies,"directional_derivative":derivative,
            "normalized_abs_derivative":abs(derivative)/normalization,
            "second_difference":energies[1]+energies[-1]-2*base,
            "second_variation_diagnostic_only":True})
    return result
