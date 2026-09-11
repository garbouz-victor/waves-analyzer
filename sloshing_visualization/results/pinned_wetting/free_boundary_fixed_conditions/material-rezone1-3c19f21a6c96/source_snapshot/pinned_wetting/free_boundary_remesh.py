"""Same-domain rezoning with conservative mixed L2 field transfer (PW2).

No free-boundary particle is moved by this numerical operation. Wall mesh
nodes may slide parametrically on the same solid; their physical velocity is
still constrained. Every projection input/output is native audit evidence.
"""
from dataclasses import replace
import numpy as np
from scipy.sparse import bmat, csr_matrix, vstack
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from skfem import Basis, ElementTriP2, BilinearForm, LinearForm, asm
from skfem.helpers import dot,grad
from .free_boundary import Operators, quadratic_minima, interface_nodes


@BilinearForm
def laplace(u,v,w): return dot(grad(u),grad(v))


@LinearForm
def transfer_rhs(v,w): return dot(v,w.old_velocity)


@LinearForm
def momentum_x(v,w): return v[0]


@LinearForm
def momentum_z(v,w): return v[1]


def shape(q):
    x,y=q
    phi=np.array([1-3*x-3*y+2*x*x+4*x*y+2*y*y,2*x*x-x,2*y*y-y,
                  4*x-4*x*x-4*x*y,4*x*y,4*y-4*x*y-4*y*y])
    dx=np.array([-3+4*x+4*y,4*x-1,0*x,4-8*x-4*y,4*y,-4*y])
    dy=np.array([-3+4*x+4*y,0*y,4*y-1,-4*x,4*x,4-4*x-8*y])
    return phi,dx,dy


def locate_and_evaluate(old,velocity,points):
    """Bounded local Newton inverse, no clipped/extrapolated outside points."""
    points=np.asarray(points); flat=points.reshape(2,-1)
    mesh=old.mesh; elements=old.scalar.element_dofs
    tree=cKDTree(mesh.p[:,mesh.t].mean(axis=1).T)
    count=flat.shape[1];cells=np.full(count,-1,dtype=int);refs=np.zeros((2,count))
    map_errors=np.zeros(count)
    for k in (8,32,128,512):
        unresolved=np.flatnonzero(cells<0)
        if not len(unresolved):break
        k=min(k,mesh.t.shape[1])
        # Bound scratch memory irrespective of total fine mesh quadrature size.
        for offset in range(0,len(unresolved),256):
            ids=unresolved[offset:offset+256]
            candidate=tree.query(flat[:,ids].T,k=k)[1].reshape(len(ids),k)
            elem=candidate.ravel();target=np.repeat(flat[:,ids],k,axis=1)
            X=mesh.p[:,elements[:,elem]]
            e1=X[:,1]-X[:,0];e2=X[:,2]-X[:,0];delta=target-X[:,0]
            determinant=e1[0]*e2[1]-e1[1]*e2[0]
            q=np.array([(delta[0]*e2[1]-delta[1]*e2[0])/determinant,
                        (e1[0]*delta[1]-e1[1]*delta[0])/determinant])
            with np.errstate(over="ignore",invalid="ignore",divide="ignore"):
                for _ in range(15):
                    phi,dx,dy=shape(q)
                    residual=target-np.einsum("ijn,jn->in",X,phi)
                    jx=np.einsum("ijn,jn->in",X,dx);jy=np.einsum("ijn,jn->in",X,dy)
                    det=jx[0]*jy[1]-jx[1]*jy[0]
                    deltaq=np.array([(residual[0]*jy[1]-residual[1]*jy[0])/det,
                                    (jx[0]*residual[1]-jx[1]*residual[0])/det])
                    q+=np.clip(deltaq,-2.,2.)
                phi,_,_=shape(q)
                error=np.linalg.norm(target-np.einsum("ijn,jn->in",X,phi),axis=0)
                valid=(q.min(axis=0)>=-1e-10)&(q.sum(axis=0)<=1+1e-10)&(error<1e-10)
            valid=valid.reshape(len(ids),k)
            found=valid.any(axis=1);choice=valid.argmax(axis=1)
            selected=np.arange(len(ids))*k+choice
            cells[ids[found]]=elem[selected[found]]
            refs[:,ids[found]]=q[:,selected[found]]
            map_errors[ids[found]]=error[selected[found]]
    if np.any(cells<0):raise ValueError(f"conservative remap has {np.sum(cells<0)} unlocated quadrature points")
    phi,_,_=shape(refs)
    vnodes=old.nodal_vector(velocity)
    values=np.einsum("ijn,jn->in",vnodes[:,elements[:,cells]],phi).reshape(points.shape)
    return values,cells,refs,float(map_errors.max())


def rezone(reference,mesh,velocity,c,orientation):
    scalar=Basis(reference,ElementTriP2(),intorder=c.intorder)
    K=asm(laplace,scalar).tocsr()
    surface=scalar.get_dofs("surface").all()
    bottom=scalar.get_dofs("bottom").all()
    walls=np.union1d(scalar.get_dofs("left").all(),scalar.get_dofs("right").all())
    displacement=np.zeros_like(mesh.p)
    for component,fixed in enumerate((np.unique(np.r_[surface,bottom,walls]),np.union1d(surface,bottom))):
        free=np.setdiff1d(np.arange(scalar.N),fixed)
        displacement[component,fixed]=mesh.p[component,fixed]-reference.p[component,fixed]
        displacement[component,free]=spsolve(K[free][:,free],-K[free][:,fixed]@displacement[component,fixed])
    new=replace(mesh,doflocs=reference.p+displacement)
    # Exact boundary coordinates, not an endpoint correction: rezoning keeps
    # every surface geometry DOF; walls remain the same straight segments.
    new.p[:,surface]=mesh.p[:,surface]
    mins,sign=quadratic_minima(new)
    if mins.min()<=0 or not np.array_equal(sign,orientation):
        raise ValueError(f"harmonic rezone invalid: cell {int(mins.argmin())}, J={float(mins.min())}")
    old=Operators(mesh,c);op=Operators(new,c)
    values,cells,refs,error=locate_and_evaluate(old,velocity,np.asarray(op.scalar.global_coordinates()))
    rhs=asm(transfer_rhs,op.vbasis,old_velocity=values)
    P=np.vstack([asm(momentum_x,op.vbasis),asm(momentum_z,op.vbasis)])
    Pold=np.vstack([asm(momentum_x,old.vbasis),asm(momentum_z,old.vbasis)])
    target_momentum=Pold@velocity
    free=op.free;C=vstack((op.B[:,free],csr_matrix(P[:,free])),format="csr")
    A=bmat([[op.M[free][:,free],C.T],[C,None]],format="csc")
    target=np.r_[rhs[free],np.zeros(op.B.shape[0]),target_momentum]
    solution=spsolve(A,target)
    v=np.zeros_like(velocity);v[free]=solution[:len(free)]
    if not np.all(np.isfinite(v)):raise ValueError("nonfinite conservative rezone projection")
    oldrow,newrow=old.measure(velocity),op.measure(v)
    volume_error=newrow["volume_m2"]-oldrow["volume_m2"]
    momentum_error=float(np.max(abs(P@v-target_momentum)))
    projection_residual=float(np.linalg.norm(A@solution-target,np.inf))
    if abs(volume_error)>1e-10 or momentum_error>1e-10 or projection_residual>1e-9:
        raise ValueError("conservative rezone constraints failed")
    evidence={"geometry":new.p.copy(),"velocity":v,"projection_multipliers":solution[len(free):],
        "source_cells":cells,"source_reference_coordinates":refs,
        "energy_work":newrow["energy"]-oldrow["energy"],"volume_change":volume_error,
        "momentum_error":momentum_error,"projection_residual":projection_residual,"inverse_map_error":error,
        "min_jacobian":float(mins.min()),"surface_geometry_change":float(np.max(abs(new.p[:,surface]-mesh.p[:,surface])))}
    return new,v,evidence
