"""Optional PW2 same-boundary mesh optimization; NOT a material constitutive law.

This module is not used by a run until explicitly selected and included in its
numerical source identity.  Its dimensionless pseudo-elastic weights affect
only mesh coordinates at a rezoning instant.  Every free-surface P2 coordinate
is held bitwise fixed; the subsequent velocity transfer is a separate step.
"""
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
from skfem import Basis, ElementTriP2

from .free_boundary import quadratic_minima


def minimum_witnesses(mesh, orientation):
    """Exact signed quadratic-J minima AND their reference-coordinate witnesses."""
    q=np.array([[0.,1.,0.,.5,0.,.5],[0.,0.,1.,0.,.5,.5]])
    mapping=mesh.mapping()
    values=(mapping.J(0,0,q)*mapping.J(1,1,q)-mapping.J(0,1,q)*mapping.J(1,0,q))*orientation[:,None]
    f,f10,f01,fh0,f0h,fhh=values.T
    aa=2*(f10+f-2*fh0);cc=2*(f01+f-2*f0h)
    dd=f10-f-aa;ee=f01-f-cc
    bb=4*(fhh-f-.25*(aa+cc)-.5*(dd+ee))
    minima=f.copy();points=np.zeros((2,len(f)))
    def take(x,y,valid):
        value=aa*x*x+bb*x*y+cc*y*y+dd*x+ee*y+f
        better=valid&(value<minima)
        minima[better]=value[better]
        points[0,better]=np.broadcast_to(x,minima.shape)[better]
        points[1,better]=np.broadcast_to(y,minima.shape)[better]
    take(1.,0.,np.ones(len(f),dtype=bool));take(0.,1.,np.ones(len(f),dtype=bool))
    with np.errstate(divide="ignore",invalid="ignore"):
        x=-dd/(2*aa);take(x,0.,(x>0)&(x<1))
        y=-ee/(2*cc);take(0.,y,(y>0)&(y<1))
        x=-(bb-2*cc+dd-ee)/(2*(aa-bb+cc));take(x,1-x,(x>0)&(x<1))
        determinant=4*aa*cc-bb*bb
        x=(bb*ee-2*cc*dd)/determinant;y=(bb*dd-2*aa*ee)/determinant
        take(x,y,(x>0)&(y>0)&(x+y<1))
    return minima,points


def feasible_lbfgs(fun, initial, max_iterations):
    """Bounded L-BFGS with explicit feasible Armijo backtracking.

    Never infer success from an unchanged objective after an infeasible trial.
    Mesh selection still uses actual valid minimum-J improvement, not this flag.
    """
    x=initial.copy();value,gradient=fun(x);history=[]
    if not np.isfinite(value) or not np.isfinite(gradient).all():
        return SimpleNamespace(x=x,fun=value,nit=0,success=False,message="nonfinite initial objective or gradient")
    if max_iterations<=0:
        return SimpleNamespace(x=x,fun=value,nit=0,success=False,message="bounded iteration limit")
    success=False;message="bounded iteration limit";iteration=0
    for iteration in range(max_iterations):
        if np.linalg.norm(gradient,np.inf)<=1e-7:
            success=True;message="scaled gradient tolerance";break
        q=gradient.copy();alphas=[]
        for s,y,rho in reversed(history):
            alpha=rho*float(s@q);alphas.append(alpha);q-=alpha*y
        scale=float(history[-1][0]@history[-1][1])/float(history[-1][1]@history[-1][1]) if history else 1.
        direction=scale*q
        for (s,y,rho),alpha in zip(history,reversed(alphas)):
            direction+=s*(alpha-rho*float(y@direction))
        direction=-direction
        descent=float(gradient@direction)
        if not np.isfinite(descent) or descent>=0:
            history=[];direction=-gradient;descent=-float(gradient@gradient)
        step=min(1.,1./max(1.,float(np.linalg.norm(direction,np.inf))))
        accepted=False
        for _ in range(60):
            candidate=x+step*direction
            next_value,next_gradient=fun(candidate)
            if (np.isfinite(next_value) and np.isfinite(next_gradient).all()
                    and next_value<=value+1e-4*step*descent and next_value<value):
                accepted=True;break
            step*=.5
        if not accepted:
            message="bounded feasible line search exhausted";break
        s=candidate-x;y=next_gradient-gradient;curvature=float(s@y)
        if curvature>1e-10*np.linalg.norm(s)*np.linalg.norm(y):
            history.append((s,y,1./curvature));history=history[-10:]
        x,value,gradient=candidate,next_value,next_gradient
    return SimpleNamespace(x=x,fun=value,nit=iteration+1,success=success,message=message)


def density_and_stress(F):
    """Dimensionless neo-Hookean mesh density and its exact derivative."""
    j = F[0, 0]*F[1, 1]-F[0, 1]*F[1, 0]
    if not np.all(np.isfinite(j)) or np.any(j <= 0):
        raise ValueError("mesh objective requires positive relative Jacobian")
    inverse_transpose = np.array([[F[1, 1], -F[1, 0]],
                                  [-F[0, 1], F[0, 0]]])/j
    logj = np.log(j)
    density = .5*(np.sum(F*F, axis=(0, 1))-2)-logj+.5*logj**2
    stress = F+(logj-1)*inverse_transpose
    return density, stress, j, inverse_transpose


class MeshObjective:
    """Uniform-cell quality weights plus a vertex Jacobian barrier.

    Initial reference triangles must be affine.  Volume quadrature alone can
    miss vertex collapse; each objective evaluation also checks the EXACT
    quadratic Jacobian minimum on the entire closed reference triangle.
    The explicit vertex barrier repels the observed surface-attached collapse.
    """
    def __init__(self, reference, current, orientation):
        self.reference, self.current = reference, current
        self.orientation = np.asarray(orientation)
        basis = Basis(reference, ElementTriP2(), intorder=4)
        self.reference_min, reference_sign = quadratic_minima(reference)
        current_min, current_sign = quadratic_minima(current)
        if (not np.all(np.isfinite(self.reference_min)) or self.reference_min.min() <= 0
                or not np.array_equal(reference_sign, self.orientation)):
            raise ValueError("mesh optimization requires a valid reference orientation")
        if (not np.all(np.isfinite(current_min)) or current_min.min() <= 0
                or not np.array_equal(current_sign, self.orientation)):
            raise ValueError("mesh optimization requires a valid current geometry")
        mids = basis.dofs.facet_dofs[0]
        expected_mids = reference.p[:, reference.facets].mean(axis=1)
        tolerance = 32*np.finfo(float).eps*max(1., float(np.max(abs(reference.p))))
        if not np.allclose(reference.p[:, mids], expected_mids, rtol=0., atol=tolerance):
            raise ValueError("mesh optimization requires an affine initial reference")
        self.nodes = basis.element_dofs
        self.grad = np.array([field[0].grad for field in basis.basis])
        self.weights = np.broadcast_to(basis.W/basis.W.sum(), basis.dx.shape)
        vertex_basis = Basis(reference, ElementTriP2(),
                             quadrature=(np.array([[0., 1., 0.], [0., 0., 1.]]),
                                         np.full(3, 1/3)))
        self.vertex_grad = np.array([field[0].grad for field in vertex_basis.basis])
        self.vertex_weight = .25/3
        self.minimum_weight = .25
        vertices=reference.p[:,reference.t]
        first,second=vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0]
        determinant=first[0]*second[1]-first[1]*second[0]
        self.reference_inverse_transpose=np.array([[second[1],-first[1]],[-second[0],first[0]]])/determinant
        # A diagonal coordinate preconditioner, not a displacement cutoff.
        local_diag = np.einsum("iaeq,eq->ie", self.grad*self.grad, self.weights)
        diagonal = np.zeros(basis.N)
        np.add.at(diagonal, self.nodes.ravel(), local_diag.ravel())
        self.scale = np.broadcast_to(1/np.sqrt(diagonal), current.p.shape).copy()
        surface = basis.get_dofs("surface").all()
        bottom = basis.get_dofs("bottom").all()
        walls = np.union1d(basis.get_dofs("left").all(), basis.get_dofs("right").all())
        allowed = np.ones(current.p.shape, dtype=bool)
        allowed[:, np.union1d(surface, bottom)] = False
        allowed[0, walls] = False
        self.allowed = allowed
        self.best = current.p.copy()
        self.best_quality = float(np.min(current_min/self.reference_min))
        self.evaluations = 0
        self.invalid_evaluations = 0

    def coordinates(self, variables):
        points = self.current.p.copy()
        points[self.allowed] += variables*self.scale[self.allowed]
        return points

    def value_gradient(self, variables):
        points = self.coordinates(variables)
        mesh = replace(self.current, doflocs=points)
        minima, witnesses = minimum_witnesses(mesh,self.orientation)
        self.evaluations += 1
        if not np.all(np.isfinite(minima)) or minima.min() <= 0:
            self.invalid_evaluations += 1
            return np.inf, np.zeros_like(variables)
        X = points[:, self.nodes]
        F = np.einsum("ine,njeq->ijeq", X, self.grad)
        value, stress, _, _ = density_and_stress(F)
        total = float(np.sum(value*self.weights))
        local_gradient = np.einsum("ijeq,njeq,eq->ine", stress, self.grad, self.weights)
        vertex_F = np.einsum("ine,njeq->ijeq", X, self.vertex_grad)
        _, _, j, inverse_transpose = density_and_stress(vertex_F)
        total += self.vertex_weight*float(np.sum(j-1-np.log(j)))
        barrier_stress = self.vertex_weight*(j-1)*inverse_transpose
        local_gradient += np.einsum("ijeq,njeq->ine", barrier_stress, self.vertex_grad)
        # The exact-minimum witness can lie inside an edge (as in actual cell
        # 3000). Envelope derivative at a unique active witness; ties can be
        # nonsmooth, hence bounded feasible descent and no convergence promise.
        x,y=witnesses
        dx=np.array([-3+4*x+4*y,4*x-1,0*x,4-8*x-4*y,4*y,-4*y])
        dy=np.array([-3+4*x+4*y,0*y,4*y-1,-4*x,4*x,4-4*x-8*y])
        witness_grad=np.einsum("ije,nje->nie",self.reference_inverse_transpose,np.array([dx,dy]).transpose(1,0,2))
        witness_F=np.einsum("ine,nje->ije",X,witness_grad)
        _,_,j,inverse_transpose=density_and_stress(witness_F)
        total+=self.minimum_weight*float(np.sum(j-1-np.log(j)))
        local_gradient+=np.einsum("ije,nje->ine",self.minimum_weight*(j-1)*inverse_transpose,witness_grad)
        gradient = np.zeros_like(points)
        for component in (0, 1):
            np.add.at(gradient[component], self.nodes.ravel(), local_gradient[component].ravel())
        quality = float(np.min(minima/self.reference_min))
        if quality > self.best_quality:
            self.best_quality, self.best = quality, points.copy()
        return total, (gradient*self.scale)[self.allowed]


def optimize_mesh(reference, current, orientation, max_iterations=80):
    """Return a valid same-boundary candidate, never force optimizer success."""
    objective = MeshObjective(reference, current, orientation)
    initial = np.zeros(int(objective.allowed.sum()))
    initial_quality = objective.best_quality
    initial_value, _ = objective.value_gradient(initial)
    result = feasible_lbfgs(objective.value_gradient, initial, max_iterations)
    new = replace(current, doflocs=objective.best)
    info = {"method":"dimensionless_mesh_neo_hookean_with_vertex_and_exact_minimum_barrier",
            "physical_constitutive_law":False, "cell_weights":"uniform",
            "mesh_mu":1., "mesh_kappa":1., "vertex_barrier_weight":.25,
            "minimum_witness_barrier_weight":.25,"line_search":"bounded_feasible_Armijo_LBFGS",
            "initial_objective":initial_value, "optimizer_last_objective":float(result.fun),
            "iterations":int(result.nit), "evaluations":objective.evaluations,
            "invalid_trials_rejected":objective.invalid_evaluations,
            "optimizer_success":bool(result.success), "optimizer_message":str(result.message),
            "initial_min_relative_jacobian":initial_quality,
            "selected_min_relative_jacobian":objective.best_quality,
            "full_free_boundary_displacement_m":0.}
    return new, info
