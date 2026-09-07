"""Exact polynomial activation; certified and diagnostic normal spacing.

Reference nodal order: (0,0),(1,0),(0,1),(.5,0),(.5,.5),(0,.5).
Bernstein bounds remain available for explicit historical comparisons only.
"""
import numpy as np
from .free_energy import transition_width
from .validation_policy import TRANSITION_PHI_LIMIT
from .p2_certification import polynomial_ranges, certified_normal_widths

REFERENCE_POINTS=np.array([[0.,0.],[1.,0.],[0.,1.],[.5,0.],[.5,.5],[0.,.5],[1/3,1/3]])


def polynomial_bounds(nodal_values, degree):
    values=np.asarray(nodal_values,dtype=float)
    if degree not in (1,2) or values.shape[-1]!=(3 if degree==1 else 6):
        raise ValueError("Expected triangle CG1 three or CG2 six reference nodal values")
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite phase values in interface activation")
    if degree==1:
        coefficients=values
    else:
        vertices=values[...,:3]
        edge=2*values[...,3:6]-.5*vertices[...,[0,1,2]]-.5*vertices[...,[1,2,0]]
        coefficients=np.concatenate((vertices,edge),axis=-1)
    return coefficients.min(axis=-1),coefficients.max(axis=-1)


def active_transition_cells(nodal_values, degree, limit=TRANSITION_PHI_LIMIT):
    result=polynomial_ranges(nodal_values,degree)
    return (result["min"]<=limit)&(result["max"]>=-limit)


def normal_widths(triangles,gradients):
    triangles=np.asarray(triangles,dtype=float)
    gradients=np.asarray(gradients,dtype=float)
    if not np.isfinite(gradients).all():
        raise ValueError("Nonfinite phase gradient in resolution measurement")
    norms=np.linalg.norm(gradients,axis=-1)
    valid=norms>1e-12
    directions=gradients/np.maximum(norms[...,None],1e-30)
    projections=np.einsum("cvi,csi->csv",triangles,directions)
    width=np.max(np.where(valid,np.ptp(projections,axis=-1),0.),axis=1)
    # A constant value within the band has no defined normal. Do not drop it.
    diameter=np.linalg.norm(triangles-np.roll(triangles,1,axis=1),axis=-1).max(axis=1)
    return np.where(valid.any(axis=1),width,diameter)


def local_resolution(triangles,values,gradients,degree,epsilon):
    active=active_transition_cells(values,degree)
    indices=np.flatnonzero(active)
    widths=certified_normal_widths(triangles[active],gradients[active,:3])["width"]
    counts=transition_width(epsilon)/widths
    return indices,widths,counts
