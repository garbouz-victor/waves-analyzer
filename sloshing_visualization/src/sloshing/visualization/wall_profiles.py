"""Actual P2 wall-normal samples; no fitting, artificial zero strip or field edits."""

import numpy as np
from ..validation.qualification.fem_evaluation import evaluate_map


def profile_points(geometry, side, z_levels, epsilon):
    if side not in ('left','right'):raise ValueError('side must be left or right')
    eps=np.asarray(epsilon,dtype=float);depth=np.asarray(z_levels,dtype=float)
    if eps.ndim!=1 or not len(eps) or eps[0]!=0 or np.any(np.diff(eps)<=0):
        raise ValueError('epsilon must start at zero and strictly increase')
    edge=geometry.x[0] if side=='left' else geometry.x[-1]
    x=edge+eps if side=='left' else edge-eps
    xx,zz=np.meshgrid(x,depth)
    return np.array([xx.ravel(),zz.ravel()]),xx.shape


def wall_normal_profiles(geometry,u,w,side,z_levels,epsilon,mapping=None):
    points,shape=profile_points(geometry,side,z_levels,epsilon)
    mapping=geometry.interpolation(points) if mapping is None else mapping
    up,wp=(evaluate_map(c,mapping).reshape(shape) for c in (u,w))
    return {'u':up,'w':wp,'speed':np.hypot(up,wp)}


def trace_speed_max(u,w):
    """Exact max of the norm of a vector P2 trace, including interior extrema.

    Input trace nodes alternate endpoint/midpoint/endpoint. No-slip zero at
    these unisolvent trace nodes proves zero on the entire wall, not just probes.
    """
    values=np.column_stack((u,w))
    maximum=float(np.linalg.norm(values,axis=1).max())
    if maximum==0:return 0.
    for c,m,e in zip(values[:-1:2],values[1::2],values[2::2]):
        a=2*c-4*m+2*e;b=-3*c+4*m-e
        # d |a*s²+b*s+c|² / ds (cubic).
        roots=np.roots([4*np.dot(a,a),6*np.dot(a,b),2*np.dot(b,b)+4*np.dot(a,c),2*np.dot(b,c)])
        for r in roots:
            if abs(r.imag)<1e-10 and 0<r.real<1:
                s=float(r.real);maximum=max(maximum,float(np.linalg.norm(a*s*s+b*s+c)))
    return maximum


def arrow_geometry(x,z,u,w,seconds,scope):
    """Conservative complete-glyph envelope, plus actual tail/tip ranges.

    Quiver head half-width is width*headwidth/2 of the 2 m x-axis span.
    With tail pivot all longitudinal vertices lie between tail and tip;
    the transverse half-width bound also covers short-arrow geometry.
    """
    dx, dz=seconds*np.asarray(u),seconds*np.asarray(w)
    head=.0024*3.4/2*(2*scope['domain_abs_x_max_m'])
    bound=float(abs(x).max()+abs(dx).max()+head)
    wall_margin=scope['domain_abs_x_max_m']-bound
    contact_margin=scope['transition_abs_x_max_m']-bound
    if contact_margin<scope['arrow_margin_extra_m']:
        raise RuntimeError('Arrow glyph envelope approaches contact strip; move the fixed seed grid inward')
    return {'arrow_seed_x_min':float(x.min()),'arrow_seed_x_max':float(x.max()),'arrow_seed_nx':x.shape[1],
            'arrow_seed_z_min':float(z.min()),'arrow_seed_z_max':float(z.max()),'fixed_arrow_seconds':float(seconds),
            'maximum_rendered_arrow_dx':float(abs(dx).max()),'maximum_rendered_arrow_length':float(np.hypot(dx,dz).max()),
            'wall_visual_margin':float(wall_margin),'contact_strip_visual_margin_m':float(contact_margin),
            'glyph_head_half_width_bound_m':float(head),'tip_x_range':[float((x+dx).min()),float((x+dx).max())]}
