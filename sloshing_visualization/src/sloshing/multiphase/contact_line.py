"""Wall crossings and apparent angles from interface fits away from wall nodes."""
import numpy as np
from .interface import fit_circle


def all_wall_crossings(polylines,wall_coordinate,axis=0,tolerance=1e-9):
    result=[]
    for component,line in enumerate(polylines):
        values=[]
        for a,b in zip(line,line[1:]):
            da,db=a[axis]-wall_coordinate,b[axis]-wall_coordinate
            if abs(da)<tolerance:
                values.append(float(a[1-axis]))
            if da*db<0:
                values.append(float((a+(-da/(db-da))*(b-a))[1-axis]))
            if abs(db)<tolerance:
                values.append(float(b[1-axis]))
        unique=[]
        for value in sorted(values):
            if not unique or abs(value-unique[-1])>tolerance:
                unique.append(value)
        result.extend({"component":component,"coordinate":value} for value in unique)
    return result


def sessile_apparent_angle(points,wall_z,epsilon,window=(2.,10.)):
    """Circle-fit static droplet angle through liquid; fit-window must be reported.

    This specialized benchmark estimator is NOT a general film/junction classifier.
    """
    points=np.asarray(points).reshape(-1,2)
    distance=points[:,1]-wall_z
    mask=(distance>=window[0]*epsilon)&(distance<=window[1]*epsilon)
    fit=fit_circle(points[mask])
    cosine=(wall_z-fit["center_z"])/fit["radius"]
    if abs(cosine)>1:
        raise ValueError("Fitted circle does not intersect wall")
    return dict(fit,theta_deg=float(np.rad2deg(np.arccos(cosine))),window=list(window))
