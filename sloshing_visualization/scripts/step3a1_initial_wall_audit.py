"""Analytic initial-data audit, NOT a new PDE or a modified initial condition."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss


def main():
    # These are the two recorded benchmark geometries; no fitting of parameters.
    R,epsilon,sigma,half_width=.3,.025,.1,.6
    nodes,weights=leggauss(512)
    x=half_width*nodes
    sample=np.linspace(-half_width,half_width,20000)  # no isolated circle-center cusp sample
    rows=[]
    for initial_angle in (60.,90.):
        cosine=0. if initial_angle==90 else np.cos(np.deg2rad(initial_angle))
        target=np.cos(np.deg2rad(60.))
        def residual(points):
            r=np.sqrt(points**2+(R*cosine)**2)
            phi=np.tanh((R-r)/(np.sqrt(2)*epsilon))
            # From direct differentiation of the declared signed-distance tanh.
            return 3*sigma/4*(R*cosine/r-target)*(1-phi**2)
        values=residual(x)
        rows.append({"theta_initial_deg":initial_angle,"theta_target_deg":60.,
            "wall_variation_L2_N_per_sqrt_m":float(np.sqrt(half_width*np.dot(weights,values**2))),
            "wall_variation_max_sampled_N_per_m":float(np.max(np.abs(residual(sample))))})
    result={"scope":"analytic natural-wall residual of the initial geometries, before FEM interpolation",
        "parameters":{"radius_m":R,"epsilon_m":epsilon,"sigma_N_per_m":sigma,"half_width_m":half_width},
        "formula":"L = (3 sigma / 4) [R cos(theta_initial)/r - cos(theta_target)] (1-phi^2), r=sqrt(x^2+R^2 cos^2(theta_initial))",
        "interpretation":"theta_initial=theta_target makes L=0 at r=R (phi=0) only, not on the entire diffuse wall trace. Geometrically compatible is not a prepared variational equilibrium.",
        "measurements":rows}
    output=Path("validation_results/step3a1/initial_wall_audit.json")
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
