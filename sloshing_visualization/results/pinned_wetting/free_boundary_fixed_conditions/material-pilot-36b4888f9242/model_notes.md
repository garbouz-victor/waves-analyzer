# PW2: full material free-boundary Navier–Stokes (2026-09-11)

This is a NEW equation-level extension, explicitly authorized by the user.
PW1 fixed-domain linear results, their failed slope screen, and all historical
STOP/COMPLETE bundles remain unchanged. Model ID:
`PW2_FULL_NS_LEFT_NAVIER_RIGHT_NOSLIP_MATERIAL_P2_V1`.

## Physical problem (fixed before the first solve)

In the actual liquid domain Ω(t), per unit density,
`∂t v + (v·∇)v = div T - g e_z`, `div v=0`,
`T=-p I+2νD(v)`. Parameters: a=1 m, d=10 m, g=9.81 m/s²,
ν=.01 m²/s, σ=0. Initial velocity is zero and the entire initial interface
is the exact straight segment z=x tan(2°).

LEFT x=-a: u=0 and `(T n)_z=-(ν/b)w`, b=.50 m only on LEFT.
Since n=-e_x, this is `ν(∂x w+∂z u)=ν w/b`.
RIGHT x=a and BOTTOM z=-d: u=w=0. At the moving free boundary,
`T n=0` with its actual normal, and material kinematics `dX/dt=v(X,t)`.
There is no right Robin term, no capillary force, no prescribed contact angle,
no detached replacement of the right branch. P2 is a fixed material particle
because both components of its velocity DOF are Dirichlet zero.

The release pressure is NOT prescribed tilted hydrostatics. The initial
pressure/acceleration are obtained from the coupled weak momentum/divergence
system with the initial geometry and actual zero-traction boundary. Initial
triangles are straight-sided; quadratic geometry edge nodes are their exact
midpoints. This also permits the exact flat hydrostatic rest benchmark.
No extra C1 compatibility condition is imposed at the pinned corner.

## Selected numerical method

Local installed scikit-fem 11 / SciPy sparse direct solver are available. Use
Taylor–Hood P2 velocity / P1 pressure on isoparametric quadratic triangles.
Every P2 geometry node is material. For one accepted step,
`X1-X0=dt*(v1+v0)/2`. Assemble mass, stress, divergence, LEFT wall friction
and gravity on `Xm=(X0+X1)/2`, solving the resulting nonlinear midpoint
system by converged geometry iteration. All pressure DOFs are retained:
free traction fixes the pressure, so no pressure constraint is discarded.

This is NOT unsteady Stokes on a fixed rectangle. The material derivative
of v∘X contains the full Eulerian convection. Relative ALE convection is
zero specifically because mesh velocity equals the full P2 fluid velocity
everywhere, not because convection has been dropped. If mesh redistribution
is introduced later, relative convection and conservative transfer must be
implemented and independently verified; it cannot reuse this zero term.

The 2D midpoint determinant identity gives exactly
`J1-J0=dt*Jm*div(vm)` pointwise. Constant pressure testing then conserves
total area, although weak divergence does not imply pointwise incompressibility.
Geometry validation uses the minimum of the quadratic Jacobian polynomial
over each reference triangle, including boundary/interior stationary points.
The full ordered quadratic interface is saved; it is never reduced to eta(x).

## Energy and independently measurable numerical work

`E=1/2∫Ω|v|² + g(∫Ω z - ∫flat z)`,
`E0=g*a³*tan²(2°)/3`. Physical losses are separately
`Db=2ν∫Ω D:D ≥0`, `Dl=(ν/b)∫LEFT w² ≥0`.
Continuous budget is `E-E0+∫(Db+Dl)dt`, normalized by E0, not background
hydrostatic energy. Signed numerical work is NOT physical heat.

For reference quadrature define δv=v1-v0, vm=(v1+v0)/2,
Ja=(J1+J0)/2, δJ=J1-J0, zm=(z1+z0)/2 and δz=z1-z0.
Measure the following explicit geometric/time defects, not a fitted residual:

`Wk=∫[(Ja-Jm) vm·δv + δJ |vm|²/2 + δJ |δv|²/8]`,
`Wp=g∫[(Ja-Jm)δz + δJ zm]`.

They follow algebraically from kinetic/potential differences and midpoint
momentum work. Check both the raw continuous budget (≤5% E0) and the split
discrete balance with independently reconstructed algebraic/pressure work.
Surface/domain quadrature is order 8 initially (degree sufficient for the
polynomial geometry-energy identities); curved inverse-map stress remains
quadrature-approximated and is part of refinement, not exact continuum proof.

## Irreversible coating: separate leading-order law

After EACH accepted step, `H_L=max(previous H_L, actual X_left,z)`;
`H_R=z2`. Store the complete lower-connected interval [-d,H] and historical
markers in every checkpoint. Rejected trials cannot update them.
This one-way subgrid coating has uncomputed thickness, volume, inertia,
separate dynamics and feedback to bulk. Navier friction does not itself
produce irreversibility. Any resolved thin liquid branch belongs to Ω(t)
and participates in mass and momentum; it is not the subgrid strip.

## Qualification and finite recovery policy

Start with a short genuine moving-domain pilot beyond .0775 s; then 1 s
and target/refinement at fixed physics. The old small-slope screen does not
apply to this nonlinear model. Nonpositive Jacobian, wall crossing, full-curve
self-intersection, loss of mass or unconverged equations reject a trial.
At most 8 dt halvings per accepted step, with a documented minimum dt.
Persistent distortion calls for targeted local refinement and a distinct
mesh/representation remedy within the inherited actual budget, not an
assertion that no continuum solution exists. No shape parameter is fitted.

Independent mathematical audit: `/root/pw2_independent` derived the midpoint
GCL and the explicit energy split independently before production. It does
not certify future numerical data. Method cross-checks:
[ALE transport discussion](https://arxiv.org/html/2003.07166),
[material free-interface weak formulation](https://pyoomph.readthedocs.io/en/latest/tutorial/ale/freesurf.html).
These references are not validation of this vessel/contact approximation.
