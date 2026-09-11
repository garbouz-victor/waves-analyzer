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

## Targeted numerical recovery (declared before its first solve)

The original material12×24 pilot reaches t=.825045166s, then a subsurface
RIGHT element loses its positive Jacobian despite bounded dt reductions.
This is a mesh/representation failure, NOT a continuum impossibility. The
first remedy is local initial refinement of the RIGHT corner at unchanged
physics. The second is same-domain harmonic rezoning with conservative field
transfer; historical original-pressure/material runs remain unchanged.

Rezoning preserves EVERY free-interface P2 coordinate exactly, including P2.
It solves reference-domain harmonic displacement for interior nodes. Wall
mesh nodes can slide tangentially as a coordinate redistribution, with their
fluid velocity still constrained no-slip on RIGHT/BOTTOM. LEFT/right x and
bottom z remain fixed. The material time steps between remaps retain full
convection exactly through moving coordinates; the instantaneous rezoning
is a separately stored, audited projection, NOT unaccounted ALE motion.

At new-mesh quadrature points, locate the old curved element and invert its
P2 map without clipping/extrapolation. Transfer by constrained L2 projection,
enforcing the new weak divergence and preserving both integrated momentum
components. Boundary geometry identity preserves liquid area; independently
measure volume/momentum defects, inverse-map errors, projection residuals and
the actual kinetic/potential change of the remap. The remap work is a signed
numerical defect, not physical dissipation or a fitted global-budget term.
Native stores the remapped starting geometry/velocity, projection multipliers
and old-cell reference coordinates before the next physical accepted step.
H and markers are not updated during remap or rejected attempts.

Two optional consistent weak-form improvements are numerical controls:
`hydrostatic_split`: approximate π=p+gz in P1, actual p=π−gz, replace bulk
gravity by free-boundary traction `Tπ n=−gz n` on the ACTUAL curved surface.
The continuous equations/LEFT tangential law are identical; pressure trial
space at finite h differs. Coupled startup still solves the tilted release.
`skew_divergence`: add `+(1/2 div vm) vm` to midpoint material momentum.
It vanishes for the exact incompressible solution; it is a consistent
geometric stabilization, not extra physical viscosity.

Audit the assembled signed work `dt vmᵀ Kg vm` separately; subtract it from
raw Wk. Audit `dt vm·(f_surface−f_body)` separately; add it to raw Wp.
This yields the remaining explicit cubic geometry terms; the continuous
energy gate still uses raw E−E0+physical losses (including measured remap
defects, never hidden by subtracting them). Independent reviewer confirmed
the algebra before enabling these controls. Numerical identity records all
three choices and their source snapshots. No b, alpha, sigma or BC is changed.

The full P2 harmonic proposal itself folded at t=.5s (RIGHT cell611,
J=−1.40781e−5); that proposal was rejected before any physical state/history
change. A further targeted mesh remedy, `straight_interior`, computes a P1
harmonic vertex target and straight internal edge-midpoints, always retaining
the complete original P2 free curve. A bounded line search also considers
straightening mid-edges without relocating vertices. Among valid proposals
choose the largest minimum Jacobian relative to initial local cell scale;
the unchanged valid mesh is included. Fractions/strategy are stored with
the same conservative transfer evidence. This changes neither liquid-domain
boundary nor BCs. It is not an assertion that every thin region is resolved;
independent full-boundary h and dt comparisons are still required.

The independent old-pilot traction diagnostic exposed a scikit-fem inverse
Newton limitation on highly curved facets. Production now integrates boundary
terms directly on polynomial edge parameters: `n ds=(-z',x')dt` on the ordered
LEFT→RIGHT free curve, and LEFT Robin by the actual arc-length measure.
No reference-plane normals or physical-point extrapolation are used. The
independent audit instead evaluates known reference-edge points and their
Jacobian-based normals; it does not import this production helper.

## Controlled accuracy set, chosen before the h solve

Time pair:12×24 plus one local RIGHT refinement, dt=.01 versus.005s.
Space pair:12×24 versus18×36, each plus one local RIGHT refinement, dt=.005s.
All use the same pressure/geometric forms and conservative straight-interior
rezone interval.05s, with fixed physical b=.50m. A global refinement factor1.5
is selected to reserve actual shared budget for recovery and complete independent
audits; it is not selected by a desired wave/record. Acceptance bounds are not
relaxed. The solver stores actual topology; the comparator verifies its increase.
These two levels provide engineering differences, not a formal convergence order.
