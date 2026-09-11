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

### Recovery after the measured pilot h failure

The12/18 full-curve pair fails the unchanged5% engineering screen at .95s:
the independently bounded Hausdorff distance is at least5.411712mm (7.75% of
the initial height span). LEFT-contact agreement alone does not qualify it.
The next pair is therefore18×36 versus24×48 at dt=.005s, with a separate
18×36 dt=.01/.005 time pair. All other physical and numerical controls remain
as above. This is refinement in response to error, not calibration to peaks.

At .95s the fine mesh's rezone fast point locator also fails on13 quadrature
points, although the last accepted full boundary and Jacobians are valid.
A bounded fallback searches all cells whose quadratic Bernstein control-point
boxes can contain a target, then uses four interior Newton seeds with a
residual-decreasing line search. The previous reference-domain and physical
inverse residual tolerances are unchanged; outside targets are not clipped,
projected or extrapolated. The failed native state is the regression witness.
This is a numerical inverse-map repair, not a change of the liquid boundary
or of its evolution law. It receives a new numerical-source identity.

The exact failed .95s remap now passes: all13 previously missed points located,
inverse residual≤3.94385e-11m, full-boundary displacement0, volume change0,
integrated momentum defect8.33e-17. Original fast search failure is reproduced
with fallback disabled; the original native SHA is unchanged. The production
record detector also requires each new local peak to equal the accumulated
global H_L, so a lower peak cannot supersede an earlier unmarked plateau.

## Numerically equivalent assembly/linear-solve acceleration (PW2 continuation)

An actual saved 24×48 step cost 20.43 s: assembly 10.86 s and five direct
SuperLU solves 5.75 s. The volume forms are now expanded into scalar-P2
component blocks using the same curved geometry, quadrature order and trial
spaces. No mass lumping, reduced integration, boundary change or viscosity
change is introduced. Scalar-block matrices are independently compared with
the original vector-form assembly. Explicit stored zeros are removed only as
sparse storage, not as physical coefficients.

A fresh-per-physical-step LU of the first current mixed matrix preconditions
bounded GMRES for later Picard iterates. Each iterate still solves its CURRENT
matrix/RHS. Krylov settings: rtol 2e-13, atol 1e-15, restart 30, at most four
restart cycles; true raw momentum and divergence block residuals must each be
at most 1e-13. Failure refactors the current matrix with direct SuperLU/COLAMD;
the existing nonlinear tolerance and final accepted PDE gates are unchanged.
Counters and true residuals are stored. A lagged factor is never treated as
the exact inverse of the new matrix. This is solver preconditioning, not new
physics or artificial dissipation.

The captured 24×48 step replay takes 2.66 s after this change, with maximum
geometry/velocity/pressure differences 1.11e-16 m / 1.28e-15 m/s / 1.68e-14
(pressure per density), compared with the previously saved direct step.
This single-step agreement is not a full-trajectory error estimate. Old native
bytes remain unchanged; new source hashes require new controlled trajectories
from the exact initial condition, not silent checkpoint migration.

### Conservative transfer: exact elimination of two global constraints

A second bounded actual profile identified the remaining cost: 38.21 of
48.65 seconds for one 24×48 remap were spent in direct solves, dominated by
the two dense integrated-momentum rows in the mixed projection. The locator
was not changed. For projection only (its pressure sign differs from NS),
let A0=[[M,B^T],[B,0]], D=[P^T;0]. Factor A0 once and solve three RHS:
y0=A0^-1[f;0], Z=A0^-1D; then S=P Z_velocity,
S lambda=P y0_velocity-target_momentum, y=y0-Z lambda.
The stored multipliers remain [projection_pressure,lambda_x,lambda_z].
This is algebraic Schur elimination of the SAME full mixed system, not the
removal of momentum or pressure constraints. The full original residual is
checked before acceptance, with fallback to its original direct solve.

On a strictly straight tilted interface, Pz-slope*Px can be dependent on weak
divergence. Numerical rank and condition of the computed 2×2 Schur matrix are
therefore recorded. A rank-revealing least-squares solution for a dependent
system is accepted only if ALL original equations still meet the same raw
residual check. There is no epsilon penalty, physical cutoff, or silently
dropped conservation law. A large condition number alone is not a physical
blocker. Both fallback and final conservation gates remain explicit.

The saved 24×48 remap replay now takes4.14s: geometry, source cells and inverse
coordinates identical, velocity difference4.98e-11m/s versus the old full
direct solve, full residual2.50e-16, zero volume change. At the previously
difficult18×36 t=.95s remap, the same check gives6.25s, unchanged full boundary,
velocity difference6.11e-12m/s, full residual2.36e-16. These are measured
algebraic replay comparisons, not full-trajectory or physical error estimates.
Numerical source identity changes; all controlled prefixes are recomputed.

### Additional initial RIGHT refinement at unchanged physics

The 24×48/local-RIGHT-level1 run loses a surface-attached cell near t=.9301s;
bounded time-step reduction and another unchanged-grid rezone do not repair it.
The next actual run uses local RIGHT level2, retaining b=.50m, the same .005s
time cap, .05s rezone interval and all equations above. This refines the initial
mesh, not a physical wall length. It reaches1s without rejected steps.
Independent whole-interface comparison against18×36/local-level1 on all201
union accepted times0–1s gives Hausdorff upper3.155922mm (4.51869% of the
initial span), maximum LEFT R/H difference4.28380µm, and P2 drift0. This is
pilot engineering evidence only;0–5s comparison and native audit remain gates.

For possible later mesh recovery a separately tested dimensionless mesh-quality
objective has been prepared, but it is NOT enabled in these runs. No physical
pseudo-elastic stress or film thickness has been added. Any future use requires
an explicit new numerical identity, same-boundary transfer and a fresh audit.

The declared normal project dependency is now SciPy>=1.12,<2, matching the
[GMRES rtol/atol API](https://docs.scipy.org/doc/scipy-1.12.0/reference/generated/scipy.sparse.linalg.gmres.html).
The actual numerical environment remains SciPy1.13.1; this dependency correction
does not change existing native numerical identities or execution libraries.

## Same-boundary quality recovery after the actual1.4583s mesh failure

The level2 run later develops a near-zero Jacobian at an INTERIOR point of a
free-surface edge of cell3000, not atP2. Its boundary tangent remains nonzero;
no extra wall contact/self-intersection is detected. Interior vertex and edge
DOFs can change the transverse map derivative while the full curve stays fixed.
This is targeted numerical recovery, not evidence of a missing physical law.

New strategy `quality_optimized` keeps the existing straight-interior/harmonic
proposals and adds a pure r-mesh optimization when their best minimum relative
Jacobian is below0.2. On each affine INITIAL reference cell, F=J_current J_ref^-1,
j=detF>0. Use dimensionless, uniformly cell-weighted
W=(F:F−2)/2−log(j)+(log(j))²/2, with P=F+(log(j)−1)F^-T.
Additional dimensionless barriers0.25*(j−1−logj) are averaged at the vertices
and evaluated at the EXACT minimum-J witness in each closed reference triangle.
The witness can be inside an edge or triangle, not only a quadrature point.
At a unique witness the envelope derivative is0.25*(j−1)F^-T; changes of active
witness can be nonsmooth. No claim of global mesh-optimality is made.

Every optimization trial requires an exact positive signed Jacobian on every
closed triangle. All free-interface P2 coordinates, wall-normal coordinates,
bottom coordinates and endpoint locations are fixed. Dimensionless displacements
use a fixed inverse-square-root initial Laplace-diagonal coordinate scaling.
Bounded feasible Armijo L-BFGS uses80iterations,10history vectors,60halvings per
line search, Armijo coefficient1e-4 and scaled-gradient tolerance1e-7. An unchanged
objective after an infeasible trial is NOT convergence. Select only a valid
candidate with a measured improvement in minimum relative J; unchanged geometry
remains allowed. These constants are mesh/solver choices, not physical elasticity,
film width, a contact angle, or any alteration to b=.50m or the wall conditions.

Rezoning occurs at the existing .05s interval OR earlier when a current accepted
mesh's exact minimum relative Jacobian drops below0.01. The independent verifier
reconstructs this trigger from geometry. This is a mesh-quality trigger normalized
by the initial local cell, not a physical length or a resolved-film cutoff.
The same conservative transfer, stored inverse maps/multipliers, exact-Dirichlet
velocity constraints and signed energy-work accounting apply. No accepted coating
state changes during a rezone or a rejected trial.

Actual read-only recovery: minJ1.06069e-10→8.12449e-7m², complete free-boundary
displacement0. Transferred volume change−3.55e-15m², momentum1.18e-16,
projection residual1.11e-16, inverse-map error1.60e-11m. The next .005s material
candidate has positive minJ5.86566e-7, momentum1.12e-14 and P2drift0. This probe
is NOT appended to the old native history or represented as an accepted run.
Independent audit REJECTS that particular .005s trial: its full quadratic
surface intersects itself at midpoint and endpoint despite positive local J
and small PDE residuals. The actual intersection roots are preserved as evidence.
Production therefore has an additional64-subdivision, edge-box-pruned early
crossing rejection on both midpoint and end curves, before accepting a step or
updating coating. This early polygon test is not a geometric certificate:
the separate verifier still checks controlled P2 curves and exact additional
wall contacts at both collocation and accepted states. New runs use earlier
adaptive mesh recovery from their initial condition, not this rejected trial.
New numerical source identity includes the mesh-quality module; controlled
trajectories must restart from the exact initial condition. The previous1s hPASS
remains historical numerical-method evidence, not refinement of this new method.

### Further spatial remedy under the same frozen equations/implementation

The fresh quality-optimized localRIGHT2 trajectory reaches1.460673828125s,
then resolved midpoint/end curve intersections persist as dt is halved. Tiny
attempts also hit the unchanged raw mixed residual check. Interior r-motion
has improved cell Jacobians but cannot fix an intersecting material surface.
The next initial mesh uses localRIGHT3, keeping24×48 global controls,dt=.005s,
the same rezone/quality controls and every physical parameter. This subdivides
the material surface edges originally at x≈.98561–.99640m, each previously
about3.599mm long in x, which feed the unresolved shoulder. It is an initial
spatial refinement, not a prescribed future film width or contact geometry.
The trajectory starts at the exact tilted line and zero velocity; no numerical
source bytes change and no old checkpoint is appended under new controls.
Complete-curve comparisons on their actual common interval remain required;
a failed coarse horizon cannot become a5s refinement PASS.

### Consistent grad-div remedy, fixed before new initial solves

The localRIGHT3 run fails at1.3953323364257728s. Its limiting free-edge
tangent has minimum norm2.05142e-7m (reference edge parameter), and the
one-sided midpoint divergence at the exact minimum-J witness is about
−9.92e4/s, despite a small pressure-tested weak-divergence residual. A bounded
500-iteration interior r-optimization improves minimum J only1.74510e-10 to
1.96549e-10m² and the next .005s candidate still inverts. No surface coordinate
was modified by that probe. This motivates controlling unresolved divergence,
not changing the physical viscosity or prescribing a surface profile.

For NEW initial-state trajectories select a fixed numerical grad-div coefficient
γ=.1m²/s (=10ν), independent of cell size, dt and desired wave peaks. Add
`Gγ(v,q)=γ∫Ωm div(v) div(q)` to midpoint momentum, keeping the entire
`B vm=0` pressure constraint. This term vanishes for exact incompressible
solutions. LEFT normal, RIGHT/BOTTOM full velocity constraints and the
LEFT-only Navier matrix are unchanged. The numerical parameter is stored in
controls/numerical identity, not the physical contract. Old γ=0 native data
remain historical method evidence, not refinement of this method.

The finite-h weak natural stress is augmented by `γ div(vm) I`. Thus report
physical free-surface traction `T_phys n`, the separate numerical normal
traction `γ div(vm)n`, and augmented traction separately; do not claim the
physical trace is pointwise zero just because the augmented weak solve passes.
On LEFT this added stress is normal, leaving its tangential Navier law unchanged.
No cancellation force, prescribed angle, surface smoothing or endpoint clamp
is introduced. The primary [grad-div consistency analysis](https://arxiv.org/abs/1610.05017)
concerns fixed-domain Dirichlet problems, not a stability/convergence proof for
this moving free-boundary/contact problem. Grad-div controls an L2 error, not
pointwise incompressibility or cusp regularity; curved quadrature and full
interface refinement still require independent checking.

Measure `Wγ=−dt γ∫Ωm (div vm)²≤0` independently, save it per accepted step
and cumulatively, and expose it without the renderer. The split identity is
`E−E0+D_bulk+D_LEFT−ΣWk−ΣWp−ΣW_rezone−ΣWγ=0`.
The raw physical continuous-energy gate remains
`|E−E0+D_bulk+D_LEFT|/E0≤.05`, with NO subtraction of Wγ to pass that gate.
Gγ is independently reconstructed in the unreduced space and checked for
symmetry, nonnegative quadratic work, action, coefficient and source identity.

New controlled set, selected before the first γ=.1 solve:12×24/global with
two local RIGHT levels, dt=.01/.005s for temporal comparison;24×48 with
the same two local RIGHT levels at .005s for spatial comparison. All use
quality_optimized rezone at .05s/unchanged quality trigger, hydrostatic split,
skew-divergence and fixed γ=.1. First pause the requested5s target at1.5s
for pilot inspection; resume the SAME source/controls native if admissible.
The local refinement targets the previously observed right region; no output
shape/record values are targets. Full5s comparisons and all unchanged gates
are required before COMPLETE; two levels imply no formal convergence order.

### Additional local h remedy within the fixed γ=.1 method

The24×48/RIGHT-local2,.005s trajectory passes an independent all-state audit
through1.5s (301states,601accepted/midpoint curves), but the full h comparison
against12×24 is not qualified: sampled distance3.741315mm, lower3.739317mm,
upper4.243162mm versus3.492077mm tolerance. The witness lies at the right
shoulder, around x=.990–.993m,z=.008–.011m, not at fixedP2. Continuation
rejects full-curve intersections and ends at1.5219157600402722s.

The next targeted numerical remedy increases ONLY initial local RIGHT
refinement2→3 on the same24×48 grid, dt=.005s, γ=.1m²/s, rezone/quality and
all physical settings unchanged. It starts from the exact tilted line/rest.
The earlier γ=0/local3 failure is historical and not a controlled test of
this γ=.1 discretization. No b, sigma, angle, contact position, future film
width or record target is changed. The actual stopping times are properties
of these finite numerical policies, not established continuum singularities.
Independent full-curve comparison is required on the explicitly stated common
accepted interval; a prefix comparison cannot qualify the missing0–5s result.

### Acceptance-guard correction after independent local3 audit

The independent all-state audit found tiny adjacent quadratic loops in the
OLD local3 accepted states294–296, first at t=1.46625s. The64-chord early
guard missed these loops inside the shared last/first chord pair. That is an
implementation defect, not a physical law and not a qualified native result.
The old native file and its FAIL audit are preserved without modification.

The NEW production guard additionally checks adjacent full quadratic edges.
Writing both edges from their shared join as a*u+b*u² and c*v+d*v²,
u,v>0, v=r*u yields the cubic cross(a-c*r,b-d*r²)=0. Exact rational native
coefficients identify collinearity and shared-tangency factors; numerical
cubic roots are tested in the original two-coordinate equations and within
both parameter intervals. Both orientations use0<r<=1. Only the algebraic
known shared parameter pair is excluded, never a neighbourhood of the join.
The independent verifier keeps its DIFFERENT full-quadratic elimination
implementation and imports no production guard. Neither finite-precision
root method claims a rigorous separation certificate.

This change affects step acceptance, not NS equations, viscosity, gamma,
wall/free-surface conditions or coating law. It changes numerical source
identity. A fresh24×48/RIGHT-local3,dt=.005,γ=.1 target starts at exact tilt
and rest; it does not resume the three invalid states or silently relabel
the old source run. The earlier dt/h data remain explicitly versioned
evidence of their own accepted prefixes. No altered historical native and
no new physical length, contact angle or interface profile is introduced.
