# Corrective PW1 — LEFT Navier / RIGHT no-slip

Model ID: `PW1_LEFT_NAVIER_RIGHT_NOSLIP_V1`. This replaces the current mission
target, not historical data. `../extension/model_decision.md` and its bundles
describe the superseded **HISTORICAL_SYMMETRIC_SIDE_SLIP_VARIANT**: left Navier /
right Navier. Its right BC did not match the corrected user contract. The original
no-slip STOP also remains unchanged. None of their P3/P4 values are new evidence.

## Physical contract, fixed before corrected runs

LEFT WALL x=−a: u=0, (Tn)_tau=−(nu/b_left)w.
RIGHT WALL x=+a: u=w=0.
BOTTOM z=−d: u=w=0.
FREE SURFACE: same linear gravity model, sigma=0.

User baseline: a=1 m,d=10 m,g=9.81 m/s²,alpha=2°,nu=.01 m²/s,
b_left=.50 m. b applies ONLY to the left wall; there is no b_right.
b_left=.25/1.00 m are separate parameter-sensitivity cases, never tuning for a
desired crest and never linked to mesh size. Density1000kg/m³ is only the inherited
illustrative energy scale. There is no continuing horizontal force.

Let v=(u,w), T=−qI+2nuD(v), where q is departure from hydrostatic pressure.
In the reference liquid domain [−a,a]×[−d,0]: v_t=div T,div v=0.
At z=0: Tn=−g eta n and eta_t=w. Initial eta=x tan(2°),v=0.
z=0 is a linearization plane, not a lid; the unbounded graph eta reconstructs the
macroscopic liquid, including eta>0. Finite d=10 m is the actual bottom, not a
claim of exact infinite depth. Nonlinear moving geometry is NOT being solved.

## FEM constraints and energy

fixed = LEFT horizontal DOFs ∪ RIGHT all velocity DOFs ∪ BOTTOM all velocity DOFs.
The shared bottom-left corner is fixed by bottom no-slip.
The reduced kinematic trace has endpoint row nonzeros [1,0].
R_L=eta(−a) evolves from the actual free vertical endpoint DOF, not an interior
crest or extrapolation. P2 is the actual current RIGHT endpoint as well as a
historical marker: eta(+a,t)=z2=a tan(2°) for every accepted state, to machine zero.
Consequently R_R=H_R=z2; no right record is created. This is a hard native/FEM
invariant, not a renderer correction.

Only LEFT contributes K_wall=(nu/b_left)∫LEFT w*v_tau ds. On LEFT n=−ex,
so ∂xw+∂zu=w/b_left. The right wall has no Navier boundary matrix or traction
diagnostic. Its two velocity traces are checked as Dirichlet constraints instead.
The full K_wall (before reduction), fixed DOFs and trace structure are stored in
native data so an erroneous right-wall matrix cannot hide behind zero DOFs.

E/rho=1/2(vᵀMv+g etaᵀS eta), dE/dt=−vᵀK_bulk v−vᵀK_left v.
Both physical losses are nonnegative and separately integrated with the accepted
SDIRK2 stages. Signed RK correction is separate and is not physical heat.
Mass and weak mixed-FEM incompressibility are checked; pointwise divergence is
not falsely required to vanish in this discretization.

Independent steady channel benchmark: nu w''+G=0,
w'(-a)=w(-a)/b_left,w(+a)=0. The formula
w=G/(2nu)(a²−x²)+Ga b_left/[nu(b_left+2a)](a−x)
has w(a)=0,−nu w''=G and
w'(-a)=2Ga²/[nu(b_left+2a)]=w(-a)/b_left.
The independent scalar FE benchmark does not import the production Navier operator.

## Independent irreversible coating law

After EVERY accepted step H_L[n]=max(H_L[n−1],R_L[n]); H_R=z2.
W_i=[−d,H_i], L_total=L_macro ∪ W_L ∪ W_R. W is model/checkpoint state available
without renderer, has no dry holes and never recedes. At RIGHT the main surface
always ends in P2: no fictional film gap is drawn above a receded R_R.
Navier slip and this memory law are separate postulates. The layer's thickness,
volume,inertia,separate dynamics,feedback to bulk and adsorption work are NOT
computed: explicit one-way leading-order subgrid approximation, not molecular proof.

P1/P2 never move. P3/P4/... arise only from actual resolved LEFT local record peaks,
confirmed by the next accepted sample, with the inherited predeclared0.1mm floor.
No continuous between-step peak finder is claimed; sampling/floor sensitivity is
reported. No P4 is a valid outcome.

## Applicability and numerical release

Unchanged hard screens: max exact|eta|/a≤.05, max exact P2|eta_x|≤.30 on ALL
surface edges, including the pinned RIGHT and sliding LEFT corners. Also report
max|eta|/b_left and omitted convection/A0,kinematic/U0 with
A0=g|tan(alpha)|,U0=sqrt(ga)|tan(alpha)|. Quadrature omitted-term indicators are
not rigorous physical-error bounds and a nonzero value alone is not failure.
Corner regularity or uniform C1 convergence is not presumed. A genuine hard-screen
failure will not be waived by omitting wall edges or changing b: preserve pilot,
check dt/mesh evidence and identify the need for nonlinear geometry.

First corrected pilot24×48,dt=.005s,T=1s and end-to-end preview. If qualified,
baseline to5s; time24×48,dt=.0025; space48×96,dt=.0025; b_left=.25/1.00 on fine
mesh/dt. Same-physics refinement compares surface,R_L,H_L,P2,records and energy.
Height norms are divided by the initial height difference;5% engineering gates
do not establish formal convergence order or5% nonlinear physical accuracy.
The declared12h/40GiB budget continues across the correction; one kernel-locked
heavy job. No paid services, push or rewriting historical bundles.
