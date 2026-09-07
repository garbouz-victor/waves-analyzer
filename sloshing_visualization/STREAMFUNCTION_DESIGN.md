# LINEAR STREAMFUNCTION SLOSHING — SF-0 formulation and numerical design

Status: design and sign-verification checkpoint.  This document does **not**
claim that a production solver or a physical solution has been obtained.  The
work is independent of STEP 3/3A CHNS; no CHNS equation, discretization, contact
model, or result is used here.

Throughout, subscripts denote partial derivatives.  The reference domain is

\[
  \Omega_0=\{(x,z):-a<x<a,\ -\infty<z<0\}.
\]

## 1. Physical problem

The fluid is two-dimensional, incompressible, Newtonian, and viscous.  The
coordinates are horizontal \(x\) and upward vertical \(z\).  The side walls are
at \(x=\pm a\), with \(a=1\ \mathrm m\), and the liquid is semi-infinite below.
The undisturbed free surface is \(z=0\).  The data are

\[
 a=1\ \mathrm m,\qquad g=9.81\ \mathrm{m\,s^{-2}},\qquad
 \nu=10^{-6}\ \mathrm{m^2\,s^{-1}},\qquad \sigma=0,
\]

with configurable density \(\rho\) and configurable dimensionless slope \(k\).
The baseline for later small-amplitude examples will be \(k_{\rm ref}=10^{-3}\),
not a finite-angle interpretation of the initial line.

The initial conditions are

\[
 \eta(x,0)=kx,\qquad \boldsymbol v(x,z,0)=\boldsymbol 0.
\]

The model is exactly linear about rest and the flat surface.  Consequently all
free-surface equations are imposed at **\(z=0\)**.  The graph \(z=kx\) is an
initial value of the unknown \(\eta\), not the boundary on which the linearized
conditions are evaluated.  Applying them on \(z=kx\) would introduce a different,
partly nonlinear geometry.

No small-\(\nu\) term is discarded in any derivation.  This is distinct from the
numerical fact that the specified \(\nu\) is very small on the gravity-wave scale.

## 2. Coordinate and sign conventions

One convention is used everywhere:

\[
 u=\psi_z,\qquad w=-\psi_x,
\]

where \(u\) is horizontal and \(w\) is upward velocity.  Hence

\[
 u_x+w_z=\psi_{zx}-\psi_{xz}=0,
 \qquad
 \omega=w_x-u_z=-\psi_{xx}-\psi_{zz}=-\Delta\psi.
\]

The perturbation pressure \(p\) is measured relative to the hydrostatic base
pressure.  The total pressure is

\[
 P=P_{\rm atm}-\rho g z+p.
\]

This convention determines the gravity sign in the normal-stress condition.

## 3. Linearized primitive equations

On the fixed reference domain \(\Omega_0\),

\[
 \begin{aligned}
 u_t&=-\frac{p_x}{\rho}+\nu\Delta u,\\
 w_t&=-\frac{p_z}{\rho}+\nu\Delta w,\\
 u_x+w_z&=0.
 \end{aligned}
\]

These are perturbation equations: gravity has been absorbed into the
hydrostatic base pressure and returns through the free-surface normal stress.

## 4. Streamfunction reduction

Take \(\partial_x\) of vertical momentum minus \(\partial_z\) of horizontal
momentum.  Mixed pressure derivatives cancel, giving

\[
 \omega_t=\nu\Delta\omega.
\]

Because \(\omega=-\Delta\psi\), multiplication by \(-1\) gives the bulk equation

\[
 \boxed{(\Delta\psi)_t=\nu\Delta^2\psi.}
\]

The componentwise symbolic regression in
`scripts/derive_streamfunction_conditions.py` also verifies that the curl of
the two primitive momentum residuals is the negative of this boxed residual.

## 5. Wall conditions

At either vertical wall physical no slip is

\[
 u=0,\quad w=0
 \quad\Longleftrightarrow\quad
 \psi_z=0,\quad\psi_x=0.
\]

Impermeability is \(u=\psi_z=0\).  It says that \(\psi\) is constant along each
connected vertical wall.  The difference between the two wall constants is

\[
 \psi(a,z)-\psi(-a,z)=\int_{-a}^{a}\psi_x\,dx
 =-\int_{-a}^{a}w\,dx,
\]

the signed vertical volume flux through a horizontal section.  Incompressibility
and impermeable side walls make that flux independent of depth; decay at
\(z\to-\infty\) makes it zero.  The constants therefore coincide.  A gauge shift
then gives the equivalent and numerically convenient statement

\[
 \boxed{\psi=0,\qquad\psi_x=0\quad\text{at }x=\pm a.}
\]

Here \(\psi=0\) is not extra physics: along a wall it is exactly the integrated
form of \(\psi_z=0\), after the zero-flux gauge has been chosen.  The second
condition is tangential no slip.

## 6. Deep-water condition

The physical condition is

\[
 \nabla\psi\longrightarrow0\quad(z\to-\infty),
\]

and the gauge \(\psi\to0\) is chosen there.  A computation on
\(-H\le z\le0\) is only a **truncated approximation to semi-infinite depth**.
It is not a finite-depth physical tank.

For truncation, the proposed far-field homogeneous surrogate is

\[
 \psi=0,\qquad\psi_z=0\quad(z=-H).
\]

It happens to imply zero velocity on the artificial line, but it must not be
called a physical no-slip bottom.  Its acceptability is conditional on an
explicit \(H\)-convergence study.

## 7. Kinematic free-surface condition

Linearizing material motion of the surface on \(z=0\) gives

\[
 \eta_t=w(x,0,t),
\]

so the primary kinematic equation is

\[
 \boxed{\eta_t=-\psi_x\quad(z=0).}
\]

## 8. Tangential stress condition

With zero gas shear and \(\sigma=0\), the linearized tangential traction is

\[
 \mu(u_z+w_x)=0.
\]

The chosen streamfunction signs give

\[
 u_z+w_x=\psi_{zz}-\psi_{xx},
\]

and therefore

\[
 \boxed{\psi_{zz}-\psi_{xx}=0\quad(z=0).}
\]

## 9. Normal stress condition

Let \(\boldsymbol T=-P\boldsymbol I+2\mu\boldsymbol e\),
\(\mu=\rho\nu\), and neglect gas viscous stress.  At the actual surface
\(z=\eta\), the normal traction balances \(-P_{\rm atm}\boldsymbol n\).  To first
order, evaluation of the hydrostatic base pressure gives

\[
 P_0(\eta)=P_{\rm atm}-\rho g\eta.
\]

Using \(\boldsymbol n=\boldsymbol e_z+O(\eta_x)\), base isotropic-traction terms
proportional to the perturbed normal cancel between the two sides.  The remaining
normal component is

\[
 -\big(P_{\rm atm}-\rho g\eta+p\big)+2\mu w_z=-P_{\rm atm},
\]

or

\[
 \boxed{p-2\rho\nu w_z=\rho g\eta.}
\]

Since \(w_z=-\psi_{xz}\), the streamfunction form is

\[
 \boxed{\frac p\rho+2\nu\psi_{xz}=g\eta\quad(z=0).}
\]

Thus the viscous streamfunction term has a **plus** sign.  This sign follows
from both the upward normal and the convention \(w=-\psi_x\).

## 10. Pressure elimination

Horizontal momentum becomes

\[
 \psi_{zt}=-\frac{p_x}{\rho}
 +\nu(\psi_{xxz}+\psi_{zzz}).
\]

Differentiating the normal-stress condition in \(x\) yields

\[
 \frac{p_x}{\rho}+2\nu\psi_{xxz}=g\eta_x.
\]

Substitution, without dropping viscosity, gives

\[
 \psi_{zt}=-g\eta_x+\nu(3\psi_{xxz}+\psi_{zzz}),
\]

or the primary dynamic condition

\[
 \boxed{
 \psi_{zt}-\nu(3\psi_{xxz}+\psi_{zzz})=-g\eta_x
 \quad(z=0).}
\]

The coefficient 3 consists of one \(\psi_{xxz}\) from \(\Delta u\) and two from
the differentiated normal viscous stress.  In the solved-for form the entire
viscous contribution has a **plus** sign.

## 11. Primary \(\psi+\eta\) formulation

The primary unknowns remain \(\psi(x,z,t)\) and \(\eta(x,t)\).  The complete
semi-infinite system is

\[
 \begin{array}{ll}
 (\Delta\psi)_t=\nu\Delta^2\psi,&(x,z)\in\Omega_0,\\[2mm]
 \eta_t=-\psi_x,&z=0,\\
 \psi_{zz}-\psi_{xx}=0,&z=0,\\
 \psi_{zt}-\nu(3\psi_{xxz}+\psi_{zzz})=-g\eta_x,&z=0,\\[1mm]
 \psi=0,\ \psi_x=0,&x=\pm a,\\
 \nabla\psi\to0,\ \psi\to0,&z\to-\infty.
 \end{array}
\]

Keeping \(\eta\) makes the tilted initial surface explicit and prevents pressure
elimination from erasing the restoring force.

## 12. Secondary \(\psi\)-only formulation

Differentiate the primary dynamic condition in time:

\[
 \psi_{ztt}-\nu(3\psi_{xxzt}+\psi_{zzzt})=-g\eta_{xt}.
\]

The differentiated kinematic condition is \(\eta_{xt}=-\psi_{xx}\).  Therefore

\[
 \boxed{
 \psi_{ztt}-g\psi_{xx}
 -\nu(3\psi_{xxzt}+\psi_{zzzt})=0
 \quad(z=0),}
\]

equivalently

\[
 \boxed{
 \psi_{ztt}=g\psi_{xx}
 +\nu(3\psi_{xxzt}+\psi_{zzzt}).}
\]

The viscous term is again **positive** in the solved-for form.  This formulation
is secondary and intended only for analytic/sign validation in SF-0/SF-1.  It is
not a production scheme until its additional initial derivative data and corner
compatibility are handled explicitly.

## 13. Initial conditions

The primary data are simply

\[
 \boxed{\eta(x,0)=kx,\qquad\psi(x,z,0)=0.}
\]

Kinematics then gives \(\eta_t(x,0)=0\).  Linearity gives exact amplitude
scaling: for a response computed at \(k_0\),

\[
 (\psi,\eta,\boldsymbol v,p)_{2k_0}=2
 (\psi,\eta,\boldsymbol v,p)_{k_0},\qquad E_{2k_0}=4E_{k_0}.
\]

A normalized unit-slope response can therefore be rescaled, subject only to the
small-amplitude validity of the linear physical model.

## 14. Initial acceleration

At \(t=0\), all spatial derivatives of \(\psi\) vanish.  At smooth interior
points of the free surface, the primary dynamic condition gives

\[
 \boxed{\psi_{zt}(x,0,0)=-gk,\qquad |x|<a.}
\]

Equivalently, normal stress initially gives \(p(x,0,0)=\rho gkx\), so horizontal
momentum gives \(u_t=-p_x/\rho=-gk\) away from the walls.

This is required initial information for the \(\psi\)-only representation.  The
homogeneous \(\psi\)-only differential equation also admits \(\psi=0\); if the
condition \(\psi_{zt}=-gk\) (or its equivalent consistent initial derivative) is
forgotten, the tilted surface is lost and a spurious permanently trivial state
results.  This is a mandatory regression test.

At a strict no-slip wall, however, \(u=\psi_z=0\) for all time, hence
\(\psi_{zt}=u_t=0\) at the corner.  It cannot simultaneously equal \(-gk\) for
\(k\ne0\).  The free-surface interior acceleration is therefore not extended by
assumption to the endpoints.

More strongly, let \(\chi=\psi_t(\cdot,\cdot,0)\) and suppose for the moment that
it is a classical smooth field.  The bulk equation at \(\psi(0)=0\) gives
\(\Delta\chi=0\).  Time differentiation of the gauged wall conditions gives
\(\chi=\chi_x=0\) on each wall, while primary surface dynamics requires
\(\chi_z=-gk\) on the open top edge.  Zero Dirichlet and normal Cauchy data for a
harmonic field on an open wall would force the smooth continuation to be zero,
contradicting the top trace.  Thus there is no globally smooth classical initial
acceleration field satisfying every trace for \(k\ne0\).  A weak evolution must
form the startup wall/corner structure instantaneously; SF-1 consistent
initialization must not replace it by a uniform acceleration or silently smooth
the tilt.

## 15. Contact-line pinning

At either contact point, wall no slip gives \(w=-\psi_x=0\).  The kinematic
condition then gives

\[
 \eta_t(\pm a,t)=0.
\]

Thus

\[
 \boxed{\eta(a,t)=ka,\qquad\eta(-a,t)=-ka.}
\]

These are pinned contact points.  The endpoints must not be described or drawn
as a moving contact line.  This is a strong physical limitation of the strict
linear no-slip model, not a post-processing choice.

## 16. Corner incompatibility

The initial free-surface pressure has the nonzero tangential gradient
\(p_x=\rho gk\).  It demands \(u_t=-gk\) on the smooth free-surface interior.
No slip demands \(u_t=w_t=0\) on the wall, including its trace at the contact
point.  In addition, free-surface stress conditions, a pinned endpoint, and the
wall's two streamfunction conditions meet at the same point.  The classical
boundary data are therefore incompatible at \((\pm a,0)\) at startup.

One should expect a singular or nonuniformly convergent contact-region response,
especially in derivatives, vorticity, and wall shear.  It is not automatically a
solver bug.  Nor may it be silently removed by smoothing \(kx\); smoothing would
define a separate initial/contact model.

Validation will report three regions separately:

* bulk proposal: \(|x|\le0.9a\);
* transition proposal: \(0.9a<|x|\le0.98a\);
* contact/corner proposal: \(|x|>0.98a\).

The thresholds are diagnostic proposals, not physical constants.  A single
global norm must not be used to claim both bulk and corner convergence.

At long time, dissipation suggests flattening in the bulk, but not uniform
convergence \(\eta\to0\) on the closed interval: the two endpoint values remain
\(\pm ka\), and increasingly sharp near-wall structure can develop.

## 17. Energy identity

For the linearized fixed domain, define

\[
 K=\frac\rho2\int_{\Omega_0}
 (\psi_x^2+\psi_z^2)\,dA,
 \qquad
 P=\frac{\rho g}{2}\int_{-a}^{a}\eta^2\,dx,
 \qquad E=K+P.
\]

Write primitive momentum in stress-divergence form,
\(\rho\boldsymbol v_t=\nabla\cdot\boldsymbol\tau\), where
\(\boldsymbol\tau=-p\boldsymbol I+2\mu\boldsymbol e\).  Dot with
\(\boldsymbol v\), integrate, and use incompressibility:

\[
 \frac{dK}{dt}
 =\int_{\partial\Omega_0}\boldsymbol v\cdot
 (\boldsymbol\tau\boldsymbol n)\,ds
 -\int_{\Omega_0}2\mu\boldsymbol e:\boldsymbol e\,dA.
\]

Wall work is zero by no slip and the deep contribution vanishes by decay.  At
the free surface, tangential traction is zero and normal traction is

\[
 \tau_{zz}=-p+2\mu w_z=-\rho g\eta.
\]

Thus its work is

\[
 \int_{-a}^{a}w\tau_{zz}\,dx
 =-\rho g\int_{-a}^{a}\eta\eta_t\,dx=-\frac{dP}{dt},
\]

where the kinematic condition was used.  Consequently

\[
 \boxed{\frac{dE}{dt}=-D.}
\]

For the selected convention,

\[
 e_{xx}=u_x=\psi_{xz},\qquad
 e_{zz}=w_z=-\psi_{xz},\qquad
 e_{xz}=\frac12(u_z+w_x)=\frac12(\psi_{zz}-\psi_{xx}).
\]

Since \(\boldsymbol e:\boldsymbol e=e_{xx}^2+e_{zz}^2+2e_{xz}^2\),

\[
 2\mu\boldsymbol e:\boldsymbol e
 =\mu\left[4\psi_{xz}^2+(\psi_{zz}-\psi_{xx})^2\right].
\]

Therefore

\[
 \boxed{
 D=\rho\nu\int_{\Omega_0}
 \left[4\psi_{xz}^2+(\psi_{zz}-\psi_{xx})^2\right]dA\ge0.}
\]

The derivation includes the surface boundary work; it is not an assumed decay
law.  It is classical for sufficiently regular solutions.  With the startup
corner incompatibility it is to be understood through the corresponding weak
energy balance (or a domain with vanishing corner cutouts and no residual corner
power).  Discrete energy balance/inequality and positive dissipation are primary
future gates.

## 18. Pressure reconstruction

Pressure is absent from the primary solve but can be reconstructed diagnostically.
Horizontal momentum gives

\[
 \boxed{\frac{p_x}{\rho}=-\psi_{zt}
 +\nu(\psi_{xxz}+\psi_{zzz}).}
\]

Vertical momentum, with \(w=-\psi_x\), independently gives

\[
 \boxed{\frac{p_z}{\rho}=\psi_{xt}
 -\nu(\psi_{xxx}+\psi_{xzz}).}
\]

Their cross-derivative mismatch is

\[
 \partial_z\!\left(\frac{p_x}{\rho}\right)
 -\partial_x\!\left(\frac{p_z}{\rho}\right)
 =-(\Delta\psi)_t+\nu\Delta^2\psi,
\]

so a single-valued pressure exists exactly when the streamfunction bulk PDE is
satisfied.  A future least-squares/Poisson integration will determine \(p\) up to
a time-dependent constant, fixed by the surface normal-stress condition.  The
independent diagnostic is then

\[
 p/\rho+2\nu\psi_{xz}-g\eta=0\quad(z=0).
\]

No production pressure reconstruction is part of SF-0.

## 19. Nondimensionalization

Use

\[
 x=a x^*,\quad z=a z^*,\quad
 t=\sqrt{a/g}\,t^*,\quad
 \eta=ka\,\eta^*,\quad
 \psi=ka\sqrt{ga}\,\psi^*,\quad
 p=\rho gka\,p^*.
\]

After division by the corresponding scales, the boxed primary equations retain
their form with \(a=g=k=1\) and

\[
 \boxed{\nu^*=\frac{\nu}{a\sqrt{ga}}.}
\]

For the specified values,

\[
 \boxed{\nu^*=3.1927542840705044\times10^{-7}.}
\]

The solver should operate in these nondimensional variables and convert
diagnostics back to SI.  This avoids mixing order-one gravity dynamics with tiny
dimensional viscous coefficients and sub-millimetre layers.  Density cancels
from velocity/surface evolution; it remains in pressure and energy scales.

## 20. Relevant scales for \(\nu=10^{-6}\)

The gravity velocity is \(U_g=\sqrt{ga}=3.1320919527\ \mathrm{m\,s^{-1}}\).
The gravity Reynolds number is

\[
 \boxed{Re_g=\frac{U_g a}{\nu}=3.1320919527\times10^6.}
\]

Thus the derivation makes no small-viscosity approximation, while the requested
physical parameter is nevertheless a high-Reynolds-number case.

Using the first inviscid frequency below, the characteristic oscillatory Stokes
layer is

\[
 \boxed{\delta_s=\sqrt{\frac{2\nu}{\omega_1}}
 =7.1378559101\times10^{-4}\ \mathrm m,}
\]

about \(0.714\ \mathrm{mm}\), versus a vessel width of \(2\ \mathrm m\).  A
uniform production grid capable of resolving this scale is impractical.  Small
\(\nu\) is not numerically easy: no-slip creates rapidly decaying high-wavenumber
modes and a stiff spatial operator.  A large stable time step may resolve bulk
surface motion while missing wall vorticity and startup transients.

The later interval \(0\le t\le5\ \mathrm s\) spans about 3.12 first-mode periods,
but no such production run belongs to SF-0.

## 21. Infinite-depth reference modes

For an inviscid infinite-depth problem with impermeable vertical walls, the odd
surface modes use

\[
 q_n=\frac{(n+1/2)\pi}{a},\qquad n=0,1,\ldots.
\]

Using the requested label \(q_1\) for the first mode (the \(n=0\) member of the
sequence above),

\[
 q_1=\frac{\pi}{2a}=1.5707963268\ \mathrm{m^{-1}},
\]

\[
 \boxed{\omega_1=\sqrt{gq_1}=3.9254951237\ \mathrm{s^{-1}},\qquad
 T_1=\frac{2\pi}{\omega_1}=1.6006096325\ \mathrm s.}
\]

This is a small-viscosity **bulk reference**, not an exact viscous no-slip mode:
an inviscid wall permits vertical tangential motion and hence does not pin the
contact points.

Orthogonality gives the initial-surface expansion

\[
 kx=\sum_{n=0}^{\infty}b_n\sin(q_nx),\qquad
 b_n=\frac1a\int_{-a}^{a}kx\sin(q_nx)\,dx
 =\boxed{\frac{2k(-1)^n}{a q_n^2}}
 =\frac{8ka(-1)^n}{(2n+1)^2\pi^2}.
\]

The coefficients decay as \(n^{-2}\), and the surface series converges to
\(\pm ka\) at the endpoints.  Each basis function has zero derivative there,
however, whereas \((kx)_x=k\).  Derivative convergence is therefore nonuniform
near the endpoints.  Time evolution of this inviscid series also does not retain
the viscous model's pinned endpoint values.  These facts delimit its use to
early-time/bulk benchmarking.

## 22. Semi-infinite numerical strategies

Three approaches were compared.

### A. Depth truncation with convergence in \(H\)

Use a nonuniform grid on \([-H,0]\) and the far-field homogeneous surrogate
\(\psi=\psi_z=0\) at \(-H\).  It preserves sparse local operators, permits strong
surface clustering, makes two conditions for the fourth-order equation explicit,
and supports direct energy quadrature.  Its approximation error is measurable by
varying \(H\).  The drawback is that semi-infinite depth is not exact at any one
\(H\).

### B. Mapped semi-infinite coordinate

A rational map such as \(z=-L s/(1-s)\), \(0\le s<1\), represents infinity
directly and clusters near the surface.  Near \(s=1\), transformed fourth-order
coefficients and quadrature weights span extreme scales; two decay conditions,
local FD closure, DAE conditioning, and accurate energy integration become more
delicate.  This is worth a small-system prototype, but not the first production
path.

### C. Modal/spectral depth representation

Laguerre/rational-Chebyshev or decaying exponential bases give high depth
accuracy for smooth modes and are excellent independent references.  The
fourth-order boundary conditions and contact singularities require many modes;
collocation is dense in depth and couples awkwardly to strongly clustered wall
resolution.  Tensor/spectral scaling and reliable corner energy are less
attractive for the eventual large system.

**Selected semi-infinite treatment:** A, a sequence of depth-truncated sparse
problems with predeclared \(H\)-convergence gates.  It is selected for controlled
error and sparse robustness, not merely ease of coding.  B and C remain
independent small-system comparisons.

## 23. Recommended discretization

The selected primary spatial method to investigate in SF-1 is a **sparse
high-order finite-difference DAE on a boundary-fitted, nonuniform tensor grid**.

Use uniform computational coordinates \(\xi\in[-1,1]\), \(r\in[0,1]\) with

\[
 x(\xi)=a\frac{\tanh(\beta_x\xi)}{\tanh\beta_x},\qquad
 z(r)=-\ell_z\sinh\!\left(r\,\operatorname{asinh}(H/\ell_z)\right).
\]

The starting feasibility parameters are \(\beta_x\simeq3.4\) and
\(\ell_z\simeq2\ \mathrm{mm}\), to be tuned by conditioning and convergence,
not frozen as physical constants.  The maps cluster at both walls and at the
free surface while coarsening smoothly into the decayed far field.

Construct first through fourth derivative matrices directly on physical nodes
with 9-point local polynomial/Fornberg weights, centered in the interior and
one-sided at boundaries.  On a generic nonuniform grid, degree-8 exactness gives
derivative-dependent formal truncation orders; no single optimistic global order
will be claimed without MMS.  Tensor/Kronecker operators keep the matrices
sparse.  The mixed \(D_{xx}D_{zz}\) term gives at most an \(O(9^2)\) local
footprint, so storage remains \(O(N_xN_z)\), although sparse-factorization fill
must be measured.

The unknown vector contains all nodal \(\psi_{ij}\), followed by all surface
\(\eta_i\).  Boundary equations are tau/collocation row replacements, not a
change in unknown physics.  The exact proposed allocation is in Section 25.

Why the other spatial candidates are secondary:

* global spectral collocation is an excellent small-grid sign/reference solver,
  but corner singularities make convergence nonuniform and dense factors scale
  poorly;
* a conforming biharmonic FEM requires robust \(C^1\) elements (for example
  Argyris-type).  Standard FEniCSx Lagrange elements are only \(C^0\) and cannot
  be assumed to provide a direct conforming biharmonic solve.  High-order
  triangular \(C^1\) traces with third-derivative surface conditions and extreme
  anisotropy need a separate library feasibility study.  They do not beat the
  tensor sparse approach at SF-0.

No CHNS discretization is reused or forced onto this problem.

## 24. DAE structure

Let \(\Psi\) contain streamfunction DOFs, \(H_\eta\) surface DOFs, and
\(y=[\Psi,H_\eta]^T\).  Before boundary row replacement, the conceptual rows are

\[
 \begin{array}{rcl}
 L\dot\Psi&=&\nu L^2\Psi,\\
 D_z\dot\Psi_s&=&\nu(3D_{xx}D_z+D_{zzz})\Psi_s-gD_xH_\eta,\\
 \dot H_\eta&=&-D_x\Psi_s,\\
 0&=&(D_{zz}-D_{xx})\Psi_s,
 \end{array}
\]

plus algebraic wall and deep conditions.  Thus

\[
 \boxed{M\dot y=Ay}
\]

has a singular mass matrix: tangential-stress, wall, and far-field rows are
algebraic.  Coefficients are constant, the initial \(kx\) remains in the state,
and no nonlinear/Newton iterations are needed.  For a fixed time step, an
implicit stage matrix can be factored/preconditioned once and reused.  The same
pencil supports later generalized eigenanalysis.

SF-1 must first check matrix rank, constraint consistency, and the differentiation
index on tiny grids.  Production assembly is forbidden until those row-count and
rank tests pass.

## 25. Corner row treatment

Blindly imposing every wall and surface equation at a top corner would
overdetermine the discretization and contradict the initial classical data.
Wall no-slip conditions take priority at corner row slots; differential
free-surface stress equations are evaluated only at interior surface nodes.

For \(N_xN_z\) streamfunction unknowns, the proposed tensor tau allocation is:

1. all \(4N_z\) row slots in the two outermost columns at each side are assigned
   to the two wall conditions evaluated at \(x=-a\) and the two at \(x=a\):
   \(\psi=0\), \(\psi_x=0\);
2. for the remaining \(N_x-4\) columns, the two bottom row slots impose
   \(\psi=0\), \(\psi_z=0\), evaluated at \(z=-H\);
3. for those same columns, the two top row slots impose tangential stress and
   primary dynamics, evaluated at \(z=0\);
4. the remaining \((N_x-4)(N_z-4)\) rows impose the bulk PDE.

The count is exact:

\[
 4N_z+4(N_x-4)+(N_x-4)(N_z-4)=N_xN_z.
\]

Concretely, rows \(i=0,1,N_x-2,N_x-1\) are wall-tau slots.  For
\(2\le i\le N_x-3\), rows \(j=0,1\) are bottom-tau slots,
\(j=N_z-2,N_z-1\) are surface-tau slots, and
\(2\le j\le N_z-3\) are bulk rows.  A row slot's grid index does not change the
location where its boundary operator is evaluated.

There are \(N_x\) additional eta equations.  At \(i=1,\ldots,N_x-2\), impose
\(\dot\eta_i=-(D_x\Psi_s)_i\).  At the two endpoints impose the equivalent wall
consequence \(\dot\eta_0=\dot\eta_{N_x-1}=0\) explicitly.  This retains endpoint
eta in the state and pins it without adding a surface-stress equation at either
corner.  The narrow set of surface collocation points consumed by wall tau rows
shrinks with refinement; its effect belongs to the separate contact-region study.

## 26. Spatial convergence strategy

The Stokes thickness requires approximately 8--12 effective points/DOFs across
\(\delta_s\) for qualification-quality wall data.  This is a proposed starting
target, not a universal guarantee.  It means a smallest wall spacing near

\[
 \delta_s/8=8.92\times10^{-5}\ \mathrm m
 \quad\hbox{to}\quad
 \delta_s/12=5.95\times10^{-5}\ \mathrm m.
\]

With \(\beta_x=3.4\), indicative x grids are:

* medium: \(N_x=321\), \(h_{x,\min}\approx9.5\times10^{-5}\ \mathrm m\);
* fine: \(N_x=385\), \(h_{x,\min}\approx7.9\times10^{-5}\ \mathrm m\);
* superfine: \(N_x=513\), \(h_{x,\min}\approx5.9\times10^{-5}\ \mathrm m\).

At \(H=6\ \mathrm m\), \(\ell_z=2\ \mathrm{mm}\) gives approximately
\(h_{z,\min}=6.8,5.4,4.5\times10^{-5}\ \mathrm m\) for
\(N_z=257,321,385\), respectively.  Candidate paired grids are therefore
321x257, 385x321, and 513x385 (about 83k, 124k, and 198k streamfunction DOFs).
These are estimates for planning, not accepted production grids.

Convergence reports must separate bulk eta/velocity/energy, transition metrics,
and contact-region vorticity/wall shear.  They will also report tangential-stress,
no-slip, endpoint, dynamic-condition, and discrete-energy residuals.  Refinement
must cover wall-normal spacing, surface-normal spacing, stencil width, and corner
row effects, rather than only total DOF.

The free surface and its viscous layer require z clustering, while both side
walls require x clustering.  Extra equal spacing in the deep far field would be
wasted.  Corner zooms need their own resolution evidence.

## 27. Depth convergence strategy

Run the same nondimensional problem at \(H=4,6,8,10\ \mathrm m\), retuning
\(N_z\) or \(\ell_z\) so the near-surface physical spacing is held fixed.  At
fixed diagnostic times compare:

* bulk \(\eta\) and bulk velocity in consistent physical subdomains;
* kinetic and total energy normalized by initial potential energy;
* fitted first-mode frequency;
* velocity/gradient magnitude in a monitor band immediately above, but excluding,
  the imposed bottom line.

The velocity on the artificial line itself is exactly zero by construction and
would be a vacuous decay metric.  Define instead
\(\epsilon_{\rm deep}=\max_{\text{bottom monitor band}}|\nabla\psi|/(kU_g)\).
Before production, require consecutive-depth changes below \(10^{-3}\) in bulk
eta and energy, below \(10^{-4}\) in fitted first frequency, and
\(\epsilon_{\rm deep}<10^{-6}\), or tighten them if spatial/time errors are
smaller.  Failure means increasing \(H\); it does not redefine a finite bottom
as physical.

For context, the inviscid first-mode depth factor is
\(e^{-q_1H}\), approximately \(1.87\times10^{-3}\),
\(8.07\times10^{-5}\), \(3.49\times10^{-6}\), and
\(1.51\times10^{-7}\) over this sequence.  Viscous and corner fields must still
be measured rather than inferred from that estimate.

## 28. Time integration strategy

Backward Euler is L-stable and is the startup/rank/reference method.  Implicit
midpoint is second order but not L-stable, so it is not primary for rapidly
decaying wall modes.

The selected production candidate is the stiffly accurate, L-stable two-stage
SDIRK2 method with \(\gamma=1-1/\sqrt2\).  It is second order, uses the same
matrix \(M-\gamma\Delta t A\) at both stages, and permits one reusable sparse
factorization/preconditioner per fixed \(\Delta t\).  BE will provide consistent
startup experiments around the incompatible corner data and an independent
first-order time-convergence baseline.

Future time studies must distinguish stability from accuracy and report bulk eta,
energy balance, wall vorticity/shear, and constraint drift.  A step resolving the
1.6 s bulk period may still fail the fast grid-scale viscous transient.  On the
fine/superfine grids, sparse assembly is \(O(N)\); direct fill and solve cost must
be prototyped.  If direct factorization is too large, retain the same DAE and use
Kronecker-aware matrix application with a sparse approximate-factorization or
multilevel preconditioner.  Reusing a stage solve over a future 5 s run is the
key cost advantage over a nonlinear method.

## 29. Validation plan

SF-0 implements only lightweight sign and identity checks.  SF-1 and later must
add the following independent gates.

1. **Symbolic/component checks.** Incompressibility, vorticity/bulk sign,
   tangential stress, normal-pressure sign, pressure elimination, psi-only sign,
   pressure-gradient curl, and strain/dissipation identity.
2. **Initial-tilt regression.** With \(\psi=0,\eta=kx\), require nonzero interior
   \(\psi_{zt}=-gk\).  A permanently zero state must fail.
3. **Pinned endpoints.** Require \(\dot\eta(\pm a)=0\) and preservation of
   \(\eta(\pm a)=\pm ka\) to time-solver tolerance.
4. **Linearity.** Runs at \(k\) and \(2k\) must give factor 2 in state, velocity,
   and pressure and factor 4 in every energy, with tolerances tied to linear-solve
   error.
5. **Energy.** Check \(E\), \(D\ge0\), and the integrated balance
   \(E(t)+\int_0^tD\,dt=E(0)\), separating time-discretization dissipation.
6. **MMS.** A future forced benchmark may use, in nondimensional form,
   \(\psi_m=e^{t+z}(1-x^2)^2(1+z)\) and
   \(\eta_m=e^t(x+x^3/5)\).  Compute
   \[
   f_\Omega=(\Delta\psi_m)_t-\nu\Delta^2\psi_m,
   \]
   \[
   f_k=\eta_{m,t}+\psi_{m,x}|_0,\quad
   f_s=(\psi_{m,zz}-\psi_{m,xx})|_0,
   \]
   \[
   f_d=[\psi_{m,zt}-\nu(3\psi_{m,xxz}+\psi_{m,zzz})
       +g\eta_{m,x}]_{z=0},
   \]
   plus evaluated wall/bottom data.  Forced boundaries deliberately avoid
   confusing MMS convergence with restrictive homogeneous corner physics.
7. **Inviscid bulk reference.** Compare small-\(\nu\) early/bulk eta frequency
   and modal content with \(q_n,\omega_n,b_n\), never as an exact no-slip answer.
8. **Independent discretizations.** Use small spectral/collocation systems and,
   at moderate viscosity and the same \(a,H,\nu,k,g\), compare the older
   velocity-pressure solver's eta, velocity, and kinetic energy.  It is an
   independent cross-formulation check, not a dependency.
9. **Never use STEP 3 CHNS as the exact reference.** CHNS has a diffuse
   interface, surface tension, slip, and different contact-line physics.

Future diagnostic APIs are planned for \(\eta\), \(u=\psi_z\), \(w=-\psi_x\),
\(\omega=-\Delta\psi\), kinetic/potential/total energy, dissipation, wall shear,
maximum no-slip error, endpoint eta error, deep-decay error, and separate
bulk/transition/contact metrics.

Future visualization, not SF-0, will show true-scale eta, bulk velocity,
vorticity/wall layers, a separate near-wall zoom, and energy versus time.  It
must never draw nonzero arrows exactly on a no-slip wall.

## 30. Known physical/model limitations

* The theory is linear about \(z=0\); it is not a finite-angle free-boundary
  solution.
* Strict no slip pins the endpoints and produces incompatible startup traces for
  the tilted surface.  No moving-contact-line physics is present.
* Surface tension is exactly zero.  It might regularize some short surface scales,
  but capillarity is outside this problem definition and no STEP 3 term is imported.
* A truncated bottom is a numerical far-field surrogate only.
* At \(\nu=10^{-6}\), bulk convergence does not imply wall-vorticity or corner
  convergence.
* Long-time bulk flattening does not imply uniform endpoint flattening.

## 31. Proposed SF-1 implementation plan

SF-1 should remain a qualification stage before any 5 s production calculation:

1. implement nondimensional configuration/scaling and the mapped tensor grid in
   a new streamfunction-only package, without touching CHNS;
2. generate and audit 9-point physical-coordinate derivative matrices on tiny
   and medium grids, including polynomial exactness and conditioning;
3. assemble the exact tau-row DAE above on tiny truncated grids and test row
   count, rank, nullspace/gauge removal, algebraic constraints, and consistent
   initialization for \(\eta=kx\);
4. verify that the corner policy pins eta and retains the nonzero interior
   acceleration; report the excluded/transition surface rows explicitly;
5. implement BE and SDIRK2 for manufactured forced problems only, with reusable
   stage solves and time-order tests;
6. implement MMS residuals, primitive-pressure curl checks, discrete energy and
   dissipation, and regional bulk/contact norms;
7. compare one tiny case against an independent spectral collocation prototype;
8. run moderate-viscosity, modest-grid spatial/time/depth feasibility studies;
9. measure CSR nonzeros, factorization fill, memory, and solve time, then decide
   whether direct reuse or a Kronecker-aware iterative solve is needed;
10. present SF-1 evidence for review before attempting \(\nu=10^{-6}\), the
    medium/fine/superfine grids, an \(H\) sweep, or any production animation.

### SF-0 decision

The single recommended path is: **semi-infinite depth approximated by a declared
\(H=4,6,8,10\) truncation study; sparse 9-point high-order nonuniform tensor-grid
FD; the explicit \(\psi+\eta\) singular-mass DAE with wall-dominant corner tau
rows; BE for startup/reference and L-stable SDIRK2 for production.**  This choice
resolves the double no-slip and free-surface derivative conditions without a
dense global discretization, targets the 0.714 mm layers with roughly 83k--198k
streamfunction unknowns, isolates the contact singularity, retains exact initial
tilt/pinning information, and permits stage-matrix reuse over a later 5 s run.
It remains conditional on SF-1 rank, MMS, energy, regional convergence, and cost
gates.
