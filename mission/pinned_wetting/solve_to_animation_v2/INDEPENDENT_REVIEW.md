# Independent numerical review — PW2

This review is separate from the production assembly and renderer. It does not
declare a five-second result or mesh-resolution PASS: those require native
audits and separate full-curve temporal/spatial comparisons. Historical linear
and symmetric-slip results are not validation data for this physical case.

## Locked physical problem

LEFT is impermeable with tangential Navier slip, b = 0.50 m; RIGHT and BOTTOM
are full no-slip. The material right surface endpoint is exactly P2. Gravity,
viscosity, sigma = 0, the initial straight inclined boundary and initial zero
velocity are the fixed contract. The retained coating is a separate one-way
lower-connected state interval on each wall. Its thickness, volume, inertia,
separate dynamics and feedback are not resolved. Navier slip alone does not
imply this irreversible law.

New left markers must be confirmed local peaks and equal the maximum of every
accepted left endpoint through their source state, including an earlier
unmarked plateau. The discrete maximum selects one of the native values, so
this global-record decision needs no interpolation or arithmetic tolerance.
Beating only the previous plotted marker is not sufficient.

## Equations and discrete geometry

The isoparametric P2 mesh follows the same full P2 velocity as the fluid during
each physical step. Thus the relative ALE convection is zero, while the
material derivative retains Eulerian nonlinear convection through the moving
map. This is not fixed-domain Stokes relabelled as nonlinear Navier–Stokes.
After a rezoning, the next physical material step starts from the explicitly
stored transferred geometry and velocity, not from an unrecorded adjustment.

For X1-X0 = dt * vm and a two-dimensional element map F, determinant
polarization gives exactly

    J1-J0 = dt * Jm * div(vm),
    Fm = (F0+F1)/2,  vm = (v0+v1)/2.

Keeping all pressure test functions includes the constant test function. The
midpoint weak incompressibility equation consequently controls global volume
through this geometric conservation identity. It does not imply pointwise
incompressibility, nor exactly B(Xn)vn = 0 at accepted endpoint states. The
independent audit reports both the collocation residual and accepted-endpoint
residual, plus strong divergence and per-step local dilation. Their numerical
effects need refinement; none is concealed by calling a weak residual a
pointwise zero.

The initial acceleration and pressure must satisfy the coupled weak momentum
and incompressibility system. A prescribed hydrostatic pressure on the tilted
release is not an admissible substitute. No pressure gauge DOF is silently
removed from this traction-boundary problem.

## Boundary signs and pressure variable

With outward normal on the left, the tangential traction is -nu/b times the
tangential velocity. Its weak contribution is positive nu/b times the LEFT
boundary mass matrix, and its physical dissipation is nonnegative. RIGHT has
no Robin contribution, including in the unreduced matrix; all its velocity
DOFs are fixed. BOTTOM also has both velocity components fixed.

The optional recovery uses pi = p + g*z, where p is the physical pressure.
Then the body gravity term is removed and the pi-stress has actual free-surface
traction -g*z*n. This is an explicit affine shift of the pressure approximation
space, not an imposed tilted equilibrium. The solid normal test velocity is
zero and the LEFT tangential law is unchanged. The native pressure variable
must identify this choice; a stored pi must not be described as physical p.

The optional +0.5*div(vm)*vm term is a consistent geometric/skew numerical
formulation: it vanishes for an exactly incompressible continuum field. It is
assembled during the solve and its signed work is independently integrated.
It is not an after-the-fact residual inserted to force an energy PASS.

## Independent energy accounting

Let ja=(abs(J1)+abs(J0))/2, jm=abs(Jm), dj=abs(J1)-abs(J0),
dv=v1-v0 and zm=(z1+z0)/2. Independently integrated raw geometric works are

    Wk = integral[(ja-jm)*vm.dot(dv)
                  + dj*|vm|^2/2 + dj*|dv|^2/8] dXref,
    Wp = g*integral[(ja-jm)*(z1-z0) + dj*zm] dXref.

The skew/geometric matrix work is subtracted from Wk. The work of the
surface-gravity load minus the original body-gravity load is added to Wp when
the hydrostatic split is used. Remap work is the independently measured energy
of the transferred field minus the energy before transfer. None of these
signed numerical terms is called physical heat. Bulk symmetric-gradient and
LEFT Navier losses remain separate nonnegative quantities.

Both the raw continuous-energy defect and the fully split algebraic budget
are checked. Passing the latter alone cannot pass the former's 5% of E0 gate.
This matters particularly because a fitted closing term could trivially make
an algebraic budget vanish without establishing numerical accuracy.

## Remapping and geometry validity

The transfer audit independently evaluates the old curved P2 map and old
velocity at the stored source reference points; points outside old cells are
rejected. It rebuilds the L2 projection system and the incompressibility and
integrated-momentum constraints, including their stored multipliers. The full
ordered surface P2 coordinates must be bitwise/unambiguously preserved to the
declared tolerance, with no wall-endpoint clamp or branch removal.

The locator recovery keeps these transfer equations and acceptance tolerances.
For each P2 edge, the quadratic Bernstein control point is
`2*midpoint-(vertex_i+vertex_j)/2`. Together with the three vertices, these six
controls enclose the complete curved triangle because the reference Bernstein
basis is nonnegative and sums to one. Their bounding boxes are therefore
complete candidate enclosures, unlike centroid proximity or nodal-coordinate
bounds. Four interior Newton seeds and residual-decreasing damping improve
inverse search, but do not constitute a theorem that every inverse will be
found. Limiting a Newton correction is not clamping the target, result or
interface. Final reference containment and physical residual tolerances remain
mandatory; unresolved targets are rejected. Independent pure tests cover an
interior image point outside the nodal box, affine-field reproduction and an
outside physical target still inside the Bernstein box. They do not replace
per-transfer native validation.

Exact quadratic determinant minima are checked over each reference triangle,
including interior and edge stationary points. A positive determinant is
necessary but not sufficient for a resolved thin or steep boundary region.
The full quadratic interface is also checked for additional solid contact.
Exact rational polynomials formed from the stored IEEE coordinates test
signed-distance extrema against LEFT, RIGHT and BOTTOM; only the known P1/P2
endpoint roots are exempt. An interior tangency can invalidate the assumed
one-contact topology without crossing a wall, and must not be missed by a
bounding-box test. A candidate separated from a wall by no more than eight
input-coordinate ULPs is labelled precision-unresolved, not asserted to be a
physical intersection. This is a floating-representation check, not a minimum
physical film thickness; a representable positive gap above that numerical
scale is not rejected merely because it is small.

The optional near-right diagnostic compares the actual inward normal ray gap
to the neighboring curved cell's exact projected normal width and retains its
geometric witness. No microscopic thickness or arbitrary film profile is
imposed. In particular, the gap intrinsically approaches zero at P2, so this
diagnostic alone is not a hard failure criterion.

Full h refinement remains essential. The comparison operates on the complete
ordered parametric P2 boundary, not sorted x coordinates or a single eta(x).
At common physical times it brackets two-sided geometric distance using
controlled P2 chord errors and a source-sampling bound. Time interpolation is
explicitly disclosed and is not claimed to have a rigorous error bound. The
method reports contacts, retained height, peak heights/times, P2 and energy
observables separately. Two refinement levels do not establish formal order.

## Independence boundary and evidence

`free_boundary_verify.py` independently reconstructs the P2/P1 operators,
physical constraints, boundary signs, midpoint residuals, GCL, coating/marker
history, remap transfer and energy. It does not import production operators,
diagnostics, summaries or the renderer. Its shared numerical substrate is
scikit-fem basis/geometry and SciPy. Boundary integration uses reference
triangle edges and transformed normals independently of production's ordered
native P2 edge-polynomial integration.

For accepted endpoints and material midpoint stages, the audit evaluates
actions on the stored fields through independently written LinearForms rather
than assembling unused sparse mass, strain and geometric-mass matrices.
Bulk power is integrated as `(nu/2)*|grad(v)+grad(v)^T|^2`; LEFT power as
`(nu/b)*w^2` on the actual LEFT wall; signed geometric power as
`0.5*div(v)*|v|^2`. No term is inferred by closing an energy residual.
The full unreduced initial LEFT-only Robin matrix witness, startup matrices
and every mixed remap projection matrix remain explicitly reconstructed.
Tiny straight and curved P2 algebra tests cover both hydrostatic-pressure and
skew-divergence variants, including nonzero values on constrained DOFs: every
velocity row, the weak-divergence vector and the split powers agree with the
independent full sparse reconstruction to roundoff. This changes only audit
evaluation cost, not its equations, tolerances, state coverage or physical
qualification. Its source digest invalidates earlier cached audit reports.

`free_boundary_compare.py` is the independent full-boundary dt/h comparator;
`free_boundary_geometry_audit.py` provides the non-gating gap witness. Negative
tests mutate native equations, physical identities, constraints, unreduced
RIGHT Robin terms, GCL, remaps, coating history, records and exports. A test
suite PASS validates those checks, not the uncomputed physical trajectory.

Each native audit JSON identifies the native digest, verifier digest, fully
audited time interval, qualification failures and whether it actually reaches
5 seconds. A short-run PASS is explicitly not release qualification. A
third-party geometry exception is a failed/incomplete audit, never a CFD PASS
or automatically a proof that the physical problem is impossible.

For the generic ALE distinction and the importance of consistent energy and
weak-divergence treatment, a primary cross-check is Fehn et al.,
[High-order arbitrary Lagrangian–Eulerian discontinuous Galerkin methods for
the incompressible Navier–Stokes equations](https://arxiv.org/html/2003.07166).
That paper is not a convergence theorem for this particular continuous-P2,
free-boundary, symmetric-gradient formulation. The identities above were
derived directly for the implemented two-dimensional midpoint map.
