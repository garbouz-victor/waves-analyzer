# STEP 3A.1 — repaired validation logic, failed startup qualification

**MODEL NOT YET VALIDATED**

Base: `15a9e07c496c757cf40d35cd10fbc721afb9805a`. Execution stopped at PHASE 6,
the geometrically compatible-angle energy baseline. This is not completion of
all STEP 3A.1 scientific gates. There was no new Laplace time series, incompatible
90→60 run, full moving-interface corridor run, tank, film or water–air calculation.

## 1. Scope and unchanged physics

The three automatic-validation defects are repaired: energy closure is mandatory,
transition cells are not selected solely by their centroid, and settling is
measured over physical time. The continuum weak form, capillary-force sign,
pressure transform, material interpolation, free-energy normalization, wetting
law and initial-field formulas were not edited. No finite wall relaxation,
clipping, smoothing, energy fitting or modified dissipation was introduced.

The new baseline retains matched rho_l=rho_g=1 kg/m³, dynamic viscosities 1 Pa s,
sigma=.1 N/m, epsilon=.025 m, M=1 m²/(Pa s), L_s=.01 m, theta_e=60°, g=a_x=0.
It is a demonstration benchmark, not water. Its initial zero-contour angle is
60° instead of the historical 90°, so these are **different initial conditions**,
not two levels of a temporal convergence series.

## 2. Exact origin of the old 456 J/m

The historical data are unchanged. Recomputing from their actual power samples:

```
D_CH(0)    = 91221.0739273854 W/m
D_CH(.01)  =     0.1288992508151518 W/m
first trapezoid = .01/2 * (D_CH(0)+D_CH(.01))
                = 456.1060141331811 J/m
full CH integral to .5 s = 456.10918847807 J/m
first interval fraction = 99.9993040383818%

E(0)  = .09424840859733032 J/m
E(.5) = .08851376023738629 J/m
actual energy decrease = .00573464835994403 J/m
```

Almost all of the unqualified integral comes from one unresolved trapezoid.
It is **not physical heat**. This calculation diagnoses invalid power integration,
but does not establish that the underlying initial PDE trajectory was resolved.
The new gate rejects this record even though its E is monotonically decreasing.

## 3. Common energy gate

`energy_validation.py` independently integrates every accepted interval,
including the initial sample, and checks

`R(t)=E(t)-E(0)+integral(D_visc+D_CH+D_slip) dt`.

Named policy values were recorded before rerunning in
[STEP3A1_DESIGN.md](STEP3A1_DESIGN.md) and `validation_policy.py`:

| Check | Limit |
|---|---:|
| Mass error / domain area | 1e-8 |
| Largest positive physical-energy step | 1e-9 J/m |
| Max abs(R) / positive initial component sum | .01 |
| Abs(R) / abs(E(t)-E(0)), when applicable | .05 |
| Minimum significant energy change | max(1e-10 J/m, 1e-6 E_scale_initial) |
| Single-interval fraction warning | >.5 of a component's total integral |

Positive scale is |K|+|G|+E_interface+|E_wall| for the nonnegative interface
energy of this model; its peak is also reported. Legacy abs(R)/abs(E(0)) is
retained. A self-review regression found that making the legacy ratio a second
mandatory gate would falsely reject an exactly balanced history when arbitrary
energy references cancel E(0). It is therefore diagnostic, not the primary gate.
This definition correction changes **none** of the measured verdicts below:
both initial contact energies equal their positive component sums. Neither the
1% nor the 5% threshold was relaxed.

Below the energy-change floor the relative-to-change metric is N/A, not a hidden
absolute defect. NaN/inf raises an explicit error. Saved cumulative integrals
are checked against independently re-integrated powers. Negative dissipation,
inconsistent energy components, monotone E with bogus D, or missing required
metrics cannot silently pass. Per-interval delta_E, power integral and defect
are saved along with the time of the worst local defect.

Contact, Laplace and flat-interface checks now require this gate. Numerical
`status=complete` is separate from `qualification_status`. A single-history
closure pass is not temporal convergence; the top-level summary requires the
corresponding measured series record, not a hardcoded eternal False or a bare
execution-complete marker.

## 4. Conservative active cells and robust normal spacing

For CG1, vertex min/max is the exact polynomial range. For CG2, with vertices
f_i and edge-midpoint values f_ij, the Bernstein coefficients are

`b_i=f_i; b_ij=2 f_ij-(f_i+f_j)/2`.

The polynomial lies in their convex hull. A cell is active when that bound
overlaps [-.9,.9]. This permits conservative false positives but eliminates the
old centroid-based false negatives. Tests include edge/vertex crossings with a
pure centroid and curved P2 contours.

grad(phi) is sampled at three vertices, three edge midpoints and the centroid.
The worst projected triangle width across these directions is used. A constant
transition cell falls back to its diameter; it is not dropped. This is a robust
sampled normal-width indicator, not a theorem about the exact directional maximum
or a measured molecular interface thickness. The nominal width remains
4.164065537 epsilon and the policy remains **at least 8 cells**.

Min/p05/p50/p95, active-cell count, worst triangle bounds/centroid and worst time
are stored. Rotation tests cover rigid covariance and continuous response when
the interface rotates relative to a fixed triangle.

### Read-only historical final-checkpoint audit

| Indicator at old t=.5 s | Minimum cells | Worst h_normal [m] |
|---|---:|---:|
| Historical centroid method | 7.488279 | .0139019442 |
| New Bernstein / seven-gradient method | **3.865905** | .0269281396 |

The worst active triangle has centroid (-.4833333,.00555556) m, bounds
x=[-.5,-.4666667], z=[0,.0166667]. There are 3,244 active cells; p05/p50/p95
counts are 8.839711 / 12.072036 / 16.546513. The previous minimum materially
understated the worst under-resolution. Only the final checkpoint was remeasured:
this is **not** a new all-time history or a rerun of the old PDE.

Configurable static `refinement_boxes` are implemented and tested. The full
spreading corridor/dry run is deferred by the PHASE 6 failure. A future corridor
must cover the entire expected droplet envelope plus several epsilon and a
geometric buffer; no box was tuned to just end at the old ±.4093 crossing.

## 5. BE startup and correct BDF2 restart

The schedule supports explicitly declared BE blocks, followed by a constant main
step. At a spacing change the first main interval is BE. After it, the current
and previous accepted states are exactly main_dt apart, so standard BDF2 can
start. Unequal history spacing is checked and rejected; old microstates are never
treated as main_dt-spaced history. The scalar y'=-y test verifies the schedule,
post-switch convergence and the history rule.

Checkpoints preserve startup phase, actual dt, history spacing, the full declared
schedule, both BDF states and cumulative accounting. Tiny CHNS restart tests
cover the middle/end of startup, the main BE interval and BDF2. Legacy v1 fixed-dt
checkpoints have explicit read-only-compatible loading; no new startup can be
silently attached to them.

The measured baseline used **pure BE**, including the main block:

| Time block [s] | dt [s] | Intervals |
|---|---:|---:|
| 0…1e-6 | 1e-8 | 100 |
| 1e-6…1e-5 | 1e-7 | 90 |
| 1e-5…1e-4 | 1e-6 | 90 |
| 1e-4…1e-3 | 1e-5 | 90 |

There are 371 power/energy samples including t=0. Trapezoids use the actual
interval length. No variable-step BDF2 formula is used. The BE/BDF2 physical-run
comparison was not executed after the first compatible candidate failed.

## 6. Startup table — measured, not an asserted convergence series

| Initial dt [s] | Case / unknowns | Min cells | D_CH(0) [W/m] | First CH integral [J/m] | Total CH integral [J/m] | Delta E [J/m] | Max abs budget defect [J/m] | Relative defect to initial scale |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| .01 | Historical 90→60 / 82,769 | 7.488 old indicator | 91221.0739 | 456.106014 | 456.109188 | -.00573464836 | 456.103467 | 4839.375789 |
| 1e-8 | New 60→60 / 58,625 | **8.833317** | 577.494578 | 2.96029596e-6 | 2.34328978e-4 | -2.30764854e-4 | **3.56418174e-6** | **5.33152001e-5** |

The old and new rows differ in geometry, liquid volume, mesh location, source and
time schedule. They are shown for audit, **not** as dt-convergence evidence.

New mesh: 6,798 triangles; velocity 27,566, pressure 3,493, phi and mu 13,783 DOFs
each. Rough CSR memory estimate 126.6 MB excludes LU fill. Runtime **2464.63 s
(41.08 min)**. Minimum nominal transition count 8.833317 over all accepted times;
worst recorded triangle centroid (.2388889,.0861111) m at 7.1e-6 s. Mass error /
area <=6.10623e-14, physical-energy growth zero, no material-coefficient failure.

## 7. Why the new baseline still fails

```
E(0)      = .06685113685786184 J/m
E(.001)   = .06662037200421221 J/m
D_CH(0)   = 577.4945782537685 W/m
D_CH(1e-8)=  14.56461356656488 W/m

first actual energy loss = 3.3436420139587586e-7 J/m
first CH trapezoid       = 2.960295959101667e-6 J/m
first local defect      = 2.6259317577060312e-6 J/m
max relative defect to actual energy change = 7.853507482988639 (785.35%)
```

The first actual loss exceeds the predeclared significance floor 6.68511e-8 J/m,
so its 785.35% defect cannot be labelled N/A. The 1% initial-scale gate passes:
the maximum is only **.00533152%**. Even the final relative-to-change defect is
only 1.544508%, but the **all-time** 5% gate fails at startup. Dropping that first
interval or weakening the policy would hide precisely the failure being tested.

The largest single interval contributes 1.7561% of the full CH integral, so the
50% sampling heuristic does not warn here. This is another useful negative
result: that heuristic alone does not prove the very early pulse is resolved.
No startup_dt is qualified, and no convergence of integrated D_CH or initial/end
fields has yet been established.

Offline uniform-block Simpson gives CH integral 2.33385418e-4 J/m, versus
2.34328978e-4 from trapezoids (difference 9.43560e-7 J/m). For the old run the
two values are 304.074376 and 456.109188 J/m. Simpson is **not** substituted into
the gate; quadrature agreement cannot replace temporal refinement.

## 8. Initial-data root-cause audit

A circular tanh with zero-contour angle theta_i uses
`r=sqrt(x²+(z+R cos(theta_i))²)` and `phi=tanh((R-r)/(sqrt(2) epsilon))`.
At the bottom wall, direct differentiation gives

`L = lambda epsilon partial_n phi + f_wall'(phi)`

`  = (3 sigma/4) [R cos(theta_i)/r - cos(theta_e)] (1-phi²)`.

Even theta_i=theta_e=60° makes L zero only at r=R, not throughout the diffuse
wall trace. This is checked symbolically. **Geometric angle compatibility is not
a prepared variational equilibrium.** The analytic max sampled |L| is .00218343
N/m for the equal-angle geometry, versus .03749999 N/m for the 90→60 mismatch;
L2 wall norms are .000859858 and .011514448 N/sqrt(m), respectively.

This explains why simply matching the apparent angles did not eliminate a stiff
initial wall variation. The consistent weak mu initialization represents that
variation, plus bulk/interpolation contributions; it was not arbitrarily capped.
These checks do not assign every part of the energy defect to quadrature alone,
nor prove the remaining weak/time discretization error is negligible. A genuinely
prepared equilibrium and/or much more resolved initial-time experiment is needed
before proceeding. Changing wall relaxation or the free-energy law is not a fix.

## 9. Settling is now in physical time

Over the last fixed .10 s window, piecewise-linear signals are fit with exact
time-weighted least-squares integrals. This avoids weighting dense microsteps more
heavily merely because there are more samples. Limits are .5 degree/s for angle,
1e-4 m/s for each contact speed, 1e-8 for peak K/E_scale, 1e-4 m/s for sampled
fluid speed and .01 for relative mu variation. Missing/short histories cannot pass.

Re-auditing the **historical** t=.4… .5 s window:

| Quantity | Measured | New gate |
|---|---:|---|
| dtheta/dt | -2.164234 degree/s | failed |
| Left contact speed | -.006482054 m/s | failed |
| Right contact speed | .006470809 m/s | failed |
| Peak K/E_scale | 2.67896e-9 | passed |
| Peak sampled fluid speed | 4.94228e-5 m/s | passed |
| Relative mu variation | absent in old history | not qualified |

The formerly acceptable .0132 degree per final step was not settling evidence.
CH contact motion need not equal fluid advection speed; both are now monitored.
There is **no new** .5 s settling result because the earlier energy gate failed.

## 10. Stopped phases and current verdict

| Phase | Result |
|---|---|
| Energy/resolution/settling repairs and cheap tests | implemented and tested |
| Compatible-angle baseline | mass/resolution/growth pass; early energy closure FAIL |
| BE vs BDF2 compatible physical series | deferred after failure |
| R=.25/.30 Laplace temporal qualification | not run in this iteration |
| Cheap incompatible startup and covering-corridor dry run | not run |
| Resolved 90→60 microstep refinement and .5 s continuation | not run |
| New physical-time settling qualification | not run |

Old Laplace pressure-only successes and shrinking BDF2 energy overshoots are not
new repaired-energy passes. A full moving-contact minimum cell count, converged
startup_dt, integrated-power convergence and settled contact speeds remain unknown.
Static box-refinement support exists, but no qualified spreading corridor is claimed.

**MODEL NOT YET VALIDATED**. Stop before subsequent physical benchmarks. The next
task is to review the variational compatibility of initial data and resolve the
initial CH pulse with measured time convergence, retaining all t=0+ energy/power
work. No threshold change, initial-energy reset, fake heat or new wall physics is
justified by the current result. Even eventual success here would still leave the
other STEP 3A angle/epsilon/slip/mobility/falling-film benchmarks outstanding.

## 11. Provenance, storage, tests and artifacts

Every new run records base commit, module-source hash, config fingerprint, Docker
image ID, MPI ranks, complete dt schedule and mesh fingerprint. For the measured
baseline, module SHA-256 is
`16eca3b39c02b21f55025c7101809b6ca1e613726bb90dd038958c9d799f67b0`, and Docker image
is `sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8`.
The exact module sources were captured in a small archive after verifying their
hash against the run. Subsequent streaming/error-handling/read-only/summary and
reference-cancellation guard fixes have a different source hash; they are not
retrospectively asserted to have generated this run. No convergence series mixes
these revisions silently.

The completed run has all 371 accepted power samples. Subsequent code also
flushes each accepted row during execution, preserving evidence on interruptions;
an injected-failure regression verifies this. Historical result trees are protected
from write/re-analysis in the updated pipeline. No old HDF5 or result was overwritten.

See [test results](validation_results/step3a1/TEST_RESULTS.md) for final counts.
Local test success is not a remote GitHub Actions result for uncommitted changes.
No commit is created because STEP 3A.1's scientific acceptance has not passed.

Artifacts:

- [Measured summary](validation_results/step3a1/summary.json) and [gate CSV](validation_results/step3a1/summary.csv)
- [Startup table](validation_results/step3a1/startup_table.csv)
- [Linear-time power plot](validation_results/step3a1/startup_power_linear.pdf)
- [Log-time power plot](validation_results/step3a1/startup_power_logtime.pdf)
- [Resolution history / worst location](validation_results/step3a1/resolution_history.pdf)
- [Historical physical-time settling audit](validation_results/step3a1/settling_history.pdf)
- [Quadrature comparison](validation_results/step3a1/quadrature_comparison.json)
- [Initial wall audit](validation_results/step3a1/initial_wall_audit.json)
- [Historical final resolution recheck](validation_results/step3a1/mesh_resolution_tests/historical_final.json)

PNG counterparts are adjacent. Raw FEM HDF5/checkpoints, nonlinear logs and caches
are ignored. No animation was generated.
