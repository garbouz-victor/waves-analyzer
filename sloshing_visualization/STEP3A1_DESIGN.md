# STEP 3A.1 — policy fixed before new calculations

Base: `15a9e07c496c757cf40d35cd10fbc721afb9805a`. The continuum equations,
constitutive parameters and wall law are unchanged. No tank/film run is authorized.
New outputs belong to `validation_results/step3a1`; STEP 3 records are historical.

## Energy and the old first trapezoid

The recorded CH powers are 91221.0739273854 and .1288992508151518 W/m at
t=0 and .01 s. Their first trapezoid is **456.1060141331811 J/m**, versus
456.10918847807 J/m over the complete .5 s run. The actual energy decrease is
.09424840859733032-.08851376023738629 = **.00573464835994403 J/m**.
This identifies an invalid time integral, not physical heat, and does not prove
that the underlying initial PDE transient was resolved. No sample will be removed.

Named policies, before rerunning:

- mass error / domain area <=1e-8;
- largest positive physical-energy step <=1e-9 J/m (unchanged);
- legacy cumulative budget defect / |E(0)| retained for continuity;
- max defect / initial positive component sum <=.01;
- defect / actual energy change <=.05 once the change exceeds
  max(1e-10 J/m, 1e-6 times the initial positive scale); report N/A below that floor;
- all energy/power/cumulative values must be finite; negative powers are rejected
  beyond 1e-12 W/m roundoff;
- trapezoidal integration uses every accepted interval, including t=0;
- single-interval dissipation fraction >.5 produces a temporal-sampling warning.

Single-history closure is not temporal convergence. The latter needs separately
measured refinement, converging integrated powers and end states. A missing
refinement record cannot authorize the repaired-core gate.

Reference-handling self-review: the legacy ratio is not a second mandatory gate
when the arbitrary initial total reference cancels to zero. A regression with
E_kin+E_interface+E_gravity+E_wall=0 but a finite positive component scale exposed
that false-failure risk. The hard 1% test uses the positive scale. This definition
correction changes none of the measured verdicts here: both recorded contact
initial states have |E(0)| equal to the positive scale. No threshold was relaxed.

## Interface resolution

CG1: exact vertex range. CG2: convert the six reference nodal values to quadratic
Bernstein coefficients: b_i=f_i, b_ij=2 f_mid-(f_i+f_j)/2. Their convex hull bounds
the polynomial; overlap with [-.9,.9] is conservative (false positives possible,
not centroid-based false negatives). No field coefficients are modified.

For active cells sample gradients at vertices, edge midpoints and centroid; use
the largest projected triangle width across these directions. This is a robust
sampled normal-spacing estimate, not an exact maximum over all directions. For
vanishing gradients use the cell diameter conservatively. Keep **8 cells** across
the nominal phi=-.9..+.9 width 4.164065537 epsilon. Save cellwise quantiles and the
worst element location. A later static corridor must include several epsilon
and an explicit geometric buffer, not just the previous crossing trajectory.

## Settling in physical time

The last **.10 s** window is fixed, independent of dt. Fit time-weighted linear
trends to theta and both crossings, with the window boundary interpolated exactly.
Require angle rate <=.5 degree/s, each contact speed <=1e-4 m/s, peak kinetic
energy / positive initial scale <=1e-8, peak sampled speed <=1e-4 m/s and peak
chemical-potential relative variation <=.01. These are strict benchmark settling
criteria, not universal wetting constants. Too short a history is not qualified.
Retain angle error <3 degrees and fit-window spread <1 degree.

## Time stepping and restart

BE is permitted on explicitly declared startup intervals. Standard BDF2 is used
only after two accepted states separated by main_dt exist. At a spacing change,
take one BE main_dt interval from the final startup state; the following BDF2
step uses this state and the final startup state, now exactly main_dt apart.
Never insert the last microstep as an unequal-spaced BDF2 history. Checkpoint the
schedule, current phase, time, dt and history spacing. Actual time intervals,
not config.dt, drive power quadrature.

For the compatible-angle baseline, declare BE blocks dt=1e-8 to 1e-6 s,
1e-7 to 1e-5 s, 1e-6 to 1e-4 s, then main_dt=1e-5 to .001 s. This is 370
accepted intervals, not variable-step BDF2. It resolves much shorter times than
the historical .01 s first trapezoid. The schedule is a candidate, not an assertion
that the pulse is resolved. Compare pure BE and BDF2; a factor-two temporal series
halves every block dt without moving its boundary. Stop at a failed compatible
gate and diagnose before moving to the deliberately incompatible case.

## Ordered gates

1. Repair energy, interface and settling logic; synthetic and FEM regressions.
2. Ordinary/container tests; compatible-angle BE/BDF energy baseline. Stop on failure.
3. Laplace radii .25/.30 temporal/energy qualification. Stop on unexplained failure.
4. Cheap incompatible startup scale study; static corridor dry run.
5. Resolved startup refinement, then .5 s contact continuation and settling.

No theta sweep, epsilon/M/L_s study, falling film or bridge is in this iteration.
Even full success permits only “CORE CONTACT MODEL VALIDATED THROUGH REPAIRED
ENERGY/RESOLUTION GATES — REMAINING STEP 3A BENCHMARKS REQUIRED”.
