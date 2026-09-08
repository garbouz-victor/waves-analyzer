# STEP 3A.8 — frozen full-CHNS phase-rate bridge protocol

Base HEAD: 1cc78a1f4108876a1610570b48138381f29a2484.
Overall MODEL NOT YET VALIDATED. No moving-contact, film, tank, BDF2,
unequal-density or body-force calculations belong to this iteration.

## Purpose and inherited inputs

STEP3A7's three COMPLETE isolated trajectories are read-only inputs, not
rerun. Their archive/history/identity hashes and qualified convergence are
verified before production execution. Full L0 minus isolated L0 measures
coupling at the SAME discretization; isolated L0 minus L2 measures the
previously qualified temporal discrepancy. Neither difference alone proves
the full continuum solution's error.

Required initial phi SHA:
6e59cd0bd3ee5c421af9694508d439e8e05a565de329ed38a14af39444d4e681.
Psi SHA: 831db9041c11ce04ab8974f34e841b62b56bd77474290e3c978c9b5f5af724d0.
Equilibrium fingerprint:
df3a6176bc442190c15ca8e4b6faa2c885dbbe89d5126196148d92b1552809fc.
Parent plan SHA:
5d291468a3520da5267a8d76b8015d0709f5e28b44fbaa5f8a7c6b264664e2b8.
Expanded L0 SHA:
21545943a51bdf364a714fde6ad07f3c6879ad48ac792f0aed901afad6776e31.
Read and verify all 1068 intervals, endpoints and common field indices from
STEP3A7. P2/P1/P2/P2, mesh, quad12, delta=1e-3, parameters and coefficients
remain unchanged. Reproduce D_CH(0)=0.11400313905309394 to relative 1e-10,
cells=8.328131075834218 to relative 1e-12 and initial CH excess
5.542888467657825e-8 J/m to absolute 1e-14. Initial u=pi=0.

## Separate algebraic bridge

Keep CHNSSolver and weak_form.residual unchanged. New validation-only BE
unknowns are (u,pi,r,mu), with current phase phi_old+dt*r everywhere.
Phase time term is r, never an absolute subtraction. Momentum retains the
existing matched-density BE time, skew convection, viscosity, pressure,
capillary and Navier terms; continuity, chemical and wall variation are
unchanged. Reject configs with unmatched density, forcing, wrong FE/quad,
or any request for a non-BE step. The Jacobian change is COLUMNS
J_rate=J_abs diag(I,I,dt I,I), not residual row scaling.

Tiny arbitrary-field residual tests cover u!=0, nonconstant mu and r!=0,
all four blocks and the wall term separately. Relative residual/Jacobian
transformation limits: 1e-11. Central FD h=1e-2,5e-3,2.5e-3,1.25e-3 checks
second-order truncation decrease; a final h=1e-5 action must have relative
defect <=1e-7. Same velocity constraints and same single pressure gauge;
no phase-rate Dirichlet condition. Tests precede source freeze.

CI root comparison uses a separate 4x4 toy, phi=0.2+1e-4*cos(2x)*cos(3z),
epsilon=1, theta=90, dt=1e-3 and strict P2 Newton controls. This is not the
production bump or its moderate dt; the coarse mass matrix amplifies
absolute-subtraction roundoff at dt=1e-6. Arbitrary-field bridge tests use
theta=60 and nonzero velocity/rate independently of this toy root test.

## Moderate root-accuracy experiments (not temporal qualification)

Independent fresh [0,1e-6] solves, both P1 and P2 mandatory:
P1 atol=rtol=1e-12; P2 atol=rtol=3e-13. Both use stol=0,max_it=30,
newtonls/bt,preonly/LU/MUMPS. Absolute guess is old physical state; rate
guess is semidiscrete CH rate with previous u/pi/mu. No solution sharing.
P1 cross-comparison failure permits the predeclared P2, but physical failure
stops downstream work. P2 must converge and independently pass physics.

Freeze phi/mu relative consistent L2 comparisons at 1e-10 (the stricter
STEP3A6 threshold), D_CH relative at 1e-10, F and mass absolute at 1e-12.
Velocity and pressure differences are scaled by U_ref sqrt(area) and
P_ref sqrt(area), U_ref=sqrt(sigma/(rho epsilon)), P_ref=sigma/epsilon;
both limits 1e-9. Report viscous/slip powers without unstable relative
division by almost-zero values. No mean subtraction of mu or pressure.

## Exact tiny proof and immutable trajectory

After moderate PASS, exactly one primary full rate step at
dt=4.657657028534679e-12, fresh historical initial state. Controls for this
and EVERY trajectory interval: atol=1e-10,rtol=1e-9,stol=0,max_it=30,
newtonls/bt,preonly/LU/MUMPS. Actual getters checked before/after solve.
Compare with the archived STEP3A6 accepted isolated tiny step. The optional
absolute tiny negative control is NOT elected; reproducing its failure is
not a transfer gate.

Then fresh full L0: first 50 intervals, diagnostic restart fork from 25
through 28, extend prefix through first REAL dt change plus five accepted
intervals (derived from plan, not a hardcoded index), cost gate, continue
the SAME path to 1068 and t=1e-4. No full CHNS L1/L2. First r guess is the
unprojected semidiscrete CH rate; thereafter previous accepted r, unscaled.
u/pi/mu guesses always previous full accepted fields, never isolated fields.
SNES reason<=0 or physical failure stops the level; no adaptive retry.

## Physical acceptance and incompressibility

Every state: finite/material, mass relative to initial AND target /area
<=1e-10, certified P2 cells>=8, exactly two wall crossings, energy growth
<=1e-9 J/m. Existing BEWorkDiagnostics every interval: absolute weak work,
decomposition and work identity <=1e-12 J/m. Physical total dissipation is
D_CH+D_visc+D_slip; B_BE remains signed numerical work, never heat.

Essential velocity and gauge coefficient errors, scaled by U_ref/P_ref,
<=1e-12. Measure both every accepted state. The qualified manufactured
weak-continuity Riesz bound <=1e-12 is retained, including the gauge row.
No qualified universal relative strong-divergence policy exists in the
project. Record strong and relative norms every step; enforce no increase
beyond the same-mesh historical stiff-bump absolute range:
strong-divergence L2 <=3.6880704671415807e-6 m/s (maximum of committed
STEP3A2 bump histories). Relative strong divergence remains diagnostic;
the old maximum 0.9301196221025294 is reported for context, not silently
treated as a previously qualified tolerance.
Tiny CI cases test weak continuity/BCs, not production range evidence.

## Transfer thresholds and denominators

At frozen STEP3A7 common field checkpoints (no time interpolation):
max ||phi_full0-phi_iso0||_M / A_phi <=1e-3, where
A_phi=||phi_initial-phi_eq||_M. Current-disturbance normalization is
diagnostic and used only above 1e-6 A_phi. Chemical difference is diagnostic
with denominator max(||mu_iso-mu_eq||,1e-6 initial chemical disturbance).

Every common scalar time (all L0 intervals):
max |F_CH_full-F_iso0|/initial_excess <=1e-3;
max significant |D_CH_full-D_CH_iso0|/D_CH_iso0 <=1e-2;
max E_kin/initial_excess <=1e-4;
max significant (D_visc+D_slip)/D_CH_full <=1e-4.
Power significance floor is 1e-8 D_CH(0); below it record absolute/initial
power differences, no meaningless relative gate. Final cumulative hydro
dissipation /initial_excess <=1e-4. Final full budget/actual nonzero
|Delta E_total| <=5%, independent of the shared applicability flag.

Additionally require coupling <= measured isolated L0->L2 uncertainty:
max common phi gap <=max common iso02 phi gap;
final CH-integral gap <=abs(I_iso2-I_iso0);
final CH-free-energy gap <=max(abs(F_iso2-F_iso0),1e-13 J/m).
The energy floor is fixed from STEP3A7 stable-work versus total-subtraction
differences <=2.47e-14 J/m and double-precision quadrature accumulation;
1e-13 is below the measured iso02 gap 6.0652038946784614e-12, so it does not
enlarge this case's transfer allowance. No thresholds will be retuned.

Report separate coupling and temporal components, their sum as a practical
triangle-style proxy (NOT a rigorous full-continuum bound), and direct
full0 versus iso2 fields. At fixed expanded audit points measure R_adv and
R_diff phase Riesz norms independently; this ratio is diagnostic only.

## Transactions, restart and archives

Private reconstructed full candidate passes all gates before physical state,
old/older, time, index, retained rate and persistent Diagnostics advance.
True r, phi_old, dt*r, rounded phi, recovered-rate DIAGNOSTIC, u/pi/mu and
independent array hashes are saved at predetermined field checkpoints.
Checkpoint includes all full state/history, cumulative CH/viscous/slip
powers and BE sums. Reuse STEP3A7 atomic/hash-chained journal architecture;
scientific FAIL forbids retry/resume. Sessions finish with explicit status
complete/failed; COMPLETE requires exact full horizon and verified archives.

Production fork limits: phi/mu/r relative L2<=1e-10; u/pi differences scaled
by U_ref/P_ref sqrt(area)<=1e-10; cumulative CH relative<=1e-10; cumulative
hydro and total energy/budget differences absolute<=1e-12 J/m. The fork
does not replace a primary field. Tiny separate-process restart crosses dt.

## Cost and source freeze

New STEP3A8 budget, explicitly replacing STEP3A4's 7200s SHORT-probe budget
for this different continuous-to-1e-4 question: preliminary scientific
moderate/tiny/pilot/fork<=1800s; full L0 including reused prefix<=14400s;
total scientific wall<=16200s. Pilot is charged once in the total and is
also included in its two applicable component limits. Build/general pytest
and report generation are separate. This change is declared BEFORE any
production solve; none of these limits may increase after measurements.
Forecast remaining L0 using only full phase-rate pilot accepted timings:
1.25*max(mean,p95), plus actual spent time/initialization. If the complete
horizon cannot fit, stop with incomplete transfer, never shorten its claim.

Freeze complete solver/runner/analysis/tested sources, policy, design,
original weak-form hash, stack, source archive, all schedules and isolated
reference fingerprints BEFORE the first production moderate solve.
Changing implementation then stops this series. Execution base and exact
dirty source snapshot remain distinct from the eventual commit.

Maximum verdict: FULL CHNS PHASE-RATE BRIDGE QUALIFIED FOR THE
MATCHED-DENSITY STIFF-BUMP CASE; HYDRODYNAMIC COUPLING IS BELOW THE QUALIFIED
ISOLATED-CH TEMPORAL UNCERTAINTY; ISOLATED-CH TEMPORAL QUALIFICATION
CONDITIONALLY TRANSFERS TO THIS FULL-CHNS BENCHMARK.
Overall MODEL NOT YET VALIDATED.
