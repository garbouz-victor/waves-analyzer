# STEP 3A.9 — divergence policy audit, then conditional transfer resumption

Overall: **MODEL NOT YET VALIDATED**. Historical STEP3A8 is read-only and its
rejected P1 candidate will never be promoted. Execution base is
`5f2d9498e1d0715f57b125ee49b4ddfd533467ae`; execution source archives, not a later
commit, identify calculations. No solver-core changes, retries, dt changes,
physics changes, pressure-space changes, or post-solve corrections are allowed.

## A. Question and predeclared decision

For continuous vector P2 / continuous scalar P1 Taylor–Hood, the enforced
equation is `(q_h, div u_h)=0` for every pressure test function. Cellwise linear
`div u_h` need not be continuous and hence need not belong to the pressure
space. Weakly divergence-free does not mean pointwise divergence-free.

Independently assemble a new pressure mass matrix and load, solve `M_Q c=b`,
and integrate the square of the resulting pressure-space Function. Compare
its norm with the existing FullObserver continuity Riesz result. Integrate
`(div u-P_Q div u)^2` independently and verify Pythagoras. Quadrature is 12.
Projection agreement tolerance is `1e-11*max(projected, Riesz, 1e-30)` plus
`64*machine_epsilon*strong_divergence`, accounting explicitly for assembly
cancellation. Pythagorean relative squared-norm tolerance is `1e-11` (floor
`1e-60`). Both the discrepancy and the tolerance are archived.

The old candidate is reassembled read-only before any spatial transient solve.
Cell contributions are exact quad12 integrals of divergence squared, not
point samples. Report percentiles and top ceil(1%,5%,10%) cell shares of the
squared norm, also their square roots. Localization is diagnostic: interface
cells have a P2 Bernstein interval intersecting [-0.9,0.9]; wall-near cells have
centroid distance <=2 epsilon from any domain wall; contact-near cells are
both interface-active and bottom-near. Report overlapping shares explicitly,
and bulk as neither interface-active nor wall-near. No visual threshold tuning.
A tiny algebraic kernel projection of the actual P1/P2 divergence matrix
demonstrates weak zero with nonzero strong divergence.

### Spatial experiment (fixed before measurement)

| mesh | nx x nz | topology | h |
| --- | --- | --- | --- |
| M0 | 96 x 48 | crossed triangles | maximum geometric cell diameter |
| M1 | 120 x 60 | crossed triangles | same definition |
| M2 | 144 x 72 | crossed triangles | same definition |

Same domain, P2/P1/P2/P2, quad12, theta60, epsilon0.025, sigma0.1,
mobility1, matched density/viscosity1, slip0.01, zero body force. Mesh allocation
guard `max_unknowns=400000` is permitted on M1/M2 only; this is not physics or
solver accuracy. M0 keeps its exact configuration and all historical arrays.
All mesh initializations and certified resolution checks are completed before
the three one-step absolute CHNS solves, in M0/M1/M2 order.

M0 loads the historical prepared equilibrium and frozen perturbation arrays.
M1/M2 each solve their own mass-constrained variational equilibrium using the
unchanged equilibrium solver and its default frozen EquilibriumPolicy.
Initial guess: the same analytic sessile cap at theta60, radius0.3, centered
on the same domain bottom. Initial mu guess is sigma/(2 radius). All meshes
use the **same** mass target `-0.6060961027495515`, not the analytic seed's
mesh-dependent mass. No branch search/retry. Each new equilibrium must pass
the existing stationarity/mass criteria and have constant mu, negligible
initial CH dissipation, and certified cells >=8. The unchanged
`zero_mass_perturbations(...)["smooth_bulk"]` defines the analytic two-bump
family and deterministic discrete zero-mass construction. Apply delta=1e-3.
Archive equilibrium, psi, initial phi and all coefficient hashes per mesh.
M0 exact hashes are independently compared with committed STEP3A4/7 policy.

Transient: one **absolute** BE interval [0,1e-6], atol=rtol=1e-12,
stol0, max_it30, newtonls/bt, preonly/LU/MUMPS, independent old-state guess.
No rate solve in spatial study. All non-strong-divergence physical gates from
STEP3A8 remain hard, including weak continuity <=1e-12. Strong divergence is
the audited output, not an acceptance veto. Record actual controls/reasons,
velocity/pressure DOFs, grad-u, symgrad-u, velocity norm/speed, kinetic energy,
viscosity dissipation, mass, topology, resolution and BE work.

Reclassification is allowed **only if all**:

1. All three weak continuity Riesz norms <=1e-12.
2. Independent projection agrees with existing Riesz, and Pythagoras passes.
3. `D_M2 < D_M1 < D_M0` strictly.
4. Least-squares slope of log strong-divergence vs log geometric h >=0.5.
5. Every other physical/material/energy/BC/gauge gate passes on all meshes.
6. All SNES and KSP solves succeed without retry.

Slope is a measured refinement slope, not a formal optimal Taylor–Hood order;
0.5 requires positive spatial convergence well above numerical noise. No
manual override. Any failure stops PART B. Optional cross-mesh velocity
transfer is not undertaken; it is unnecessary for this decision.

## B. Predeclared production policy, activated ONLY after audit PASS

Version `step3a9-v1`: weak continuity Riesz <=1e-12 is hard. Strong divergence
must be finite and is measured every step, but its historical envelope
3.6880704671415807e-6 is stored only in historical metadata, inactive. For
strong divergence >=1e-12, pressure-projection fraction <=1e-6 is hard; below
this significance floor the absolute weak gate suffices. Relative strong
divergence is diagnostic. No optional five-times alert is enabled.

All STEP3A8 physical, coupling, energy and cost thresholds are unchanged.
Independent P1/P2 abs/rate solves at dt1e-6 are mandatory (P1 atol=rtol1e-12;
P2 3e-13). At **both** levels phi/mu/D relative <=1e-10, F/mass absolute
<=1e-12, scaled u/pi <=1e-9, and relative strong-divergence abs/rate difference
<=1e-4. Within-formulation P1/P2 strong-divergence change is diagnostic.

Only then: one fresh full-rate target dt4.657657028534679e-12 with production
atol1e-10, rtol1e-9, stol0, max_it30, newtonls/bt/preonly/LU/MUMPS, no retry.
Then exact STEP3A7 L0 continuous prefix50, checkpoint25 fork through28,
first actual block boundary+5, cost gate, same trajectory through1068 and
t1e-4. First rate guess semidiscrete; later previous accepted rate unscaled.
Scalar/weak/projection/strong/physical/work gates every interval. True retained
rate, selected common fields, atomic checkpoints and COMPLETE required.
No full CHNS L1/L2, new optimizer, isolated reruns or STEP3B.

Coupling against same-clock isolated L0: phi/initial disturbance<=1e-3,
F_CH/initial excess<=1e-3, significant D_CH relative<=1e-2,
kinetic/initial excess<=1e-4, hydro power/CH<=1e-4 and cumulative
hydro/initial excess<=1e-4. Full energy closure<=5% of actual energy change.
Phi, integrated CH and final CH-free-energy coupling must also not exceed
measured isolated L0–L2 discrepancies (energy precision floor1e-13).
Verify these references from complete archives; never substitute L0–L1.

## Two freezes, provenance, costs and stop conditions

Audit freeze stores decision/hierarchy hashes, source archive, core file hashes,
input artifacts and stack before spatial runs. Core files full_phase_rate,
full_rate_bridge, weak_form, free_energy, material, boundary must match base.
After audit PASS, production runner/diagnostics/policy/source/stack/schedule
are frozen separately before new production solves. No implementation changes
within either frozen scientific series. A discovered solver bug stops A9.

Scientific wall budgets: audit<=3600s; preliminary moderate/tiny/pilot<=1800s;
continuous L0<=14400s; total<=19800s. No post-hoc increase. Initialization,
equilibrium preparation and JIT inside scientific sessions count; general
pytest/container build do not. Record SNES, diagnostics, projection/resolution,
BE work and I/O separately. Forecast uses 1.25*max(mean,p95) of full accepted
step wall times. Timeouts create failed/incomplete records, not acceptance.

Historical directories remain read-only including JIT caches (new cache path).
New reports explicitly mark NOT RUN, never use zeros for absent results, and
retain any failure. No remote CI claim without observation, no automatic push.

## Tests

Pure decision cases (increasing divergence, shallow slope, failed weak/physics/
solver/projection all veto), same-root comparison, inactive historical cap,
finite/fraction/weak gates, budget logic. Pinned tiny FEM: independent projection,
kernel demonstration, unchanged residual/Jacobian/BC bridge, new policy
candidate rejection, transactional history/restart regression. Production
meshes and trajectories are excluded from CI.
