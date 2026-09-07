# STEP 3A — measured model-development checkpoint

**MODEL NOT YET VALIDATED**

This report records a stopped validation sequence, not completion of STEP 3A's
acceptance criteria. No tank bridge, film demo, water–air calculation, or film
animation was run. No commit named “add validated … model” is justified yet.
The base remains `c6dcf3faf28ff2cf9f303e5d5d8832d9ed2f7517`.

## 1. Selected model and scientific audit

The independent source branch is `src/sloshing/multiphase/`. It uses the
volume-averaged AGG mass/momentum structure, not the former linear equations.
Primary references: [Abels–Garcke–Grün](https://arxiv.org/abs/1104.1336) for mass
transport, and [Qian–Wang–Sheng](https://arxiv.org/abs/cond-mat/0602293) for the
variational wall/slip setting. The conservative-potential extension and pressure
convention below are explicitly derived in [STEP3_DESIGN.md](STEP3_DESIGN.md),
which was written and self-reviewed before the large solver implementation.

Coordinates `(x,z)` are dimensional metres, z upwards. Velocity u is m/s,
mechanical p and chemical potential mu_ch are Pa, phi is dimensionless and is
positive in liquid. With dynamic viscosity eta in Pa s,

```
rho = rho_bar + rho_prime phi
rho_bar=(rho_l+rho_g)/2; rho_prime=(rho_l-rho_g)/2
Psi=g z-a_x x; b=-grad Psi
f=lambda/(4 epsilon)(phi²-1)² + lambda epsilon/2 |grad phi|²
mu_0=lambda/epsilon(phi³-phi)-lambda epsilon Delta phi
mu_ch=mu_0+rho_prime Psi
J=-rho_prime M grad mu_ch

div u=0
phi_t+div(phi u)=div(M grad mu_ch)
(rho u)_t+div(u tensor (rho u+J))
    = -grad p + div(2 eta D(u)) + mu_0 grad phi + rho b
```

Thus `rho_t+div(rho u+J)=0`: solenoidal volume-averaged velocity does **not** imply
zero diffusive mass flux. The actual weak pressure is
`pi=p-phi mu_0+rho_bar Psi`; force becomes `-grad pi-phi grad mu_ch`.
Gravity must not be added a second time. Laplace measurements recover mechanical
p; testing the pressure jump of pi alone would be incorrect.

The tensor convention is `partial_j(u_i J_j)`. Full nonlinear transport and
`2 eta D:D` are retained. Affine rho/eta are never clipped. Quadrature samples
are checked for positivity/nonfinite values; polynomial phi can overshoot ±1.
High-density-ratio admissibility and gravity-driven unequal-density dynamics
have **not** been validated. Time-evolved benchmarks below are matched-density
demonstration fluids; the independent operator test also exercises rho_l/rho_g=1.2.

## 2. Free energy, surface tension and dissipation

Planar `phi=tanh(s/(sqrt(2) epsilon))` gives
`sigma=2 sqrt(2) lambda/3`, hence `lambda=3 sigma/(2 sqrt(2))`.
SymPy verifies the stationary chemical equation and the analytic energy integral;
independent numerical integration reproduces sigma within the test tolerance
2e-10 relative. This normalization is also checked by measured pressure jumps.

Per unit out-of-plane length,

```
E = integral [rho |u|²/2 + rho Psi + f] dA + integral f_wall ds  [J/m]
dE/dt = -D_viscous-D_CH-D_slip                                [W/m]
D_viscous=integral 2 eta D:D dA
D_CH=integral M |grad mu_ch|² dA
D_slip=integral (eta/L_s)|u_tau|² ds
```

This continuum identity assumes fixed conservative forcing and impermeable,
stationary, chemically no-flux boundaries. The gravity contribution to mu_ch
is necessary for the stated gravitational-energy budget with diffusive mass
transport. No production gravity benchmark is inferred from the algebra alone.

## 3. Contact regularization and wall law

```
f_wall=-sigma cos(theta_e)(3 phi-phi³)/4
lambda epsilon partial_n phi+f_wall'(phi)=0
n dot M grad mu_ch=0
u_n=0
(2 eta D n)_tau+(eta/L_s)u_tau=0
```

The angle is measured through liquid. Young's sign follows from
`f_wall(-1)-f_wall(+1)=sigma cos(theta_e)` and is symbolically checked at
60°, 90°, 120°. This first version uses instantaneous variational wetting,
not finite-rate wall relaxation; its relaxation dissipation is zero. A future
finite-rate law would require the matching uncompensated Young stress/GNBC.

Diffuse epsilon and positive Navier slip length remove the *sharp-interface +
strict-no-slip + moving-line* combination. They are explicit regularizations,
not molecular measurements. Hysteresis is absent. The sessile benchmark uses
an impermeable slipping/wetting **horizontal bottom wall**, free-slip other
walls, and no chemical flux. It tests contact motion but is not a high-wall tank
release. All crossings are retained; no “first crossing” is silently selected.

## 4. Parameters, groups and environment

All three benchmark configs have rho_l=rho_g=1 kg/m³, eta_l=eta_g=1 Pa s,
sigma=0.1 N/m, L_s=0.01 m, g=a_x=0. They are **not water**.

| Case | epsilon [m] | M [m²/(Pa s)] | theta_e | Reference L/U | dt [s] |
|---|---:|---:|---:|---|---:|
| Flat | .05 | .1 | 90° | 1 m / .1 m/s | .002 |
| Laplace | .03 | 1 | 90° neutral walls | .25 m / .1 m/s | .02, .01 |
| Contact | .025 | 1 | 60° | .3 m / .1 m/s | .01 |

For example contact Re=.03, We=.03, Ca=1, Cn=.08333, L_s/L=.03333,
Pe_CH=.09. Pe_CH uses chemical scale sigma/L: `U L²/(M sigma)`.
Density/viscosity ratios are one. Bo=0 and l_c is infinite for these matched,
zero-gravity benchmarks; no finite capillary length is invented. These reference
groups do not substitute for measured contact-line speeds or a mobility study.

The immutable Docker base is
`dolfinx/dolfinx@sha256:f7cce2a2271bf838c080751348c471064acb41fef0330e2c08178a688f71890d`.
Measured stack: Python 3.12.3, DOLFINx/Basix/FFCx 0.10.0, UFL 2025.2.0,
PETSc/petsc4py 3.24.0, mpi4py 4.1.1, MPICH 4.3.1, NumPy 2.2.6,
SciPy 1.16.2, h5py 3.12.1, Matplotlib 3.9.4. MUMPS is available.
See [stack.json](docker/step3/stack.json). The old `.venv` was not modified.

## 5. FEM, time scheme, solver and cost

Velocity CG2 vector / pressure CG1 on triangles; phi/mu both CG2, with CG1 phase
comparison. Degree-12 quadrature, automatic UFL Jacobian, monolithic SNES Newton
with backtracking and LU/MUMPS. One pressure DOF fixes the closed-container
constant; no pressure penalty/compressibility stabilization is added.

Backward-Euler startup followed by fixed-step BDF2. Skew momentum transport uses
rho u+J and temporal counterpart `D_t(rho u)-0.5 D_t(rho)u`. Conservative CH
advection preserves its constant test. This is continuum-consistent; it is
**not a theorem of monotonic physical energy for BDF2**.

Each run logs Newton/linear iterations and residuals. No adaptive dt or silent
dt reduction is active. Rank-local checkpoint functions preserve state, both
BDF histories, config, time, energy accounting and MPI partition; exact restart
is unit-tested. Later runs now write periodic `running` checkpoints as well.
Loading running data requires explicit permission in the API; failed data are
not complete. Automatic interrupted-benchmark CLI continuation is still absent.

Interface roots use actual CG1/CG2 edge polynomials. Arc chords/connectivity are
an approximation to the FEM zero contour, not raster thresholds. The apparent
angle is fitted over z-distance windows 2–6, 2–10, 3–8 epsilon, never one wall node.
The width convention is phi=-.9…+.9: `4.1640655 epsilon`. The reported cell count
is a local projected-spacing indicator using interface-cell centroids, not a
proof that every geometric or molecular scale is resolved.

Flat runs cost 3.15/17.04/57.26 s (CG2 n=16/32/48), CG1 n=48 cost 50.98 s.
The resolved R=.25 Laplace case had 59,325 unknowns and took 64.11 s. The
contact run had 82,769 unknowns and took 964.57 s including per-step contour
analysis. Docker has 8 CPUs/8.06 GB; observed contact container RAM was ~401 MiB.
No large depth-10 mesh was allocated. Static circle-band refinement was used;
its failure to cover the **moving** interface is a measured limitation below.

## 6. Flat-interface spatial result — passed for this benchmark

Analytic tanh initial data, not a discrete equilibrium built from FEM matrices.
“Max speed” here is the maximum over velocity DOF samples and saved time steps.

| Phase | n | Max spurious speed [m/s] | Surface-energy error [J/m] | Max mass error / area | Min transition cells |
|---|---:|---:|---:|---:|---:|
| CG2 | 16 | 1.24171e-5 | 6.12885e-5 | 2.10e-15 | 3.30 |
| CG2 | 32 | 1.29036e-6 | 5.77580e-6 | 3.17e-15 | 6.65 |
| CG2 | 48 | 2.65084e-7 | 1.16911e-6 | 4.60e-15 | 9.99 |
| CG1 | 48 | 4.35239e-5 | 2.87239e-4 | 5.76e-15 | 9.88 |

Spurious currents and surface-energy errors decrease. CG2 is materially better
than CG1 at the same spacing. Coarse under-resolved cases are diagnostics, not
qualified interfaces. Physical energy did not grow. For fixed n=16, trapezoidal
energy-budget error decreases with dt=.002/.001/.0005:
0.0051662 / 0.0025745 / 0.0012793 relative. It is not an exact discrete identity.

## 7. Independent manufactured operator check

SymPy generates analytic solenoidal velocity, nonzero pressure, phi, mu and
coupled sources. The latter are confined to the benchmark; the production
residual has no source callback enabled. Non-matched density ratio 1.2 and
viscosity ratio 2 exercise variable-coefficient terms. Errors are L2 Riesz
representatives of weak defects, **not raw algebraic residuals**.

| n | Momentum defect [Pa] | CH defect [m/s] | Chemical defect [Pa m] | Strong div L2 [m/s] |
|---:|---:|---:|---:|---:|
| 6 | .0410458 | 3.12766e-4 | 1.84359e-4 | .00607481 |
| 12 | .0126204 | 8.61768e-5 | 6.32111e-5 | .00157836 |
| 24 | .00344652 | 2.47824e-5 | 2.15368e-5 | .000398412 |

Weak divergence is 3.78e-17…1.20e-16. An initial test incorrectly demanded its
monotone decrease at roundoff; it was replaced by a roundoff bound **and the
additional strong norm**, not by hiding divergence. Deliberately reversing the
pressure sign gives momentum defect .719395 at n=24 (versus .00344652).
This is an independent **operator-consistency** test, not full solution-error
MMS or proof of all transient unequal-density physics.

## 8. Laplace pressure — pressure passes; complete gate does not

Pressure averages exclude a 4-epsilon band. R is fitted to the computed phi=0
contour; the initial radius changes slightly through finite-epsilon relaxation.

| Initial R [m] | dt [s] | Measured R [m] | Delta p [Pa] | sigma/R [Pa] | Relative pressure error | Largest energy increase [J/m] |
|---:|---:|---:|---:|---:|---:|---:|
| .25 | .02 | .2413938 | .4138424 | .4142609 | .10101% | none |
| .30 | .02 | .2944377 | .3393050 | .3396304 | .09581% | 1.98013e-7 |
| .30 | .01 | .2944377 | .3393050 | .3396304 | .09581% | 1.64874e-8 |

All have >8 estimated transition cells. Relative chemical-potential standard
deviation is 2.89e-6 for the first case, ~6.33e-7 for the second radius.
The BDF2 physical-energy overshoot shrinks substantially with dt/2 but still
exceeds the pre-existing 1e-9 J/m gate. That gate was **not relaxed**. The earlier
coarser .25 run had excellent pressure but failed resolution and was not accepted.
No epsilon→0 Laplace convergence claim is made from this table.

## 9. Contact-angle / spreading attempt — gate failed

Starting apparent angle 90°, imposed Young angle 60°, t=0…0.5 s.
Measured final angle **60.13451°**. Contact crossings move from [-.3,.3] m to
[-.4093156,.4093111] m: the more-wetting drop spreads in the correct direction.
Half-width change .1093133 m; fit-window spread <1°; last angle change .01320°.
Mass error / area <=7.39e-14, impermeability norm zero. Final strong divergence
is 1.698e-7 m/s, relative to velocity gradient .002251. This is not a statement
that weak divergence alone qualifies the motion.

Two reasons prevent qualification:

1. From **t=.21 s** the transition enters coarser elements at the edge of the
   initial circular refinement band. Minimum indicator falls to **7.4883 cells**,
   with h_normal=.013902 m versus the finest edge spacing .008333 m. The 8-cell
   policy is retained. A static band around the initial interface was insufficient.
2. The initial 90° geometry is incompatible with the *instantaneous* 60° wetting
   law. Weak chemical-potential initialization represents that wall variation as
   a very fast discrete CH transient. Initial D_CH is **91,221.1 W/m**; trapezoidal
   integration with dt=.01 yields **456.109 J/m**, although E(0)=**.0942484 J/m**
   and E(.5)=**.0885138 J/m**. This is emphatically **not physical dissipated heat**.
   The raw relative budget defect is 4,839.38, not an acceptable energy closure.

The likely mechanism is under-resolved stiff startup / power quadrature, rather
than a sign violation: physical energy decreases and the following fixed-space
temporal experiment supports this diagnosis. It does **not** yet prove acceptable
energy closure on the resolved moving-interface run.

| dt [s] | Max relative energy defect | Final E [J/m] | Max mass error / area |
|---:|---:|---:|---:|
| .0004 | .417810 | .0937567713 | 2.00e-15 |
| .0002 | .205987 | .0937556823 | 2.31e-15 |
| .0001 | .100656 | .0937554112 | 3.70e-15 |

This last experiment deliberately uses the **same coarse 2,397-unknown system**
to isolate time error; its spatial resolution is unsuitable for contact physics.
Runtime was 1.48/1.96/3.08 s. The defect decreases, but remains too large. No
curve fitting, energy correction, omitted initial power, smoothing, added
surface tension or slip/mobility tuning was used to disguise it.

## 10. MPI, restart, tests and software corrections

Identical tiny n=1/n=2 MPI runs differ in E_total by 5.55e-17 J/m, phase mass by
7.47e-17 m², sampled speed by 8.13e-20 m/s, strong divergence by 2.28e-18 m/s.
Restart reproduces the uninterrupted BDF state and cumulative accounting in its
unit test. It requires the same rank count/partition, not arbitrary repartition.

Development issues are recorded rather than counted as passed physics:
official-image PYTHONPATH preservation; Docker Desktop UID mapping; missing
edge-to-cell connectivity before local refinement; NumPy bool JSON conversion;
root-only atomic MPI status creation; and near-tangent endpoint root conditioning.
The latter came from phi≈3.14e-15 at an exact mesh vertex; a declared 1e-12 root
tolerance recognizes the endpoint without modifying any FEM coefficient.
Contour, CLI analysis-only and restart behavior have regression tests; the
environment and MPI corrections were checked by the repeated benchmark runs.
Failed output directories were not silently reused.

Final ordinary suite: **148 passed**, one DOLFINx module skipped outside its
container, two expensive old validation tests deselected; 15 old warnings visible.
Both expensive old validation tests were also run separately: **2 passed**.
Container STEP 3 tests: **26 passed** (20 overlap with ordinary pure tests, six
additional FEM tests). Runtimes were 28.71 s, 57.65 s and 6.19 s respectively.
See [test results](validation_results/step3/TEST_RESULTS.md). A dedicated
cheap-container CI workflow is added, but no remote GitHub run is claimed for
these uncommitted changes.

## 11. Missing studies and physical limits

Not yet run: PDE angle cases 90/120°, angle epsilon/h convergence, complete
spreading slip/mobility sensitivity, falling-film flux/profile benchmark,
draining patch, interface-thickness study, and transient unmatched-density
gravity validation. The Nusselt-with-slip helper and resolution classifier are
unit-tested formulas only; **no computed wall film is claimed**.

Epsilon is numerical diffuse thickness, not molecular thickness. L_s=.01 m is a
declared benchmark regularization, not a measured microscopic slip length. M is
large for rapid static relaxation and cannot be transferred to a physical tank
without sensitivity analysis. No contact-angle hysteresis, wetting experiment,
wall-film thickness, high-wall recession or drainage rate has been validated.
An h>=4 epsilon and >=6 cells classification would still require convergence;
the helper deliberately calls it `resolved_candidate`, not a physical proof.

## 12. Readiness and next action

**MODEL NOT YET VALIDATED**

The sequence stops at the failed contact/time-energy gates. It is not legitimate
to advance to a tank or publish a film animation. Next work is to qualify the
stiff startup with explicitly controlled temporal refinement (potentially
adaptive startup), cover the entire moving interface with adequate refinement,
and repeat the failed tests without changing thresholds. Only then proceed
through 60/90/120°, epsilon, M and L_s studies and the independent falling-film
benchmark. Adding finite wall relaxation would be a separately reviewed wall-law
change, not a shortcut to passing energy tests.

Old fine HDF5 SHA-256 remains
`fa299d74c9487a40a5fac55966c1288c950bc0a46116f80e7da0063f794cf09a`.
Historical STEP1_6 report and summary hashes are unchanged. No old PDE was rerun.

Artifacts: [gate summary](validation_results/step3/summary.json),
[gate CSV](validation_results/step3/summary.csv),
[gate diagnostics](validation_results/step3/gate_diagnostics.pdf),
[flat convergence](validation_results/step3/flat/flat_convergence.pdf),
[contact diagnostics](validation_results/step3/contact/theta60_resolved/contact.pdf).
Raw HDF5/checkpoints and compiler caches remain ignored. See README for exact
build/run commands. Benchmark `complete` means execution finished, **not** that
its scientific gate passed. Commit A is intentionally not created.

## STEP 3A.1 follow-up (base 15a9e07)

The automatic energy, interface-activation and settling gates were repaired;
see [STEP3A1_REPORT.md](STEP3A1_REPORT.md). Historical numbers above are retained,
not retroactively converted into passes. The new Bernstein-based final-checkpoint
audit finds only **3.8659** minimum nominal transition cells where the older
centroid indicator reported 7.4883. Physical-time settling also rejects the old
per-step apparent-angle success.

A new 60°/60° geometrically compatible, spatially resolved BE startup baseline
has a small overall energy defect (0.00533% of initial scale), but a **785%**
first-interval defect relative to actual energy change. Direct curved tanh data
are not a variational wall equilibrium merely because their zero-contour angle
equals theta_e. The sequence stopped at that baseline; no new Laplace series,
90→60 long run, tank or film followed. Current verdict remains
**MODEL NOT YET VALIDATED**. The repaired gates must not be weakened to advance.
