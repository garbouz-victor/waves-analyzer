# STEP 1.6 — animation qualification

## Goal

Квалифицировать конкретный dataset текущей linear pinned-contact модели,
не создавать production animation. База: e3b4ca1b93af22be6ade99cfcf56048a77805875.
Реально рассчитаны оба полных dataset 0…5 s, выполнены field/path gates.
Итог: **NOT QUALIFIED FOR PHYSICAL ANIMATION** для alpha=0.02°:
на fine локальный slope достигает 0.373408. Это объяснимое нарушение
small-slope qualification policy у pinned contact, не instability solver.
Bulk fields и контрольные trajectories согласуются; подробности ниже.

## Initial audit and config cleanup

HEAD соответствует заданному commit, исходный worktree чистый. Проверены
interfaces FEMSystem, SDIRK2, Diagnostics и SnapshotWriter. Математическая
модель и assembled operators STEP 1.5 не меняются. Давление по-прежнему
восстанавливается в момент snapshot, его gauge задаётся traction.
Для нового dataset нужен streaming writer без visualization grid с полными
FEM coefficients, точными regional P2 slopes и restart accounting.

Исторический configs/linear_safe.json переименован git mv в
configs/validation_alpha_0_2.json без изменения параметров. Добавлены
animation_candidate.json и animation_candidate_fine.json. Default config
и production defaults сохранены. Большие raw artifacts нового этапа ignored.

## Dataset parameters

Для обоих кандидатов a=1 m, d=10 m, g=9.81 m/s², alpha=0.02°,
nu=0.01 m²/s, SDIRK2, dt=0.00125 s, snapshot_dt=0.005 s, t_end=5 s.
Medium: 80×140; fine: 120×210. Grading: x=1.8, z=4.0.
Linearized Navier–Stokes, P2/P1, full strain tensor, no-slip на трёх solid
boundaries, sigma=0, fixed domain и pinned endpoints остаются прежними.
Snapshots содержат FEM u,w,q,eta,omega, а не визуальные samples.

## Why alpha=0.02° was selected

STEP 1.5 предложил этот угол только по линейности; здесь создаётся настоящий
run с такими initial data. Независимая short regression production path
(8×16, SDIRK2, dt=0.00125, t_end=0.025) измерила:

| Regression metric | Measured value |
|---|---:|
| tan(0.02°)/tan(0.2°) | 0.09999959790465966 |
| Max absolute velocity scaling difference, m/s | 1.7973e-17 |
| Max absolute eta scaling difference, m | 2.1684e-19 |
| Max relative energy scaling difference | 6.6443e-16 |

Реальный 0.02° dataset не заменяется масштабированной копией 0.2°.
Линейность позволяет предложить следующий угол после измерения максимального
slope, но сама по себе не квалифицирует ещё не рассчитанный dataset.

## Medium 0...1 s gate

Реальный run alpha=0.02°, medium, SDIRK2, dt=0.00125, snapshot_dt=0.005
завершён за 1172.53 s; 201 snapshot. No-slip/contact errors=0,
volume error<=3.838e-19 m², weak div<=4.951e-16 m/s,
relative energy budget<=1.151e-13, step energy growth=0.

| Metric | Maximum | Time, s |
|---|---:|---:|
| Global/contact slope | 0.20584348 | 1.000 |
| Bulk slope, abs(x)<=0.9 | 0.000453650 | 0.855 |
| Intermediate slope, abs(x)<=0.98 | 0.00411992 | 0.515 |
| Strong div L2, m/s | 1.73437e-4 | 0.590 |
| Relative strong div | 0.0842806 | 0.775 |
| Bulk strong div L2, m/s | 2.66313e-6 | 0.830 |
| Bulk relative strong div | 0.0360716 | 0.005 |

Slope first exceeds 0.1 at the saved time 0.545 s, never exceeds 0.3 by
1 s. Thus the conservative small-slope gate is **conditional**, not passed
unconditionally. Global relative div exceeds 5% on 42% of active snapshots
(active means grad norm above 1% of its peak): medium is questionable for
particle advection until mesh/path comparison. Bulk and global divergence
are explicitly separated; the large global contribution is predominantly
outside the bulk strip. No-slip, volume, energy and symmetry are intact,
so this is not an unexplained numerical failure blocking the 5 s investigation.
Eta antisymmetry error<=4.581e-18 m; wall omega norm asymmetry<=1.757e-17 m/s.

Pressure performance note: stiff accuracy makes the second SDIRK stage
pressure identical (up to linear-solve roundoff) to instantaneous pressure
recovery: its stage derivative satisfies B k_v=0. This identity is tested
against independent recovery at both nu=0.01 and 1. Later runs reuse that
already calculated q and do not retain a redundant large pressure LU.
The PDE and all stage equations are unchanged. Streaming writes are flushed
in bounded batches; incomplete files remain invalid caches.

## Medium 0...5 s result

The completed 201-snapshot prefix was continued, not overwritten. Extension
1...5 s took 2845.29 s; prefix plus extension took 4017.82 s. The resulting
file contains exactly 1001 snapshots, including both endpoints. Every saved
state was checked; internal-step energy accounting also spans the restart.

| Metric over 0...5 s | Measured value | Time of maximum, s |
|---|---:|---:|
| Global/contact slope | 0.269734707 | 4.295 |
| Bulk slope | 0.000453650108 | 0.855 |
| Intermediate slope | 0.00411992169 | 0.515 |
| Strong divergence L2, m/s | 0.000173436672 | 0.590 |
| Relative strong divergence | 0.0842805615 | 0.775 |
| Bulk divergence L2, m/s | 0.00000266313254 | 0.830 |
| Omega L2, m/s | 0.00369655598 | 0.350 |
| Left wall omega L2, m/s | 0.00256520731 | 0.340 |
| Pointwise omega maximum, 1/s | 0.295283328 | 0.565 |

At 5 s, E=8.74319352e-8 m^4/s^2, from E(0)=3.98439617e-7;
cumulative viscous dissipation=3.11008161e-7, and signed RK correction
=-4.78737434e-13 in the same units. Maximum relative budget residual is
3.83854e-13; no positive internal-step energy increment occurred. All six
no-slip residuals and contact error are zero. Volume error<=4.86621e-19 m²;
weak divergence<=4.95065e-16 m/s. Eta antisymmetry<=5.41424e-18 m,
left/right omega norm difference<=2.16841e-17 m/s.

Relative strong divergence exceeds 5% on 22.2% of active saved states,
not merely an isolated initial sample. Its bulk maximum is 3.60716% at
0.005 s and its bulk final value is 0.170471%. This distinction motivates
the independent fine-grid and trajectory comparisons; Bv≈0 is not used
to dismiss the strong-divergence issue.

The medium run is numerically stable but **conditional** by slope policy.
The first saved violation of 0.1 occurs at 0.545 s; 0.3 is never crossed.
Scaling this medium result alone would suggest alpha=0.00741469310° for
a maximum slope of 0.1, or 0.00370734657° for 0.05. These angles have NOT
been run and are not final recommendations: fine may impose a stronger
constraint. The absence of unexplained invariant failure permits fine-short.

## Fine 0...1 s gate

The actual fine run completed in 2117.88 s with 201 snapshots. Initial
assembly, factorization and pressure recovery took 699.5 s before t=0 was
saved; this startup cost is included, not hidden as background work.

| Maximum over 0...1 s | Medium | Fine |
|---|---:|---:|
| Global/contact slope | 0.205843481 | 0.251406582 |
| Bulk slope | 0.000453650108 | 0.000453718610 |
| Intermediate slope | 0.00411992169 | 0.00408994856 |
| Strong divergence L2, m/s | 1.73436672e-4 | 1.35739326e-4 |
| Bulk divergence L2, m/s | 2.66313254e-6 | 9.80723154e-7 |
| Relative strong divergence | 0.0842806 | 0.0644279 |
| Bulk relative strong divergence | 0.0360716 | 0.0162287 |
| Omega L2, m/s | 0.003696556 | 0.003716665 |
| Pointwise omega maximum, 1/s | 0.2952833 | 0.3579024 |

Fine no-slip and contact errors are zero; volume error<=1.10300e-19 m²,
weak divergence<=7.79562e-16 m/s, relative energy budget<=1.38448e-13,
and maximum positive internal-step energy increment=0. Eta antisymmetry
error<=1.08692e-17 m; wall omega norm asymmetry<=2.53704e-17 m/s.
Slope first exceeds 0.1 at saved time 0.490 s; it remains below 0.3 through
1 s. Fine relative divergence exceeds 5% on 24.5% of active short-run states,
versus 42% on medium over the same short window. These percentages must not
be compared to the medium **five-second** fraction 22.2% without noting the
different windows.

The numerical field comparison uses 237114 common triangles and 1422684
degree-four quadrature points. Its gate is reported separately from the
slope policy: more localized contact gradients do not erase the measured
improvement of bulk and integral derivative norms.

| Short-window difference | Max absolute L2 difference | Relative to peak fine norm | Time at maximum, s |
|---|---:|---:|---:|
| eta full, m^(3/2) | 2.40418398e-6 | 0.84354% | 1.000 |
| eta bulk, m^(3/2) | 4.98013751e-8 | 0.020465% | 1.000 |
| eta intermediate, m^(3/2) | 6.08002026e-8 | 0.021989% | 1.000 |
| eta contact, m^(3/2) | 2.40341506e-6 | 3.47735% | 1.000 |
| velocity, m²/s | 3.52779128e-7 | 0.042318% | 0.575 |
| omega full, m/s | 2.09041077e-4 | 5.62443% | 0.550 |
| omega left wall, m/s | 1.47804427e-4 | 5.73043% | 0.550 |
| omega right wall, m/s | 1.47804427e-4 | 5.73043% | 0.550 |
| omega surface layer, m/s | 2.09040574e-4 | 5.67534% | 0.550 |

The comparison took 176.44 s; integrated overlay area=20.000000000000036 m².
Global and bulk strong divergence, absolute and relative, are lower on fine
at **every active common snapshot**. Velocity reflection errors are at most
5.94143e-17 (medium) and 1.16010e-16 m/s (fine). The predeclared numerical
gate passes without changing limits. This authorizes the long fine run,
not unconditional animation qualification. Wall-vorticity field error is
materially larger than surface/velocity error and must remain visible.

## Fine 0...5 s result

Fine completed all 1001 snapshots. Extension 1...5 s took 5283.23 s;
subsequent complete-payload verification brought the writer invocation to
5300.8 s. Prefix plus extension calculation time is 7401.11 s (about
123.35 minutes), excluding the separately timed field comparisons.
The original short prefix is preserved. No-slip and contact errors are
zero; volume error<=1.10300e-19 m², weak divergence<=7.79562e-16 m/s,
relative energy budget<=4.52015e-13, and positive internal-step energy
increment=0. Eta antisymmetry<=1.08692e-17 m and wall omega norm
asymmetry<=2.55872e-17 m/s.

The actual maximum slope is **0.3734076287 at 4.320 s**. The first saved
violation of 0.3 is **2.365 s** (snapshot resolution 0.005 s); no exact
continuous-time crossing is claimed. Final slope=0.318061816. This excludes
alpha=0.02° from global linear-small-slope animation qualification, without
identifying a broken solver. The later maximum demonstrates why a short
validation could not qualify the five-second contact behavior.

At 5 s K=2.68147843e-8, P=6.03743569e-8, E=8.71891411e-8 m^4/s².
Cumulative viscous dissipation=3.11250955e-7 and signed RK correction
=-4.78749481e-13 in the same units. Strong divergence maxima are unchanged
from fine-short: 1.35739326e-4 m/s globally and 9.80723154e-7 in bulk.
Final global relative divergence=0.009395376, bulk=0.000622347.
Global relative divergence exceeds 5% on 9.5% of active five-second states.

Exact tan scaling of the measured full-run maximum suggests
alpha=0.005356077171° for target slope 0.1, or 0.002678038591° for 0.05.
Neither angle has been run. These are finite-mesh, finite-time extrapolations,
not a guarantee of a bounded continuum corner slope. No automated sequence
of decreasing angles has been launched.

## Small-slope history

| Maximum over 0...5 s | Medium | Fine |
|---|---:|---:|
| Global | 0.269734707 at 4.295 s | 0.373407629 at 4.320 s |
| Bulk, abs(x)<=0.9 | 0.000453650108 | 0.000453718610 |
| Intermediate, abs(x)<=0.98 | 0.00411992169 | 0.00408994856 |
| Contact, 0.98<abs(x)<=1 | 0.269734707 | 0.373407629 |
| First saved slope>0.1 | 0.545 s | 0.490 s |
| First saved slope>0.3 | none | 2.365 s |
| Final global slope | 0.217267226 | 0.318061816 |

Both meshes resolve a small bulk slope; the global maximum is entirely
in the contact strips. Medium alone is conditional by the project policy,
but the fine reference fails its hard 0.3 criterion. Selecting only medium
or only the final frame to call the candidate safe would be misleading.

## Strong divergence history

Taylor–Hood is weakly, not pointwise, divergence-free. Over the full interval,
maximum global L2 divergence decreases from 1.73436672e-4 to 1.35739326e-4
m/s; bulk L2 decreases from 2.66313254e-6 to 9.80723154e-7 m/s. Maximum
relative values are 8.42806% and 6.44279% globally, 3.60716% and 1.62287%
in bulk. Both bulk-relative maxima occur at the first nonzero snapshot,
t=0.005 s, and are not masked by excluding startup.

At t=5 s, global relative divergence is 1.85187% (medium) vs 0.939538%
(fine), and bulk relative divergence is 0.170471% vs 0.0622347%.
Medium exceeds 5% on 22.2% of active snapshots, fine on 9.5%; the active
criterion uses gradient norm>1% of its own peak. Thus medium carries the
specified particle-advection caution even though Bv is near roundoff.
Actual trajectory comparisons, not weak-divergence arguments alone, are
required to assess the seeded bulk paths.

## Medium/fine convergence

All 1001 common snapshots were compared, not just the illustrated times.
The full comparison took 1126.81 s. Relative errors use the **peak fine norm
over 0...5 s**, including the initial surface norm; no division by instantaneous
zero velocity is performed. Fine is a reference, not an exact solution.

| Field/region | Max absolute difference | Max relative | Final relative | Time at max, s |
|---|---:|---:|---:|---:|
| eta full, m^(3/2) | 5.61803819e-6 | 1.97116% | 1.91941% | 4.495 |
| eta bulk, m^(3/2) | 1.84848238e-7 | 0.075961% | 0.038100% | 4.470 |
| eta intermediate, m^(3/2) | 2.00767519e-7 | 0.072609% | 0.036554% | 4.470 |
| eta contact, m^(3/2) | 5.61448348e-6 | 8.12325% | 7.91361% | 4.495 |
| velocity, m²/s | 6.61011786e-7 | 0.079293% | 0.066940% | 4.850 |

Velocity difference grows from the short-window maximum 3.52779e-7 to
6.61012e-7 m²/s, and bulk eta difference from 4.98014e-8 to 1.84848e-7
m^(3/2). This late change is real and reported; both remain below 0.08%
of the respective peak reference norms and vary with the oscillation phase.
It is not a large/non-convergent bulk discrepancy. In contrast, contact
surface differences dominate the full surface error and grow substantially.
Two production meshes do not establish an asymptotic convergence order.

Maximum energy difference=2.49690791e-10 m^4/s² at 4.810 s, final difference
=2.42794081e-10. The maximum is about 0.06267% of initial energy. Global
and bulk divergence, absolute and relative, decrease on fine at **every
active common time**, where this comparison uses fine velocity norm>1%
of its peak. Velocity reflection errors over 5 s remain <=5.94143e-17
(medium) and 1.16010e-16 m/s (fine).

| Selected probe difference over 5 s, m/s | Maximum | Time, s |
|---|---:|---:|
| u(0,-0.1) | 6.71511538e-7 | 4.880 |
| w(-0.75,-0.2), mirrored right equivalent | 6.21963197e-7 | 4.870 |
| u(0,-1) | 1.21011875e-7 | 4.860 |
| u(0,-3) | 1.47745241e-8 | 0.030 |
| u(-0.95,-0.05), mirrored right equivalent | 6.24300880e-7 | 0.390 |
| w(-0.95,-0.05), mirrored right equivalent | 8.64915134e-7 | 2.215 |

The compact CSV includes 0.1, 0.5, 1, 2, 3, 4, 5 s plus automatically
selected times of maximum kinetic energy, contact slope, wall omega norm,
and velocity difference. Full raw curves are reproducible and ignored.
The predeclared long-window field gate passes without changing thresholds.

## Vorticity comparison

Vorticity is the DG1 representation of exact P2 FEM derivatives, not finite
differences of a visualization grid. The following are **norms of the field
difference**, not merely differences of norms.

| Region | Max L2 difference, m/s | Max relative | Final relative | Time at max, s |
|---|---:|---:|---:|---:|
| Full domain | 2.09041077e-4 | 5.62443% | 0.78833% | 0.550 |
| Left wall | 1.47804427e-4 | 5.73043% | 0.80235% | 0.550 |
| Right wall | 1.47804427e-4 | 5.73043% | 0.80235% | 0.550 |
| Surface layer z>=-1 | 2.09040574e-4 | 5.67534% | 0.79537% | 0.550 |

No new late-time omega discrepancy exceeds the short-window maximum.
Peak full omega norm is 0.00369656 vs 0.00371666 m/s; peak left-wall norm
is 0.00256521 vs 0.00257929 m/s. Those norm differences are much smaller
than the approximately 5.7% spatial field difference and must not substitute
for it. Pointwise omega maxima are 0.2952833 vs 0.3579024 1/s, about 21.2%
apart; these corner maxima are not declared resolved. Integral wall structure
is coherent under the stated 10% screening criterion, but less accurate
than bulk velocity. No temporal or spatial smoothing has been used.

## Energy/volume/no-slip diagnostics

All saved snapshots and all internal-step energy increments were checked.

| Diagnostic | Medium | Fine |
|---|---:|---:|
| All six no-slip residual maxima, m/s | 0 | 0 |
| Pinned endpoint error, m | 0 | 0 |
| Volume error maximum, m² | 4.86620428e-19 | 1.10299572e-19 |
| Weak divergence maximum, m/s | 4.95064695e-16 | 7.79561908e-16 |
| Relative energy budget residual maximum | 3.83853362e-13 | 4.52014239e-13 |
| Positive internal-step energy increment maximum | 0 | 0 |
| Eta antisymmetry maximum, m | 5.41423460e-18 | 1.08691268e-17 |
| Left/right omega norm asymmetry, m/s | 2.16840434e-17 | 2.55871713e-17 |
| E(0), m^4/s² | 3.98439617e-7 | 3.98439617e-7 |
| E(5), m^4/s² | 8.74319352e-8 | 8.71891411e-8 |
| Cumulative viscous dissipation, m^4/s² | 3.11008161e-7 | 3.11250955e-7 |
| Signed cumulative RK correction, m^4/s² | -4.78737434e-13 | -4.78749481e-13 |

The initial energy also agrees with the analytic tilted-surface integral
g*a³*tan(alpha)²/3. No unmodeled energy gain is used to sustain oscillations.
The signed RK term is tracked separately from positive viscous dissipation.

## Measured sloshing period and damping

The diagnostic is a(t)=(1/a) integral eta(x,t) sin(pi*x/(2a)) dx, evaluated
by quadrature of the calculated P2 trace. This is a projection, not a
potential-flow replacement of the viscous PDE. Event times use linear
zero-crossing interpolation and a local three-sample quadratic extremum
estimate; the underlying saved fields are not smoothed or changed.

| Medium signal event | Time, s | Amplitude, m |
|---|---:|---:|
| First zero crossing | 0.416116170 | 0 |
| First opposite extremum | 0.809724177 | -0.000248137142 |
| Next zero crossing | 1.222267262 | 0 |
| Next positive extremum | 1.615824218 | 0.000220813449 |

Twice the median separation of successive zero crossings gives
T=1.612640536 s. Individual half-cycle estimates span 1.610907...1.617142 s.
The inviscid deep-tank reference 2*pi/sqrt(g*k1*tanh(k1*d)), k1=pi/(2*a),
is 1.600609632 s, about 0.75% shorter. It is used only for interpretation.
Successive same-sign extrema have amplitude ratios 0.77781...0.77975,
with measured log decrements 0.24878...0.25127 per cycle. These are extracted
from actual viscous evolution, not imposed exponential damping.

The first kinetic-energy maximum at 0.395 s is close to the first modal
zero at 0.416 s; a multimode viscous surface need not become exactly flat
at a modal zero, especially because endpoints remain pinned.

Fine gives T=1.613247280 s, first zero=0.416097553 s, first opposite
extremum at 0.809794008 s with amplitude -2.48242972e-4 m, and next
zero=1.222442416 s. Its same-sign amplitude ratios are 0.77756...0.77992.
The two period estimators differ by 0.000606744 s (about 0.0376%);
the last zero crossings differ by 0.000583913 s. These are interpolated
event estimates from 0.005 s samples, not claims of exact event timing.
Cycle-to-cycle variation of the multimode signal is larger than that
medium/fine period difference; the appropriate characteristic period to
quote is approximately **1.613 s**, not a single exact eigenfrequency.

## Particle trajectory qualification

Forty deterministic particles were integrated over the actual medium and fine
0…5 s fields: 28 bulk seeds (x=-0.75,-0.5,-0.25,0,0.25,0.5,0.75 at
z=-0.1,-0.3,-1,-3 m), and 12 near-wall seeds (x=±0.90,±0.95 at
z=-0.05,-0.2,-1 m). None is placed in the contact corner. Spatial evaluation
uses the actual triangular P2 coefficients, not a regular visualization grid.
Temporal interpolation is linear between saved states. RK4 uses a maximum
step of 0.0025 s and ends steps exactly at snapshot knots.

| Mesh comparison, 0…5 s | Bulk | Near-wall |
|---|---:|---:|
| Maximum medium/fine path separation, m | 4.02706092e-7 | 7.60592773e-7 |
| Maximum final separation, m | 1.62703942e-7 | 4.21564438e-7 |
| Time of maximum separation, s | 4.460 | 4.320 |
| Maximum fine displacement from seed, m | 5.27814768e-4 | 6.01134850e-4 |
| Separation / maximum fine displacement | 0.07630% | 0.12653% |
| Invalid trajectories, either mesh | 0 | 0 |

The largest bulk difference is for seed (-0.75,-0.1); the largest near-wall
difference is for (-0.95,-0.05). These are comfortably below the declared
0.005 m bulk criterion. Path maxima are measured at common 0.005 s output
times, not certified continuous-time suprema. First-invalid checks are also
made at RK stages: crossing the reference box or reconstructed P2 surface
permanently invalidates a path, without reflection or silent clipping. No
such event occurred in any of the mesh, sampling or RK-step studies.

The wall-limit diagnostic evaluates both sides, z=-0.05,-0.2,-1 m and
t=0.1,0.5,1,2,5 s on both meshes. Epsilon decreases through
0.02,0.01,0.005,0.001,0.0001,0.00001 to exactly zero. Across these 60 cases,
speed(epsilon=1e-5)/speed(epsilon=0.02) ranges from 0.00042417 to
0.00065660. The maximum directly evaluated speed **on** the wall is
2.27692e-21 m/s (evaluation roundoff); constrained DOF residuals are zero.
No interpolation-induced finite wall velocity appears. No power law was
assumed or fitted.

The entire particle phase, including the extra dense-snapshot PDE run,
took 546.2 s. Actual maximum particle excursions are only about 0.5–0.6 mm
at this angle. A future renderer must not silently magnify those excursions.
Reproducible paths do not validate nonlinear net transport; see the model
limitation below.

## Snapshot temporal interpolation test

A separate real medium run to 0.5 s saved every 0.0025 s, while retaining
PDE dt=0.00125 s. It contains 201 snapshots and took about 503.8 s of
calculation. At all common times its u,w,eta coefficients equal the
0.005 s dataset exactly: maximum difference is 0.0. Thus the comparison
isolates saved-time interpolation, not a different PDE trajectory.

| Comparison over 0…0.5 s | Bulk maximum, m | Near-wall maximum, m |
|---|---:|---:|
| snapshot_dt=0.005 vs 0.0025 s | 9.27683540e-9 | 1.29330557e-8 |
| RK maximum step=0.0025 vs 0.00125 s, same saved field | 8.44426612e-15 | 2.92054080e-15 |

Sampling differences are at most 12.94 nm, below the predeclared 10 µm
criterion. RK step differences are roundoff-level consistency measurements,
not a certified physical ODE error bound. Independent synthetic polynomial
interpolation and RK-order unit tests supplement this actual FEM experiment.
There is no measured reason here to reduce future snapshot_dt below 0.005 s.
The dense-sampling/RK refinement window is **0…0.5 s**, not the full 5 s;
the medium/fine path comparison itself covers all 5 s.

## Qualification methods and units

Regional surface slopes are maxima of the exact linear derivative on each
P2 edge, including the points where a region cuts through an edge. Bulk is
abs(x)<=0.9, intermediate abs(x)<=0.98, contact 0.98<abs(x)<=1. No finite
differences of a visualization grid are used. The policy 0.1/0.3 is an
engineering criterion of this project, not a universal mathematical theorem.

Cross-mesh volume integration overlays the actual triangles of both graded
meshes, splits along both diagonals and the region boundaries, then applies
degree-four quadrature on every intersection triangle. Squared differences
of P2 velocity and discontinuous P1 FEM vorticity are integrated exactly up
to roundoff. Surface differences likewise use the union of both P2 partitions.
Each relative field error has denominator max_t ||fine|| over the comparison
window; it is not divided by a vanishing instantaneous initial velocity.
Fine is a reference, not an exact solution; two production meshes alone do
not measure an asymptotic convergence order.

Predeclared comparison screening limits are 1% for bulk eta, 2% for velocity,
and 10% for wall omega L2 differences relative to peak reference norms.
Fine must reduce maximum strong-divergence L2 both globally and in bulk.
These are transparent screening limits, separate from solver tolerances.
No threshold is changed to fit the subsequently measured data.

| Quantity | Units |
|---|---|
| eta, particle displacement | m |
| u,w | m/s |
| q (pressure divided by density) | m²/s² |
| omega | 1/s |
| Surface L2 norm | m^(3/2) |
| Velocity domain L2 norm | m²/s |
| Strong/weak divergence L2, omega L2 | m/s |
| Energy per density and unit out-of-plane span | m^4/s² |
| Viscous dissipation rate in this normalization | m^4/s³ |
| Surface integral, volume error per span | m² |
| Slope, relative divergence/errors | dimensionless |

The SDIRK budget is E(t)+D_cumulative(t)+Q_RK(t)-W(t)-E(0)=residual.
Production W=0. Q_RK is signed and is not physical heat; SDIRK2 is L-stable
but not algebraically stable. Physical energy monotonicity is separately
checked at every internal step, rather than inferred from the budget identity.

At t=0.1 s the actual medium probes give u(0,-0.1)=-3.01708e-4 m/s,
w(-0.75,-0.2)=+3.00272e-4 and w(+0.75,-0.2)=-3.00272e-4 m/s.
At t=1 s these directions reverse. The central deep probe u(0,-3) is
-3.38150e-6 at 0.1 s: penetration is much weaker at depth, without any
prescribed depth attenuation. These probes do not establish a universal
simultaneous deep return flow; such a claim must be read from the actual
field, not imposed from a schematic circulation picture.

## Contact-region limitations

There is a useful algebraic explanation, separate from the independent MMS
benchmark and the actual transient runs. Let eta_static have the prescribed
antisymmetric endpoint values and satisfy (S eta_static)_interior=0, where
S is the exact P2 surface mass matrix. Then C eta_static=0, so v=q=0 with
this trace is a discrete stationary state. Its cross-energy with any trace
having zero endpoints vanishes. Pinned evolution therefore retains this
localized component while the dynamic component dissipates. The implementation
constructs S from exact one-dimensional P2 integrals and tests C eta_static=0
against the actual FEM trace operator.

| Mesh | Static discrete max slope | Static potential energy, m^4/s² |
|---|---:|---:|
| medium | 0.254558024 | 7.25559930e-10 |
| fine | 0.387722573 | 4.76418619e-10 |

These numbers are algebraic diagnostics, **not** substituted 5 s solutions
and not proof that equilibrium has been reached. The actual medium energy
at 5 s is still far above its static floor. A decreasing energy and an
increasing, more localized corner gradient are therefore not contradictory.
This also cautions against claiming a mesh-independent limiting contact
slope merely by reducing alpha on one finite mesh.

The measured field/path gates support numerical interpretation of bulk surface
motion, bulk velocities, viscous energy exchange and damping, depth
penetration, and integral wall-vorticity structure away from pointwise corner
claims. Not resolved by this model: pointwise limiting corner omega or slope,
wall-film thickness/deposition, moving contact lines, wetting/dewetting and
wall drainage. No smoothing, clipping, artificial surface tension, slip or
fictitious film has been added.

Particle qualification is also limited by the linear model, independently
of interpolation accuracy. Expanding a path as X=X0+epsilon*X1+epsilon²*X2
shows that advection through epsilon*v1 produces an order-epsilon² term
grad(v1)*X1, while the order-epsilon² Eulerian velocity v2 is absent from
this PDE. Thus tiny net drift or apparent mixing over a cycle is not a
validated second-order transport prediction, even if the numerical RK path
itself is highly reproducible.

## Final qualification decision

**NOT QUALIFIED FOR PHYSICAL ANIMATION** for the requested alpha=0.02°
dataset under this project's global linear-small-slope policy.

The machine-readable [summary](validation_results/animation_qualification/summary.json)
has status `failed`. Its only failure is local slope >0.3 on fine: the first
saved violation is at 2.365 s, and max slope=0.37340763 at 4.320 s. Medium
alone would have suggested a conditional slope result (max 0.26973471);
the reference run changes the decision. Selecting medium merely to avoid
the higher fine-mesh contact slope would conceal unresolved structure.

This is **not** an unexplained numerical failure: both completed datasets
preserve no-slip, pinned points, volume and the SDIRK energy budget; physical
energy never increases at an internal step. Fine reduces strong divergence,
bulk field differences remain below 0.08% of peak fine norms, integral wall
vorticity differs by at most 5.731%, symmetry is maintained and all seeded
paths are stable. The separate numerical field gates therefore pass.

A further caution remains even apart from slope: medium relative strong
divergence exceeds 5% for 22.2% of active snapshots. Fine reduces this to
9.5%; neither Taylor–Hood field is pointwise divergence-free. The successful
particle experiment is evidence for these seeds and durations, not a proof
of exact local volume-preserving advection everywhere.

No smoothing, changed thresholds or physical-model modifications were used
to change this negative qualification decision. Completing this iteration
means reporting that decision, not granting an animation release.

## Recommended dataset for STEP 2 animation

Neither current alpha=0.02° dataset is globally qualified for release as a
linear-small-slope physical animation. **Fine** is the preferred numerical
reference for bulk interpretation and for the next confirmation, because
its strong divergence is smaller; medium's corner slope is not evidence
that it is physically safer.

Using the measured fine maximum and exact amplitude linearity,

    tan(alpha_new) = tan(0.02°) * target_slope / 0.37340762866770183

gives alpha_new=0.00535607717° for target 0.1, or 0.00267803859° for
target 0.05. These are **extrapolations, not additional calculated runs**.
A rounded conservative next candidate is:

| Parameter | Next candidate, NOT YET RUN / NOT YET QUALIFIED |
|---|---:|
| alpha | 0.0025° |
| nu | 0.01 m²/s |
| mesh | fine, 120×210 intervals |
| x_grading / z_grading | 1.8 / 4.0 |
| integrator | SDIRK2 |
| dt | 0.00125 s |
| snapshot_dt | 0.005 s |
| t_end | 5.0 s |
| a / d / g | 1 m / 10 m / 9.81 m/s² |
| Predicted maximum slope by tan scaling | 0.0466759517 |

A separate reduced-angle confirmatory dataset and qualification decision
are required before STEP 2. This iteration does not launch an automatic
angle-reduction loop or implement any renderer. Reducing alpha reduces
absolute gradients/displacements, but does **not** reduce relative spatial
errors, relative strong divergence, or establish a continuum corner limit.

Subject to that confirmation, bulk free-surface motion, bulk velocity,
energy exchange/damping, penetration depth and the tested away-from-corner
paths are appropriate numerical interpretation targets. Integral wall
vorticity is supported at the measured approximately 5.7% comparison level;
pointwise corner vorticity/slope remains unqualified. Wall films, moving
contact lines, wetting/dewetting and nonlinear net particle transport cannot
be inferred from this model. Full 5 s PDE time refinement was not repeated
here: the chosen integrator/dt inherit STEP 1.5 temporal validation, while
this iteration measures long-time spatial and saved-time path sensitivity.

## Reproducibility and artifacts

From `sloshing_visualization/`, with the installed requirements:

```bash
env OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    python scripts/qualify_animation_dataset.py --phase all
python -m pytest -q
python -m pytest -q -m validation
```

Individual phase commands and cache rules are documented in the
[README](README.md#step-16-commands-and-gates). A completed qualification
command may exit successfully while summary.status is `failed`: execution
success and scientific qualification are deliberately different outcomes.

Medium/fine full HDF5 files contain 1001 snapshots each and occupy about
1.073 / 2.365 GB. They are local ignored artifacts, not Git payloads.
Complete short prefixes are preserved in separate files. Failed, running,
incompatible or incomplete HDF5 files are rejected, and every FEM/probe
snapshot is checked for finite values on first access to a file revision.
Cache fingerprints include config, path, size and modification time, not
a cryptographic content hash. Recomputed summaries reuse compatible complete
data; dependency versions are recorded in HDF5. Across different library
versions compare tolerances, not a promise of bitwise identical LU results.

Compact tables: [dataset metrics](validation_results/animation_qualification/dataset_metrics.csv),
[full field comparison](validation_results/animation_qualification/comparison_full.csv),
[particle comparison](validation_results/animation_qualification/particle_mesh.csv),
[particle sampling](validation_results/animation_qualification/particle_sampling.csv).
Scientific figures (also available as PNG):
[surface](validation_results/animation_qualification/surface_profiles.pdf),
[slopes](validation_results/animation_qualification/slope_history.pdf),
[energy](validation_results/animation_qualification/energy_history.pdf),
[vorticity](validation_results/animation_qualification/vorticity_history.pdf),
[divergence](validation_results/animation_qualification/divergence_history.pdf),
[mesh comparison](validation_results/animation_qualification/medium_vs_fine.pdf),
[particles](validation_results/animation_qualification/particle_qualification.pdf),
[depth](validation_results/animation_qualification/depth_penetration.pdf),
[modal signal](validation_results/animation_qualification/modal_signal.pdf).
The complete local test inventory and measured outcomes are in
[TEST_RESULTS.md](validation_results/animation_qualification/TEST_RESULTS.md).
Final runs: 94 default tests passed (25.06 s), 2 opt-in validation tests
passed (91.17 s); a separate invocation of all 23 new tests passed (24.13 s).
Thus all 96 distinct tests, including the 73 existing tests, pass. Expected
warnings from deliberately under-resolved test meshes remain visible.
Rebuilding the summary from the same completed datasets produced exactly
equal parsed JSON values. All required PDF/PNG pairs and local document
links were verified. CI is configured, but has not been run on GitHub.
No MP4, Plotly, PyVista renderer or production animation was added.
