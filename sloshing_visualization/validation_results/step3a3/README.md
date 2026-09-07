# STEP3A3 measured checkpoint

Verdict: **MODEL NOT YET VALIDATED**. See [the report](../../STEP3A3_REPORT.md)
and [predeclared design](../../STEP3A3_DESIGN.md).

| Artifact | Scientific status |
| --- | --- |
| linear_operator | PASS, exact frozen historical input hashes |
| spectrum | Fast seeded Ritz evidence PASS; single polynomial exponential FAIL |
| spectrum/full_ritz_residual_audit.json | Independent actual operator residual PASS for all 33 converged significant vectors |
| krylov/restarted | Interrupted at finite cost guard; NOT qualified |
| krylov/rational | Unbalanced factorization prototype interrupted; NOT qualified |
| krylov/rational_balanced | Qualified exponential reference, M-Krylov 32/48/64 |
| low_mode | Isolated CH temporal/energy/BE-work PASS; NOT full CHNS |
| stiff_bump/planner | Reduced accuracy PASS, production cost STOP |
| stiff_bump/nonlinear_level0/1/2 | NOT RUN; full sparse verification also NOT RUN after cost STOP |
| linear_operator/initial_audit | t=0 only, exactly identical nonlinear initial field; ZERO accepted intervals |

The fixed candidate is 1796 intervals, estimated 4.65 hours for level0 and
32.54 hours for plan/half/quarter. Limits fixed before PDE runs were 3/10 hours.
No runtime tuning, threshold relaxation, physical change or full nonlinear
continuum qualification is implied by the linear candidate.

## Reproduction and storage

From the repository root, the tiny CI-equivalent checks use the pinned image:

```sh
bash sloshing_visualization/docker/step3/run.sh python3 -m pytest -q tests/test_step3a3_linear_ch.py tests/test_step3a3_krylov.py tests/test_step3a3_timestep_plan.py tests/test_step3a3_be_work.py
bash sloshing_visualization/docker/step3/run.sh python3 scripts/step3a3_report.py
```

The report command only reads measured inputs and regenerates STEP3A3 derived
summary/PDF/report artifacts. It does not execute a PDE or modify a plan.
Scientific stage scripts refuse to overwrite occupied result directories and
reject output outside STEP3A3. A new scientific replay needs a separate output
namespace/checkout, not deletion or overwriting of this checkpoint.

Tracked: recovered exact phi_eq/psi/initial-phi arrays, sparse M/K/H/mass,
qualified rational basis/reference (~41 MiB together), low-mode coefficients,
all accepted low-mode histories/BE work/SNES/resolution diagnostics, source
archives, all failed/interrupted-stage summaries, immutable candidate plan,
JUnit XML and figures. The large failed m=1024 basis/reference are local-only;
their hashes and all level summaries remain tracked. Replaying the separate
full Ritz audit requires that archived local basis (or regeneration of it).

The historical perturbation lacked a standalone coefficient file: its
unchanged initializer was replayed once, BOTH existing SHA256 values matched,
and the recovered arrays were frozen. Later stages only read those frozen
arrays. No new bump was designed.

Matrix/reference powers use q=psi, normalized by delta. Multiply those
quadratic energies/powers/integrals by delta²=1e-6 when comparing with the
physical delta=1e-3 nonlinear bump. The report performs that conversion.

Each completed calculation records actual source hash plus source archive,
base git SHA, pinned image, equilibrium/input/mesh fingerprints and quadrature.
Changes in analysis implementation are versioned in separate result folders;
old failed records were not relabelled PASS. Interrupted early prototypes have
archives and terminal status but incomplete final run provenance.

`tests/docker_invocation_error.xml` is an early local command with one incorrect
filename and zero tests executed, not a hidden test failure. Final full suites
are `tests/pure.xml` and `tests/docker.xml`. No GitHub-hosted workflow was
dispatched and nothing was pushed.

Physical heat is D_CH+D_visc+D_slip only. Signed B_BE is diagnostic numerical
work and is never added to physical dissipation. New stiff-bump BE work fields
are null/NOT RUN; measured low-mode work is explicitly labelled separately.
