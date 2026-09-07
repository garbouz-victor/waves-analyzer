# STEP 3A.1 local test results

Base commit: `15a9e07c496c757cf40d35cd10fbc721afb9805a`, with the uncommitted
STEP 3A.1 repairs. These are code checks, not scientific qualification.

## Ordinary suite

From `sloshing_visualization/`:

```bash
.venv/bin/python -m pytest -q
```

**198 passed, 2 skipped, 2 deselected, 15 warnings in 24.94 s.**

The two skipped modules require DOLFINx and are exercised in the container below.
The two deselected tests are the existing expensive `validation` marker tests;
they were not rerun in this iteration. Warnings are the existing FontTools
deprecation and deliberately under-resolved/small-slope linear-model test cases.
No legacy production dataset was recalculated.

## Pinned FEniCSx suite

```bash
bash docker/step3/run.sh python3 -m pytest -q \
  tests/test_step3_model.py tests/test_step3_interface.py tests/test_step3_fenicsx.py \
  tests/test_step3a1_energy.py tests/test_step3a1_resolution.py \
  tests/test_step3a1_settling.py tests/test_step3a1_schedule.py \
  tests/test_step3a1_summary.py tests/test_step3a1_initial_wall.py \
  tests/test_step3a1_fenicsx.py
```

**85 passed, 1 warning in 9.52 s.** Some pure tests also run in the ordinary
suite; the counts are not disjoint.

The warning is intentional: the injected-failure/history-flush regression uses
an under-resolved tiny interface (1.2492 nominal transition cells) to test evidence
preservation, not physical convergence. The test checks that accepted samples
survive a solver exception and that execution is marked failed.

Container image ID:
`sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8`.
DOLFINx/Basix/FFCx 0.10.0, PETSc 3.24, Python 3.12.3; exact pinned stack is
documented in `docker/step3/stack.json`.

## Covered new regressions

- Monotone energy with bogus dissipation fails; NaN/inf raises; saved power
  integrals and energy components must agree; arbitrary initial-energy reference
  cancellation does not replace the positive characteristic scale.
- CG1 vertex/edge crossing and conservative CG2 Bernstein activation; curved
  interface, centroid false negatives, normal-width refinement and rotation.
- Physical-time settling is independent of sampling dt; small degrees per step
  cannot conceal a large physical rate; missing equilibrium metrics cannot pass.
- Scalar BE startup / constant-step BDF2 transition; no unequal-spacing BDF2
  history; actual CHNS checkpoint/restart during four scheme phases.
- Actual P1/P2 DOLFINx resolution, box refinement, read-only historical output,
  accepted-step history preservation, and tiny coupled solver smoke tests.
- Measured summary dependencies distinguish execution complete, single-history
  closure and scientific convergence-series qualification.
- Symbolic initial-wall variation shows why matching the zero-contour angle does
  not make a circular diffuse tanh profile a variational equilibrium.

## CI and scientific verdict

The ordinary and STEP 3A Docker workflow test commands pass locally. The Docker
workflow now includes the new cheap tests. **Remote GitHub Actions was not run for
these uncommitted changes**; no remote-green claim is made.

Re-running `scripts/step3_analyze.py` produced byte-identical `summary.json`,
compatible-series qualification and historical energy audit JSON. `git diff
--check` passed; the continuum weak form, free energy, material/wall laws,
initial-field formulas and historical `validation_results/step3/` have no diff.
The legacy fine HDF5 was also rehashed without modification:
`fa299d74c9487a40a5fac55966c1288c950bc0a46116f80e7da0063f794cf09a`.

The compatible-angle run completed numerically but failed the predeclared
all-time energy-change-relative budget gate. **MODEL NOT YET VALIDATED.**
Following the requested phase order, no subsequent Laplace/contact/film run or
success-conditioned commit was made. See `STEP3A1_REPORT.md` for measurements.
