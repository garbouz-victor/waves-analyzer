# STEP 3A.2 test results

Scientific verdict: **MODEL NOT YET VALIDATED**. Green software tests do not
override the failed measured perturbation energy/dissipation gates.

Final local checks, 2026-09-07, against the STEP 3A.2 working implementation
based on HEAD `04904c00c90b54da0ea03525399407ef2ce98b6c`:

| suite | passed | skipped | deselected | failures/errors | elapsed |
|---|---:|---:|---:|---:|---:|
| all ordinary pure tests | 230 | 8 | 2 | 0 | 21.46 s |
| STEP 3 Docker CI-equivalent | 123 | 0 | 0 | 0 | 9.16 s |

Machine-readable final logs: [tests_pure_final.xml](tests_pure_final.xml) and
[tests_docker_final.xml](tests_docker_final.xml). Earlier XML snapshots are
retained; the top-level summary uses only these final files.

## Ordinary pure suite

Run from `sloshing_visualization/`:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg .venv/bin/python -m pytest -q --junitxml=validation_results/step3a2/tests_pure_final.xml
```

Host Python 3.9. The eight skips are unavailable FEniCSx modules/cases in the
host environment; the two slow validation cases are deselected by the existing
ordinary-suite configuration. The 15 warnings concern existing small-slope
validation examples and fontTools, not STEP 3A.2 failures.

New coverage includes all requested quadratic extrema cases, labelled
near-singular fallback, stationary lines, a true Bernstein false positive,
350 seeded random polynomials checked at 5151 dense points each, and 250 random
gradient hulls checked at 4000 points each. Sampling is used only in tests.
Constant-gradient width, rotation covariance and origin-hull diameter fallback
are checked explicitly. Scalar-root bracketing, failure safety, cache identity
and exact backward-energy polynomial remainders are also covered.

## Pinned STEP 3 Docker suite

Run from the repository root:

```bash
bash sloshing_visualization/docker/step3/run.sh bash -lc 'python3 -m pytest -q tests/test_step3_model.py tests/test_step3_interface.py tests/test_step3_fenicsx.py tests/test_step3a1_*.py tests/test_step3a2_*.py --junitxml=validation_results/step3a2/tests_docker_final.xml --tb=short'
```

Python 3.12, DOLFINx 0.10.0, PETSc 3.24.0. Docker image:
`sha256:b1cd670d366e54d00cf53541dd69bdf260265cb293f34bc9c2c69669181c6ec8`.
The one warning is the deliberately under-resolved tiny scientific smoke mesh.

The new FEM cases cover neutral-wall uniform equilibria (including phi=1.02,
without clipping), constrained mass and Riesz residual, independent zero-mass
first variations, incomplete/tampered/incompatible prepared caches, constant-mu
initialization, stationary/transient chemical-block agreement, one and several
BE steps, standard-history BDF2 preservation, and actual FE P2 coefficients and
gradient reconstruction. These are tiny cases, not the expensive scientific
equilibrium or perturbation series.

The STEP 3 GitHub Actions workflow includes `tests/test_step3a2_*.py` and the
new config path. Its equivalent command passed locally. **Remote GitHub Actions
was not run; no push was made.** No remote CI status is inferred from local tests.

## Scientific measurements are separate

Two certified 60-degree equilibria and the 100-step BE/BDF2 preservation tests
pass. Three prescribed nonzero BE levels show decreasing endpoint and budget
errors, but fail integrated-dissipation convergence and final actual-change
closure. All scientific attempts, including failed preparations and the first
run's superseded qualification interpretation, remain recorded in
[STEP3A2_REPORT.md](../../STEP3A2_REPORT.md) and the per-run artifacts.

## Final read-only artifact checks

All 53 JSON documents parsed with nonstandard NaN/Infinity values rejected.
The two prepared H5 coefficient hashes, embedded/external metadata, complete
status, integrity fingerprints and current equilibrium-algorithm hashes agree.
All seven successful preparation/transient source archives reproduce the
recorded runtime source SHA256. Three failed preparation archives are retained,
but their original failure summaries lack complete runtime provenance; no
missing runtime metadata was invented.

All five transient histories include t=0 and exactly 100, 100, 10, 20 and 40
accepted steps, with matching geometry, resolution and nonlinear records. Every
state has at least 8.3281310758 certified cells and zero active-cell diameter or
polynomial-range fallbacks. The historical t=0.5 checkpoint still has SHA256
`14c5be2cb7d3d0c05fb8b902924a9a49c8c2efd4d9c99981747951d3d9c4a948`.
The three PDFs were regenerated from measurements and visually inspected.
`git diff --check` passed; the historical result trees and physical weak-form/
free-energy/material/boundary source files have no diff.
