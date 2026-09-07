# STEP 3A software verification checkpoint

Base `c6dcf3faf28ff2cf9f303e5d5d8832d9ed2f7517`, with the uncommitted STEP 3A
source documented in STEP3_MODEL_REPORT.md. These are software tests, **not** a
positive model-qualification verdict.

| Command, from sloshing_visualization | Measured result | Runtime |
|---|---|---:|
| `.venv/bin/python -m pytest -q` | 148 passed, 1 skipped module, 2 deselected, 15 warnings | 28.71 s |
| `.venv/bin/python -m pytest -q -m validation` | 2 passed, 1 skipped module, 1 warning | 57.65 s |
| `bash docker/step3/run.sh python3 -m pytest -q tests/test_step3_model.py tests/test_step3_interface.py tests/test_step3_fenicsx.py` | 26 passed | 6.19 s |

The module skipped outside Docker requires DOLFINx. Its six coupled FEM tests
run in the pinned container. Twenty pure tests overlap between the ordinary and
container commands; the counts must not simply be added as distinct tests.

New checks cover analytic surface-energy normalization, CH/density mass identity,
momentum flux orientation, gravity/pressure transform, Young-angle signs, material
admissibility, tilt sign, BDF derivative, film-resolution helper, polynomial
interface roots/connectivity and angle fits, analysis-only CLI safety, tiny
coupled mass/energy behavior, exact BDF checkpoint restart and wetting smoke tests.

Old small-slope/strong-divergence warnings and a fontTools deprecation warning
remain visible; they were not suppressed. Historical tests and solver files
were not edited. The two expensive validation tests were collected before the
final CLI-only tests were added; their underlying solver/tests did not change.

Local container success is **not** a remote GitHub Actions result. No changes
have been committed or pushed because the scientific gate is still
**MODEL NOT YET VALIDATED**. See summary.json for failed and missing benchmarks.
