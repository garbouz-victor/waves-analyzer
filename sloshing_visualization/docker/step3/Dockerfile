# Official release, immutable multi-architecture manifest. Benchmark platform: amd64.
FROM dolfinx/dolfinx@sha256:f7cce2a2271bf838c080751348c471064acb41fef0330e2c08178a688f71890d
USER root
RUN python3 -m pip install --no-cache-dir \
    h5py==3.12.1 matplotlib==3.9.4 sympy==1.13.3 pytest==8.3.5 scikit-fem==11.0.0
# scikit-fem is required by the historical sloshing package's eager __init__ only.
# The multiphase discretization imports DOLFINx, not that solver.
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg \
    PYTHONPATH=/work/src:${PYTHONPATH} XDG_CACHE_HOME=/work/validation_results/step3/cache \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR /work
CMD ["python3", "scripts/step3_benchmarks.py", "--help"]
