"""Record actual container libraries; never infer them from a floating tag."""
import importlib
import json
import platform
from pathlib import Path

from mpi4py import MPI
from petsc4py import PETSc


def main():
    names = ("dolfinx", "petsc4py", "mpi4py", "ufl", "basix", "ffcx", "numpy",
             "scipy", "h5py", "matplotlib", "pytest", "sympy", "skfem")
    result = {name: importlib.import_module(name).__version__ for name in names}
    result.update(python=platform.python_version(), platform=platform.machine(),
                  petsc=list(PETSc.Sys.getVersion()), mpi=MPI.Get_library_version().strip(),
                  mumps=PETSc.Sys.hasExternalPackage("mumps"),
                  base_image="dolfinx/dolfinx@sha256:"
                  "f7cce2a2271bf838c080751348c471064acb41fef0330e2c08178a688f71890d")
    print(json.dumps(result, indent=2))
    if MPI.COMM_WORLD.rank == 0:
        Path("docker/step3/stack.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
