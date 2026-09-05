"""Streaming HDF5 snapshots and a separate per-run diagnostics.csv."""

import csv
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform

import h5py
import numpy as np

from .diagnostics import PINNED_MESSAGE
from .postprocess.regular_grid import RegularGrid
from .postprocess.vorticity import VorticityProjector


class SnapshotWriter:
    def __init__(self, path, solver, visualization_grid=True):
        self.path = Path(path)
        self.solver = solver
        self.fem = solver.fem
        self.projector = VorticityProjector(self.fem)
        self.grid = RegularGrid(self.fem, self.projector) if visualization_grid else None
        self.file = None
        self.csv_file = None
        self.count = 0

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Never silently overwrite a reproducible run or its diagnostics.
        diagnostics_dir = self.path.with_suffix("")
        if self.path.exists() or (diagnostics_dir / "diagnostics.csv").exists():
            raise FileExistsError(f"Run already exists: {self.path}; choose another output path")
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        self.file = h5py.File(self.path, "x")
        try:
            self.csv_file = (diagnostics_dir / "diagnostics.csv").open("x", newline="", encoding="utf-8")
            self._metadata()
        except BaseException:
            if self.csv_file is not None:
                self.csv_file.close()
            self.file.close()
            raise
        return self

    def _metadata(self):
        f, h = self.fem, self.file
        h.attrs.update({
            "schema_version": "1.5", "status": "running", "config_json": self.solver.config.to_json(),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "python_version": platform.python_version(),
            "dependencies_json": json.dumps({name: importlib.metadata.version(name)
                                              for name in ("numpy", "scipy", "scikit-fem", "h5py")}),
            "model": "linearized incompressible Navier–Stokes; fixed domain; sigma=0",
            "integrator": self.solver.config.integrator+", monolithic velocity/surface",
            "pressure_convention": "q at snapshot time; constant determined by free-surface traction",
            "vorticity_convention": "omega = dw/dx - du/dz; elementwise DG P1, exact derivative of P2",
            "energy_units": "m^4/s^2; per unit density and out-of-plane span",
            "contact_line": PINNED_MESSAGE,
        })
        for key, data in {"coordinates": f.mesh.p, "triangles": f.mesh.t,
                          "velocity_dof_coordinates": f.scalar.doflocs,
                          "velocity_element_dofs": f.scalar.element_dofs,
                          "pressure_dof_coordinates": f.pressure.doflocs,
                          "pressure_element_dofs": f.pressure.element_dofs,
                          "surface_x": f.surface_x,
                          "omega_dof_coordinates": self.projector.basis.doflocs,
                          "omega_element_dofs": self.projector.basis.element_dofs}.items():
            h.create_dataset(f"mesh/{key}", data=data)
        for name, facets in f.mesh.boundaries.items():
            h.create_dataset(f"mesh/boundaries/{name}", data=facets)
        if self.grid is not None:
            h.create_dataset("visualization/xv", data=self.grid.xv)
            h.create_dataset("visualization/zv", data=self.grid.zv)
        a, d = f.config.a, f.config.d
        self.probe_points = np.array([[0, -.75*a, .75*a, 0, 0, 0],
                                      [-min(.1, d*.05), -min(.2, d*.1), -min(.2, d*.1),
                                       -min(1., d*.25), -min(3., d*.5), -.9*d]])
        self.probes = f.scalar.probes(self.probe_points)
        h.create_dataset("probes/coordinates", data=self.probe_points)

    def _append(self, key, value):
        data = np.asarray(value)
        if key not in self.file:
            options = {"compression": "gzip", "compression_opts": 1, "shuffle": True} if data.ndim else {}
            self.file.create_dataset(key, shape=(0,) + data.shape, maxshape=(None,) + data.shape,
                                     dtype=data.dtype, **options)
        ds = self.file[key]
        ds.resize(self.count + 1, axis=0)
        ds[self.count] = data

    def append(self, state, row):
        f = self.fem
        velocity = f.expand_velocity(state.velocity)
        pressure = self.solver.integrator.instantaneous_pressure(state)
        omega = self.projector.evaluate(velocity)
        self._append("times", state.time)
        self._append("fields/u", velocity[f.component_dofs[0]])
        self._append("fields/w", velocity[f.component_dofs[1]])
        self._append("fields/q", pressure)
        self._append("fields/eta", state.eta)
        self._append("fields/omega", omega)
        self._append("probes/u", self.probes @ velocity[f.component_dofs[0]])
        self._append("probes/w", self.probes @ velocity[f.component_dofs[1]])
        if self.grid is not None:
            for name, data in self.grid.evaluate(velocity, pressure, omega).items():
                self._append(f"visualization/{name}", data)
        for name, value in row.items():
            self._append(f"diagnostics/{name}", value)
        if self.count == 0:
            self.csv_writer = csv.DictWriter(self.csv_file, fieldnames=list(row))
            self.csv_writer.writeheader()
        self.csv_writer.writerow(row)
        self.csv_file.flush()
        self.count += 1
        self.file.attrs["snapshots_written"] = self.count
        self.file.flush()

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.attrs["status"] = "complete" if exc_type is None else "failed"
        if exc_type is not None:
            self.file.attrs["failure"] = f"{exc_type.__name__}: {exc_value}"
        self.csv_file.close()
        self.file.close()


def save_run(solver, path, visualization_grid=True, progress=None):
    with SnapshotWriter(path, solver, visualization_grid) as writer:
        for state, row in solver.snapshots():
            writer.append(state, row)
            if progress is not None:
                progress(state, row)
    return Path(path)
