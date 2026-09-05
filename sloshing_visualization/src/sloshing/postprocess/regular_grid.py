"""Reusable FEM probe operators; no finite differences of plotted samples."""

import numpy as np


class RegularGrid:
    def __init__(self, fem, vorticity):
        self.fem = fem
        self.xv = np.linspace(-fem.config.a, fem.config.a, fem.config.visualization_nx)
        self.zv = np.linspace(-fem.config.d, 0, fem.config.visualization_nz)
        xx, zz = np.meshgrid(self.xv, self.zv)
        self.shape = xx.shape
        points = np.vstack((xx.ravel(), zz.ravel()))
        # Chunk construction limits temporary element-finder memory on large meshes.
        from scipy.sparse import vstack
        def probes(basis):
            return vstack([basis.probes(points[:, i:i + 1024])
                           for i in range(0, points.shape[1], 1024)], format="csr")
        self.velocity_probes = probes(fem.scalar)
        self.pressure_probes = probes(fem.pressure)
        self.vorticity_probes = probes(vorticity.basis)

    def evaluate(self, velocity, pressure, omega):
        u, w = [np.asarray(self.velocity_probes @ velocity[ids]).reshape(self.shape)
                for ids in self.fem.component_dofs]
        return {"u": u, "w": w,
                "q": np.asarray(self.pressure_probes @ pressure).reshape(self.shape),
                "omega": np.asarray(self.vorticity_probes @ omega).reshape(self.shape)}
