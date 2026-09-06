"""Direct P2/DG1 evaluation on the actual graded triangular FEM mesh."""

import numpy as np


class FEMGeometry:
    def __init__(self, coordinates, triangles, velocity_dofs, omega_dofs=None):
        self.p = np.asarray(coordinates)
        self.t = np.asarray(triangles)
        self.velocity_dofs = np.asarray(velocity_dofs)
        self.omega_dofs = np.asarray(omega_dofs) if omega_dofs is not None else None
        self.vertices = self.p[:, self.t].transpose(2, 1, 0)
        self.origin = self.vertices[:, 0]
        self.inverse = np.linalg.inv((self.vertices[:, 1:]-self.vertices[:, :1]).transpose(0, 2, 1))
        self.x = np.unique(self.p[0, self.p[1] == self.p[1].min()])
        self.z = np.unique(self.p[1, self.p[0] == self.p[0].min()])
        self.nz = len(self.z)-1
        center = self.vertices.mean(axis=1)
        cells = self.cell_indices(center.T)
        self.cell_triangles = np.full(((len(self.x)-1)*self.nz, 4), -1, dtype=np.int32)
        occupied = np.zeros(len(self.cell_triangles), dtype=int)
        for tri, cell in enumerate(cells):
            if occupied[cell] == 4:
                raise ValueError("Unsupported non-tensor triangulation")
            self.cell_triangles[cell, occupied[cell]] = tri
            occupied[cell] += 1
        if np.any(occupied == 0):
            raise ValueError("Mesh has an empty tensor cell")

    @classmethod
    def from_hdf5(cls, h):
        m = h["mesh"]
        return cls(m["coordinates"][:], m["triangles"][:], m["velocity_element_dofs"][:], m["omega_element_dofs"][:])

    def cell_indices(self, points):
        x, z = np.asarray(points)
        if np.any((x < self.x[0]-1e-13) | (x > self.x[-1]+1e-13) | (z < self.z[0]-1e-13) | (z > self.z[-1]+1e-13)):
            raise ValueError("FEM query outside reference domain")
        # Restrict cell INDEX at closed-domain endpoints; never change coordinates.
        i = np.minimum(np.maximum(np.searchsorted(self.x, x, side="right")-1, 0), len(self.x)-2)
        j = np.minimum(np.maximum(np.searchsorted(self.z, z, side="right")-1, 0), len(self.z)-2)
        return i*self.nz+j

    def locate(self, points):
        points = np.asarray(points)
        if points.ndim != 2 or points.shape[0] != 2:
            raise ValueError("Expected points shape (2,n)")
        candidates = self.cell_triangles[self.cell_indices(points)]
        uv = np.einsum("ncij,ncj->nci", self.inverse[candidates], points.T[:, None]-self.origin[candidates])
        bary = np.concatenate((1-uv.sum(axis=2, keepdims=True), uv), axis=2)
        valid = (candidates >= 0) & np.all(bary >= -2e-10, axis=2)
        if not np.all(np.any(valid, axis=1)):
            raise ValueError("No FEM element contains query")
        selected = np.argmax(valid, axis=1)
        n = np.arange(points.shape[1])
        return candidates[n, selected], bary[n, selected]

    def interpolation(self, points, degree=2):
        elements, bary = self.locate(points)
        if degree == 1:
            if self.omega_dofs is None:
                raise ValueError("DG1 dof mapping is missing")
            return self.omega_dofs[:, elements].T, bary
        if degree != 2:
            raise ValueError("Supported degrees: P2 velocity and DG1 curl")
        l0, l1, l2 = bary.T
        weights = np.column_stack((l0*(2*l0-1), l1*(2*l1-1), l2*(2*l2-1),
                                   4*l0*l1, 4*l1*l2, 4*l0*l2))
        return self.velocity_dofs[:, elements].T, weights

    def evaluate(self, coefficients, points, degree=2):
        ids, weights = self.interpolation(points, degree)
        return np.sum(np.asarray(coefficients)[ids]*weights, axis=1)


def evaluate_map(coefficients, mapping):
    ids, weights = mapping
    return np.sum(coefficients[ids]*weights, axis=1)
