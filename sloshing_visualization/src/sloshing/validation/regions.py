"""Exact integrals of DG1 vorticity on geometrically clipped triangles."""

import numpy as np


def clipped_mass(vertices, axis, bound, direction):
    polygon = [(vertices[i], np.eye(3)[i]) for i in range(3)]
    result = []
    for i, (point, bary) in enumerate(polygon):
        previous, old_bary = polygon[i-1]
        inside = direction*(point[axis]-bound) >= 0
        was_inside = direction*(previous[axis]-bound) >= 0
        if inside != was_inside:
            s = (bound-previous[axis])/(point[axis]-previous[axis])
            result.append((previous+s*(point-previous), old_bary+s*(bary-old_bary)))
        if inside:
            result.append((point, bary))
    matrix = np.zeros((3, 3))
    for j in range(1, len(result)-1):
        tri = [result[0], result[j], result[j+1]]
        p = np.array([v[0] for v in tri])
        bary = np.array([v[1] for v in tri])
        area = abs(np.linalg.det((p[1:]-p[0]).T))/2
        matrix += bary.T @ (area/12*(np.ones((3, 3))+np.eye(3))) @ bary
    return matrix


class VorticityRegions:
    def __init__(self, fem, projector):
        self.element_dofs = projector.basis.element_dofs.T
        vertices = fem.mesh.p[:, fem.mesh.t].transpose(2, 1, 0)
        areas = abs(np.linalg.det((vertices[:, 1:]-vertices[:, :1]).transpose(0, 2, 1)))/2
        full = areas[:, None, None]/12*(np.ones((3, 3))+np.eye(3))
        self.matrices = {"full": full}
        cuts = {"wall_left": (0, -.9*fem.config.a, -1),
                "wall_right": (0, .9*fem.config.a, 1),
                "surface_layer": (1, -min(1., fem.config.d), 1)}
        for name, (axis, bound, direction) in cuts.items():
            inside = direction*(vertices[:, :, axis]-bound) >= 0
            matrix = np.zeros_like(full)
            matrix[np.all(inside, axis=1)] = full[np.all(inside, axis=1)]
            partial = np.flatnonzero(np.any(inside, axis=1) & ~np.all(inside, axis=1))
            for i in partial:
                matrix[i] = clipped_mass(vertices[i], axis, bound, direction)
            self.matrices[name] = matrix

    def norms(self, omega):
        values = omega[self.element_dofs]
        return {key: float(np.sqrt(max(0., np.einsum("ei,eij,ej->", values, matrix, values))))
                for key, matrix in self.matrices.items()}
