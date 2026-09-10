"""Tensor graded coordinates, triangulated without modifying the domain."""

import numpy as np
from skfem import MeshTri


def generate_mesh(config):
    nx, nz = config.resolution
    s = np.linspace(-1, 1, nx + 1)
    x = config.a * (np.tanh(config.x_grading * s) / np.tanh(config.x_grading)
                    if config.x_grading else s)
    r = np.linspace(0, 1, nz + 1)
    depth = config.d * (np.expm1(config.z_grading * r) / np.expm1(config.z_grading)
                        if config.z_grading else r)
    z = -depth[::-1]
    # Exact endpoint coordinates avoid ambiguous boundary selection.
    x[[0, -1]] = [-config.a, config.a]
    z[[0, -1]] = [-config.d, 0.0]
    # Reflect diagonals across x=0 to avoid a directional bias at the contacts.
    xx, zz = np.meshgrid(x, z, indexing="ij")
    points = np.vstack((xx.ravel(), zz.ravel()))
    ids = np.arange(points.shape[1]).reshape(nx + 1, nz + 1)
    triangles = []
    extra_points = []
    for i in range(nx):
        for j in range(nz):
            ll, lr = ids[i, j], ids[i + 1, j]
            ul, ur = ids[i, j + 1], ids[i + 1, j + 1]
            if nx % 2 and i == nx // 2:
                # An odd Nx needs four symmetric triangles in the center column.
                center = points.shape[1] + len(extra_points)
                extra_points.append(((x[i] + x[i+1])/2, (z[j] + z[j+1])/2))
                triangles.extend(((ll, lr, center), (lr, ur, center),
                                  (ur, ul, center), (ul, ll, center)))
            elif i < nx // 2:
                triangles.extend(((ll, lr, ur), (ll, ur, ul)))
            else:
                triangles.extend(((ll, lr, ul), (lr, ur, ul)))
    if extra_points:
        points = np.hstack((points, np.asarray(extra_points).T))
    mesh = MeshTri(points, np.asarray(triangles, dtype=np.int32).T)
    return mesh.with_boundaries({
        "left": lambda p: p[0] == -config.a,
        "right": lambda p: p[0] == config.a,
        "bottom": lambda p: p[1] == -config.d,
        "surface": lambda p: p[1] == 0.0,
    })
