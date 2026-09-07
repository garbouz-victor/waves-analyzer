"""DOLFINx rectangle with optional grading and pre-solve cost guard."""
import numpy as np
from mpi4py import MPI
from dolfinx import mesh as dmesh
from .free_energy import transition_width

WALL_IDS = {"left": 1, "right": 2, "bottom": 3, "top": 4}


def estimate_unknowns(config):
    vertices = (config.nx+1)*(config.nz+1)
    p2 = (2*config.nx+1)*(2*config.nz+1)
    phase = vertices if config.phase_degree == 1 else p2
    return {"cells": 2*config.nx*config.nz, "velocity_dofs": 2*p2,
            "pressure_dofs": vertices, "phi_dofs": phase, "mu_dofs": phase,
            "total_unknowns": 2*p2+vertices+2*phase,
            "matrix_memory_estimate_bytes": (2*p2+vertices+2*phase)*180*12,
            "memory_note": "rough 180 nnz/row CSR estimate; LU fill NOT included"}


def make_mesh(config, comm=MPI.COMM_WORLD):
    estimate = estimate_unknowns(config)
    if estimate["total_unknowns"] > config.max_unknowns:
        raise ValueError(f"Cost guard: {estimate}; raise max_unknowns explicitly after resource review")
    domain = dmesh.create_rectangle(comm, [[config.x_min, config.z_min],
                                           [config.x_max, config.z_max]],
                                    [config.nx, config.nz], dmesh.CellType.triangle)
    xyz = domain.geometry.x
    if config.x_grading:
        center, half = (config.x_min+config.x_max)/2, (config.x_max-config.x_min)/2
        s = (xyz[:, 0]-center)/half
        xyz[:, 0] = center+half*np.tanh(config.x_grading*s)/np.tanh(config.x_grading)
    if config.z_grading:
        # Cluster near z=0 when it is internal; otherwise near the interval center.
        center = 0.0 if config.z_min < 0 < config.z_max else (config.z_min+config.z_max)/2
        scale = np.where(xyz[:, 1] >= center, config.z_max-center, center-config.z_min)
        s = (xyz[:, 1]-center)/scale
        xyz[:, 1] = center+scale*np.sinh(config.z_grading*s)/np.sinh(config.z_grading)
    for level in range(config.refinement_levels):
        domain.topology.create_entities(1)
        domain.topology.create_connectivity(1,2)
        domain.topology.create_connectivity(1,0)
        edges=np.arange(domain.topology.index_map(1).size_local,dtype=np.int32)
        mid=dmesh.compute_midpoints(domain,1,edges)
        radius=np.hypot(mid[:,0]-config.refinement_circle_x,mid[:,1]-config.refinement_circle_z)
        selected=edges[np.abs(radius-config.refinement_circle_radius)<config.refinement_band_m]
        domain,_,_=dmesh.refine(domain,selected)
    xyz=domain.geometry.x
    # Recompute actual FE counts after grading/local refinement, before assembly.
    domain.topology.create_entities(1)
    nv=domain.topology.index_map(0).size_global
    ne=domain.topology.index_map(1).size_global
    np2=nv+ne
    nq=nv if config.phase_degree==1 else np2
    actual_total=2*np2+nv+2*nq
    estimate.update(cells=domain.topology.index_map(2).size_global,
                    velocity_dofs=2*np2,pressure_dofs=nv,phi_dofs=nq,mu_dofs=nq,
                    total_unknowns=actual_total,matrix_memory_estimate_bytes=actual_total*180*12)
    if actual_total>config.max_unknowns:
        raise ValueError(f"Refined mesh cost guard before FEM assembly: {estimate}")
    fdim = domain.topology.dim-1
    locators = {"left": lambda x: np.isclose(x[0], config.x_min),
                "right": lambda x: np.isclose(x[0], config.x_max),
                "bottom": lambda x: np.isclose(x[1], config.z_min),
                "top": lambda x: np.isclose(x[1], config.z_max)}
    facets = {side: dmesh.locate_entities_boundary(domain, fdim, fun)
              for side, fun in locators.items()}
    all_facets = np.concatenate(list(facets.values()))
    values = np.concatenate([np.full(len(facets[s]), WALL_IDS[s], dtype=np.int32)
                             for s in facets])
    order = np.argsort(all_facets)
    tags = dmesh.meshtags(domain, fdim, all_facets[order], values[order])
    ncells = domain.topology.index_map(2).size_local
    triangles = xyz[domain.geometry.dofmap[:ncells], :2]
    edges = np.linalg.norm(triangles-np.roll(triangles, 1, axis=1), axis=2)
    v1, v2 = triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]
    twice_area = np.abs(v1[:, 0]*v2[:, 1]-v1[:, 1]*v2[:, 0])
    min_alt = twice_area/edges.max(axis=1)
    hmax = comm.allreduce(float(edges.max()), op=MPI.MAX)
    hmin = comm.allreduce(float(edges.min()), op=MPI.MIN)
    aspect = comm.allreduce(float(np.max(edges.max(axis=1)/min_alt)), op=MPI.MAX)
    estimate.update(h_min_m=hmin, h_max_m=hmax, aspect_ratio_max=aspect,
                    epsilon_over_h_max=config.epsilon/hmax,
                    transition_cells_conservative=transition_width(config.epsilon)/hmax,
                    interface_resolution_policy="phi=-0.9..0.9 width / h_normal >=8; h_max is conservative")
    if comm.rank == 0:
        print(f"STEP3 mesh/cost: {estimate}", flush=True)
    return domain, tags, facets, estimate
