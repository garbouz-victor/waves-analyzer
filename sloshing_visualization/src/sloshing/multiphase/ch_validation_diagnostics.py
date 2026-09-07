"""Deterministic observation selection and cached *mesh-only* wall topology."""
import numpy as np
from .interface import quadratic_edge_roots
from .ch_timestep_planner import FrozenSchedule


class ValidationSchedule(FrozenSchedule):
    def __init__(self, plan, fingerprints, sparse_result, role, authorization, factor=1):
        if (sparse_result.get("qualification_status") != "passed" or
                sparse_result.get("plan_sha256") != plan["plan_sha256"] or
                sparse_result.get("input_fingerprints") != fingerprints):
            raise ValueError("Matching full-sparse verified plan required")
        if role not in ("isolated_CH_cost_pilot", "isolated_CH", "full_CHNS_probe"):
            raise ValueError("Unknown execution role")
        if not authorization.get(role+"_execution_authorized", False):
            raise ValueError("This scientific role is not authorized")
        super().__init__(plan, fingerprints, factor=factor, analysis_only=True)
        self.role = role


def observation_indices(blocks, tau_fast, tau_energy):
    ends = np.cumsum([b["steps"] for b in blocks]).tolist()
    times = np.concatenate([b["t_start"]+b["dt"]*np.arange(1, b["steps"]+1) for b in blocks])
    targets = [5*tau_fast, tau_energy, 5*tau_energy, 1e-5]
    indices = [int(np.searchsorted(times, t, side="left"))+1 for t in targets[:3]]
    indices.append(int(np.searchsorted(times, targets[3], side="right")))
    if not (all(0 < i <= len(times) for i in indices) and all(a < b for a, b in zip(indices, indices[1:]))):
        raise ValueError("Probe times must be distinct increasing existing checkpoints")
    n = len(times)
    audits = set(range(1, min(20, n)+1)) | set(ends)
    if n >= 21:
        audits.update(np.rint(np.geomspace(21, n, 40)).astype(int).tolist())
    return {"BE_work": sorted(audits), "geometry": sorted(set([0]+ends+indices)),
        "probe_indices": indices, "probe_times": [float(times[i-1]) for i in indices],
        "probe_target_times": targets, "block_end_indices": ends}


class WallTrace:
    """Exact P2 bottom-facet roots using cached scalar DOF/geometry maps."""
    def __init__(self, solver):
        from dolfinx import fem
        if solver.comm.size != 1 or solver.config.phase_degree != 2:
            raise ValueError("This STEP3A4 topology cache requires serial P2")
        self.solver = solver
        space, self.mapping = solver.space.sub(2).collapse()
        coords = space.tabulate_dof_coordinates()
        self.edges = []
        for facet in solver.facets["bottom"]:
            dofs = fem.locate_dofs_topological(space, 1, np.asarray([facet], dtype=np.int32))
            if len(dofs) != 3:
                raise ValueError("P2 facet must have exactly three DOFs")
            dofs = dofs[np.argsort(coords[dofs, 0])]
            x = coords[dofs, 0]
            if (not np.allclose(coords[dofs, 1], solver.config.z_min, atol=1e-12, rtol=0) or
                    not np.isclose(x[1], .5*(x[0]+x[2]), atol=1e-12, rtol=0)):
                raise ValueError("Unexpected bottom P2 edge geometry")
            self.edges.append((dofs, x[0], x[2]))

    def crossings(self):
        phi = self.solver.state.x.array[self.mapping]
        values = []
        for dofs, left, right in self.edges:
            a, mid, b = phi[dofs]
            if max(abs(a), abs(mid), abs(b)) < 1e-12:
                raise RuntimeError("Ambiguous identically-zero wall trace")
            values.extend(left+(right-left)*r for r in quadratic_edge_roots(a, mid, b))
        unique = []
        for x in sorted(values):
            if not unique or x-unique[-1] > 1e-9:
                unique.append(float(x))
        return unique
