"""Independent read-only audit of LEFT Navier / RIGHT full no-slip.

No production Navier/LeftNavier FEM is imported.  Full wall operators, including
constrained rows, are checked: a spurious right Robin term cannot disappear
behind zero right velocity or behind the reduced variational residual.
"""
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from skfem import asm

from ..config import SimulationConfig
from ..diagnostics import quadratic_edge_max, surface_slopes
from ..fem_spaces import FEMSystem, scalar_mass
from ..validation.contact import surface_value
from .extension_verify import _csv, _energy, _first_peak, _hash, _json, _record_markers, _surface_export
from .mission import sha256
from .reference import load_native

MODEL_ID = "PW1_LEFT_NAVIER_RIGHT_NOSLIP_V1"
LABEL = "DECLARED_EXTENSION — LEFT NAVIER / RIGHT NO-SLIP"
SOURCE_KIND = "declared_left_Navier_right_noslip_linear_FEM"
CONTACT = "actual_FE_trace_left_free_right_fixed"
GAMMA = 1.0 - 1.0 / np.sqrt(2.0)


def _coordinate_check(actual, expected, axis_scale):
    """Only tolerate float64 elementary-function/averaging roundoff.

    The graded mesh uses tanh/expm1 and the surface uses P2 midpoints. NumPy
    1.21.5 versus 2.0.2 changes saved x coordinates by at most 2.22e-16 m in
    this case. Eight binary64 epsilons times each physical axis extent cover
    these few operations; this is NOT a mesh-accuracy or PDE tolerance.
    Shapes, finiteness and (elsewhere) topology/DOF identities remain exact.
    There is deliberately no allclose default relative tolerance.
    """
    actual, expected = np.asarray(actual), np.asarray(expected)
    tolerance = 8. * np.finfo(np.float64).eps * np.asarray(axis_scale, dtype=float)
    shape_exact = actual.shape == expected.shape
    finite = bool(np.isfinite(actual).all() and np.isfinite(expected).all())
    result = {"passed": False, "shape_exact": shape_exact, "finite": finite,
              "absolute_tolerance_m": tolerance.tolist(), "max_absolute_difference_m": None}
    if shape_exact and finite and np.isfinite(tolerance).all() and np.all(tolerance > 0):
        delta = abs(actual - expected)
        result.update(passed=bool(np.all(delta <= tolerance)),
                      max_absolute_difference_m=float(np.max(delta)))
    return result


def _physical(c):
    if "slip_length_m" in c or "right_slip_length_m" in c or "b_right" in c:
        raise ValueError("forbidden_symmetric_or_right_slip_parameter")
    return {"model_id": MODEL_ID, "a_m": c["a"], "d_m": c["d"], "g_m_s2": c["g"],
            "nu_m2_s": c["nu"], "alpha_deg": c["alpha_deg"], "left_slip_length_m": c["left_slip_length_m"],
            "sigma_N_m": 0., "bottom": "no_slip", "left_wall": "impermeable_Navier_tangential_slip",
            "right_wall": "no_slip", "geometry": "linearized_free_surface",
            "layer": "irreversible_lower_connected_zero_volume_one_way"}


def _operators(c):
    _physical(c)
    slip = float(c["left_slip_length_m"])
    if not np.isfinite(slip) or slip <= 0 or c["integrator"] != "sdirk2":
        raise ValueError("Finite positive LEFT slip and SDIRK2 required")
    f = FEMSystem(SimulationConfig(**{k: v for k, v in c.items() if k != "left_slip_length_m"}))
    left_normal = np.intersect1d(f.wall_dofs["left"], f.component_dofs[0])
    fixed = np.unique(np.r_[left_normal, f.wall_dofs["right"], f.wall_dofs["bottom"]])
    free = np.setdiff1d(np.arange(f.velocity.N), fixed)
    # Independent scalar assembly, unlike the production vector bilinear form.
    scalar = asm(scalar_mass, f.scalar.boundary("left", intorder=6)).tocoo()
    vertical = f.component_dofs[1]
    wall_full = coo_matrix((c["nu"] / slip * scalar.data,
        (vertical[scalar.row], vertical[scalar.col])), shape=(f.velocity.N, f.velocity.N)).tocsr()
    result = {"fem": f, "free": free, "fixed": fixed, "wall_full": wall_full,
              "M": f.M_full[free][:, free].tocsr(), "bulk": f.K_full[free][:, free].tocsr(),
              "wall": wall_full[free][:, free].tocsr(), "B": f.B_full[:, free].tocsr(),
              "R": f.R_full[:, free].tocsr()}
    result["C"] = (result["R"].T @ f.S).tocsr()
    return result


def _read_wall_matrix(h):
    group = h["operators/K_wall_full"]
    shape = tuple(int(value) for value in group.attrs["shape"])
    return csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)


def _operator_contract(h, op):
    """Return structural failures independently of the computed trajectory."""
    errors = []
    if not np.array_equal(h["mesh/free_velocity_dofs"][:], op["free"]):
        errors.append("wrong_asymmetric_free_DOFs_or_right_vertical_unconstrained")
    if not np.array_equal(h["operators/fixed_velocity_dofs"][:], op["fixed"]):
        errors.append("wrong_asymmetric_fixed_DOFs")
    if (list(h["operators"].attrs["endpoint_trace_nnz"]) != [1, 0]
            or h["operators"].attrs["wall_support"] != "left_only"):
        errors.append("wrong_asymmetric_operator_declaration")
    wall = _read_wall_matrix(h)
    expected = op["wall_full"]
    if wall.shape != expected.shape or not np.isfinite(wall.data).all():
        errors.append("invalid_full_wall_matrix")
    else:
        difference = (wall - expected).tocsr()
        scale = max(float(np.max(abs(expected.data))), 1e-30)
        defect = float(np.max(abs(difference.data))) if difference.nnz else 0.
        if defect > 1e-12 * scale:
            errors.append("full_wall_operator_is_not_positive_LEFT_ONLY")
        right = op["fem"].wall_dofs["right"]
        if (wall[right].nnz and np.max(abs(wall[right].data)) > 1e-14 * scale
                or wall[:, right].nnz and np.max(abs(wall[:, right].data)) > 1e-14 * scale):
            errors.append("forbidden_RIGHT_Navier_term_even_if_constrained")
    return errors


def _slope_witness(times, ids, x, eta):
    """Exact P2 edge derivative maxima, with reproducible accepted-state IDs."""
    left,middle,right=eta[:,:-2:2],eta[:,1::2],eta[:,2::2]
    widths=x[2::2]-x[:-2:2]
    values=np.stack([(-3*left+4*middle-right)/widths,(left-4*middle+3*right)/widths],axis=2)
    maximum_index=np.unravel_index(int(np.argmax(abs(values))),values.shape)
    violating=np.flatnonzero(np.any(abs(values)>.30,axis=(1,2)))
    def describe(index):
        step,edge,side=map(int,index)
        return {"time_s":float(times[step]),"state_id":ids[step],"accepted_step":step,
                "edge_index":edge,"edge_bounds_m":[float(x[2*edge]),float(x[2*edge+2])],
                "x_m":float(x[2*edge+2*side]),"one_sided_edge_endpoint":"left" if side==0 else "right",
                "eta_x":float(values[index]),"absolute_slope":float(abs(values[index])),"hard_limit":.30}
    first=None
    if len(violating):
        step=int(violating[0]); edge,side=np.unravel_index(int(np.argmax(abs(values[step]))),values[step].shape)
        first=describe((step,edge,side))
    return {"maximum_slope_witness":describe(maximum_index),"first_slope_screen_failure":first}


def surface_applicability(native):
    """Cheap raw-geometry diagnostic only: never claims an equation/BC audit."""
    with h5py.File(native,"r") as h:
        keys=sorted(h["accepted"],key=int)
        groups=[h["accepted"][key] for key in keys]
        x=h["mesh/surface_x"][:]
        times=np.array([group.attrs["time_s"] for group in groups])
        ids=[group.attrs["state_id"] for group in groups]
        eta=np.array([group["eta"][:] for group in groups])
        c=json.loads(h.attrs["config"])
        identity=json.loads(h.attrs["identity"])
    witness=_slope_witness(times,ids,x,eta)
    amplitude=quadratic_edge_max(eta[:,:-2:2],eta[:,1::2],eta[:,2::2])/c["a"]
    return {"scope":"RAW_GEOMETRY_SCREEN_ONLY_NOT_AN_OPERATOR_AUDIT","run_id":identity["run_id"],
            "max_eta_over_a":amplitude,"linear_applicability_passed":amplitude<=.05 and witness["first_slope_screen_failure"] is None,
            **witness}


def verify_native(native, *, heartbeat=lambda: None):
    """Audit a complete asymmetric trajectory; pilot horizon may be below 5 s."""
    path = Path(native)
    failures = []
    def require(condition, reason):
        if not condition and reason not in failures:
            failures.append(reason)
    report = {"passed": False, "numerical_checks_passed": False, "passed_equations_and_BCs": False,
              "linear_applicability_passed": False, "scope": "LEFT_NAVIER_RIGHT_NOSLIP_NATIVE", "model_id": MODEL_ID,
              "native": str(path)}
    try:
        d = load_native(path)
        c, identity = d["config"], d["identity"]
        physical = _physical(c)
        require(identity.get("model_id") == MODEL_ID and identity.get("physical_model_label") == LABEL
                and identity.get("source_kind") == SOURCE_KIND and identity.get("synthetic") is False,
                "not_explicit_asymmetric_source_identity")
        require(d["contact_method"] == CONTACT, "contact_is_not_asymmetric_actual_FE_trace")
        require(identity.get("physical") == physical and identity.get("physical_config_hash") == _hash(physical),
                "asymmetric_physical_identity_mismatch")
        require(all(c.get(key) == value for key, value in
                    {"a": 1., "d": 10., "g": 9.81, "nu": .01, "alpha_deg": 2.}.items()),
                "unapproved_change_to_corrective_baseline")
        require(_json(path.parent / "identity.json") == identity, "native_identity_file_mismatch")
        resolved = _json(path.parent / "resolved_case.json")
        require(resolved["controls"] == c and resolved["physical"] == physical
                and resolved["physical_model_label"] == LABEL, "resolved_case_mismatch")
        require(resolved.get("film_thickness_m") is None
                and resolved.get("retained_layer_volume_feedback") == "neglected_at_leading_order", "layer_model_changed")
        floor = resolved["extension_contract"]["record_height_tolerance_m"]
        require(floor == .0001, "undeclared_record_floor")
        source = identity["source_files"]
        require(identity["source_hash"] == _hash(source), "source_aggregate_hash")
        require({"pinned_wetting/asymmetric.py", "pinned_wetting/asymmetric_run.py"} <= set(source),
                "symmetric_source_reused_as_asymmetric")
        for name, digest in source.items():
            relative = Path(name)
            require(not relative.is_absolute() and ".." not in relative.parts, "invalid_source_path")
            if not relative.is_absolute() and ".." not in relative.parts:
                require(sha256(path.parent / "source_snapshot" / relative) == digest, "source_snapshot_integrity")
        op = _operators(c)
        f, R, dt = op["fem"], op["R"], c["dt"]
        times, eta, H = d["times"], d["eta"], d["H"]
        n = len(times)
        geometry_witness=_slope_witness(times,d["ids"],d["x"],eta)
        require(n == round(c["t_end"] / dt) + 1 and np.allclose(times, np.arange(n)*dt, atol=1e-12, rtol=0),
                "native_time_cadence_or_horizon")
        require(n > 1 and times[0] == 0., "missing_initial_or_advanced_state")
        surface_coordinates = _coordinate_check(d["x"], f.surface_x, c["a"])
        report["mesh_coordinate_roundoff"] = {"surface": surface_coordinates}
        require(surface_coordinates["passed"], "native_surface_mesh_mismatch")
        require(eta.shape == (n, len(f.surface_x)) and H.shape == (n, 2), "native_surface_shape")
        require(np.isfinite(eta).all() and np.isfinite(H).all(), "nonfinite_surface_or_memory")
        require(np.max(abs(eta[0] - f.config.slope*f.surface_x)) < 1e-12, "initial_not_straight")
        require(eta[0, -1] > eta[0, 0] and d["initial_velocity_max"] == 0., "initial_geometry_or_velocity")
        require(d["ids"] == [f"{identity['run_id']}:{i}" for i in range(n)] and all(d["accepted"]), "accepted_state_identity")
        require([R.getrow(i).nnz for i in (0, R.shape[0]-1)] == [1, 0], "endpoint_trace_not_left_free_right_fixed")
        require(np.max(abs(eta[:, -1] - eta[0, -1])) <= 1e-12, "RIGHT_native_contact_moves_despite_P2_label")
        require(np.max(abs(H[:, 1] - eta[0, -1])) <= 1e-12, "RIGHT_H_not_initial_z2")
        expected_H = np.maximum.accumulate(eta[:, [0, -1]], axis=0)
        require(np.allclose(H, expected_H, rtol=0, atol=1e-12), "memory_without_accepted_contact")
        intervals = np.stack([np.full_like(H, -c["d"]), H], axis=-1)
        require(np.allclose(d["intervals"], intervals, rtol=0, atol=1e-12), "dry_hole_or_wrong_W")
        catalogs = _record_markers(identity["run_id"], times, eta, H, floor)
        for i, coating in enumerate(d["coatings"]):
            require(coating["run_id"] == identity["run_id"] and coating["accepted_step"] == i
                    and abs(coating["time"]-times[i]) < 1e-12 and coating["bottom"] == -c["d"], "checkpoint_coating_identity")
            require(np.allclose(coating["H"], H[i], rtol=0, atol=1e-12), "checkpoint_coating_H")
            require(coating["markers"] == catalogs[i], "marker_not_native_resolved_record")
            require([m for m in coating["markers"] if m["side"] == "R"] == [catalogs[0][1]], "forbidden_new_RIGHT_record")
        maxima = {key: 0. for key in ("stage_momentum_relative", "weak_divergence", "mass_relative",
            "endpoint_kinematic_absolute", "split_energy_relative", "energy_budget_relative", "max_eta_over_a",
            "max_slope", "max_eta_over_left_b", "convective_indicator", "kinematic_indicator", "left_wall_traction_L2",
            "right_velocity_max", "left_normal_velocity_max", "bottom_velocity_max")}
        bulk_loss = wall_loss = numerical = 0.
        weights = (1-GAMMA, GAMMA)
        surface_basis = f.velocity.boundary("surface", intorder=6)
        scalar_surface = f.scalar.boundary("surface", intorder=6)
        left_basis = f.velocity.boundary("left", intorder=6)
        with h5py.File(path, "r") as h:
            require(bool(h.attrs.get("completed", False)), "native_not_complete")
            require(h.attrs.get("model_id") == MODEL_ID, "native_header_asymmetric_model_id")
            require(h.attrs.get("physical_model_label") == LABEL, "native_header_model")
            require(h.attrs.get("wetting_sampling") == "every_accepted_step_endpoint", "unknown_contact_sampling")
            require(sorted(h["accepted"], key=int) == [str(i) for i in range(n)], "noncontiguous_native_groups")
            bulk_coordinates = _coordinate_check(h["mesh/points"][:], f.mesh.p, [[c["a"]], [c["d"]]])
            topology_exact = np.array_equal(h["mesh/triangles"][:], f.mesh.t)
            report["mesh_coordinate_roundoff"].update(bulk=bulk_coordinates, topology_exact=topology_exact)
            require(bulk_coordinates["passed"] and topology_exact, "native_bulk_mesh_mismatch")
            for reason in _operator_contract(h, op):
                require(False, reason)
            if failures:
                report["failures"] = failures
                return report
            previous = h["accepted/0/velocity"][:]
            E0 = sum(_energy(op, previous, eta[0]))
            scale = max(E0, 1e-30)
            for i in range(n):
                if i % 25 == 0:
                    heartbeat()
                group = h["accepted"][str(i)]
                v = group["velocity"][:]
                require(v.shape == (len(op["free"]),) and np.isfinite(v).all(), "invalid_native_velocity")
                scalars = json.loads(group.attrs["integrator_scalars"])
                require(scalars["step"] == i and abs(scalars["time"]-times[i]) < 1e-12, "integrator_state_time")
                require(all(np.isfinite(value) for value in scalars.values()), "nonfinite_integrator_scalars")
                checked_velocities = [v]
                if i:
                    v1, p1, p2 = (group[key][:] for key in ("stage1_velocity", "stage1_pressure", "stage2_pressure"))
                    require(v1.shape == v.shape and p1.shape == p2.shape == (f.pressure.N,)
                            and all(np.isfinite(q).all() for q in (v1,p1,p2)), "invalid_stage_fields")
                    k1 = (v1-previous)/(GAMMA*dt)
                    k2 = (v-previous-dt*weights[0]*k1)/(GAMMA*dt)
                    e1 = eta[i-1]+GAMMA*dt*(R@v1)
                    predicted = eta[i-1]+dt*(weights[0]*(R@v1)+weights[1]*(R@v))
                    maxima["endpoint_kinematic_absolute"] = max(maxima["endpoint_kinematic_absolute"], float(np.max(abs(predicted-eta[i]))))
                    for vv,pp,kk,ee in ((v1,p1,k1,e1),(v,p2,k2,eta[i])):
                        terms=(op["M"]@kk,op["bulk"]@vv,op["wall"]@vv,c["g"]*(op["C"]@ee),-(op["B"].T@pp))
                        residual=float(np.linalg.norm(sum(terms)))/max(sum(np.linalg.norm(term) for term in terms),1e-30)
                        maxima["stage_momentum_relative"]=max(maxima["stage_momentum_relative"],residual)
                        weak=op["B"]@vv
                        maxima["weak_divergence"]=max(maxima["weak_divergence"],float(np.sqrt(max(0.,weak@f.pressure_mass_lu.solve(weak)))))
                    bulk_loss+=dt*(weights[0]*float(v1@(op["bulk"]@v1))+weights[1]*float(v@(op["bulk"]@v)))
                    wall_loss+=dt*(weights[0]*float(v1@(op["wall"]@v1))+weights[1]*float(v@(op["wall"]@v)))
                    numerical+=dt**2*GAMMA**2*(sum(_energy(op,k2,R@v))-sum(_energy(op,k1,R@v1)))
                    checked_velocities.append(v1)
                for vv in checked_velocities:
                    full=np.zeros(f.velocity.N)
                    full[op["free"]]=vv
                    for key,dofs in (("right_velocity_max",f.wall_dofs["right"]),
                        ("bottom_velocity_max",f.wall_dofs["bottom"]),
                        ("left_normal_velocity_max",np.intersect1d(f.wall_dofs["left"],f.component_dofs[0]))):
                        maxima[key]=max(maxima[key],float(np.max(abs(full[dofs]))))
                kinetic,potential=_energy(op,v,eta[i])
                energy=kinetic+potential
                remainder=energy-E0+bulk_loss+wall_loss+numerical
                maxima["energy_budget_relative"]=max(maxima["energy_budget_relative"],abs(remainder)/max(E0,abs(energy-E0)))
                for key,value in {"bulk_dissipated_energy":bulk_loss,"wall_dissipated_energy":wall_loss,
                    "dissipated_energy":bulk_loss+wall_loss,"rk_energy_correction":numerical,"work_energy":0.}.items():
                    maxima["split_energy_relative"]=max(maxima["split_energy_relative"],abs(scalars[key]-value)/scale)
                for key,value in {"kinetic_energy":kinetic,"gravitational_potential_energy":potential,
                    "cumulative_bulk_dissipation":bulk_loss,"cumulative_wall_dissipation":wall_loss,
                    "cumulative_viscous_dissipation":bulk_loss+wall_loss,"cumulative_work":0.,
                    "rk_energy_correction":numerical,"energy_balance_residual":remainder}.items():
                    require(abs(d["diagnostics"][i][key]-value)<=1e-8*scale+1e-14,"diagnostic_energy_not_native")
                maxima["mass_relative"]=max(maxima["mass_relative"],abs(float(f.surface_weights@(eta[i]-eta[0])))/(2*c["a"]*c["d"]))
                amplitude=quadratic_edge_max(eta[i,:-2:2],eta[i,1::2],eta[i,2::2])
                maxima["max_eta_over_a"]=max(maxima["max_eta_over_a"],amplitude/c["a"])
                maxima["max_eta_over_left_b"]=max(maxima["max_eta_over_left_b"],amplitude/c["left_slip_length_m"])
                maxima["max_slope"]=max(maxima["max_slope"],max(float(np.max(abs(s))) for s in surface_slopes(f.surface_x,eta[i])))
                full=np.zeros(f.velocity.N); full[op["free"]]=v
                field=f.velocity.interpolate(full)
                adv=np.einsum("ij...,j...->i...",field.grad,np.asarray(field))
                maxima["convective_indicator"]=max(maxima["convective_indicator"],float(np.max(np.linalg.norm(adv,axis=0)))/(c["g"]*max(abs(f.config.slope),1e-12)))
                trace=surface_basis.interpolate(full)
                scalar=np.zeros(f.scalar.N); scalar[f.surface_dofs]=eta[i]
                graph=scalar_surface.interpolate(scalar)
                omitted=np.asarray(graph)*trace.grad[1,1]-np.asarray(trace)[0]*graph.grad[0]
                maxima["kinematic_indicator"]=max(maxima["kinematic_indicator"],float(np.max(abs(omitted)))/(np.sqrt(c["g"]*c["a"])*max(abs(f.config.slope),1e-12)))
                wall_field=left_basis.interpolate(full)
                defect=c["nu"]*(-wall_field.grad[1,0]-wall_field.grad[0,1]+np.asarray(wall_field)[1]/c["left_slip_length_m"])
                maxima["left_wall_traction_L2"]=max(maxima["left_wall_traction_L2"],float(np.sqrt(np.sum(defect**2*left_basis.dx))))
                previous=v
        limits={"stage_momentum_relative":1e-9,"weak_divergence":1e-9,"mass_relative":1e-6,
                "endpoint_kinematic_absolute":1e-11,"split_energy_relative":1e-8,"energy_budget_relative":.05,
                "max_eta_over_a":.05,"max_slope":.30,"right_velocity_max":1e-12,
                "left_normal_velocity_max":1e-12,"bottom_velocity_max":1e-12}
        for key,limit in limits.items():
            require(maxima[key]<=limit,"limit:"+key)
        require(bulk_loss>=0 and wall_loss>0,"nonpositive_or_omitted_physical_losses")
        applicability_failures={"limit:max_eta_over_a","limit:max_slope"}
        numerical_passed=not any(reason not in applicability_failures for reason in failures)
        report.update(run_id=identity["run_id"],source_hash=identity["source_hash"],physical_config_hash=identity["physical_config_hash"],
            physical=physical,native_sha256=sha256(path),time_range=[float(times[0]),float(times[-1])],accepted_states=n,
            maxima=maxima,endpoint_trace_rows_nonzeros=[1,0],full_wall_support="left_only",right_eta_and_H_initial_z2=True,
            final_integrals_per_density={"bulk":bulk_loss,"LEFT_wall":wall_loss,"signed_RK":numerical},
            H_final_m=H[-1].tolist(),markers=catalogs[-1],
            numerical_checks_passed=numerical_passed,passed_equations_and_BCs=numerical_passed,
            linear_applicability_passed=not any(reason in applicability_failures for reason in failures),
            numerical_check_scope="All source/structural/accepted-field/equation/BC/coating/energy checks; only amplitude/slope screens excluded.",
            **geometry_witness,
            applicability_note="Hard amplitude/slope screens include the pinned RIGHT corner; no uniform physical-error bound is claimed.")
    except (KeyError,IndexError,ValueError,TypeError,OSError,FloatingPointError) as error:
        failures.append("native_schema_or_read_error:"+str(error))
    report["passed"],report["failures"]=not failures,failures
    return report


def comparison_metrics(first, second, *, same_physics=True):
    """Same exact P2/common-time norm, guarded by the asymmetric physical ID."""
    a,b=load_native(first),load_native(second)
    for d in (a,b):
        if d["identity"].get("model_id")!=MODEL_ID or d["identity"].get("physical")!=_physical(d["config"]):
            raise ValueError("symmetric_or_undeclared_native_in_asymmetric_comparison")
    if same_physics and _physical(a["config"])!=_physical(b["config"]):
        raise ValueError("mixed_asymmetric_physics_or_left_slip_in_refinement")
    ta,tb=np.round(a["times"],12),np.round(b["times"],12)
    common,ia,ib=np.intersect1d(ta,tb,return_indices=True)
    if len(common)<2 or common[0]!=0 or common[-1]!=ta[-1] or common[-1]!=tb[-1]:
        raise ValueError("refinement_common_horizon")
    edges=np.unique(np.r_[a["x"][::2],b["x"][::2]])
    query=np.column_stack([edges[:-1],.5*(edges[:-1]+edges[1:]),edges[1:]])
    normalization=abs(float(a["eta"][0,-1]-a["eta"][0,0]))
    maximum=0.
    for i,j in zip(ia,ib):
        error=surface_value(a["x"],a["eta"][i],query)-surface_value(b["x"],b["eta"][j],query)
        maximum=max(maximum,quadratic_edge_max(error[:,0],error[:,1],error[:,2]))
    pa,pb=_first_peak(a),_first_peak(b)
    if (pa is None)!=(pb is None):
        raise ValueError("first_peak_identity_unresolved")
    return {"normalized_surface_Linf":maximum/normalization,
        "normalized_macro_Linf":float(np.max(abs(a["eta"][ia][:,[0,-1]]-b["eta"][ib][:,[0,-1]])))/normalization,
        "normalized_H_Linf":float(np.max(abs(a["H"][ia]-b["H"][ib])))/normalization,
        "first_peak_time_fraction":None if pa is None else abs(pa[0]-pb[0])/pb[0],
        "first_peak_height_normalized":None if pa is None else abs(pa[1]-pb[1])/normalization}


def channel_reference(config):
    """Analytical LEFT-Robin/RIGHT-Dirichlet channel; no vessel bottom/top BC."""
    c={**config,"nx":8,"nz":12}
    op=_operators(c); f=op["fem"]
    vertical=f.component_dofs[1]
    free=np.setdiff1d(vertical,f.wall_dofs["right"])
    G=.001
    matrix=f.K_full+op["wall_full"]
    rhs=G*(f.M_full[free][:,vertical]@np.ones(len(vertical)))
    velocity=np.zeros(f.velocity.N); velocity[free]=spsolve(matrix[free][:,free],rhs)
    x=f.velocity.doflocs[0,vertical]; a=c["a"]; b=c["left_slip_length_m"]; nu=c["nu"]
    exact=G/(2*nu)*(a*a-x*x)+G*a*b/(nu*(b+2*a))*(a-x)
    maximum=float(np.max(abs(velocity[vertical]-exact)))
    left_value=2*G*a*a*b/(nu*(b+2*a))
    left_derivative=2*G*a*a/(nu*(b+2*a))
    return {"passed":maximum/max(float(np.max(abs(exact))),1e-30)<1e-9,
        "model_id":MODEL_ID,"benchmark":"steady_channel_LEFT_Navier_RIGHT_no_slip",
        "relative_error":maximum/max(float(np.max(abs(exact))),1e-30),"max_absolute_error_m_s":maximum,
        "left_Robin_residual_m_s":float(left_value-b*left_derivative),
        "right_Dirichlet_value_m_s":0.,"exact_formula":"G/(2nu)(a^2-x^2)+Ga*b/[nu(b+2a)]*(a-x)"}
