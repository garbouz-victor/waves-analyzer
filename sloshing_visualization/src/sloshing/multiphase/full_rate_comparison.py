"""Read-only full L0 minus isolated L0, separately from isolated L0-minus-L2."""
from pathlib import Path
import numpy as np
from .phase_rate_history import read_json,load_archive,validate_complete,atomic_json
from .coupling_transfer import scalar_coupling,field_coupling,transfer_gates


def references(root,schedules):
    result=[]
    for i,schedule in enumerate(schedules):
        folder=Path(root)/f"level{i}"; marker=read_json(folder/"COMPLETE.json")
        if marker["identity"]["schedule_sha256"]!=schedule.sha256:
            raise ValueError("Historical COMPLETE schedule mismatch")
        rows=validate_complete(folder,marker["identity"],schedule)
        result.append({"identity":marker["identity"],"marker":marker,"rows":rows})
    return result


def compare(full_root,iso_root,schedules,M,phi_eq,mu_eq,policy):
    full_root=Path(full_root); iso_root=Path(iso_root)
    refs=references(iso_root,schedules); mark=read_json(full_root/"COMPLETE.json")
    full=validate_complete(full_root,mark["identity"],schedules[0])
    if mark["identity"]["schedule_sha256"]!=schedules[0].sha256 or (full_root/"FAILED.json").exists():
        raise ValueError("Continuous full L0 completion required")
    if not all(all(r["full_physical_checks"].values()) for r in full):
        raise ValueError("Independent full physical checks failed")
    norm=lambda x:float(np.sqrt(max(0.,x@(M@x))))
    initial,_=load_archive(iso_root/"level0/fields/000000")
    initial_full,_=load_archive(full_root/"fields/000000")
    if not np.array_equal(initial["phi"],initial_full["phi"]): raise ValueError("Initial phi differs")
    A_phi=norm(initial["phi"]-phi_eq); A_mu=norm(initial["mu"]-mu_eq)
    excess=policy["initial_excess"]; initial_D=refs[0]["rows"][0]["CH_dissipation"]
    scalars=[]
    for r,iso in zip(full,refs[0]["rows"]):
        if r["time"]!=iso["time"] or r["dt"]!=iso["dt"]: raise ValueError("No interpolation; identical L0 clock required")
        c=scalar_coupling(dict(D=r["CH_dissipation"],F=r["F_CH"],kinetic=r["E_kin"],
            visc=r["viscous_dissipation"],slip=r["slip_dissipation"],mass=r["phase_mass"]),
            dict(D=iso["CH_dissipation"],F=iso["E_total"],mass=iso["phase_mass"]),initial_D,excess,policy)
        scalars.append(dict(step=r["step"],time=r["time"],**c))
    fields=[]
    for parent in schedules[0].common_parent_indices:
        f,fm=load_archive(full_root/"fields"/f"{parent:06d}",mark["identity"])
        a,am=load_archive(iso_root/"level0/fields"/f"{parent:06d}",refs[0]["identity"])
        b,bm=load_archive(iso_root/"level2/fields"/f"{4*parent:06d}",refs[2]["identity"])
        if not fm["time"]==am["time"]==bm["time"]: raise ValueError("Common field time mismatch")
        metrics=field_coupling(norm(f["phi"]-a["phi"]),norm(b["phi"]-a["phi"]),norm(f["phi"]-b["phi"]),
            norm(a["phi"]-phi_eq),A_phi,norm(f["mu"]-a["mu"]),norm(a["mu"]-mu_eq),A_mu,policy)
        fields.append(dict(parent_step=parent,time=fm["time"],**metrics))
    last=full[-1]; iso0=refs[0]["rows"][-1]; iso2=refs[2]["rows"][-1]
    final={"I_hydro":last["cumulative_viscous_dissipation"]+last["cumulative_slip_dissipation"],
        "excess":excess,"I_CH_full":last["cumulative_CH_dissipation"],"I_iso0":iso0["cumulative_CH_dissipation"],
        "I_iso2":iso2["cumulative_CH_dissipation"],"F_CH_full":last["F_CH"],
        "F_iso0":iso0["E_total"],"F_iso2":iso2["E_total"],"budget":last["energy_budget_defect"],
        "DeltaE":last["E_total"]-full[0]["E_total"]}
    gates=transfer_gates(scalars,fields,final,policy)
    energy={"Delta_E_total":final["DeltaE"],"integral_D_CH":final["I_CH_full"],
        "integral_D_visc":last["cumulative_viscous_dissipation"],"integral_D_slip":last["cumulative_slip_dissipation"],
        "integral_D_total":final["I_CH_full"]+final["I_hydro"],"budget_defect":final["budget"],
        "budget_over_change":gates["measured"]["budget_over_change"],"sum_B_BE":last["work_sums"]["B_BE"],
        "max_weak_work_defect":max(abs(r["work"]["weak_work_defect"]) for r in full[1:])}
    return {"coupling":gates,"fields":fields,"scalars":scalars,"final":final,"energy":energy,
        "A_phi":A_phi,"A_mu":A_mu,"initial_excess":excess,"no_interpolation":True,
        "isolated_reference_markers":[r["marker"] for r in refs],"full_marker":mark}
