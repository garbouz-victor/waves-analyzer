import numpy as np
import pytest
from sloshing.config import SimulationConfig
from sloshing.fem_spaces import FEMSystem
from sloshing.validation.qualification.fem_evaluation import FEMGeometry,evaluate_map
from sloshing.visualization.data import grid_map,surface_extrema,prepare_data,HISTORY_FIELDS


def test_precomputed_actual_p2_map_matches_polynomial_without_relocation(monkeypatch):
    f=FEMSystem(SimulationConfig(nx=4,nz=8))
    g=FEMGeometry(f.mesh.p,f.mesh.t,f.scalar.element_dofs)
    x,z=np.linspace(-1,1,11),np.linspace(-1.5,0,9)
    mp=grid_map(g,x,z)
    monkeypatch.setattr(g,"locate",lambda *_: (_ for _ in ()).throw(AssertionError("unnecessary relocation")))
    xx,zz=np.meshgrid(x,z)
    for amp in (1.,.01):
        field=amp*(f.scalar.doflocs[0]**2+f.scalar.doflocs[1])
        np.testing.assert_allclose(evaluate_map(field,mp).reshape(xx.shape),amp*(xx**2+zz),atol=2e-13)


def test_surface_range_includes_non_nodal_extremum():
    x=np.array([-1.,0.,1.]);eta=(x-.25)**2
    lo,hi=surface_extrema(x,eta)
    np.testing.assert_allclose([lo,hi],[0,1.5625],atol=1e-15)


def test_missing_dataset_never_runs_solver(tmp_path):
    with pytest.raises(FileNotFoundError,match="renderer never runs PDE"):
        prepare_data(tmp_path/"absent.h5",tmp_path/"absent2.h5",tmp_path/"out")
    assert not (tmp_path/"out").exists()


def test_energy_history_uses_actual_solver_diagnostic_name(short_run):
    solver, states = short_run
    row=states[0][1]  # snapshots yield (State, measured diagnostics)
    assert all(key in row for key in (HISTORY_FIELDS["potential_energy"],HISTORY_FIELDS["kinetic_energy"],HISTORY_FIELDS["total_energy"]))


def test_render_interpolation_exact_wall_velocity():
    f=FEMSystem(SimulationConfig(nx=6,nz=12))
    g=FEMGeometry(f.mesh.p,f.mesh.t,f.scalar.element_dofs)
    rng=np.random.default_rng(25)
    full=f.expand_velocity(rng.normal(size=len(f.free)))
    p=np.array([(s,z) for s in (-1,1) for z in np.linspace(-10,0,73)]).T
    mp=g.interpolation(p)
    for ids in f.component_dofs:
        assert abs(evaluate_map(full[ids],mp)).max()<1e-12
