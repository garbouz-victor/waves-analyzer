"""Pure root/cache identity checks and tiny serial FEniCSx equilibrium smoke."""
from pathlib import Path
import numpy as np
import pytest
from sloshing.multiphase.equilibrium import safeguarded_scalar_root, equilibrium_config_fingerprint
from sloshing.multiphase.config import ModelConfig


def config(**changes):
    path=Path(__file__).resolve().parents[1]/"configs/step3/benchmark_flat.json"
    return ModelConfig.from_json(path).changed(nx=8,nz=8,dt=1e-5,t_end=4e-5,
        snes_atol=1e-10,snes_rtol=1e-9,**changes)


def test_scalar_root_synthetic_monotone():
    root,history=safeguarded_scalar_root(lambda x:np.exp(x)-2,(0.,2.))
    assert root==pytest.approx(np.log(2),abs=1e-12)
    assert all(0<=r["mu_star"]<=2 for r in history)


def test_failed_bracket_does_not_continue():
    evaluated=[]
    def f(x):
        evaluated.append(x)
        return x*x+1
    with pytest.raises(RuntimeError,match="not bracketed"):
        safeguarded_scalar_root(f,(-1,1))
    assert evaluated==[-1.,1.]


def test_nonfinite_or_unconverged_root_cannot_pass():
    with pytest.raises(RuntimeError,match="Nonfinite"):
        safeguarded_scalar_root(lambda x:np.nan,(0,1))
    with pytest.raises(RuntimeError,match="did not converge"):
        safeguarded_scalar_root(lambda x:x**3-.123,(0,1),max_iterations=1)


def test_equilibrium_identity_ignores_only_time_and_execution_controls():
    c=config()
    assert equilibrium_config_fingerprint(c)==equilibrium_config_fingerprint(c.changed(dt=2e-5,time_scheme="be"))
    for update in ({"sigma":.2},{"nx":10},{"quadrature_degree":14},{"theta_equilibrium_deg":60.},
                   {"phase_degree":1},{"wetting_walls":()}):
        assert equilibrium_config_fingerprint(c)!=equilibrium_config_fingerprint(c.changed(**update))


def fem_solver(c):
    pytest.importorskip("dolfinx",reason="tiny equilibrium test requires pinned FEniCSx container")
    from sloshing.multiphase.solver import CHNSSolver
    return CHNSSolver(c)


@pytest.mark.parametrize("value",[1.,1.02])
def test_uniform_neutral_equilibrium_mass_constant_mu_and_no_clipping(value,tmp_path):
    from sloshing.multiphase.equilibrium import solve_equilibrium,save_prepared,load_prepared
    from sloshing.multiphase.free_energy import bulk_derivative
    c=config(theta_equilibrium_deg=90.,g=0.,a_x=0.)
    s=fem_solver(c)
    from sloshing.multiphase.equilibrium_diagnostics import transient_term_residuals
    s.initialize(lambda x:np.full(x.shape[1],value))
    p=solve_equilibrium(s)
    assert p.metadata["mass_residual_relative_to_domain"]<1e-10
    assert p.metadata["stationarity"]["normalized"]<1e-8
    assert p.phi_coefficients==pytest.approx(value,abs=1e-11)
    assert p.metadata["mu_star"]==pytest.approx(bulk_derivative(value,c.sigma,c.epsilon),abs=1e-11)
    s.initialize_prepared_equilibrium(p)
    assert np.ptp(s.state.sub(3).collapse().x.array)==0.
    assert s.old.x.array==pytest.approx(s.state.x.array,abs=0)
    assert s.older.x.array==pytest.approx(s.state.x.array,abs=0)
    assert np.max(abs(s.state.sub(0).collapse().x.array))==0.
    terms=transient_term_residuals(s,p.metadata["mu_star"])
    assert terms["chemical_plus_stationary_algebraic_norm"]<1e-12
    assert terms["phase_algebraic_norm_diagnostic"]<1e-11
    assert terms["momentum_algebraic_norm_diagnostic"]<1e-12
    save_prepared(tmp_path,p)
    loaded=load_prepared(tmp_path,s,p.metadata["target_mass"])
    assert loaded.phi_coefficients==pytest.approx(p.phi_coefficients,abs=0)
    with pytest.raises(ValueError,match="target mass"):
        load_prepared(tmp_path,s,p.metadata["target_mass"]+1e-8)
    loaded.phi_coefficients[0]+=1e-7
    with pytest.raises(ValueError,match="integrity"):
        s.initialize_prepared_equilibrium(loaded)
    loaded=load_prepared(tmp_path,s,p.metadata["target_mass"])
    loaded.metadata["status"]="running"
    with pytest.raises(ValueError,match="Incomplete"):
        s.initialize_prepared_equilibrium(loaded)
    s.config=c.changed(quadrature_degree=14)
    with pytest.raises(ValueError,match="config mismatch"):
        s.initialize_prepared_equilibrium(p)
    s.config=c


@pytest.mark.parametrize("scheme",["be","bdf2"])
def test_nonuniform_equilibrium_preservation_one_and_several_steps(scheme):
    from sloshing.multiphase.equilibrium import solve_equilibrium
    from sloshing.multiphase.initialization import flat_interface
    s=fem_solver(config(theta_equilibrium_deg=90.,g=0.,a_x=0.,time_scheme=scheme))
    from sloshing.multiphase.diagnostics import Diagnostics
    from sloshing.multiphase.equilibrium_diagnostics import PreservationDiagnostics,transient_term_residuals
    s.initialize(flat_interface(s.config.epsilon))
    p=solve_equilibrium(s)
    assert np.ptp(p.phi_coefficients)>.5
    s.initialize_prepared_equilibrium(p)
    base=s.state.sub(2).collapse()
    d=Diagnostics(s);drift=PreservationDiagnostics(s,base)
    initial=d.measure()
    assert initial["CH_dissipation"]<1e-16
    terms=transient_term_residuals(s,p.metadata["mu_star"])
    assert terms["chemical_plus_stationary_algebraic_norm"]<1e-12
    for _ in range(s.schedule.nsteps):
        s.advance();row=d.measure();delta=drift.measure()
        assert delta["relative_phi_L2_drift"]<1e-8
        assert row["speed_max_dof_sample"]<1e-8
        assert abs(row["mass_error_relative_to_domain"])<1e-10
        assert row["CH_dissipation"]<1e-16


def test_zero_mass_energy_first_variations():
    from sloshing.multiphase.equilibrium import solve_equilibrium,first_variations
    s=fem_solver(config(theta_equilibrium_deg=90.,g=0.,a_x=0.,z_min=0,z_max=.6,x_min=-.6,x_max=.6).changed(nx=24,nz=12))
    s.initialize(lambda x:np.ones(x.shape[1]))
    p=solve_equilibrium(s)
    s.initialize_prepared_equilibrium(p)
    variations=first_variations(s,s.state.sub(2).collapse())
    assert len(variations)==3
    for row in variations:
        assert abs(row["mass"])<1e-14
        assert row["normalized_abs_derivative"]<1e-8
        assert row["second_difference"]>=-1e-14


def test_fem_actual_p2_extrema_and_constant_gradient_fallback():
    s=fem_solver(config(mesh_diagonal="crossed"))
    from sloshing.multiphase.diagnostics import interface_resolution,resolution_cell_data
    s.state.sub(2).interpolate(lambda x:100*((x[0]-.06)**2+(x[1]-.06)**2)+1)
    s.state.x.scatter_forward()
    r=interface_resolution(s)
    assert not r["interface_present"]
    s.state.sub(2).interpolate(lambda x:np.full(x.shape[1],.345))
    s.state.x.scatter_forward()
    data=resolution_cell_data(s)
    assert (data[5]["diameter_fallback"]).all()
    assert (data[2]==0).all()
