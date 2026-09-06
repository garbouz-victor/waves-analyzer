import numpy as np
import pytest
import matplotlib.pyplot as plt

from sloshing.config import SimulationConfig
from sloshing.fem_spaces import FEMSystem
from sloshing.validation.qualification.fem_evaluation import FEMGeometry,evaluate_map
from sloshing.visualization.wall_profiles import profile_points,wall_normal_profiles,trace_speed_max,arrow_geometry
from sloshing.visualization.panels import FlowPanel,contact_corner_mask,contact_mask,solid_walls
from sloshing.visualization.no_slip_panels import WallProfilePanel
from sloshing.visualization.scope import load_scope


@pytest.fixture
def system():
    f=FEMSystem(SimulationConfig(nx=6,nz=12))
    return f,FEMGeometry(f.mesh.p,f.mesh.t,f.scalar.element_dofs)


@pytest.mark.parametrize('side,edge',[('left',-1.),('right',1.)])
def test_wall_profiles_match_p2_polynomial_and_independent_skfem_probes(system,side,edge):
    f,g=system;x,z=f.scalar.doflocs
    u=(x-edge)*(z+2);w=(x-edge)**2
    eps=np.array([0.,1e-5,1e-4,.001,.01,.05,.15]);depth=np.array([-.05,-.2,-1.])
    points,shape=profile_points(g,side,depth,eps)
    p=wall_normal_profiles(g,u,w,side,depth,eps)
    np.testing.assert_allclose(p['u'],((points[0]-edge)*(points[1]+2)).reshape(shape),atol=5e-15,rtol=2e-12)
    np.testing.assert_allclose(p['w'],((points[0]-edge)**2).reshape(shape),atol=5e-15,rtol=2e-12)
    for key,c in (('u',u),('w',w)):
        np.testing.assert_allclose(p[key],(f.scalar.probes(points)@c).reshape(shape),atol=5e-15,rtol=2e-12)
    assert abs(p['speed'][:,0]).max()<1e-14


def test_exact_wall_trace_and_bottom_for_constrained_fem(system):
    f,g=system;rng=np.random.default_rng(29)
    full=f.expand_velocity(rng.normal(size=len(f.free)))
    wallz=np.sort(np.r_[g.z,(g.z[:-1]+g.z[1:])/2])
    for x in (-1,1):
        mp=g.interpolation(np.array([np.full_like(wallz,x),wallz]))
        uv=[evaluate_map(full[ids],mp) for ids in f.component_dofs]
        assert trace_speed_max(*uv)<=f.config.wall_tolerance
    points=np.array([np.linspace(-1,1,61),np.full(61,-10.)])
    for ids in f.component_dofs:assert abs(g.evaluate(full[ids],points)).max()<=f.config.wall_tolerance


def test_trace_norm_checks_between_nodes_not_only_nodal_max():
    s=np.array([0.,.5,1.]);u=1-(s-.25)**2
    assert np.isclose(trace_speed_max(u,2*u),np.sqrt(5),atol=1e-14)


def test_production_reflection_of_wall_normal_profiles(short_run):
    solver,states=short_run;f=solver.fem
    g=FEMGeometry(f.mesh.p,f.mesh.t,f.scalar.element_dofs)
    for state,_ in states:
        full=f.expand_velocity(state.velocity);u,w=(full[ids] for ids in f.component_dofs)
        left=wall_normal_profiles(g,u,w,'left',[-.05,-.2,-1.],[0,1e-4,.001,.01,.15])
        right=wall_normal_profiles(g,u,w,'right',[-.05,-.2,-1.],[0,1e-4,.001,.01,.15])
        np.testing.assert_allclose(left['speed'],right['speed'],atol=5e-13,rtol=1e-9)
        np.testing.assert_allclose(left['u'],right['u'],atol=5e-13,rtol=1e-9)
        np.testing.assert_allclose(left['w'],-right['w'],atol=5e-13,rtol=1e-9)


def test_arrow_envelope_has_fixed_margin_and_does_not_shorten_vectors():
    scope=load_scope();x,z=np.meshgrid(np.linspace(-.85,.85,25),[-1.,-.15])
    u=np.ones((3,2,25))*.0008;u[1]*=-1;w=u*.2;original=u.copy()
    m=arrow_geometry(x,z,u,w,102.,scope)
    assert 1-abs(x).max()>=.15-1e-14
    assert m['wall_visual_margin']>.02 and m['contact_strip_visual_margin_m']>.02
    np.testing.assert_array_equal(u,original)
    with pytest.raises(RuntimeError):arrow_geometry(x,z,u*3,w,102.,scope)


def test_flow_walls_above_quiver_tail_pivot_and_corner_only_mask():
    scope=load_scope();x=np.linspace(-1,1,9);z=np.linspace(-1.5,0,7)
    qx,qz=np.meshgrid(np.linspace(-.85,.85,5),[-1.,-.15])
    h={'x':x,'z':z,'eta_x':x,'eta':np.array([x*.0003]),'arrow_x':qx,'arrow_z':qz,
       'seeds':np.array([[0.,-.3]]),'u':np.zeros((1,7,9)),'w':np.zeros((1,7,9)),
       'arrow_u':np.ones((1,2,5))*.0008,'arrow_w':np.zeros((1,2,5))}
    m={'scope':scope,'velocity_color_range_m_per_s':[0,.002],'fixed_arrow_seconds':102.,'arrow_key_speed_m_per_s':.0005}
    fig,ax=plt.subplots();p=FlowPanel(ax,h,m,tracers=False);p.update(0);fig.canvas.draw()
    assert p.q.pivot=='tail'
    assert 'no-slip' in p.wall_label.get_text()
    assert len(p.walls)==2 and min(w.get_zorder() for w in p.walls)>p.q.get_zorder()
    for w in p.walls:assert w.get_linewidth()>=2
    for patch in p.masks:
        assert patch.get_y()==-.12 and np.isclose(patch.get_height(),.16)
        assert np.isclose(patch.get_width(),.02)
    plt.close(fig)


def test_public_patch_geometry_and_full_depth_bottom_wall():
    fig,ax=plt.subplots()
    surface=contact_mask(ax,load_scope(),warning=False)
    assert all(a.get_y()==0 and a.get_height()==1 for a in surface)
    corners=contact_corner_mask(ax,load_scope(),warning=False)
    assert all(a.get_y()>-1.5 for a in corners)
    walls=solid_walls(ax,-10,bottom=True)
    np.testing.assert_array_equal(walls[-1].get_ydata(),[-10,-10])
    plt.close(fig)


@pytest.mark.parametrize('component',['speed','u','w'])
def test_profile_axes_fixed_and_zero_visible(component):
    eps=np.linspace(0,.15,11);v=np.array([np.tile(eps,(3,1))*.001,np.tile(eps,(3,1))*.00001])
    h={'profile_epsilon':eps,'profile_depths':np.array([-.05,-.2,-1.]),'profile_left_'+component:v}
    m={'no_slip_visualization':{'profile_epsilon_range_m':[0,.15],'profile_speed_ylim_m_per_s':[0,.001],
                              'profile_component_ylim_m_per_s':{'u':[-.001,.001],'w':[-.001,.001]}}}
    fig,ax=plt.subplots();p=WallProfilePanel(ax,h,m,component)
    p.update(0);limits=ax.get_ylim();p.update(1)
    assert ax.get_ylim()==limits
    assert all(line.get_ydata()[0]==0 for line in p.lines)
    plt.close(fig)
