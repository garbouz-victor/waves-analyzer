import numpy as np
import matplotlib.pyplot as plt
import pytest
from sloshing.visualization.panels import FlowPanel,SurfacePanel
from sloshing.visualization.scope import load_scope,CORNER_WARNING


@pytest.mark.parametrize("field",["speed","omega"])
def test_two_actual_frames_keep_clim_and_arrow_scale(field):
    x=np.linspace(-1,1,5);z=np.linspace(-1.5,0,4)
    xx,zz=np.meshgrid(x,z-.15)
    h={"x":x,"z":z,"eta_x":x,"eta":np.array([x*.0003,x*.00003]),
       "arrow_x":xx,"arrow_z":zz,"seeds":np.array([[0.,-.3]]),
       "u":np.array([np.ones_like(xx)*.001,np.ones_like(xx)*.00001]),
       "w":np.zeros((2,4,5)),"omega":np.array([np.ones_like(xx)*.02,np.ones_like(xx)*-.0002])}
    h["arrow_u"],h["arrow_w"]=h["u"],h["w"]
    m={"scope":load_scope(),"velocity_color_range_m_per_s":[0,.002],"vorticity_color_range_per_s":[-.03,.03],
       "fixed_arrow_seconds":100,"arrow_key_speed_m_per_s":.0005,"surface_ylim_mm":[-.4,.4]}
    fig,(ax,sax)=plt.subplots(1,2)
    p=FlowPanel(ax,h,m,field=field,tracers=False);s=SurfacePanel(sax,h,m)
    vertices=p.image.get_coordinates()
    assert vertices[:,:,1].max()==0 and vertices[:,:,1].min()==-1.5
    assert vertices[:,:,0].min()==-1 and vertices[:,:,0].max()==1
    p.update(0);s.update(0)
    clim=p.image.get_clim();scale=p.q.scale;limits=sax.get_ylim()
    p.update(1);s.update(1)
    assert p.image.get_clim()==clim and p.q.scale==scale and sax.get_ylim()==limits
    assert any(t.get_text()==CORNER_WARNING for t in ax.texts)
    np.testing.assert_allclose(s.line.get_ydata(),h["eta"][1]*1000)
    plt.close(fig)
