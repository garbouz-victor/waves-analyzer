"""Synchronized physical wall-normal profiles and a geometrically honest wall zoom."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from .panels import display_cell_edges,contact_corner_mask,solid_walls
from .scope import CORNER_WARNING,FIXED_DOMAIN

PROFILE_COLORS=['#bb5032','#16799a','#775598']


class WallProfilePanel:
    def __init__(self,ax,h,meta,component='speed',compact=False):
        self.ax,self.h,self.component=ax,h,component
        m=meta['no_slip_visualization'];self.lines=[]
        for depth,color in zip(h['profile_depths'][:],PROFILE_COLORS):
            line,=ax.plot(h['profile_epsilon'][:],np.zeros(len(h['profile_epsilon'])),lw=1.6,color=color,label=f'z={depth:g} m')
            self.lines.append(line)
        limits=m['profile_speed_ylim_m_per_s'] if component=='speed' else m['profile_component_ylim_m_per_s'][component]
        ax.set(xlim=m['profile_epsilon_range_m'],ylim=np.array(limits)*1000,
               xlabel='distance from left wall, epsilon (m)',ylabel=('|v|' if component=='speed' else component)+' (mm/s)')
        ax.axhline(0,color='.5',lw=.6)
        ax.plot(0,0,'o',color='black',ms=5,zorder=10,clip_on=False)
        ax.set_title('No-slip: zero at wall, flow inside' if component=='speed' else component+' component, fixed physical scale',fontsize=11 if compact else 13)
        ax.grid(alpha=.2)
        ax.legend(ncol=3,fontsize=7.5 if compact else 10,loc='upper right')
        ax.tick_params(labelsize=9 if compact else 11)
        if compact:ax.xaxis.label.set_size(9);ax.yaxis.label.set_size(9)

    def update(self,i):
        values=self.h['profile_left_'+self.component][i]
        for line,v in zip(self.lines,values):line.set_ydata(v*1000)


class NoSlipFigure:
    def __init__(self,h,meta):
        self.h,self.meta=h,meta;self.fig=plt.figure(figsize=(19.2,10.8),dpi=100)
        self.time=self.fig.text(.05,.955,'',fontsize=19,weight='bold')
        self.fig.text(.05,.91,'No-slip means zero exactly on the wall — not zero throughout a finite strip.',fontsize=13)
        ax=self.fig.add_axes([.06,.16,.34,.68]);self.ax=ax
        # True aspect ratio: this is a narrow 0.2 m wall-normal cut, not a wider tank.
        ax.set(xlim=(-1,-.8),ylim=(-1.2,.04),xlabel='x (m)',ylabel='z (m)',title='Left wall zoom; equal x/z aspect')
        ax.set_aspect('equal',adjustable='box')
        self.field=ax.pcolormesh(display_cell_edges(h['wall_zoom_x'][:]),display_cell_edges(h['wall_zoom_z'][:]),
            np.zeros((len(h['wall_zoom_z']),len(h['wall_zoom_x']))),cmap='viridis',norm=Normalize(*meta['velocity_color_range_m_per_s']),shading='flat')
        cb=self.fig.colorbar(self.field,ax=ax,fraction=.035,pad=.07,extend='max');cb.set_label('speed (m/s)')
        self.q=ax.quiver(h['zoom_arrow_x'][:],h['zoom_arrow_z'][:],np.zeros_like(h['zoom_arrow_x'][:]),np.zeros_like(h['zoom_arrow_z'][:]),
            pivot='tail',angles='xy',scale_units='xy',scale=1/meta['no_slip_visualization']['zoom_fixed_arrow_seconds'],
            width=.012,headwidth=3.4,color='white',minlength=0,zorder=4)
        ax.quiverkey(self.q,.5,-.095,.0005,'0.5 mm/s; fixed zoom arrow scale',labelpos='S',coordinates='axes',fontproperties={'size':9})
        contact_corner_mask(ax,meta['scope'],warning=False);self.walls=solid_walls(ax,-1.2,sides=(-1,))
        ax.annotate('solid wall\nu = w = 0',xy=(-1,-.65),xytext=(-1.13,-.5),fontsize=11,
                    annotation_clip=False,bbox={'facecolor':'white','edgecolor':'.6'},arrowprops={'arrowstyle':'-','lw':1.5},zorder=31)
        for depth,color in zip(h['profile_depths'][:],PROFILE_COLORS):
            ax.plot([-1,-.85],[depth,depth],color=color,lw=.8,ls=':',zorder=6)
            ax.plot(-1,depth,'o',color='black',ms=4,zorder=32,clip_on=False)
        self.panels=[WallProfilePanel(self.fig.add_axes([.51,y,.43,.17]),h,meta,key) for y,key in ((.68,'speed'),(.415,'u'),(.15,'w'))]
        self.wall_text=self.fig.text(.51,.885,'',fontsize=11)
        self.fig.text(.05,.057,CORNER_WARNING+'; patch depth is a drawing convention.',fontsize=11)
        self.fig.text(.05,.023,FIXED_DOMAIN,fontsize=11)

    def update(self,i):
        self.time.set_text(f'No-slip and wall-normal velocity gradients  |  t={self.h["times"][i]:.3f} s  |  nu=0.01 m²/s')
        self.field.set_array(np.hypot(self.h['wall_zoom_u'][i],self.h['wall_zoom_w'][i]).ravel())
        self.q.set_UVC(self.h['zoom_arrow_u'][i],self.h['zoom_arrow_w'][i])
        for panel in self.panels:panel.update(i)
        self.wall_text.set_text('Actual FEM wall max |v|: '+', '.join(f'{s}={self.h["wall_"+s+"_speed_max"][i]:.2g} m/s' for s in ('left','right')))
        return self.fig
