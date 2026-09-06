"""Static scientific wall diagnostics; exact saved snapshot times and fixed scales."""

import json
from pathlib import Path
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from .panels import display_cell_edges,contact_corner_mask,solid_walls
from .no_slip_panels import WallProfilePanel,NoSlipFigure


def no_slip_previews(cache,meta,output):
    dest=Path(output)/'frames_preview';dest.mkdir(exist_ok=True)
    with h5py.File(cache,'r') as h:
        for t in (.1,.395,.810,1.615):
            i=int(np.argmin(abs(h['times'][:]-t)))
            fig,ax=plt.subplots(figsize=(8,10),dpi=120)
            speed=np.hypot(h['wall_zoom_u'][i],h['wall_zoom_w'][i])
            field=ax.pcolormesh(display_cell_edges(h['wall_zoom_x'][:]),display_cell_edges(h['wall_zoom_z'][:]),speed,
                cmap='viridis',norm=Normalize(*meta['velocity_color_range_m_per_s']),shading='flat')
            fig.colorbar(field,ax=ax,fraction=.04,pad=.07,label='speed (m/s)',extend='max')
            ax.quiver(h['zoom_arrow_x'][:],h['zoom_arrow_z'][:],h['zoom_arrow_u'][i],h['zoom_arrow_w'][i],
                pivot='tail',scale_units='xy',angles='xy',scale=1/meta['no_slip_visualization']['zoom_fixed_arrow_seconds'],
                color='white',width=.009,headwidth=3.4,minlength=0,zorder=4)
            ax.scatter(h['zoom_arrow_x'][:],h['zoom_arrow_z'][:],s=9,facecolor='black',edgecolor='white',lw=.5,zorder=5,label='FEM arrow tails')
            contact_corner_mask(ax,meta['scope'],warning=False);solid_walls(ax,-.35,sides=(-1,))
            ax.set(xlim=(-1.015,-.79),ylim=(-.35,.04),xlabel='x (m)',ylabel='z (m)')
            ax.set_aspect('equal',adjustable='box')
            ax.annotate('solid wall\nu = w = 0',xy=(-1,-.28),xytext=(-.94,-.31),zorder=31,
                bbox={'facecolor':'white','edgecolor':'none'},arrowprops={'arrowstyle':'-','color':'black'})
            ax.text(-.965,.02,'unresolved\ncontact corner',fontsize=9,color='.2')
            ax.legend(loc='lower right',fontsize=8)
            fig.suptitle(f'Actual P2 wall velocity = {h["wall_left_speed_max"][i]:.2g} m/s\nt = {h["times"][i]:.3f} s; velocity inside the liquid is not zero',fontsize=13)
            fig.text(.08,.025,'Zoom arrows: fixed 20 s visual multiplier; no zero-speed strip added.',fontsize=10)
            fig.savefig(dest/f'no_slip_t{t:g}.png');plt.close(fig)
        f=NoSlipFigure(h,meta);f.update(79).savefig(dest/'no_slip_diagnostic.png',dpi=100);plt.close(f.fig)


def static_profiles(cache,meta,output):
    with h5py.File(cache,'r') as h:
        fig,axes=plt.subplots(2,3,figsize=(17,9),constrained_layout=True)
        rows=[]
        for ax,t in zip(axes.flat,(.1,.395,.810,1.220,1.615,4.040)):
            i=int(np.argmin(abs(h['times'][:]-t)))
            p=WallProfilePanel(ax,h,meta);p.update(i)
            ax.set_title(f't = {h["times"][i]:.3f} s; left wall FEM speed = {h["wall_left_speed_max"][i]:.2g}')
            for eps in (0,1e-5,1e-4,.001,.01,.05,.15):
                j=int(np.argmin(abs(h['profile_epsilon'][:]-eps)))
                rows.append({'time_s':float(h['times'][i]),'epsilon_m':float(h['profile_epsilon'][j]),
                    'depths_m':h['profile_depths'][:].tolist(),'left_speed_m_per_s':h['profile_left_speed'][i,:,j].tolist(),
                    'right_speed_m_per_s':h['profile_right_speed'][i,:,j].tolist()})
        fig.suptitle('No-slip: v = 0 exactly at the wall, not in a finite-width strip. Actual P2 FEM profiles; fixed axes.',fontsize=15)
        for ext in ('pdf','png'):fig.savefig(Path(output)/('no_slip_profiles.'+ext),dpi=140)
        plt.close(fig)
    (Path(output)/'no_slip_profile_samples.json').write_text(json.dumps(rows,indent=2)+'\n')
