#!/usr/bin/env python3
"""Annotate actual encoded frame using deterministic Matplotlib glyph coordinates."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sloshing.visualization.panels import MainFigure


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['before','after'],required=True)
    p.add_argument('--output',type=Path,default=Path('output/step2'))
    a=p.parse_args();out=a.output
    source=out/'archive_step2_4172787' if a.stage=='before' and (out/'archive_step2_4172787').exists() else out
    meta=json.loads((source/'animation_metadata.json').read_text())
    if a.stage=='before' and meta['schema']!='step2-display-v1':
        raise RuntimeError('Before requires the archived STEP 2 v1 cache and actual old MP4; never relabel new output as before')
    with h5py.File(source/'display_data.h5','r') as h:
        times=h['times'][:];i=int(np.argmin(abs(times-.395)))
        x,z=h['arrow_x'][:],h['arrow_z'][:]
        u,w=h['arrow_u'][:],h['arrow_w'][:];s=meta['fixed_arrow_seconds']
        moving=np.any(np.hypot(u,w)>0,axis=0)
        near=np.min(x[moving]);col=int(np.argmin(abs(x[0]-near)))
        row=int(np.argmax(np.hypot(u[i,:,col],w[i,:,col])))
        tail=np.array([x[row,col],z[row,col]])
        tip=tail+s*np.array([u[i,row,col],w[i,row,col]])
        metric={'stage':a.stage,'time_s':float(times[i]),'index':i,'seed_x_range':[float(x.min()),float(x.max())],
                'nearest_nonzero_left_seed_x':float(near),'nearest_seed_distance_m':float(near+1),
                'fixed_arrow_seconds':s,'maximum_rendered_arrow_dx_m':float(abs(s*u).max()),
                'maximum_rendered_arrow_length_m':float((s*np.hypot(u,w)).max()),
                'all_tip_x_range':[float((x+s*u).min()),float((x+s*u).max())],
                'annotated_tail_xz_m':tail.tolist(),'annotated_tip_xz_m':tip.tolist()}
        if a.stage=='before':
            # Coordinate transform of the historical 4172787 layout, including
            # colorbar space. The image below is the OLD ENCODED MP4, not a redraw.
            fig=plt.figure(figsize=(19.2,10.8),dpi=100);ax=fig.add_axes([.055,.20,.555,.66])
            fig.colorbar(ScalarMappable(),ax=ax,pad=.025,fraction=.034)
            ax.set(xlim=(-1,1),ylim=(-1.5,.04));ax.set_aspect('equal',adjustable='box')
        else:
            model=MainFigure(h,meta);model.update(i);fig,ax=model.fig,model.flow.ax
        fig.canvas.draw()
        locations=ax.transData.transform(np.array([[-1,tail[1]],tail,tip]))
        plt.close(fig)
    dest=out/'geometry_audit';dest.mkdir(exist_ok=True)
    raw=dest/(a.stage+'_decoded.png')
    subprocess.run(['ffmpeg','-v','error','-i',str(source/'bulk_flow_explained.mp4'),'-vf',f'select=eq(n\\,{i})',
                    '-frames:v','1','-y',str(raw)],check=True)
    pixels=plt.imread(raw);height=pixels.shape[0]
    fig,ax=plt.subplots(figsize=(12,8),dpi=130)
    ax.imshow(pixels);xs=locations[:,0];ys=height-locations[:,1]
    ax.set(xlim=(xs.min()-75,xs.max()+350),ylim=(ys.max()+200,ys.min()-120));ax.axis('off')
    for label,point,offset,color in zip(['solid wall x = -1',f'arrow tail = FEM seed\nx = {tail[0]:.7f} m',
                                       f'arrow tip\nx = {tip[0]:.7f} m'],zip(xs,ys),[(100,-100),(140,-40),(160,30)],['#ffcc00','#00ffff','#ff7070']):
        ax.annotate(label,xy=point,xytext=offset,textcoords='offset points',fontsize=12,color='black',
                    bbox={'facecolor':'white','alpha':.95,'pad':5},arrowprops={'arrowstyle':'->','color':color,'lw':2})
        ax.plot(*point,'o',color=color,ms=6)
    fig.suptitle(f'{a.stage.upper()}: actual encoded frame t={metric["time_s"]:.3f} s; coordinates from renderer',fontsize=14)
    fig.savefig(out/(a.stage+'_arrow_geometry.png'),bbox_inches='tight');plt.close(fig)
    (out/(a.stage+'_arrow_geometry.json')).write_text(json.dumps(metric,indent=2)+'\n')
    print(json.dumps(metric,indent=2))


if __name__=='__main__':main()
