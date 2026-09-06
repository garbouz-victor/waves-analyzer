#!/usr/bin/env python3
"""All-snapshot FEM and glyph checks, with public Matplotlib artist geometry."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sloshing.visualization.panels import FlowPanel
from sloshing.visualization.no_slip_panels import NoSlipFigure
from sloshing.visualization.wall_profiles import arrow_geometry,wall_normal_profiles,trace_speed_max
from sloshing.validation.qualification.fem_evaluation import FEMGeometry,evaluate_map

FINE_HASH='fa299d74c9487a40a5fac55966c1288c950bc0a46116f80e7da0063f794cf09a'


def glyph_bounds(q,ax):
    """Transform public Path/offset collections from display back to data coordinates."""
    offsets=q.get_offset_transform().transform(q.get_offsets())
    points=np.concatenate([q.get_transform().transform(path.vertices)+offset for path,offset in zip(q.get_paths(),offsets)])
    values=ax.transData.inverted().transform(points)
    return np.min(values,axis=0),np.max(values,axis=0)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=Path('output/step2'))
    args=p.parse_args();out=args.output;start=time.perf_counter()
    meta=json.loads((out/'animation_metadata.json').read_text());source=Path(meta['sources']['fine']['path'])
    digest=hashlib.sha256()
    with source.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):digest.update(block)
    assert digest.hexdigest()==FINE_HASH,'Production FEM dataset changed'
    with h5py.File(source,'r') as src,h5py.File(out/'display_data.h5','r') as h:
        g=FEMGeometry.from_hdf5(src);z=np.sort(np.r_[g.z,(g.z[:-1]+g.z[1:])/2])
        maps={s:g.interpolation(np.array([np.full_like(z,x),z])) for s,x in (('left',-1.),('right',1.))}
        wall={s:0. for s in maps}
        for i in range(len(h['times'])):
            u,w=src['fields/u'][i],src['fields/w'][i]
            for side,mp in maps.items():
                value=trace_speed_max(evaluate_map(u,mp),evaluate_map(w,mp))
                assert value==h['wall_'+side+'_speed_max'][i]
                wall[side]=max(wall[side],value)
            if i%200==0:print(f'FEM wall check {i+1}/1001; elapsed={time.perf_counter()-start:.1f}s',flush=True)
        assert max(wall.values())<=meta['parameters']['wall_tolerance']
        for i in (0,20,79,162,244,323,808,1000):
            for side in maps:
                prof=wall_normal_profiles(g,src['fields/u'][i],src['fields/w'][i],side,h['profile_depths'][:],h['profile_epsilon'][:])
                for key,val in prof.items():np.testing.assert_array_equal(val,h['profile_'+side+'_'+key][i])
        symmetry=float(abs(h['profile_left_speed'][:]-h['profile_right_speed'][:]).max())
        assert symmetry<=1e-12
        envelope=arrow_geometry(h['arrow_x'][:],h['arrow_z'][:],h['arrow_u'][:],h['arrow_w'][:],meta['fixed_arrow_seconds'],meta['scope'])
        # Public glyph vertices across ALL frames. Only quiver is drawn in this
        # small audit canvas; equal aspect and physical x span match FlowPanel.
        fig,ax=plt.subplots(figsize=(6,5),dpi=50)
        ax.set(xlim=(-1,1),ylim=(-1.5,.04));ax.set_aspect('equal',adjustable='box')
        q=ax.quiver(h['arrow_x'][:],h['arrow_z'][:],h['arrow_u'][0],h['arrow_w'][0],
            pivot='tail',angles='xy',scale_units='xy',scale=1/meta['fixed_arrow_seconds'],width=.0024,headwidth=3.4,minlength=0)
        low=np.array([np.inf,np.inf]);high=-low
        for i in range(len(h['times'])):
            q.set_UVC(h['arrow_u'][i],h['arrow_w'][i]);fig.canvas.draw()
            lo,hi=glyph_bounds(q,ax);low=np.minimum(low,lo);high=np.maximum(high,hi)
            if i%200==0:print(f'Public glyph check {i+1}/1001; elapsed={time.perf_counter()-start:.1f}s',flush=True)
        plt.close(fig)
        assert low[0]>-.98+.02 and high[0]<.98-.02 and high[1]<0
        fig,ax=plt.subplots();panel=FlowPanel(ax,h,meta,tracers=False)
        assert panel.q.pivot=='tail'
        assert all(line.get_zorder()>panel.q.get_zorder() for line in panel.walls)
        assert all(patch.get_y()==-.12 and np.isclose(patch.get_height(),.16) for patch in panel.masks)
        plt.close(fig)
        zoom=NoSlipFigure(h,meta)
        for i in (0,20,79,162,244,323,808,1000):
            zoom.update(i);zoom.fig.canvas.draw();lo,hi=glyph_bounds(zoom.q,zoom.ax)
            assert lo[0]>-.98 and hi[0]<-.8 and hi[1]<-.12
        plt.close(zoom.fig)
    result={'status':'passed','fine_content_sha256':digest.hexdigest(),'wall_speed_max_m_per_s':wall,
            'full_frame_count':1001,'all_frame_glyph_data_bounds':[low.tolist(),high.tolist()],
            'actual_glyph_wall_margin_m':float(min(1+low[0],1-high[0])),
            'conservative_geometry':envelope,'speed_profile_symmetry_max_m_per_s':symmetry,
            'profile_source_exact_match':True,'public_artist_geometry':True,'runtime_s':time.perf_counter()-start}
    (out/'no_slip_inspection.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
