#!/usr/bin/env python3
"""Final read-only source/video checks; write small provenance and test inventories."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import xml.etree.ElementTree as ET

import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sloshing.validation.qualification.dataset import dataset_fingerprint
from sloshing.visualization.render import signature
from sloshing.visualization.panels import MainFigure
import matplotlib.pyplot as plt


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'output/step2')
    out=parser.parse_args().output.resolve()
    meta=json.loads((out/'animation_metadata.json').read_text())
    history={'STEP1_6_REPORT.md':'85da88bb6726ba486c2048d9a9048d6bfc88e65b0a9ee672f20ff3c0fcd456ef',
             'validation_results/animation_qualification/summary.json':'997266bcecc6a2a21947119dba5c9d132571b423d27e832547365b91e4157858'}
    for path,expected in history.items():
        if sha256(ROOT/path)!=expected:raise RuntimeError('Historical file changed: '+path)
    sources={}
    for key,fp in meta['sources'].items():
        if dataset_fingerprint(fp['path'])!=fp:raise RuntimeError('Source revision changed: '+key)
        with h5py.File(fp['path'],'r') as h:
            p=h['mesh/coordinates'][:]
            mesh={'nx_intervals':len(np.unique(p[0]))-1,'nz_intervals':len(np.unique(p[1]))-1,
                  'triangles':h['mesh/triangles'].shape[1]}
        sources[key]={**fp,'content_sha256':sha256(fp['path']),'resolved_mesh':mesh}
    inspected=json.loads((out/'inspection_summary.json').read_text())
    for name in ('bulk_flow_explained.mp4','bulk_flow_clean.mp4','bulk_vorticity.mp4','free_surface_true_scale.mp4','tracer_model_comparison.mp4','no_slip_boundary_layer.mp4'):
        check=inspected['videos'][name]
        if not check['decode_passed'] or sha256(out/name)!=check['sha256']:raise RuntimeError('Video inspection stale: '+name)
    browser=json.loads((out/'explorer_inspection.json').read_text())
    if browser['status']!='passed' or browser['html_sha256']!=sha256(out/'explorer.html'):raise RuntimeError('Browser inspection stale')
    no_slip=json.loads((out/'no_slip_inspection.json').read_text())
    if no_slip['status']!='passed' or no_slip['fine_content_sha256']!=sources['fine']['content_sha256']:
        raise RuntimeError('No-slip inspection stale')
    if no_slip['full_frame_count']!=meta['frame_count']:
        raise RuntimeError('No-slip inspection must cover every source snapshot')
    if no_slip['signature']!=signature(meta):raise RuntimeError('No-slip renderer inspection stale')
    ci=json.loads((out/'ci_summary.json').read_text())
    if ci['status']!='passed':raise RuntimeError('CI was not accepted')
    subprocess.run(['git','diff','--exit-code',ci['tested_commit'],'--','src','tests','configs',
                    'scripts/render_step2.py',str(ROOT.parent/'.github/workflows/test.yml')],cwd=ROOT,check=True)
    tests=[];runs=[]
    for name in ('unit_final.xml','validation_final.xml'):
        root=ET.parse(out/'test_runs'/name).getroot()
        for suite in root.iter('testsuite'):
            runs.append({'xml':name,**suite.attrib})
        for case in root.iter('testcase'):
            if any(case.find(k) is not None for k in ('failure','error','skipped')):raise RuntimeError('Non-passing final test: '+case.attrib['name'])
            tests.append({'id':case.attrib['classname']+'::'+case.attrib['name'],'status':'passed'})
    (out/'test_summary.json').write_text(json.dumps({'runs':runs,'passed':len(tests),'tests':tests},indent=2)+'\n')
    with h5py.File(out/'display_data.h5','r') as h:
        figure=MainFigure(h,meta);figure.fig.canvas.draw();ax=figure.flow.ax
        # Complete glyph bounds were checked in all 1001 frames. Also allow
        # for the finite white/black wall stroke in the actual 1920px layout.
        stroke=max(line.get_linewidth() for line in ax.lines if line.get_zorder()>=29)
        half_stroke_m=stroke*figure.fig.dpi/72/2*2/ax.get_window_extent().width
        ink_margin=no_slip['actual_glyph_wall_margin_m']-half_stroke_m
        conservative_ink_margin=meta['wall_visual_margin']-half_stroke_m
        assert conservative_ink_margin>meta['scope']['arrow_margin_extra_m']
        plt.close(figure.fig)
        i=int(np.argmax(h['kinetic_energy'][:]));depth=h['depths'][:];q=h['depth_q'][i]
        ids=[int(np.argmin(abs(depth-target))) for target in (-1.,-3.,-5.,-10.)]
        sequence=[]
        seeds=h['seeds'][:]
        for t in (.1,1.):
            k=int(np.argmin(abs(h['times'][:]-t)))
            row={'time_s':float(h['times'][k]),'samples':[]}
            for point in ((0.,-.1),(-.75,-.1),(.75,-.1)):
                j=int(np.argmin(np.linalg.norm(seeds-point,axis=1)))
                row['samples'].append({'x_z_m':seeds[j].tolist(),'u_w_m_per_s':h['seed_v'][k,j].tolist()})
            sequence.append(row)
        interpretation={'energy_initial':float(h['total_energy'][0]),'energy_final':float(h['total_energy'][-1]),
                        'first_global_K_max_time_s':float(h['times'][i]),'K_max':float(h['kinetic_energy'][i]),
                        'depth_Q_relative_to_surface_at_K_max':[{'z_m':float(depth[j]),'ratio':float(q[j]/q[-1])} for j in ids],
                        'calculated_direction_samples':sequence}
    paths=sorted((ROOT/'src/sloshing/visualization').glob('*.py'))+sorted((ROOT/'src/sloshing/visualization').glob('*.html'))
    paths += [ROOT/'configs/animation_scope.json',ROOT/'scripts/render_step2.py']
    code={str(p.relative_to(ROOT)):sha256(p) for p in paths}
    output={'status':'passed','source_commit_sha':meta['source_commit_sha'],'sources':sources,
            'historical_sha256_unchanged':history,'render_code_sha256':code,'all_six_mp4_verified':True,
            'no_slip_verified':True,'no_slip_inspection_sha256':sha256(out/'no_slip_inspection.json'),
            'main_wall_outer_stroke_halfwidth_m':half_stroke_m,'actual_arrow_to_wall_ink_margin_m':ink_margin,
            'conservative_arrow_to_wall_ink_margin_m':conservative_ink_margin,
            'ci_code_matches_tested_snapshot':ci['tested_commit'],'ci_run_url':ci['url'],
            'browser_verified':True,'test_count':len(tests),'interpretation_samples':interpretation}
    (out/'artifact_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
