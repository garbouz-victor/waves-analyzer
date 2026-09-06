#!/usr/bin/env python3
"""Explicit, recoverable archive before changing the display schema/renderer.

Never touches FEM runs. Refuses an existing archive, and never deletes files.
"""
import argparse
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,default=Path('output/step2'))
p.add_argument('--archive-name',default='archive_step2_4172787')
a=p.parse_args();out=a.output.resolve()
if Path(a.archive_name).name!=a.archive_name:raise ValueError('Archive name must be a single directory name')
dest=out/a.archive_name
targets=[v for v in out.iterdir() if v.is_file() and not v.name.startswith('before_arrow_geometry') and
         (v.suffix in ('.mp4','.json') or v.name in ('display_data.h5','key_phases.pdf','key_phases.png'))]
if (out/'frames_preview').exists():targets.append(out/'frames_preview')
dest.mkdir(exist_ok=False)
for v in targets:v.rename(dest/v.name)
print(f'Archived {len(targets)} artifacts to {dest}; recoverable, no FEM data touched.')
