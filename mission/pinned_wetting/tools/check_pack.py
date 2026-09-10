#!/usr/bin/env python3
"""Check this handoff package, never mark a scientific task complete."""
from __future__ import annotations
import ast
import json
from pathlib import Path
import sys

REQUIRED = (
 'AGENTS.md','AGENTS.fragment.md','MISSION.md','PHYSICS.md','CONFIG.yaml',
 'ACCEPTANCE.yaml','RENDERING.md','EXECUTION.md','REPO_MAP.md','EXECPLAN.md',
 'state.json','OUTPUT_SCHEMA.md','SOURCES.md','TOOLS.md',
 'tools/check_wetting_csv.py','tests/test_wetting_checker.py'
)

def main() -> int:
    root=Path(__file__).resolve().parents[1]
    errors=[]
    for rel in REQUIRED:
        path=root/rel
        if not path.is_file() or not path.stat().st_size:
            errors.append(f'Missing/empty: {rel}')
    for path in root.rglob('*.py'):
        try: ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
        except (SyntaxError,UnicodeError) as exc: errors.append(str(exc))
    try:
        state=json.loads((root/'state.json').read_text(encoding='utf-8'))
        if state.get('mission_id') != 'pinned-wetting-PW1': errors.append('Wrong mission ID')
    except (OSError,ValueError) as exc: errors.append(str(exc))
    yaml_check='SKIPPED: PyYAML not installed; runtime tools do not require it'
    try:
        import yaml
    except ImportError:
        pass
    else:
        try:
            config=yaml.safe_load((root/'CONFIG.yaml').read_text(encoding='utf-8'))
            acceptance=yaml.safe_load((root/'ACCEPTANCE.yaml').read_text(encoding='utf-8'))
            if config['mission_id'] != acceptance['mission_id']: errors.append('Mission IDs differ')
            if config['simulation']['t_end_s'] != acceptance['simulation']['t_end_s']:
                errors.append('Horizon mismatch')
            if not config['physics']['irreversible_wetting']: errors.append('Irreversibility disabled')
            if not config['physics']['connected_retained_layer']: errors.append('Connected layer disabled')
            yaml_check='PASS'
        except Exception as exc: errors.append(f'YAML: {exc}')
    print(json.dumps({'status':'PACKAGE_CHECK_FAILED' if errors else 'PACKAGE_OK',
                      'scientific_status':'NOT_RUN','yaml_check':yaml_check,'errors':errors},
                     ensure_ascii=False,indent=2))
    return bool(errors)
if __name__=='__main__': sys.exit(main())
