#!/usr/bin/env python3
"""Check DECLARED model fields. This is not an independent PDE verifier."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys

DEFAULT_CONTRACT = Path(__file__).resolve().parents[1] / 'CONTRACT.json'

def read_json(path: Path) -> dict:
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f'duplicate JSON key: {key}')
            result[key] = value
        return result
    def bad(value):
        raise ValueError(f'nonfinite JSON constant: {value}')
    value = json.loads(Path(path).read_text(encoding='utf-8'),
                       object_pairs_hook=pairs, parse_constant=bad)
    if not isinstance(value, dict):
        raise ValueError('JSON root must be an object')
    return value

def _compare(expected, actual, path, errors):
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            errors.append(f'{path}: expected object'); return
        for key in sorted(set(expected) - set(actual)):
            errors.append(f'{path}.{key}: missing')
        for key in sorted(set(actual) - set(expected)):
            errors.append(f'{path}.{key}: undeclared physical key')
        for key in sorted(set(actual) & set(expected)):
            _compare(expected[key], actual[key], f'{path}.{key}', errors)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            errors.append(f'{path}: list length/type mismatch'); return
        for i, (e, a) in enumerate(zip(expected, actual)):
            _compare(e, a, f'{path}[{i}]', errors)
    elif type(expected) is bool:
        if type(actual) is not bool or actual != expected:
            errors.append(f'{path}: expected boolean {expected!r}, got {actual!r}')
    elif isinstance(expected, (int, float)):
        if (type(actual) not in (int, float) or not math.isfinite(actual)
                or actual != expected):
            errors.append(f'{path}: expected {expected!r}, got {actual!r}')
    elif type(actual) is not type(expected) or actual != expected:
        errors.append(f'{path}: expected {expected!r}, got {actual!r}')

def contract_errors(contract: dict, resolved: dict) -> list[str]:
    errors = []
    if 'physical_contract' not in contract:
        return ['baseline: missing physical_contract']
    _compare(contract['physical_contract'], resolved.get('physical_contract'),
             'physical_contract', errors)
    return errors

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('resolved', nargs='?', type=Path)
    parser.add_argument('--contract', type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument('--example', action='store_true', help='print header template; no files written')
    args = parser.parse_args(argv)
    try:
        expected = read_json(args.contract)
        if args.example:
            print(json.dumps({'physical_contract': expected['physical_contract'],
                 'numerics': {'method': 'TO_BE_IMPLEMENTED', 'source_hash': 'TO_BE_FILLED'},
                 'provenance': {'actual_execution_HEAD': 'TO_BE_RECORDED'}},
                 ensure_ascii=False, indent=2))
            return 0
        if args.resolved is None:
            parser.error('resolved_case.json required unless --example')
        errors = contract_errors(expected, read_json(args.resolved))
        report = {'scope': 'DECLARED_CONTRACT_ONLY_NOT_PDE', 'passed': not errors,
                  'cfd_verified': False, 'errors': errors}
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 2 if errors else 0
    except (ValueError, OSError, KeyError) as exc:
        print(json.dumps({'scope': 'DECLARED_CONTRACT_ONLY_NOT_PDE', 'passed': False,
                          'cfd_verified': False, 'errors': [str(exc)]}, ensure_ascii=False))
        return 2

if __name__ == '__main__':
    sys.exit(main())
