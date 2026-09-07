"""Capture the exact uncommitted module source of an active/finished run, only on hash match."""
import argparse
import json
from pathlib import Path
import sys
import tarfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sloshing.multiphase.provenance import source_hash


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run",type=Path)
    args=parser.parse_args()
    provenance=json.loads((args.run/"provenance.json").read_text())
    if source_hash()!=provenance["multiphase_source_sha256"]:
        raise ValueError("Current source does not match this run; refusing a misleading source archive")
    source=Path(__file__).resolve().parents[1]/"src/sloshing/multiphase"
    target=args.run/"multiphase_source.tar.gz"
    with target.open("xb") as stream:
        with tarfile.open(fileobj=stream,mode="w:gz") as archive:
            for path in sorted(source.rglob("*.py")):
                archive.add(path,arcname=str(Path("multiphase")/path.relative_to(source)))
    print(str(target),source_hash())


if __name__=="__main__":
    main()
