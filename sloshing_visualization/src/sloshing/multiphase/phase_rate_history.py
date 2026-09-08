"""Small serial scientific journal and atomic accepted-state archives.

No rejected state is a restart source. Torn/orphan suffixes are preserved,
not silently interpreted as an accepted complete trajectory.
"""
import hashlib
import copy
import json
import os
from pathlib import Path
import shutil
import tempfile
import fcntl
import numpy as np
from .linearized_ch import array_hash, json_hash

VERSION = "accepted-phase-rate-trajectory-v1"


def read_json(path):
    return json.loads(Path(path).read_text())


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n").encode()
    descriptor, name = tempfile.mkstemp(prefix=".pending-", dir=path.parent)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(data); stream.flush(); os.fsync(stream.fileno())
    os.replace(name, path)


def archive(path, arrays, metadata):
    """Publish a new directory atomically; refuse all archive overwrites."""
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    arrays = {k: np.array(v, copy=True) for k,v in arrays.items()}
    if metadata.get("accepted") is not False and not all(np.isfinite(v).all() for v in arrays.values()):
        raise ValueError("Nonfinite archive coefficients")
    record = {**metadata, "archive_version": VERSION,
        "array_sha256": {k: array_hash(v) for k,v in arrays.items()}}
    record["semantic_sha256"] = json_hash(record)
    pending = Path(tempfile.mkdtemp(prefix=".pending-", dir=path.parent))
    with (pending/"arrays.npz").open("wb") as stream:
        np.savez_compressed(stream, **arrays); stream.flush(); os.fsync(stream.fileno())
    atomic_json(pending/"metadata.json", record)
    os.rename(pending, path)
    return record


def load_archive(path, identity=None):
    path = Path(path); meta = read_json(path/"metadata.json")
    if meta["semantic_sha256"] != json_hash({k:v for k,v in meta.items() if k!="semantic_sha256"}):
        raise ValueError("Archive metadata SHA mismatch")
    if identity is not None and meta["identity"] != identity:
        raise ValueError("Restart source/schedule/physics/stack identity mismatch")
    with np.load(path/"arrays.npz", allow_pickle=False) as z:
        arrays = {k:z[k].copy() for k in z.files}
    if {k:array_hash(v) for k,v in arrays.items()} != meta["array_sha256"]:
        raise ValueError("Archive array SHA mismatch")
    return arrays, meta


class Journal:
    def __init__(self, path, *, load=True):
        self.path=Path(path); self.path.parent.mkdir(parents=True,exist_ok=True)
        self.rows=self.read(self.path) if load and self.path.exists() else []

    @staticmethod
    def read(path, data=None):
        rows=[]; previous=None
        raw=Path(path).read_bytes() if data is None else data
        if raw and not raw.endswith(b"\n"):
            raise ValueError("Torn journal record")
        for line in raw.splitlines():
            row=json.loads(line)
            if row["previous_record_sha256"]!=previous or row["record_sha256"]!=json_hash(
                    {k:v for k,v in row.items() if k!="record_sha256"}):
                raise ValueError("Journal hash chain corrupted")
            rows.append(row); previous=row["record_sha256"]
        return rows

    def append(self, value):
        row={**copy.deepcopy(value),"previous_record_sha256":self.rows[-1]["record_sha256"] if self.rows else None}
        row["record_sha256"]=json_hash(row)
        with self.path.open("ab") as stream:
            stream.write((json.dumps(row,sort_keys=True,allow_nan=False)+"\n").encode())
            stream.flush(); os.fsync(stream.fileno())
        self.rows.append(row)
        return row

    def prefix(self):
        raw=self.path.read_bytes()
        return {"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),"records":len(self.rows)}


def recover_prefix(folder, checkpoint):
    """Preserve any crash suffix before returning a verified restart prefix."""
    folder=Path(folder)
    if (folder/"FAILED.json").exists():
        raise ValueError("Scientific failure forbids retry/resume")
    path=folder/"scalar_history.jsonl"; raw=path.read_bytes(); p=checkpoint["journal_prefix"]
    prefix=raw[:p["bytes"]]
    if hashlib.sha256(prefix).hexdigest()!=p["sha256"]:
        raise ValueError("Checkpoint journal prefix SHA mismatch")
    rows=Journal.read(path,prefix)
    if len(rows)!=p["records"] or rows[-1]["step"]!=checkpoint["accepted_steps"]:
        raise ValueError("Checkpoint journal count mismatch")
    if raw!=prefix:
        orphan=Path(tempfile.mkdtemp(prefix="orphan-",dir=folder))
        shutil.move(path,orphan/path.name)
        with path.open("xb") as stream:
            stream.write(prefix); stream.flush(); os.fsync(stream.fileno())
        for group in ("fields","checkpoints","residual_audits"):
            directory=folder/group
            if not directory.exists(): continue
            for child in directory.iterdir():
                if child.name.isdigit() and int(child.name)>checkpoint["accepted_steps"]:
                    (orphan/group).mkdir(exist_ok=True)
                    shutil.move(child,orphan/group/child.name)
        timing=folder/"timing_history.jsonl"
        if timing.exists():
            old=Journal.read(timing); shutil.move(timing,orphan/timing.name)
            new=Journal(timing)
            for r in old:
                if r["step"]<=checkpoint["accepted_steps"]:
                    new.append({k:v for k,v in r.items() if k not in ("record_sha256","previous_record_sha256")})
    return Journal(path)


class WriterLock:
    def __init__(self, folder):
        Path(folder).mkdir(parents=True,exist_ok=True)
        self.stream=(Path(folder)/".writer.lock").open("a")
        try: fcntl.flock(self.stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            self.stream.close(); raise RuntimeError("Concurrent trajectory writer forbidden")

    def close(self):
        self.stream.close()


def validate_complete(folder, identity, schedule, marker=None):
    folder=Path(folder); marker=read_json(folder/"COMPLETE.json") if marker is None else marker
    if marker["identity"]!=identity or marker["accepted_steps"]!=schedule.nsteps:
        raise ValueError("Invalid completion identity/count")
    rows=Journal.read(folder/"scalar_history.jsonl")
    if (len(rows)!=schedule.nsteps+1 or [r["step"] for r in rows]!=list(range(schedule.nsteps+1)) or
            marker["scalar_history_sha256"]!=file_hash(folder/"scalar_history.jsonl") or
            rows[-1]["time"]!=schedule.intervals[-1].end or marker["final_time"]!=rows[-1]["time"]):
        raise ValueError("Incomplete/corrupt completed history")
    from .nonlinear_accuracy import NonlinearControls
    for row,interval in zip(rows[1:],schedule.intervals):
        if (row["time"]!=interval.end or row["dt"]!=interval.dt or row["block"]!=interval.block or
                row["solver"]["reason"]<=0 or row["solver"]["dt_reductions"]!=0 or
                row["solver"]["actual_controls"]!=NonlinearControls().record() or not all(row["physical_gates"].values())):
            raise ValueError("History does not implement the frozen accepted BE schedule/controls/gates")
    for i in schedule.field_indices:
        _,meta=load_archive(folder/"fields"/f"{i:06d}",identity)
        if meta["step"]!=i or meta["time"]!=(schedule.intervals[i-1].end if i else 0.):
            raise ValueError("Wrong common field checkpoint time")
    arrays,meta=load_archive(folder/"checkpoints"/f"{schedule.nsteps:06d}",identity)
    if (not meta["accepted"] or meta["accepted_steps"]!=schedule.nsteps or
            meta["time"]!=marker["final_time"] or meta["array_sha256"]["state"]!=marker["final_state_sha256"]):
        raise ValueError("Completion not backed by accepted final checkpoint")
    return rows
