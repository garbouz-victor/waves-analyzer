from copy import deepcopy
import pytest
from sloshing.visualization.scope import compatible_cache
from sloshing.visualization.render import require_gate,write_json,signature


def test_changed_dataset_invalidates_display_metadata():
    old={"sources":{"fine":{"bytes":100,"mtime_ns":5,"config_sha256":"abc"}},"schema":"step2"}
    compatible_cache(old,deepcopy(old))
    for key,value in (("bytes",101),("mtime_ns",6),("config_sha256","def")):
        changed=deepcopy(old);changed["sources"]["fine"][key]=value
        with pytest.raises(RuntimeError,match="fingerprint"):compatible_cache(old,changed)


def test_full_render_requires_inspected_preview(tmp_path):
    meta={"source":"test"}
    with pytest.raises(RuntimeError):require_gate(tmp_path,meta,"static_preview")
    path=tmp_path/"static_preview_gate.json"
    write_json(path,{"status":"awaiting_inspection","signature":signature(meta)})
    with pytest.raises(RuntimeError):require_gate(tmp_path,meta,"static_preview")
    write_json(path,{"status":"accepted_after_inspection","signature":signature(meta)})
    require_gate(tmp_path,meta,"static_preview")
