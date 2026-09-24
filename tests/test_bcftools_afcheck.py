"""Offline fixture tests for the bcftools recheck comparison (tool 42).

No network: fixtures are written to tmp_path and the module's file paths are
monkeypatched. Exercises _hits schema handling, compare_records, and
run_compare's vintage split / legacy-filter explanation logic.
"""
import json
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "bcftools_afcheck", Path(__file__).parent.parent / "scripts" / "bcftools_afcheck.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def _hit(pos, gt="1/1", targets=("dsx-v3-1_proto",), gq=90, dp=20):
    return {"pos": pos, "ref": "C", "alt": "A,T,G", "gt": gt,
            "gq": gq, "dp": dp, "targets": list(targets)}


def test_hits_handles_both_schemas():
    assert m._hits({"target_hits": [1]}) == [1]
    assert m._hits({"hits": [2]}) == [2]
    assert m._hits({}) == []


def test_compare_records_concordant_and_phased():
    a = {"target_hits": [_hit(48711510, gt="1|1")]}
    b = {"target_hits": [_hit(48711510, gt="1/1")], "n_rows": 10}
    c = m.compare_records(a, b)
    assert c["concordant"] and c["n_pysam"] == 1 and c["n_bcftools"] == 1


def test_compare_records_discordant():
    a = {"target_hits": [_hit(48711510)]}
    b = {"target_hits": [_hit(48711510), _hit(48711500, targets=("dsx-v3-1_PAM",))]}
    c = m.compare_records(a, b)
    assert not c["concordant"] and len(c["only_bcftools"]) == 1


def test_run_compare_vintage_split(tmp_path, monkeypatch):
    scan_in = tmp_path / "on_target.jsonl"
    recheck = tmp_path / "recheck.jsonl"
    out_json = tmp_path / "out.json"
    out_md = tmp_path / "out.md"
    monkeypatch.setattr(m, "SCAN_IN", str(scan_in))
    monkeypatch.setattr(m, "RECHECK_OUT", str(recheck))
    monkeypatch.setattr(m, "RESULT_JSON", str(out_json))
    monkeypatch.setattr(m, "RESULT_MD", str(out_md))

    hit_c = _hit(48711510)
    # new-schema row: has n_rows -> vintage 'new'
    scan_rows = [
        {"sample": "S1", "n_rows": 3000, "n_nonref": 5, "target_hits": [hit_c]},
        # legacy row: "hits" key, no n_rows, quality-filtered (missing a GQ 18 call)
        {"sample": "S2", "n_nonref": 4, "hits": [hit_c]},
    ]
    scan_in.write_text("\n".join(json.dumps(r) for r in scan_rows) + "\n")
    lowq = _hit(48711500, targets=("dsx-v3-1_PAM",), gq=18, dp=6)
    bcf_rows = [
        {"sample": "S1", "n_rows": 3001, "n_nonref": 5, "target_hits": [hit_c]},
        {"sample": "S2", "n_rows": 3001, "n_nonref": 6, "target_hits": [hit_c, lowq]},
    ]
    recheck.write_text("\n".join(json.dumps(r) for r in bcf_rows) + "\n")

    m.run_compare()
    out = json.loads(out_json.read_text())
    assert out["n_compared"] == 2 and out["n_errors"] == 0
    assert out["n_scan_vintage_new"] == 1 and out["n_scan_vintage_legacy"] == 1
    assert out["n_new_schema_target_concordant"] == 1
    assert out["new_schema_target_concordance"] == 1.0
    assert out["n_legacy_discordant"] == 1
    assert out["n_legacy_discordant_explained_by_legacy_filter"] == 1
    # the one-row edge delta is the documented coordinate convention
    assert out["n_window_rowcount_offby1"] == 1
    assert out["window_rowcount_offby1_fraction"] == 1.0
