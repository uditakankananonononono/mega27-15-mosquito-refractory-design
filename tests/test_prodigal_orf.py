"""Prodigal coding-potency probe (tool 41): results/prodigal_orf.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "prodigal_orf.json")
M = os.path.join(BASE, "results", "prodigal_orf.md")


def load():
    return json.load(open(J))


def test_two_real_genes_are_top_signals():
    d = load()
    ctrl_max = d["control_window"]["max_orf_score"]
    assert d["nested_gene"]["orf"]["score"] > 3 * ctrl_max
    assert d["second_orf"]["score"] > 2 * ctrl_max
    assert d["nested_gene"]["orf"]["start"] == 26981  # LOC11175624, plus strand
    assert d["second_orf"]["start"] == 8265           # exact dsx exon stop boundary


def test_lead_exon_map_matches_miniprot_audit():
    m = load()["lead_exon_map"]
    assert m["dsx-v3-2"] == "8265-9460" and m["dsx-v3-3"] == "11821-11955"
    assert m["dsx-v3-1"] is None and m["kyrou"] is None


def test_caveat_and_md():
    d = load()
    assert "prokaryotic" in d["caveat"]
    assert os.path.exists(M) and "NESTED" in open(M).read().upper() or "nested" in open(M).read()
