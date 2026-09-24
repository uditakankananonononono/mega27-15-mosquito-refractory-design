"""HMMER ortholog verification (tool 32/40): results/hmmer_dsx_ortholog.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "hmmer_dsx_ortholog.json")
M = os.path.join(BASE, "results", "hmmer_dsx_ortholog.md")


def load():
    return json.load(open(J))


def test_dsx_locus_confirmed():
    d = load()
    assert d["dsx_locus_hits"]["isoforms_hit"] == 9
    assert d["dsx_locus_hits"]["best_evalue"] == 3.6e-51
    assert "LOC1270904" in d["dsx_locus_hits"]["gene"]
    assert "47,610,877-47,700,920" in d["dsx_locus_hits"]["locus"]


def test_family_census():
    c = load()["dm_domain_family_census"]
    assert c["genes_with_DM_domain_homology"] == 3
    assert len(c["DmrtA2_hits"]) == 2 and len(c["dmd4_hits"]) == 1
    assert c["best_non_dsx_evalue"] == 4.1e-21
    assert c["separation_orders_of_magnitude"] == 30.1


def test_drift_caveat_and_jackhmmer():
    d = load()
    assert d["jackhmmer_2_iter"]["dsx_best_evalue"] == 4.1e-104
    assert "1.1e-62" in d["jackhmmer_2_iter"]["contamination_caveat"]
    assert d["phmmer_hits_E_lt_1e-3"] == 25


def test_md_mirrors_json():
    d = load()
    md = open(M).read()
    assert "3.6e-51" in md and "4.1e-104" in md and "30.1" in md
    assert "LOC1270904" in md
