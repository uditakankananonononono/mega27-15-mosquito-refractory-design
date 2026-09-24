"""MUSCLE+FastTree DM-family phylogeny (tools 39-40): results/dm_phylogeny.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "dm_phylogeny.json")
M = os.path.join(BASE, "results", "dm_phylogeny.md")


def load():
    return json.load(open(J))


def test_ortholog_clade():
    d = load()
    assert d["local_support_ortholog_clade"] >= 0.99
    assert "dsx_AGAM" in d["newick"] and "Dsx_DROME" in d["newick"]


def test_divergence_separation():
    d = load()
    assert d["dsx_nearest_paralog_distance"] > d["dsx_ortholog_distance"]
    assert d["divergence_ratio"] > 1.5


def test_all_four_taxa_distances():
    d = load()
    assert len(d["patristic_distances"]) == 6


def test_md_exists():
    assert os.path.exists(M) and "MUSCLE" in open(M).read()
