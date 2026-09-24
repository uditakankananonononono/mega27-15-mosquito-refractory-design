"""BLAT placement cross-check (tool 38): results/blat_placement.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "blat_placement.json")
M = os.path.join(BASE, "results", "blat_placement.md")


def load():
    return json.load(open(J))


def test_single_marginal_locus_block_both_configs():
    d = load()
    for run in d["runs"]:
        assert run["hits_with_block_in_dsx_locus"] == 1
        bh = run["best_locus_hit"]
        assert bh["score"] == 66 and bh["strand"] == "++"
        assert bh["blocks"][0][0] == 47619289  # the single in-locus block
        assert run["n_hits_outscoring_best_locus_hit"] == 13
        assert run["top_score_overall"] == 102


def test_parsing_artifact_documented():
    d = load()
    assert "reverse-complement" in d["parsing_caveat"]
    joined = " ".join(d["findings"])
    assert "70.5 Mb" in joined and "114-116 Mb" in joined


def test_sensitivity_ordering_stated():
    joined = " ".join(load()["findings"])
    assert "miniprot" in joined and "HMMER" in joined


def test_md_exists():
    assert os.path.exists(M) and "BLAT" in open(M).read()
