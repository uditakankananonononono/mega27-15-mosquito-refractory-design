"""DIAMOND census cross-check (tool 37): results/diamond_census.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "diamond_census.json")
M = os.path.join(BASE, "results", "diamond_census.md")


def load():
    return json.load(open(J))


def test_three_gene_census_reproduced():
    c = load()["census"]
    assert c["genes_with_DM_domain_homology"] == 3
    assert c["dsx_isoforms_hit"] == 9 and c["dsx_accessions_missing"] == []
    assert len(c["dmrta2_hits"]) == 2 and len(c["dmd4_hits"]) == 1
    assert c["background_hits"] == 0


def test_neighbour_accession_not_hit():
    targets = [h["target"] for h in load()["hits"]]
    assert "XP_061505154.1" not in targets  # LOC1272222, adjacent accession block
    assert "XP_061505155.1" in targets


def test_paralog_separation():
    hits = {h["target"]: h["evalue"] for h in load()["hits"]}
    best_dsx = min(h["evalue"] for h in load()["hits"] if h["gene"] == "LOC1270904")
    best_paralog = min(hits["XP_061504876.1"], hits["XP_310668.5"])
    assert best_dsx < 1e-30 and best_paralog > 1e-20  # >= 10 orders separation


def test_md_exists():
    assert os.path.exists(M) and "3/3" in open(M).read()
