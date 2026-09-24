"""miniprot DM-domain mapping (tool 33/40): results/miniprot_dm_domain.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "miniprot_dm_domain.json")
M = os.path.join(BASE, "results", "miniprot_dm_domain.md")


def load():
    return json.load(open(J))


def test_dm_mapping():
    d = load()["dm_domain_mapping"]
    assert d["genomic_span"] == "NC_064601.1:47692251-47692565, minus strand"
    assert d["identity"] == 0.713 and d["positive"] == 0.8241
    assert len(d["refseq_exons_overlapping"]) == 3


def test_lead_exon_audit():
    la = load()["lead_sites_vs_refseq_annotation"]
    assert la["dsx-v3-2"]["in_refseq_exon"] and la["dsx-v3-3"]["in_refseq_exon"]
    assert not la["dsx-v3-1"]["in_refseq_exon"] and not la["kyrou"]["in_refseq_exon"]
    assert la["kyrou"]["distance_to_nearest_exon_boundary_bp"] == 523
    assert la["dsx-v3-1"]["distance_to_nearest_exon_boundary_bp"] == 30


def test_limitations_recorded():
    d = load()
    assert len(d["honest_limitations"]) == 3
    assert "blocked" in d["exonerate_status"]


def test_md_mirrors_json():
    md = open(M).read()
    for s in ("47,692,251-47,692,565", "71.3%", "523 bp", "30 bp"):
        assert s in md
