"""miniprot genome-wide census (tool 33 extension): results/miniprot_census.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "miniprot_census.json")
M = os.path.join(BASE, "results", "miniprot_census.md")


def load():
    return json.load(open(J))


def test_dsx_recovered_chr2_full_structure():
    scans = load()["per_contig_scans"]
    chr2 = scans["NC_064601.1"]["mrna_features"]
    assert len(chr2) == 1
    f = chr2[0]
    assert (f["start"], f["end"], f["strand"]) == (47545556, 47692565, "-")
    for acc in ("NC_064602.1", "NC_064600.1", "NC_083487.1", "NW_scaffolds_187"):
        assert scans[acc]["mrna_features"] == []


def test_memory_engineering_documented():
    me = load()["memory_engineering"]
    assert "OOM" in me["whole_genome"] and "OOM" in me["chr2_default_M1"]
    assert "fits" in me["M2"] and "loses ALL hits" in me["M3"]
    chr2_peak = load()["per_contig_scans"]["NC_064601.1"]["peak_rss_gb"]
    assert 0.5 < chr2_peak < 1.5


def test_paralogs_missed_but_controls_map():
    pl = load()["paralog_locus_checks"]
    for gene in ("dmd4", "dmrta2"):
        assert pl[gene]["dsx_query_default_sensitivity_mrnas"] == []
        assert pl[gene]["dsx_query_loose_chaining_paf_lines"] == 0
        ctrl = pl[gene]["native_protein_control_mrnas"]
        assert len(ctrl) == 1 and "Identity=1.0000" in ctrl[0]["attrs"]


def test_verdict_and_md():
    v = load()["verdict"]
    assert "1/3" in v and "HMMER" in v
    md = open(M).read()
    assert "sensitivity limit" in md and "3/3" in md
