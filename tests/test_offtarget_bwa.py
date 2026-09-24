"""BWA census + correction (tool 31/40): results/offtarget_bwa.json integrity."""
import json, os, re

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "offtarget_bwa.json")
M = os.path.join(BASE, "results", "offtarget_bwa.md")


def load():
    return json.load(open(J))


def test_census_tiers_and_total():
    d = load()
    assert d["census_per_mm_tier"] == {"0": 427, "1": 227, "2": 2886, "3": 28765}
    assert d["census_total"] == 32305
    assert d["guides"] == 421


def test_superset_and_correction():
    d = load()
    s = d["superset_check"]
    assert s["bowtie2_alignments_missing_from_bwa"] == 0
    assert s["bwa_only_alignments"] == 8027
    assert s["bwa_only_by_tier"] == {"2": 142, "3": 7885}
    assert d["bowtie2_census_total"] == 24278
    assert "32,305" in d["correction"]


def test_sequence_verification():
    v = load()["bwa_only_sequence_verification"]
    assert v["sampled"] == 50 and v["verified_hamming_le_3"] == 50 and v["seed"] == 7


def test_leads_unchanged_and_pam_less():
    leads = load()["leads"]
    assert leads["dsx-v3-1"]["cds_overlap_sites"] == []
    assert leads["dsx-v3-2"]["cds_overlap_sites"] == []
    # no MM1/MM2 for any lead (superset adds no new 1-mismatch site anywhere)
    for name in ("dsx-v3-1", "dsx-v3-2", "dsx-v3-3", "kyrou"):
        tiers = leads[name]["mm_tiers_total_inCDS"]
        assert "1" not in tiers and "2" not in tiers
    v33 = leads["dsx-v3-3"]["cds_overlap_sites"]
    kyr = leads["kyrou"]["cds_overlap_sites"]
    assert v33 == [["NC_064601.1", 38008162, 3, "LOC1270942"]]
    assert kyr == [["NC_064601.1", 41285163, 3, "LOC1274366"]]


def test_cds_context_totals():
    ctx = load()["cds_context_on_corrected_set"]
    assert ctx["per_mm_tier_total_inCDS"] == {
        "0": [427, 244], "1": [227, 29], "2": [2886, 412], "3": [28765, 3452]}
    assert ctx["library_guides_with_at_least_one_CDS_overlapping_near_cognate"] == 308


def test_md_mirrors_json():
    d = load()
    md = open(M).read()
    for n in (d["census_total"], d["census_per_mm_tier"]["3"],
              d["superset_check"]["bwa_only_alignments"],
              d["cds_context_on_corrected_set"]
              ["library_guides_with_at_least_one_CDS_overlapping_near_cognate"]):
        assert f"{n:,}" in md
