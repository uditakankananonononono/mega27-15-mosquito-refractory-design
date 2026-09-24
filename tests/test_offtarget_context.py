"""Consistency tests for results/offtarget_context_bedtools.json (tool 30/40).

Hermetic: reads only committed results files. Verifies the bedtools
functional-context audit against the committed bowtie2 census histogram.
"""
import json
import os

ROOT = os.path.join(os.path.dirname(__file__), "..")


def load():
    with open(os.path.join(ROOT, "results", "offtarget_context_bedtools.json")) as fh:
        return json.load(fh)


def test_tier_totals_match_bowtie2_census():
    d = load()
    tiers = d["per_mm_tier_total_inCDS"]
    # committed bowtie2 library histogram (results/offtarget_bowtie2.json):
    # MM0 427 (incl. 421 self), MM1 227, MM2 2744, MM3 20880
    assert tiers["0"][0] == 427
    assert tiers["1"][0] == 227
    assert tiers["2"][0] == 2744
    assert tiers["3"][0] == 20880


def test_in_cds_never_exceeds_total():
    d = load()
    for total, in_cds in d["per_mm_tier_total_inCDS"].values():
        assert 0 <= in_cds <= total


def test_leads_have_no_cleavable_cds_near_cognate():
    d = load()
    leads = d["leads"]
    # no lead has any MM1 or MM2 near-cognate at all (bowtie2 census)
    for name in ("dsx-v3-1", "dsx-v3-2", "dsx-v3-3", "kyrou"):
        tiers = leads[name]["mm_tiers_total_inCDS"]
        assert "1" not in tiers and "2" not in tiers
    # the only CDS-overlapping MM3 sites (v3-3, kyrou) are resolved PAM-less
    res = d["lead_cds_hits_pam_resolution"]
    assert set(res) == {"dsx-v3-3", "kyrou"}
    for rec in res.values():
        assert rec["pam"] is None
        assert "Cas9-inert" in rec["verdict"]


def test_md_mirrors_json_counts():
    d = load()
    with open(os.path.join(ROOT, "results", "offtarget_context_bedtools.md")) as fh:
        md = fh.read()
    assert str(d["library_guides_with_at_least_one_CDS_overlapping_near_cognate"]) in md
    for _, (total, in_cds) in d["per_mm_tier_total_inCDS"].items():
        assert f"{total:,}" in md and f"{in_cds:,}" in md
