"""Genome-completeness correction (BWA rerun): results/offtarget_bwa_complete.json integrity."""
import json, os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(BASE, "results", "offtarget_bwa_complete.json")
M = os.path.join(BASE, "results", "offtarget_bwa_complete.md")


def load():
    return json.load(open(J))


def test_complete_census_tiers_and_total():
    d = load()
    assert d["census_per_mm_tier"] == {"0": 427, "1": 227, "2": 2974, "3": 29308}
    assert d["census_total"] == 32936
    assert d["prior_census_4record_total"] == 32305
    assert d["guides"] == 421


def test_scaffold_delta_no_mm0_mm1():
    d = load()
    assert d["scaffolds_per_mm_tier"] == {"2": 88, "3": 543}
    assert "0" not in d["scaffolds_per_mm_tier"] and "1" not in d["scaffolds_per_mm_tier"]
    assert d["chromosomes_per_mm_tier"] == {"0": 427, "1": 227, "2": 2886, "3": 28765}


def test_forward_map_and_verification():
    d = load()
    assert d["forward_map_check"]["old_alignments_missing_from_complete"] == 0
    v = d["scaffold_hit_verification"]
    assert v["sampled"] == 50 and v["verified_hamming_le_3"] == 50 and v["seed"] == 11
    assert v["scaffold_alignments_total"] == 631


def test_correction_decomposition():
    d = load()["correction_decomposition"]
    assert d["previously_reported_guides_with_cds_near_cognate"] == 308
    assert d["same_32305_sites_accession_fixed"] == 367
    assert d["complete_genome"] == 367
    assert d["scaffold_only_new_guides"] == 0
    assert "accession" in d["explanation"]


def test_v32_new_cds_site_pam_less():
    v = load()["v32_new_cds_site_verification"]
    assert v["site"] == "NC_064602.1:2111186-2111205"
    assert v["hamming"] == 3 and v["mismatch_positions_1based"] == [4, 5, 7]
    assert v["is_ngg"] is False and v["pam_3prime_3bp"] == "CCC"
    assert "aconitate hydratase" in v["gene"]
    leads = load()["leads"]
    assert leads["dsx-v3-2"]["cds_overlap_sites"] == [["NC_064602.1", 2111186, 3, "LOC1278105"]]
    assert leads["dsx-v3-1"]["cds_overlap_sites"] == []
    # no MM1/MM2 for any lead on the complete genome either
    for name in ("dsx-v3-1", "dsx-v3-2", "dsx-v3-3", "kyrou"):
        tiers = leads[name]["mm_tiers_total_inCDS"]
        assert "1" not in tiers and "2" not in tiers


def test_md_companion():
    md = open(M).read()
    assert "32,936" in md and "367/421" in md and "LOC1278105" in md
