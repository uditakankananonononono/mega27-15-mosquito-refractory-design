"""Pins the AgamP4 re-annotation of the Kyrou site and the v3 leads."""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = json.load(open(os.path.join(ROOT, "results", "agamp4_reannotation.json")))


def test_kyrou_site_junction_spanning_under_agamp4():
    g = RES["guides"]["kyrou_2018"]
    assert g["exact_match"] and g["mismatches_vs_AgamP4"] == 0
    # 9 exonic bp (48,714,640-48,714,648) inside RB exon 48,712,957-48,714,648
    ov = g["exon_overlaps"]["AGAP004050-RB"][0]
    assert ov[0:2] == [48_712_957, 48_714_648]
    assert ov[2:] == [48_714_640, 48_714_648]
    # and those 9 bp are coding (CDS 48,714,557-48,714,648)
    cov = g["cds_overlaps"]["AGAP004050-RB"][0]
    assert cov[2:] == [48_714_640, 48_714_648]
    assert RES["kyrou_junction_spanning"] is True


def test_female_exon_135bp_present():
    assert [48_715_161, 48_715_295] in RES["female_exon_135bp"]


def test_leads_are_assembly_sensitive():
    g2 = RES["guides"]["dsx-v3-2"]
    assert g2["mismatches_vs_AgamP4"] == 1
    # the single mismatch is the PAM-proximal seed base (guide position 20)
    assert g2["aligned_sequence"][:19] == g2["guide"][:19]
    assert g2["aligned_sequence"][19] != g2["guide"][19]
    g1 = RES["guides"]["dsx-v3-1"]
    assert g1["mismatches_vs_AgamP4"] == 3
