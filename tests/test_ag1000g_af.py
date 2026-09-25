"""Tests for scripts/ag1000g_allele_freq.py (Ag1000G on-target AF analysis)."""
import importlib.util
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    "agaf", os.path.join(ROOT, "scripts", "ag1000g_allele_freq.py"))
agaf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agaf)


def hit(pos, ref, alt, gt, gq=60, dp=20):
    return {"pos": pos, "ref": ref, "alt": alt, "gt": gt,
            "gq": gq, "dp": dp, "targets": []}


def rec(sid, hits):
    return {"sample": sid, "n_rows": 3250, "n_nonref": len(hits),
            "target_hits": hits}


def test_position_map_seed_and_roles():
    pm = agaf.position_map("dsx-v3-1")
    assert pm[48711501]["role"] == "seed"          # PAM-proximal base
    assert pm[48711501]["dist_from_pam"] == 0
    assert pm[48711508]["role"] == "seed"          # dist 7, last seed base
    assert pm[48711509]["role"] == "proto"         # dist 8
    assert pm[48711520]["dist_from_pam"] == 19
    assert pm[48711500]["role"] == "PAM-N"
    assert pm[48711499]["role"] == "PAM-G"
    assert pm[48711498]["role"] == "PAM-G"


def test_guide_allele_is_intact_at_assembly_mismatch():
    # v3-2 pos 48,712,765: AgamP4 plus ref A mismatches the guide; guide
    # allele is G. A sample homozygous G is fully intact (pilot: 12/12).
    pm = agaf.position_map("dsx-v3-2")
    assert pm[48712765]["ref_plus"] == "A"
    assert pm[48712765]["guide_plus"] == "G"
    assert not agaf.allele_disrupts("G", pm[48712765])
    assert agaf.allele_disrupts("A", pm[48712765])   # reference = mismatch!
    assert agaf.allele_disrupts("T", pm[48712765])


def test_multiallelic_split_and_hom_alt():
    records = [rec("S1", [hit(48712765, "A", "G,C", "2/2")]),
               rec("S2", [])]
    res = agaf.compute(records)
    site = res["sites"][48712765]
    assert site["alt_ac"] == {"C": 2}
    assert site["an"] == 4                     # both samples callable
    assert site["ref_ac"] == 2                 # S2 hom-ref
    assert site["alt_af"]["C"] == 0.5
    assert res["sample_status"]["S1"]["dsx-v3-2"] == 2   # hom compromised
    assert res["sample_status"]["S2"]["dsx-v3-2"] == 2  # absent row = AgamP4 hom-ref, mismatches guide


def test_lowqual_excluded_and_denominator_shrinks():
    records = [rec("S1", [hit(48712765, "A", "G", "1/1", gq=9)]),
               rec("S2", [hit(48712765, "A", "G", "0/1", gq=99)]),
               rec("S3", [])]
    res = agaf.compute(records)
    site = res["sites"][48712765]
    assert site["n_lowqual"] == 1
    assert site["n_called"] == 1
    assert site["an"] == 4                     # S1's 2 alleles removed
    assert site["alt_ac"] == {"G": 1}
    assert site["ref_ac"] == 2 + 1             # S3 hom-ref + S2 ref allele
    assert res["sample_status"]["S1"]["dsx-v3-2"] is None  # failed quality leaves uncertain genotype


def test_pam_g_disrupts_pam_n_tolerated():
    records = [rec("S1", [hit(48711499, "C", "T", "0/1"),   # G2: het disrupt
                          hit(48711500, "C", "A", "1/1")]), # N: tolerated
               rec("S2", [])]
    res = agaf.compute(records)
    assert res["sample_status"]["S1"]["dsx-v3-1"] == 2  # another AgamP4 reference position mismatches guide
    pm = agaf.position_map("dsx-v3-1")
    assert agaf.allele_disrupts("T", pm[48711499])
    assert not agaf.allele_disrupts("A", pm[48711500])


def test_unphased_het_at_two_sites_called_het():
    # both sites are guide==ref positions: 0/1 gt -> exactly one disrupted allele
    records = [rec("S1", [hit(48712770, "T", "C", "0/1"),   # plus ref T (pos 6)
                          hit(48712775, "G", "A", "0/1")])] # plus ref G (pos 11)
    res = agaf.compute(records)
    assert res["sample_status"]["S1"]["dsx-v3-2"] == 2       # reference at assembly-mismatch seed site is disruptive


def test_guide_summary_counts():
    records = [rec("S1", [hit(48712765, "A", "G,C", "2/2")]),   # hom comp
               rec("S2", [hit(48712765, "A", "G,C", "1/2")]),   # G/C het
               rec("S3", [hit(48712765, "A", "G", "1/1")]),     # intact
               rec("S4", [])]                                   # intact
    res = agaf.compute(records)
    s = agaf.guide_summary(res)["dsx-v3-2"]
    assert s["hom_compromised"] == 2 and s["het"] == 1 and s["intact"] == 1
    assert s["carrier_fraction"] == 0.75
    assert s["compromised_allele_fraction"] == 5 / 8


def test_windows_match_reannotation_verified_bases():
    # pins the verified AgamP4 reference/aligned bases used downstream
    pm = agaf.position_map("kyrou")
    assert "".join(pm[p]["ref_plus"] for p in range(48714640, 48714660)) == \
        "CCGCTTGACCTGTGTTAAAC"
    pm2 = agaf.position_map("dsx-v3-2")
    assert pm2[48712765]["dist_from_pam"] == 0   # the assembly-sensitive base


def test_fetch_error_is_not_a_homozygous_reference_observation():
    records = [rec("S1", [hit(48711498, "C", "T", "0/1")]),
               {"sample": "S2", "error": "upstream timeout"}]
    res = agaf.compute(records)
    assert res["n_samples"] == 1
    assert res["n_fetch_errors"] == 1
    assert res["error_sample_ids"] == ["S2"]
    assert "S2" not in res["sample_status"]
    assert res["sites"][48711498]["an"] == 2
    assert agaf.guide_summary(res)["dsx-v3-1"]["compromised_allele_fraction"] == 1.0


def test_legacy_hits_are_scored_not_silently_ignored():
    old = {"sample": "S1", "hits": [hit(48711498, "C", "T", "0/1")]}
    new = rec("S2", [])
    res = agaf.compute([old, new])
    assert res["n_samples"] == 2
    assert res["sites"][48711498]["alt_ac"] == {"T": 1}
    assert agaf.guide_summary(res)["dsx-v3-1"]["hom_compromised"] == 2


def test_hom_reference_at_assembly_mismatch_is_compromised():
    res = agaf.compute([rec("S1", [])])
    assert res["sample_status"]["S1"]["dsx-v3-1"] == 2
    assert res["sample_status"]["S1"]["dsx-v3-2"] == 2
    assert res["sample_status"]["S1"]["kyrou"] == 0

def test_uncallable_excluded_from_guide_denominator():
    res = agaf.compute([rec("S1", [hit(48712765, "A", "G", "1/1", gq=10)]),
                        rec("S2", [hit(48712765, "A", "G", "1/1")])])
    z = agaf.guide_summary(res)["dsx-v3-2"]
    assert z["n_callable"] == 1 and z["intact"] == 1
