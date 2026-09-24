"""Hermetic tests for scripts/ag1000g_near_cognate.py (synthetic records)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from ag1000g_allele_freq import GUIDES, position_map  # noqa: E402
from ag1000g_near_cognate import (  # noqa: E402
    classify_variant, enumerate_alleles, guide_index, summarize)


def hit(pos, ref, alt, gt, gq=90, dp=25):
    return {"pos": pos, "ref": ref, "alt": alt, "gt": gt, "gq": gq,
            "dp": dp, "targets": ["kyrou_proto"]}


def meta_at(pos):
    return position_map("kyrou")[pos]


def other_base(*avoid):
    for b in "ACGT":
        if b not in avoid:
            return b
    raise AssertionError


def test_guide_index_orientation():
    pa, pb = GUIDES["kyrou"]["proto"]
    assert guide_index(pb, "kyrou") == 1   # 5' guide end at highest coord
    assert guide_index(pa, "kyrou") == 20  # PAM-proximal base = index 20


def test_classify_seed_mismatch():
    pos = 48714642  # dist_from_pam = 2 -> seed
    m = meta_at(pos)
    assert m["role"] == "seed"
    alt = other_base(m["ref_plus"], m["guide_plus"])
    assert classify_variant(pos, alt, m) == "seed_mismatch"


def test_classify_distal_mismatch():
    pos = 48714655  # dist_from_pam = 15 -> distal proto
    m = meta_at(pos)
    assert m["role"] == "proto"
    alt = other_base(m["ref_plus"], m["guide_plus"])
    assert classify_variant(pos, alt, m) == "distal_mismatch"


def test_classify_pam_positions():
    # PAM-N at highest coordinate tolerated
    m_n = meta_at(48714639)
    assert m_n["role"] == "PAM-N"
    assert classify_variant(48714639, "A", m_n) == "pam_n_variant"
    # PAM-G: any non-C (plus) disrupts
    m_g = meta_at(48714637)
    assert m_g["role"] == "PAM-G"
    assert classify_variant(48714637, "A", m_g) == "pam_disrupted"
    assert classify_variant(48714637, "C", m_g) == "wild_type"


def test_single_variant_seed_allele_enumerated():
    pos = 48714642
    m = meta_at(pos)
    alt = other_base(m["ref_plus"], m["guide_plus"])
    recs = [{"sample": "S1", "target_hits": [hit(pos, m["ref_plus"], alt,
                                                 "1/1")]}]
    res = enumerate_alleles(recs)
    sv = res["guides"]["kyrou"]["single_variant_alleles"]
    assert len(sv) == 1
    a = sv[0]
    assert a["class"] == "seed_mismatch"
    assert a["ac"] == 2                      # homozygous -> 2 chromosomes
    assert a["mismatch_guide_indices"] == [guide_index(pos, "kyrou")]
    assert a["seed_mismatches"] == 1
    # allele sequence differs from guide at exactly one position
    diffs = [i for i, (x, y) in enumerate(
        zip(a["proto_allele_minus"], GUIDES["kyrou"]["guide"])) if x != y]
    assert len(diffs) == 1


def test_multi_variant_phase_certain_and_ambiguous():
    p1, p2 = 48714642, 48714655
    m1, m2 = meta_at(p1), meta_at(p2)
    a1 = other_base(m1["ref_plus"], m1["guide_plus"])
    a2 = other_base(m2["ref_plus"], m2["guide_plus"])
    hom = {"sample": "HOM", "target_hits":
           [hit(p1, m1["ref_plus"], a1, "1/1"),
            hit(p2, m2["ref_plus"], a2, "1/1")]}
    het = {"sample": "HET", "target_hits":
           [hit(p1, m1["ref_plus"], a1, "0/1"),
            hit(p2, m2["ref_plus"], a2, "0/1")]}
    res = enumerate_alleles([hom, het])
    g = res["guides"]["kyrou"]
    # phase-certain multi allele from HOM only
    assert len(g["multi_variant_phase_certain"]) == 1
    multi = g["multi_variant_phase_certain"][0]
    assert multi["ac"] == 2
    assert multi["n_mismatch_vs_guide"] == 2
    assert multi["seed_mismatches"] == 1    # p1 seed, p2 distal
    assert multi["class"] == "seed_mismatch"
    # HET unphased: counted as ambiguous, singles still counted
    assert g["n_ambiguous_multisite_samples"] == 1
    sv_ac = {a["pos"]: a["ac"] for a in g["single_variant_alleles"]}
    assert sv_ac[p1] == 3   # 2 (HOM) + 1 (HET)
    assert sv_ac[p2] == 3


def test_lowqual_calls_excluded():
    pos = 48714642
    m = meta_at(pos)
    alt = other_base(m["ref_plus"], m["guide_plus"])
    recs = [{"sample": "LQ", "target_hits": [hit(pos, m["ref_plus"], alt,
                                                 "1/1", gq=10)]}]
    res = enumerate_alleles(recs)
    assert res["n_lowqual_calls"] == 1
    assert res["guides"]["kyrou"]["single_variant_alleles"] == []


def test_summary_counts():
    pos = 48714639  # PAM-N
    recs = [{"sample": "S1", "target_hits": [hit(pos, "C", "A", "0/1")]}]
    res = enumerate_alleles(recs)
    s = summarize(res)
    assert s["kyrou"]["by_class_counts"].get("pam_n_variant") == 1
    assert s["kyrou"]["by_class_chromosomes"].get("pam_n_variant") == 1
