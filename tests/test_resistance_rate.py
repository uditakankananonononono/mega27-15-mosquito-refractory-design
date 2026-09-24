"""Hermetic tests for mosqdesign/resistance_rate.py (synthetic sequence)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mosqdesign.resistance_rate import (  # noqa: E402
    cut_site_plus, microhomology_deletions, mutational_target_size,
    predict, window_pos_sets, GUIDE_WINDOWS, REGION_OFFSET)


def test_cut_site_three_nt_from_pam():
    # minus-strand guides: PAM upstream, cut 3 nt into proto from PAM end
    for g, w in GUIDE_WINDOWS.items():
        pa, pb = w["proto"]
        sa, sb = w["pam"]
        assert sb + 1 == pa                      # PAM adjacent to proto
        assert cut_site_plus(g) == pa + 2        # cut between +2|+3


def test_mutational_target_size():
    # 2 PAM-G positions x 3 + 8 seed positions x 3 = 30
    assert mutational_target_size("kyrou") == 30


def test_microhomology_detection_synthetic():
    # plant a 4-bp homology: positions 10-13 == positions 19-22; a 9-bp
    # deletion starting at 14 removes 14..22 with 4-bp breakpoint homology
    seq = list("ACGT" * 40)          # 160 bp
    seq[10:14] = list("GGCC")
    seq[19:23] = list("GGCC")
    seq = "".join(seq)
    pam_pos, seed_pos = set(), set()
    dels = microhomology_deletions(seq, cut=40, pam_pos=pam_pos,
                                   seed_pos=seed_pos, mh_min=2,
                                   max_del=40, flank=60)
    found = [d for d in dels if d["start"] == 14 and d["del_len"] == 9]
    assert found and found[0]["mh_len"] >= 4


def test_resistance_classification_overlaps_seed():
    seq = "ACGT" * 60
    pam_pos = set()
    seed_pos = set(range(100, 108))  # seed far from any MH
    dels = microhomology_deletions(seq, cut=40, pam_pos=pam_pos,
                                   seed_pos=seed_pos, mh_min=2,
                                   max_del=20, flank=40)
    # any deletion whose interval misses 100..107 is non-resistance
    assert all(d["resistance_generating"] ==
               bool(set(range(d["start"], d["start"] + d["del_len"]))
                    & seed_pos)
               for d in dels)


def test_predict_runs_on_region_fasta():
    path = os.path.join(os.path.dirname(__file__), "..",
                        "data", "agamp4", "AgamP4_2R_dsx_region.fasta")
    if not os.path.exists(path):
        return  # data file not present in minimal CI checkout
    from mosqdesign.resistance_rate import load_region
    seq = load_region(path)
    p = predict("kyrou", seq)
    assert p["mutational_target_size_L"] == 30
    assert p["mh_deletions_total"] > 0
    assert 0.0 <= p["mh_resistance_fraction_w"] <= 1.0
    lo, hi = p["e_interval"]
    assert 0 < lo <= hi


def test_window_pos_sets_bounds():
    pam, seed = window_pos_sets("kyrou")
    sa, sb = GUIDE_WINDOWS["kyrou"]["pam"]
    pa, _ = GUIDE_WINDOWS["kyrou"]["proto"]
    assert pam == {sa - REGION_OFFSET, sa + 1 - REGION_OFFSET}
    assert len(seed) == 8 and min(seed) == pa - REGION_OFFSET
