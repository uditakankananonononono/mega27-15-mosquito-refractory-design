from mosqdesign.grna_scan import GrnaSite
from mosqdesign.rank import composite_score, rank_sites
import numpy as np


def test_composite_prefers_high_efficacy_no_offtarget():
    good = composite_score(0.8, 0, 0.5, 1)
    bad = composite_score(0.8, 2, 0.5, 1)
    assert good > bad


def test_gc_penalty():
    assert composite_score(0.5, 0, 0.1, 1) < composite_score(0.5, 0, 0.5, 1)


def test_rank_orders_by_score():
    sites = [GrnaSite("g", "+", 0, "ATGCCGTAACGTTAGCCTGA", "TGG", "N" * 30),
             GrnaSite("g", "+", 10, "GATCCGTTAGCAAGTCCTAG", "AGG", "N" * 30)]
    ranked = rank_sites(sites, np.array([0.9, 0.1]), [0, 0])
    assert ranked[0].efficacy == 0.9
    assert ranked[0].score >= ranked[1].score
