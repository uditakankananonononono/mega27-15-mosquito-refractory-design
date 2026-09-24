from mosqdesign.v3_rank import (passes_specificity_filter, score_v3,
                                specificity_penalty)


def test_filter_rejects_any_exact_or_single_mismatch_hit():
    assert passes_specificity_filter(0, 0)
    assert not passes_specificity_filter(1, 0)
    assert not passes_specificity_filter(0, 1)


def test_penalty_scales_with_mismatch_closeness():
    assert specificity_penalty(1, 0) > specificity_penalty(0, 1) > 0
    assert specificity_penalty(2, 10) == 2 * 0.50 + 10 * 0.05


def test_score_v3_combines_terms():
    assert score_v3(1.2, 0.75, 1, 4) == 1.2 + 0.75 - (0.50 + 0.20)
    assert score_v3(1.0, 0.0, 0, 0) == 1.0
