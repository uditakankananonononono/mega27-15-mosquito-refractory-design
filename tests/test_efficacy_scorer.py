import os
import numpy as np
import pytest

from mosqdesign.models.efficacy_cnn import EfficacyScorer, one_hot

CKPT = os.path.join(os.path.dirname(__file__), "..", "assets", "doench2016_efficacy_cnn.pt")


@pytest.mark.skipif(not os.path.exists(CKPT), reason="checkpoint not bundled")
def test_scorer_outputs_finite_scores():
    scorer = EfficacyScorer(CKPT)
    scores = scorer.score(["N" * 4 + "ACGT" * 5 + "TGG" + "N" * 3,
                           "N" * 4 + "AAAA" * 5 + "AGG" + "N" * 3])
    assert scores.shape == (2,)
    assert np.isfinite(scores).all()
    assert scores[0] != scores[1]


def test_one_hot_basic():
    arr = one_hot("ACGT", 6)
    assert arr.shape == (4, 6)
    assert arr[:, :4].sum() == 4
