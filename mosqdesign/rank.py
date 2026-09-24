"""Design ranking: combine CNN efficacy with sequence-quality and off-target penalties."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mosqdesign.grna_scan import GrnaSite, homopolymer_run
from mosqdesign.models.efficacy_cnn import gc_content


@dataclass
class RankedDesign:
    site: GrnaSite
    efficacy: float
    offtarget_hits: int
    gc: float
    homopolymer: int
    score: float


def composite_score(efficacy: float, offtarget_hits: int, gc: float, homopolymer: int) -> float:
    """Higher is better. Penalties chosen a priori and documented in the paper:
    - each intra-panel off-target (<=3 mismatches): -0.5
    - GC outside [0.3, 0.8]: -0.2
    - homopolymer run >= 4: -0.1 per base over 3
    """
    score = float(efficacy)
    score -= 0.5 * offtarget_hits
    if not 0.3 <= gc <= 0.8:
        score -= 0.2
    if homopolymer >= 4:
        score -= 0.1 * (homopolymer - 3)
    return score


def rank_sites(sites: list[GrnaSite], efficacies: np.ndarray, offtarget_hits: list[int]) -> list[RankedDesign]:
    ranked = []
    for site, eff, ot in zip(sites, efficacies, offtarget_hits):
        gc = gc_content(site.protospacer)
        hp = homopolymer_run(site.protospacer)
        ranked.append(RankedDesign(site, float(eff), ot, gc, hp,
                                   composite_score(float(eff), ot, gc, hp)))
    return sorted(ranked, key=lambda r: r.score, reverse=True)
