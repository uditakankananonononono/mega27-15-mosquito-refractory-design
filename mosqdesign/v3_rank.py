"""v3 ranking: constraint-aware score fused with chromosome-wide specificity.

v2 ranks by CNN efficacy + exon/splice constraint bonus but its off-target
term came from the 110 kb locus only. The chromosome-2 scan (118 Mb) supplies
the real specificity burden per guide. v3 applies:
  - hard filter: any exact or 1-mismatch genomic off-target disqualifies
    (mm0 + mm1 == 0 required);
  - penalty: 0.50 per 2-mismatch hit, 0.05 per 3-mismatch hit
    (CFD-style decay: distal mismatches tolerate better).
The Kyrou 2018 guide is carried as a published control row.
"""
from __future__ import annotations

MM2_PENALTY = 0.50
MM3_PENALTY = 0.05


def specificity_penalty(mm2: int, mm3: int) -> float:
    return MM2_PENALTY * mm2 + MM3_PENALTY * mm3


def passes_specificity_filter(mm0: int, mm1: int) -> bool:
    return mm0 == 0 and mm1 == 0


def score_v3(efficacy: float, constraint_bonus: float, mm2: int, mm3: int) -> float:
    return efficacy + constraint_bonus - specificity_penalty(mm2, mm3)
