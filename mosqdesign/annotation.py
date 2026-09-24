"""Exon-structure annotation for the dsx locus (from the NCBI gene table).

The gene table lists genomic exon intervals for the representative transcript
(XM_061649168). Used to score functional constraint: exonic sites near splice
junctions are where indels most reliably destroy gene function - the property that
makes end-joining resistance alleles nonfunctional (Kyrou 2018's design insight).
"""
from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class Exon:
    start: int  # inclusive, genomic
    end: int    # inclusive, genomic (start <= end after normalization)
    length: int


def parse_gene_table(path: str) -> dict[str, list[Exon]]:
    """Parse every per-transcript exon table. Returns {transcript_id: [Exon, ...]}."""
    tables: dict[str, list[Exon]] = {}
    current: str | None = None
    for line in open(path):
        m = re.search(r"Exon table for\s+mRNA\s+(XM_\d+(?:\.\d+)?)", line)
        if m:
            current = m.group(1)
            tables[current] = []
            continue
        m = re.match(r"^(\d+)-(\d+)\s", line)
        if m and current is not None:
            a, b = int(m.group(1)), int(m.group(2))
            tables[current].append(Exon(min(a, b), max(a, b), abs(b - a) + 1))
    return {k: sorted(v, key=lambda e: e.start) for k, v in tables.items()}


def exon_union(tables: dict[str, list[Exon]]) -> list[Exon]:
    """Merged exonic intervals across all transcript variants."""
    all_exons = sorted((e for exons in tables.values() for e in exons), key=lambda e: e.start)
    merged: list[Exon] = []
    for e in all_exons:
        if merged and e.start <= merged[-1].end + 1:
            prev = merged[-1]
            merged[-1] = Exon(prev.start, max(prev.end, e.end), 0)
        else:
            merged.append(Exon(e.start, e.end, 0))
    for e in merged:
        e.length = e.end - e.start + 1
    return merged


def exon_at(pos: int, exons: list[Exon]) -> Exon | None:
    for e in exons:
        if e.start <= pos <= e.end:
            return e
    return None


def dist_to_splice_junction(pos: int, exon: Exon) -> int:
    """Distance to the nearest end of the containing exon (0 = at the junction base)."""
    return min(pos - exon.start, exon.end - pos)


def dist_to_nearest_exon_edge(pos: int, exons: list[Exon]) -> int:
    """Distance to the nearest exon boundary, inside or outside exons."""
    best = None
    for e in exons:
        if e.start <= pos <= e.end:
            d = min(pos - e.start, e.end - pos)
        else:
            d = min(abs(pos - e.start), abs(pos - e.end))
        best = d if best is None else min(best, d)
    return best if best is not None else 10**9


def constraint_bonus(site_start_genomic: int, site_len: int, exons: list[Exon],
                     junction_window: int = 15) -> tuple[float, str]:
    """A priori bonus: majority-exonic site +0.25; any base of the site within
    `junction_window` of an exon boundary (inside or outside the exon) +0.5 more.
    Returns (bonus, label)."""
    positions = list(range(site_start_genomic, site_start_genomic + site_len))
    n_exonic = sum(1 for p in positions if exon_at(p, exons) is not None)
    best_jdist = min(dist_to_nearest_exon_edge(p, exons) for p in positions)
    bonus = 0.0
    labels = []
    if n_exonic > site_len / 2:
        bonus += 0.25
        labels.append(f"exonic({n_exonic}/{site_len})")
    elif n_exonic:
        labels.append(f"partially_exonic({n_exonic}/{site_len})")
    else:
        labels.append("intronic/intergenic")
    if best_jdist < junction_window:
        bonus += 0.5
        labels.append(f"splice_proximal(d={best_jdist})")
    return bonus, ",".join(labels)
