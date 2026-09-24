"""Scan gene sequences for SpCas9 (NGG) target sites and build CNN-ready contexts."""
from __future__ import annotations

from dataclasses import dataclass

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]


@dataclass
class GrnaSite:
    gene: str
    strand: str          # '+' or '-'
    position: int        # 0-based protospacer start on the scanned (plus) sequence
    protospacer: str     # 20-nt, 5'->3' as the guide would bind
    pam: str             # 3-nt PAM
    context30: str       # 4bp upstream + protospacer + PAM + 3bp downstream (padded with N at edges)


def _context(seq: str, start: int) -> str:
    """30-mer context for a protospacer starting at `start` on the plus strand."""
    s = seq.upper()
    core = s[start:start + 23]          # 20-nt protospacer + 3-nt PAM
    up = s[max(0, start - 4):start]
    down = s[start + 23:start + 26]
    return ("N" * (4 - len(up)) + up + core + down + "N" * (3 - len(down)))


def scan_ngg_sites(seq: str, gene: str = "seq") -> list[GrnaSite]:
    """All NGG PAM sites on both strands of `seq`."""
    seq = seq.upper()
    sites: list[GrnaSite] = []
    for i in range(len(seq) - 23 + 1):
        window = seq[i:i + 23]
        pam = window[20:23]
        if pam[1:] == "GG" and pam[0] in "ACGT" and all(b in "ACGT" for b in window):
            sites.append(GrnaSite(gene, "+", i, window[:20], pam, _context(seq, i)))
        rc_window = reverse_complement(window)
        # minus-strand site: PAM (CCN) at the 5' end of the plus-strand window
        if window[:2] == "CC" and window[2] in "ACGT" and all(b in "ACGT" for b in window):
            rc = reverse_complement(window)
            rc_context = reverse_complement(_context(seq, i))
            # reverse-complementing the padded context preserves the 4+20+3+3 layout
            sites.append(GrnaSite(gene, "-", i, rc[:20], rc[20:23], rc_context))
    return sites


def count_offtargets(site: GrnaSite, panel: list[tuple[str, str]], max_mismatches: int = 3) -> int:
    """Count panel sequences (name, seq) containing a <=max_mismatches near-match
    to the site's protospacer on either strand. Excludes the site's own locus."""
    proto = site.protospacer
    rc_proto = reverse_complement(proto)
    hits = 0
    L = len(proto)
    for name, s in panel:
        s = s.upper()
        for k in range(len(s) - L + 1):
            w = s[k:k + L]
            if sum(1 for a, b in zip(w, proto) if a != b) <= max_mismatches:
                if not (name == site.gene and site.position <= k + L and k <= site.position + L):
                    hits += 1
            if sum(1 for a, b in zip(w, rc_proto) if a != b) <= max_mismatches:
                hits += 1
    return hits


def homopolymer_run(seq: str) -> int:
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best
