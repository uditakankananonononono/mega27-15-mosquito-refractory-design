"""Genome-scale off-target counting via block-sharded pigeonhole indexing.

The 110 kb dsx-locus scan keeps one k-mer index in memory. A chromosome-scale
scan (AgamP5 chromosome 2R is ~61 Mb; the full genome ~280 Mb) cannot, so the
reference is streamed in overlapping blocks: each block gets its own 5-mer
index, candidate positions are resolved to global coordinates, and a per-guide
set of global start positions deduplicates hits found in neighbouring block
overlaps. Block size is a memory knob, not an accuracy knob: with an overlap
of L-1 bases every 20-mer window is fully contained in at least one block.
"""
from __future__ import annotations

from dataclasses import dataclass

from mosqdesign.grna_scan import build_kmer_index, reverse_complement

try:
    import numpy as _np
    _LUT = _np.full(256, 4, dtype=_np.uint16)   # non-ACGT -> sentinel base 4
    for _i, _c in enumerate(b"ACGT"):
        _LUT[_c] = _i
except ImportError:                              # pragma: no cover
    _np = None


@dataclass
class OfftargetHit:
    guide: str          # protospacer sequence (5'->3')
    chrom: str          # chromosome / record name
    position: int       # 0-based start on the plus strand
    strand: str         # '+' if the plus-strand window matches the guide, '-' for its RC
    mismatches: int


def iter_fasta_blocks(path: str, block_bases: int = 5_000_000, overlap: int = 19):
    """Yield (record_name, global_start, sequence) blocks from a FASTA file.

    Blocks for the same record overlap by `overlap` bases so windows spanning a
    boundary are captured whole. Parsing is line-streamed; only one block of
    sequence is held at a time.
    """
    name = None
    pending = ""        # un-emitted tail of the current record
    emitted = 0         # bases of the current record already emitted
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None and pending:
                    yield name, emitted, pending
                name = line[1:].split()[0]
                pending, emitted = "", 0
                continue
            pending += line.upper()
            while len(pending) >= block_bases:
                keep = pending[:block_bases]
                yield name, emitted, keep
                emitted += block_bases - overlap
                pending = pending[block_bases - overlap:]
    if name is not None and pending:
        yield name, emitted, pending


def _codes5(seq: str) -> "_np.ndarray":
    """uint16 5-mer codes for every window; windows touching a non-ACGT base
    are set to 65535 so they never match a query code (0..1023)."""
    v = _LUT[_np.frombuffer(seq.encode(), dtype=_np.uint8)]
    n = len(v) - 4
    codes = v[:n].astype(_np.uint32)
    invalid = v[:n] > 3
    for shift in range(1, 5):
        codes = codes * 4 + v[shift:shift + n]
        invalid |= v[shift:shift + n] > 3
    out = codes.astype(_np.uint16)
    out[invalid] = 65535
    return out


def _code_of(mer: str) -> int:
    c = 0
    for b in mer:
        c = c * 4 + "ACGT".index(b)
    return c


def _candidate_starts(proto: str, codes: "_np.ndarray", L: int,
                      max_start: int) -> "_np.ndarray":
    """All block-relative window starts sharing an exact 5-mer block with the
    protospacer or its reverse complement (pigeonhole candidate generation)."""
    bl = L // 4
    parts = []
    for q in (proto, reverse_complement(proto)):
        for i in range(4):
            code = _code_of(q[i * bl:(i + 1) * bl])
            pos = _np.flatnonzero(codes == code)
            if pos.size:
                base = pos.astype(_np.int64) - i * bl
                parts.append((base[:, None] + _np.arange(bl)).ravel())
    if not parts:
        return _np.empty(0, dtype=_np.int64)
    starts = _np.unique(_np.concatenate(parts))
    return starts[(starts >= 0) & (starts <= max_start)]


def _verify(proto: str, s_arr: "_np.ndarray", starts: "_np.ndarray",
            L: int, max_mismatches: int):
    """Vectorised mismatch verification of candidate windows; yields
    (offset, strand, mismatches), one row per offset ('+' wins strand ties)."""
    if starts.size == 0:
        return
    rc = reverse_complement(proto)
    win = s_arr[starts[:, None] + _np.arange(L)]
    p_arr = _np.frombuffer(proto.encode(), dtype=_np.uint8)
    rc_arr = _np.frombuffer(rc.encode(), dtype=_np.uint8)
    mm_p = (win != p_arr).sum(axis=1)
    mm_m = (win != rc_arr).sum(axis=1)
    keep = (mm_p <= max_mismatches) | (mm_m <= max_mismatches)
    for st, a, b in zip(starts[keep], mm_p[keep], mm_m[keep]):
        yield int(st), ("+" if a <= b else "-"), int(min(a, b))


def _block_hits(proto: str, seq: str, max_mismatches: int):
    """Pigeonhole candidates in one block; yield (offset, strand, mismatches).

    numpy path: vectorised 5-mer codes + flatnonzero lookup; falls back to the
    pure-python k-mer index when numpy is unavailable. Semantics are identical
    (pinned by tests comparing both against the naive counter).
    """
    if _np is None:                              # pragma: no cover
        yield from _block_hits_slow(proto, seq, max_mismatches)
        return
    L = len(proto)
    codes = _codes5(seq) if len(seq) >= 5 else _np.array([], dtype=_np.uint16)
    s_arr = _np.frombuffer(seq.encode(), dtype=_np.uint8)
    starts = _candidate_starts(proto, codes, L, len(seq) - L)
    yield from _verify(proto, s_arr, starts, L, max_mismatches)


class BlockIndex:
    """Codes + byte view for one sequence block, built once and shared across
    every guide scanned against that block."""

    def __init__(self, seq: str):
        self.seq = seq
        if _np is not None and len(seq) >= 5:
            self.codes = _codes5(seq)
            self.s_arr = _np.frombuffer(seq.encode(), dtype=_np.uint8)
        else:                                    # pragma: no cover
            self.codes = None
            self.s_arr = None

    def hits(self, proto: str, max_mismatches: int):
        if self.codes is None:                   # pragma: no cover
            yield from _block_hits_slow(proto, self.seq, max_mismatches)
            return
        starts = _candidate_starts(proto, self.codes, len(proto),
                                   len(self.seq) - len(proto))
        yield from _verify(proto, self.s_arr, starts, len(proto), max_mismatches)


def _block_hits_slow(proto: str, seq: str, max_mismatches: int):
    """Pure-python fallback: k-mer index per block (also the test reference)."""
    index = build_kmer_index(seq, 5)
    L = len(proto)
    bl = L // 4
    rc_proto = reverse_complement(proto)
    seen_offsets: set[int] = set()
    for qseq, strand in ((proto, "+"), (rc_proto, "-")):
        blocks = [qseq[i * bl:(i + 1) * bl] for i in range(4)]
        candidates: set[int] = set()
        for b in blocks:
            for c in index.get(b, []):
                for off in range(0, L - bl + 1):
                    start = c - off
                    if 0 <= start <= len(seq) - L:
                        candidates.add(start)
        for start in candidates:
            if start in seen_offsets:
                continue
            w = seq[start:start + L]
            ref = proto if strand == "+" else rc_proto
            mm = sum(1 for a, b in zip(w, ref) if a != b)
            if mm <= max_mismatches:
                seen_offsets.add(start)
                yield start, strand, mm


def scan_genome_offtargets(protospacers: list[str], fasta_path: str,
                           block_bases: int = 5_000_000, max_mismatches: int = 3,
                           exclude: dict[str, tuple[str, int]] | None = None
                           ) -> list[OfftargetHit]:
    """Count every <=max_mismatches genomic near-match for each protospacer.

    `exclude` maps protospacer -> (chrom, position) of its own on-target site;
    that single locus (and its RC shadow) is dropped from the counts.
    Returns one row per hit so callers can aggregate, histogram, or map them.
    """
    exclude = exclude or {}
    global_seen: dict[str, set[tuple[str, int]]] = {p: set() for p in protospacers}
    hits: list[OfftargetHit] = []
    proto_set = set(protospacers)
    for chrom, gstart, seq in iter_fasta_blocks(fasta_path, block_bases):
        block = BlockIndex(seq)
        for proto in protospacers:
            L = len(proto)
            for off, strand, mm in block.hits(proto, max_mismatches):
                gpos = gstart + off
                key = (chrom, gpos)
                if key in global_seen[proto]:
                    continue
                excl = exclude.get(proto)
                if excl and excl[0] == chrom and abs(excl[1] - gpos) < L:
                    global_seen[proto].add(key)
                    continue
                global_seen[proto].add(key)
                hits.append(OfftargetHit(proto, chrom, gpos, strand, mm))
    return hits
