import random
import textwrap

from mosqdesign.genome_scan import iter_fasta_blocks, scan_genome_offtargets
from mosqdesign.grna_scan import (GrnaSite, count_offtargets, reverse_complement)


def _write_fasta(path, name, seq, wrap=60):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        fh.write("\n".join(textwrap.wrap(seq, wrap)) + "\n")


def _random_seq(n, seed=7):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_iter_fasta_blocks_covers_every_window():
    seq = _random_seq(1000)
    _write_fasta("/tmp/gs_blocks.fa", "chrT", seq)
    blocks = list(iter_fasta_blocks("/tmp/gs_blocks.fa", block_bases=120, overlap=19))
    assert blocks[0][0] == "chrT"
    # every 20-mer start position must lie inside some block
    covered = set()
    for _, start, s in blocks:
        for i in range(len(s) - 19):
            covered.add(start + i)
    assert covered == set(range(len(seq) - 19))


def test_iter_fasta_blocks_flushes_final_partial_block():
    seq = _random_seq(250)
    _write_fasta("/tmp/gs_flush.fa", "chrT", seq)
    blocks = list(iter_fasta_blocks("/tmp/gs_flush.fa", block_bases=100, overlap=19))
    joined = "".join(s if i == 0 else s[19:] for i, (_, _, s) in enumerate(blocks))
    assert joined == seq


def test_scan_matches_naive_counter_with_boundary_hit():
    rng = random.Random(11)
    chrom = _random_seq(600)
    proto = "ATGCCGTAACGTTAGCCTGA"
    # plant the on-target at 50 and a 2-mismatch near-match straddling the block edge
    near = list(proto)
    near[3] = {"A": "C", "C": "A", "G": "A", "T": "A"}[near[3]]
    near[17] = {"A": "T", "C": "T", "G": "T", "T": "G"}[near[17]]
    near = "".join(near)
    assert sum(1 for a, b in zip(near, proto) if a != b) == 2
    chrom = chrom[:50] + proto + chrom[70:]
    chrom = chrom[:195] + near + chrom[215:]
    # also plant the reverse complement of the guide (minus-strand hit)
    chrom = chrom[:400] + reverse_complement(proto) + chrom[420:]
    _write_fasta("/tmp/gs_naive.fa", "chrT", chrom)

    hits = scan_genome_offtargets([proto], "/tmp/gs_naive.fa",
                                  block_bases=200, max_mismatches=3,
                                  exclude={proto: ("chrT", 50)})
    positions = sorted(h.position for h in hits)
    # naive whole-sequence counter as ground truth
    site = GrnaSite("chrT", "+", 50, proto, "TGG", "N" * 30)
    naive = count_offtargets(site, [("chrT", chrom)], max_mismatches=3)
    assert len(hits) == naive == 2          # near-match at 195 and RC at 400; on-target excluded
    assert 195 in positions and 400 in positions
    strands = {h.position: h.strand for h in hits}
    assert strands[400] == "-"              # found via the reverse-complement query


def test_no_double_count_across_overlap():
    chrom = _random_seq(500)
    proto = "GGACTTGCACGATCGTACCA"
    chrom = chrom[:208] + proto + chrom[228:]   # inside every overlap zone of block 200
    _write_fasta("/tmp/gs_dup.fa", "chrT", chrom)
    hits = scan_genome_offtargets([proto], "/tmp/gs_dup.fa",
                                  block_bases=200, max_mismatches=0,
                                  exclude={proto: ("chrT", 208)})
    assert hits == []                        # on-target excluded exactly once, no duplicates
    hits2 = scan_genome_offtargets([proto], "/tmp/gs_dup.fa",
                                   block_bases=200, max_mismatches=0)
    assert len(hits2) == 1 and hits2[0].position == 208


def test_multiple_guides_independent_counts():
    chrom = _random_seq(400)
    g1 = "ACGTCAGTTGGCACTAGGCA"
    g2 = "TTGCGACCATGGAACCTGTA"
    chrom = chrom[:100] + g1 + chrom[120:]
    chrom = chrom[:300] + g2 + chrom[320:]
    _write_fasta("/tmp/gs_multi.fa", "chrT", chrom)
    hits = scan_genome_offtargets([g1, g2], "/tmp/gs_multi.fa",
                                  block_bases=150, max_mismatches=0)
    by_guide = {}
    for h in hits:
        by_guide.setdefault(h.guide, []).append(h.position)
    assert by_guide == {g1: [100], g2: [300]}


def test_numpy_path_matches_pure_python_path():
    import random as _r
    from mosqdesign.genome_scan import _block_hits, _block_hits_slow
    rng = _r.Random(23)
    seq = "".join(rng.choice("ACGT") for _ in range(3000))
    proto = "CAGTGCATCGATCGGTTACC"
    seq = seq[:1500] + proto + seq[1520:]          # exact match
    mut = proto[:7] + ("A" if proto[7] != "A" else "C") + proto[8:]  # 1 mm
    seq = seq[:200] + mut + seq[220:]
    fast = sorted(_block_hits(proto, seq, 2))
    slow = sorted(_block_hits_slow(proto, seq, 2))
    assert fast == slow
    assert any(pos == 1500 and mm == 0 for pos, _, mm in fast)
    assert any(pos == 200 and mm == 1 for pos, _, mm in fast)


def test_n_bases_never_match():
    from mosqdesign.genome_scan import _block_hits
    proto = "ACGTACGTACGTACGTACGT"
    seq = "N" * 100 + proto + "N" * 100
    hits = list(_block_hits(proto, seq, 3))
    assert [h[0] for h in hits] == [100]
    seq_n = "N" * 100 + proto[:10] + "N" + proto[10:] + "N" * 100
    assert list(_block_hits(proto, seq_n, 0)) == []


def test_hit_at_sequence_edge():
    from mosqdesign.genome_scan import _block_hits
    proto = "TTGACCGGTACCGTAACGGA"
    seq = "GGGG" + proto                      # exact match starting at 4, ends at the very edge
    hits = list(_block_hits(proto, seq, 0))
    assert (4, "+", 0) in hits
    seq2 = proto + "TTTT"                     # match starting at 0
    assert (0, "+", 0) in list(_block_hits(proto, seq2, 0))
