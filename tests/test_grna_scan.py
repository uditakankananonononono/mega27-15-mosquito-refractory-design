from mosqdesign.grna_scan import (scan_ngg_sites, reverse_complement,
                                  homopolymer_run, count_offtargets, GrnaSite)


def test_finds_plus_strand_site():
    seq = "AAAA" + "ACGTACGTACGTACGTACGT" + "TGG" + "AAAA"
    sites = scan_ngg_sites(seq, "g")
    plus = [s for s in sites if s.strand == "+" and s.pam == "TGG"]
    assert len(plus) == 1
    s = plus[0]
    assert s.protospacer == "ACGTACGTACGTACGTACGT"
    assert s.position == 4
    assert len(s.context30) == 30


def test_finds_minus_strand_site():
    proto = "ACGTACGTACGTACGTACGT"
    pam = "TGG"
    target = reverse_complement(proto + pam)
    seq = "AAAA" + target + "AAAA"
    sites = scan_ngg_sites(seq, "g")
    minus = [s for s in sites if s.strand == "-"]
    assert len(minus) == 1
    assert minus[0].protospacer == proto
    assert minus[0].pam == pam
    assert len(minus[0].context30) == 30


def test_context_padding_at_edges():
    seq = "ACGTACGTACGTACGTACGTTGG"
    (s,) = [x for x in scan_ngg_sites(seq, "g") if x.strand == "+"]
    assert len(s.context30) == 30
    assert s.context30[:4] == "NNNN"


def test_no_ambiguous_bases_in_sites():
    seq = "ACGTNACGTACGTACGTACGTACGTTGG"
    for s in scan_ngg_sites(seq, "g"):
        assert all(b in "ACGT" for b in s.protospacer + s.pam)


def test_homopolymer_run():
    assert homopolymer_run("AAATG") == 3
    assert homopolymer_run("ACGT") == 1


def test_count_offtargets_excludes_self():
    proto = "ATGCCGTAACGTTAGCCTGA"  # not a reverse-complement palindrome
    site = GrnaSite("g1", "+", 4, proto, "TGG", "N" * 30)
    panel = [("g1", "AAAA" + proto + "TGG" + "AAAA"),
             ("g2", "TTTT" + "ATGCCGTAACGTTAGCCTAA" + "TGG")]  # 1 mismatch in g2
    assert count_offtargets(site, panel, max_mismatches=1) == 1
