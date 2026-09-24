"""Hermetic tests for scripts/cross_assembly_verify.py (synthetic genome)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from cross_assembly_verify import (  # noqa: E402
    blocks, find_all, hamming, rc, search_unit)


def test_rc():
    assert rc("AACG") == "CGTT"


def test_blocks_cover_23():
    bs = blocks("A" * 23)
    assert bs[0][0] == 0 and bs[-1][1] == 23
    # contiguous coverage
    assert [b for a, b in bs[:-1]] == [a for a, b in bs[1:]]


def test_hamming():
    assert hamming("ACGT", "ACGA") == 1
    assert hamming("ACGT", "ACGT") == 0


def test_find_all_overlapping():
    assert find_all("AAAA", "AA") == [0, 1, 2]


def test_search_exact_and_near():
    unit = "AGTCCTGAACTGGCATTCGATGC"  # 23, non-periodic
    near = unit[:-1] + ("A" if unit[-1] != "A" else "C")
    genome = "TTTTTT" + unit + "CCCCC" + near + "GGGGGG"
    res = search_unit(genome, unit)
    assert res["exact"] == [6]
    # the second copy has one substitution at its last base
    assert res["near"][1] == [6 + len(unit) + 5]
    assert res["near"][2] == [] and res["near"][3] == []


def test_search_no_hits():
    res = search_unit("A" * 100, "C" * 23)
    assert res["exact"] == [] and all(res["near"][k] == [] for k in (1, 2, 3))


def test_plus_units_orientation():
    # pinned to the twice-verified window constants: minus-strand kyrou
    # guide GTTTAACACAGGTCAAGCGG + TGG PAM => plus 5'->3' unit is
    # rc(minus unit) = CCA PAM first, then protospacer.
    from cross_assembly_verify import plus_units
    u = plus_units()
    assert u["kyrou"] == "CCACCGCTTGACCTGTGTTAAAC"
    assert all(len(v) == 23 for v in u.values())


def test_guide_units_orientation():
    from cross_assembly_verify import guide_units
    u = guide_units()
    # rc(guide + GGG) for v3-1: CCA PAM-plus first, then proto
    assert u["dsx-v3-1"] == "CCCACCCTAACGCATACTGCCCA"
    assert u["kyrou"] == "CCACCGCTTGACCTGTGTTAAAC"  # guide == AgamP4 ref
    assert all(len(v) == 23 for v in u.values())
