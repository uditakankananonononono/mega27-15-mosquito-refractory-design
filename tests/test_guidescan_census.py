"""Offline fixture tests for the GuideScan census check (tool 43).

Exercises classify_hit orientation/PAM logic, resolve_own_sites, and
compare on synthetic genomes - no genome download, no guidescan binary.
"""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "guidescan_census", Path(__file__).parent.parent / "scripts" / "guidescan_census.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

# non-palindromic 20-mer (its reverse complement differs), so orientation
# is decidable in fixtures
G = "AACCGGTTAACCGGTTAACC"


def test_classify_forward_ngg():
    # genome[16:36] == G; 3' flank genome[36:39] = 'AGG'? build explicitly
    gen = {"c": ("A" * 10 + G + "AGG" + "T" * 20)}
    assert m.classify_hit(gen, {"g1": G}, "c", 10, 30, "g1") == "ngg"


def test_classify_forward_pamless():
    gen = {"c": ("A" * 10 + G + "AAA" + "T" * 20)}
    assert m.classify_hit(gen, {"g1": G}, "c", 10, 30, "g1") == "pamless"


def test_classify_reverse_ngg():
    rg = m.rc(G)
    # protospacer rc at 10..30; guide strand PAM = rc(genome[7:10])
    gen = {"c": ("AAC" + rg[:0] + "A" * 7 + rg + "T" * 20)}
    # place rc of NGG ('CCN') upstream: genome[7:10] = 'CCN' -> rc = 'NGG'
    gen = {"c": ("A" * 7 + "CCT" + rg + "T" * 20)}
    assert m.classify_hit(gen, {"g1": G}, "c", 10, 30, "g1") == "ngg"


def test_classify_ambiguous_tie():
    # a guide whose forward and reverse-complement Hamming distances tie
    gen = {"c": "A" * 40}
    g = "A" * 10 + "T" * 10
    assert m.classify_hit(gen, {"g1": g}, "c", 0, 20, "g1") == "ambiguous"


def test_resolve_own_sites_and_kmers(tmp_path):
    gen = {"c": ("N" * 5 + G + "AGG" + "N" * 30)}
    guides = {"g1": G}
    bed = tmp_path / "hits.bed"
    bed.write_text("c\t5\t25\tg1\t0\n")
    own, problems = m.resolve_own_sites(gen, guides, str(bed))
    assert not problems and own["g1"] == ("c", 28, "+")
    kpath = tmp_path / "kmers.csv"
    m.write_kmers(own, guides, str(kpath))
    lines = kpath.read_text().splitlines()
    assert lines[0] == "id,sequence,pam,chromosome,position,sense"
    assert lines[1] == f"g1,{G},NGG,c,28,+"


def test_compare_exact_and_discordant():
    import collections
    guides = {"g1": "X" * 20, "g2": "Y" * 20}
    gs = collections.defaultdict(lambda: collections.defaultdict(int))
    cen = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"ngg": 0, "pamless": 0}))
    gs["g1"][0] = 1; cen["g1"][0]["ngg"] = 1
    gs["g2"][3] = 5; cen["g2"][3]["ngg"] = 4  # one discordant cell
    disc = m.compare(gs, cen, guides)
    assert len(disc) == 1 and disc[0]["guide"] == "g2" and disc[0]["tier"] == 3
