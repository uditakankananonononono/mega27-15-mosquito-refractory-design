import os
import pytest

from mosqdesign.annotation import (parse_gene_table, exon_at,
                                   dist_to_splice_junction, constraint_bonus)

TABLE = os.path.join(os.path.dirname(__file__), "..", "data", "dsx_gene_table.txt")


@pytest.mark.skipif(not os.path.exists(TABLE), reason="gene table not fetched")
def test_gene_table_parses_all_variant_tables():
    tables = parse_gene_table(TABLE)
    assert len(tables) == 9                     # 9 RefSeq transcript variants
    assert len(tables["XM_061649168.1"]) == 8   # representative: 8 exons
    lengths = [e.length for e in tables["XM_061649168.1"]]
    assert 135 in lengths  # the small exon containing the validated Kyrou site


@pytest.mark.skipif(not os.path.exists(TABLE), reason="gene table not fetched")
def test_kyrou_site_annotation_discrepancy_documented():
    """The validated Kyrou site is INTRONIC in the current RefSeq model (~0.5 kb
    from the nearest junction), although the paper describes it as spanning the
    intron4-exon5 boundary under the older AgamP4/VectorBase annotation. This test
    pins the discrepancy so it is reported, not smoothed over."""
    from mosqdesign.annotation import exon_union, dist_to_nearest_exon_edge
    exons = exon_union(parse_gene_table(TABLE))
    kyrou_start = 47622174  # genomic start of the published gRNA site (plus-strand coords)
    assert exon_at(kyrou_start, exons) is None  # intronic in current RefSeq model
    d = dist_to_nearest_exon_edge(kyrou_start, exons)
    assert 100 < d < 1000  # near a junction but not on it, in this annotation


def test_constraint_bonus_logic():
    exons = [type("E", (), {"start": 1000, "end": 1100, "length": 101})]
    from mosqdesign.annotation import Exon
    exons = [Exon(1000, 1100, 101)]
    b_far, lab_far = constraint_bonus(1040, 23, exons)          # mid-exon
    b_junc, lab_junc = constraint_bonus(1090, 23, exons)        # crosses junction region
    b_out, lab_out = constraint_bonus(2000, 23, exons)          # outside
    assert b_junc > b_far > b_out == 0.0
    assert "splice_proximal" in lab_junc and "exonic" in lab_far
