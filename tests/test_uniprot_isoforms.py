"""Hermetic check of the archived UniProt dsx isoform evidence.

The drive strategy targets the female-specific dsx exon; the archived
UniProt records must document sex-specific isoforms in A. gambiae
(NCBI taxonomy 7165), corroborating the isoform biology the constraint
model relies on.
"""
import csv
import os


def _rows():
    path = os.path.join(os.path.dirname(__file__), "..", "data", "uniprot_dsx_anoga.tsv")
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_two_sex_specific_isoforms_archived():
    rows = _rows()
    assert len(rows) == 2
    names = {r["Protein names"] for r in rows}
    assert any("Female-specific doublesex" in n for n in names)
    assert any("Male-specific doublesex" in n for n in names)


def test_isoform_lengths_positive_and_distinct():
    rows = _rows()
    lengths = [int(r["Length"]) for r in rows]
    assert all(n > 100 for n in lengths)
    assert len(set(lengths)) == 2  # splice isoforms differ in length
