"""Independent re-derivation of the three AF-scan windows with Biopython.

The AF analysis (scripts/ag1000g_allele_freq.py) hard-codes per-guide
constants (ref_minus, PAM coordinates) that were originally derived by
manual slicing. This script re-derives every constant from the archived
AgamP4 region FASTA using Biopython's own parser and reverse-complement,
so a transcription error in the constants cannot survive silently.

Hermetic: reads only data/agamp4/AgamP4_2R_dsx_region.fasta.
Writes results/window_verification_biopython.json.
"""
from __future__ import annotations

import json
import os

from Bio import SeqIO
from Bio.Seq import Seq

import importlib.util

_spec = importlib.util.spec_from_file_location(
    'ag1000g_allele_freq', os.path.join(os.path.dirname(__file__), 'ag1000g_allele_freq.py'))
_agaf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_agaf)
GUIDES = _agaf.GUIDES

FASTA = "data/agamp4/AgamP4_2R_dsx_region.fasta"
OUT = "results/window_verification_biopython.json"
REGION_START = 48700000  # 1-based first base of the archived region


def derive(fasta_path: str = FASTA) -> dict:
    """Re-derive each guide's plus-strand protospacer, ref_minus and PAM."""
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    seq = str(rec.seq).upper()

    def window(a: int, b: int) -> str:  # 1-based inclusive genomic -> seq
        return seq[a - REGION_START : b - REGION_START + 1]

    report = {}
    for name, g in GUIDES.items():
        pa, pb = g["proto"]
        xa, xb = g["pam"]
        proto_plus = window(pa, pb)
        ref_minus = str(Seq(proto_plus).reverse_complement())
        pam_minus = str(Seq(window(xa, xb)).reverse_complement())
        expected_pam = {"dsx-v3-1": "GGG", "dsx-v3-2": "TGG", "kyrou": "TGG"}[name]
        checks = {
            "proto_plus": proto_plus,
            "ref_minus_derived": ref_minus,
            "ref_minus_constant": g["ref_minus"],
            "ref_minus_match": ref_minus == g["ref_minus"],
            "pam_minus_derived": pam_minus,
            "pam_expected": expected_pam,
            "pam_match": pam_minus == expected_pam,
        }
        checks["pass"] = checks["ref_minus_match"] and checks["pam_match"]
        report[name] = checks
    return report


def main() -> None:
    report = derive()
    os.makedirs("results", exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=2)
    for name, r in report.items():
        print(f"{name}: {'PASS' if r['pass'] else 'FAIL'} "
              f"(ref_minus {r['ref_minus_derived']}, PAM {r['pam_minus_derived']})")
    if not all(r["pass"] for r in report.values()):
        raise SystemExit("window verification failed")


if __name__ == "__main__":
    main()
