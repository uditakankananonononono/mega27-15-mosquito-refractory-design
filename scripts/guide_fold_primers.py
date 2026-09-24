"""Guide-RNA folding cross-check (ViennaRNA) and validation-amplicon
primer design (Primer3) for the three named dsx candidates.

Folding rationale: stable secondary structure within the spacer competes
with R-loop formation and is associated with reduced Cas9 activity in
published guide-design studies; a spacer-only MFE is a heuristic screen,
not an efficacy model, and is reported as such.

Primer rationale: section 4.6 stages genomic-DNA amplicon validation of
each candidate window; Primer3 designs the amplicon primer pairs against
the archived AgamP4 region so the assay targets exactly the scanned
window.

Hermetic: archived FASTA only, no network.
Writes results/guide_fold_primers.{json,md}.
"""
from __future__ import annotations

import importlib.util
import json
import os

import primer3
import RNA
from Bio import SeqIO

_spec = importlib.util.spec_from_file_location(
    "ag1000g_allele_freq", os.path.join(os.path.dirname(__file__), "ag1000g_allele_freq.py"))
_agaf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_agaf)
GUIDES = _agaf.GUIDES

FASTA = "data/agamp4/AgamP4_2R_dsx_region.fasta"
REGION_START = 48700000
HAIRPIN_MFE_THRESHOLD = -3.0  # kcal/mol; spacer-only heuristic screen


def fold_spacer(guide_dna: str) -> tuple[str, float]:
    """ViennaRNA MFE structure of the 20-nt spacer (DNA -> RNA)."""
    structure, mfe = RNA.fold(guide_dna.replace("T", "U"))
    return structure, float(mfe)


def design_amplicon(name: str, fasta_path: str = FASTA) -> dict:
    """Primer3 pair flanking the guide window on the AgamP4 region."""
    g = GUIDES[name]
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    template = str(rec.seq).upper()
    lo = min(g["pam"][0], g["proto"][0])
    hi = max(g["pam"][1], g["proto"][1])
    res = primer3.bindings.design_primers(
        {
            "SEQUENCE_ID": name,
            "SEQUENCE_TEMPLATE": template,
            "SEQUENCE_TARGET": [lo - REGION_START, hi - lo + 1],
        },
        {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_NUM_RETURN": 1,
            "PRIMER_PRODUCT_SIZE_RANGE": [[250, 450]],
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_OPT_TM": 60.0,
        },
    )
    if res.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        return {"name": name, "designed": False, "explain": res.get("PRIMER_LEFT_EXPLAIN", "")}
    return {
        "name": name,
        "designed": True,
        "left": res["PRIMER_LEFT_0_SEQUENCE"],
        "right": res["PRIMER_RIGHT_0_SEQUENCE"],
        "left_tm": round(res["PRIMER_LEFT_0_TM"], 2),
        "right_tm": round(res["PRIMER_RIGHT_0_TM"], 2),
        "product_size": res["PRIMER_PAIR_0_PRODUCT_SIZE"],
        "amplicon_spans": [lo, hi],
    }


def run() -> dict:
    report = {}
    for name, g in GUIDES.items():
        structure, mfe = fold_spacer(g["guide"])
        report[name] = {
            "spacer": g["guide"],
            "mfe_structure": structure,
            "mfe_kcal_mol": round(mfe, 2),
            "hairpin_flag": mfe <= HAIRPIN_MFE_THRESHOLD,
            "amplicon": design_amplicon(name),
        }
    return report


def main() -> None:
    report = run()
    os.makedirs("results", exist_ok=True)
    with open("results/guide_fold_primers.json", "w") as fh:
        json.dump(report, fh, indent=2)
    lines = [
        "| candidate | spacer MFE (kcal/mol) | structure | hairpin flag | amplicon (bp) | left primer | right primer |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, r in report.items():
        a = r["amplicon"]
        lines.append(
            f"| {name} | {r['mfe_kcal_mol']} | {r['mfe_structure']} | "
            f"{'yes' if r['hairpin_flag'] else 'no'} | "
            f"{a.get('product_size', '-')} | {a.get('left', '-')} | {a.get('right', '-')} |"
        )
    with open("results/guide_fold_primers.md", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
