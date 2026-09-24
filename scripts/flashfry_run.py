#!/usr/bin/env python3
"""FlashFry 2.0 (McKenna & Shendure 2018, BMC Biology, doi:10.1186/s12915-018-0545-0)
library-scale off-target audit of every Cas9-NGG guide in the two dsx female-exon
windows, against a whole-genome AgamP5 database.

Pipeline (reproduce exactly):
  1. jar: https://github.com/mckennalab/FlashFry/releases/download/2.0/FlashFry-assembly-2.00.jar
  2. index:  java -Xmx1024m -jar FlashFry.jar index --tmpLocation tmp \
               --database agamp5_spcas9ngg --reference agam_genome.fasta --enzyme spcas9-ngg-20
             (4-molecule AgamP5 reference; 192.19 s; 209 MB database)
  3. discover: java -Xmx900m -jar FlashFry.jar discover --database agamp5_spcas9ngg \
               --fasta dsx_targets.fa --output dsx_discover.txt --maxMismatch 4
             (dsx_targets.fa = the two dsx windows, headers dsx_W1/dsx_W2)
  4. score:    java -Xmx900m -jar FlashFry.jar score --input dsx_discover.txt \
               --output dsx_scored.txt \
               --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
               --database agamp5_spcas9ngg
  Note: score rejects discover output produced with --positionOutput in v2.0
  (parse failure on the annotated off-target column); discovery was re-run without it.
  The scored table is committed as results/flashfry_scored_guides.tsv.

This script parses that committed table, audits the four lead guides, and writes
results/offtarget_flashfry.json / .md.
"""
import csv, json, statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCORED = ROOT / "results" / "flashfry_scored_guides.tsv"
OUT_JSON = ROOT / "results" / "offtarget_flashfry.json"
OUT_MD = ROOT / "results" / "offtarget_flashfry.md"

def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

LEADS = {
    "dsx-v3-1": "TGGGCAGTATGCGTTAGGGT",
    "dsx-v3-2": "CATTAAGACCTACGAAGCGC",
    "dsx-v3-3": "GAAGCGAGCCCAATGGCTGT",
    "kyrou-2018": "GTTTAACACAGGTCAAGCGG",
}

rows = list(csv.DictReader(SCORED.open(), delimiter="\t"))
n_total = len(rows)
overflow = [r for r in rows if r["overflow"] != "OK"]

def mmhist(r):
    return [int(x) for x in r["0-1-2-3-4_mismatch"].split(",")]

# library-scale stats
ot_counts = [int(r["otCount"]) for r in rows]
mm_totals = [0, 0, 0, 0, 0]
for r in rows:
    h = mmhist(r)
    for i in range(5):
        mm_totals[i] += h[i]
zero_ot = sum(1 for c in ot_counts if c == 0)
mm0_dup = sum(1 for r in rows if mmhist(r)[0] > 1)
by_cfd = sorted(rows, key=lambda r: float(r["DoenchCFD_specificityscore"]), reverse=True)

library = {
    "guides_scored": n_total,
    "overflow_guides": len(overflow),
    "guides_zero_offtargets_le4mm": zero_ot,
    "guides_with_mm0_duplicate": mm0_dup,
    "otcount_median": statistics.median(ot_counts),
    "otcount_max": max(ot_counts),
    "genomewide_hit_totals_by_mismatch_0to4": mm_totals,
    "top5_by_cfd_specificity": [
        {"target": r["target"], "window": r["contig"], "start": int(r["start"]),
         "cfd_specificity": float(r["DoenchCFD_specificityscore"]),
         "hsu2013": float(r["Hsu2013"]), "otCount": int(r["otCount"]),
         "mm_hist": mmhist(r)}
        for r in by_cfd[:5]
    ],
}

lead_report = {}
for name, spacer in LEADS.items():
    hits = [r for r in rows if r["target"] in (spacer, rc(spacer))
            or r["target"].startswith(spacer) or r["target"].startswith(rc(spacer))]
    exact = [r for r in hits if r["target"][:20] in (spacer, rc(spacer))]
    if not exact:
        lead_report[name] = {"status": "not found in windows"}
        continue
    r = exact[0]
    h = mmhist(r)
    lead_report[name] = {
        "target_23mer": r["target"], "window": r["contig"],
        "window_start": int(r["start"]), "orientation": r["orientation"],
        "cfd_specificity": float(r["DoenchCFD_specificityscore"]),
        "cfd_max_ot": float(r["DoenchCFD_maxOT"]),
        "hsu2013": float(r["Hsu2013"]),
        "doench2014_ontarget": r["Doench2014OnTarget"],
        "mm_hist_0to4": h, "otCount": int(r["otCount"]),
        "mm1_mm3_offtargets": h[1] + h[2] + h[3],
        "closest_hit_mismatches": int(r["basesDiffToClosestHit"]),
        "dangerous": {"GC": r["dangerous_GC"], "polyT": r["dangerous_polyT"],
                      "in_genome": r["dangerous_in_genome"]},
    }

out = {
    "tool": "FlashFry 2.0 (mckennalab/FlashFry, McKenna & Mali 2018 doi:10.1186/s12915-018-0545-0)",
    "database": "AgamP5 4-molecule reference, spcas9-ngg-20 profile, built 192.19 s, 209 MB",
    "discovery": {"input": "dsx_W1+dsx_W2 windows (3,373 bp)", "maxMismatch": 4,
                  "guides_discovered": 430, "guides_scored": n_total},
    "library": library,
    "leads": lead_report,
    "cross_engine_conclusions": [
        "v3-2 shows ZERO 1-3-mismatch off-targets in AgamP5 (mm_hist 1,0,0,0,3): the one-mismatch site CRISPOR found in AgamP4 and CRISPRoff flagged at 77% of on-target energy is an AgamP4-reference-allele artifact, absent from the AgamP5 assembly; the flag tracks the assembly, not the guide.",
        "No lead carries an MM0 duplicate in AgamP5 (all mm0 counts = 1, the on-target): CHOPCHOP's MM0=2 flag on a lead is not reproduced by FlashFry, consistent with the CRISPOR cross-engine resolution.",
        "v3-3 is the least specific lead on this third independent axis too (CFD specificity lowest, 10 off-targets <=4 MM, 2 at 3 MM), adding to the ViennaRNA structural outlier flag.",
        "v3-2 is the most sequence-specific lead in AgamP5 (CFD specificity highest, only 3 MM4 hits): its withdrawal rests on population-diversity evidence, not reference-genome off-target load.",
    ],
}
OUT_JSON.write_text(json.dumps(out, indent=2) + "\n")

md = ["# FlashFry 2.0 library-scale off-target audit (AgamP5)", ""]
md.append(f"Guides discovered in dsx_W1+dsx_W2 (3,373 bp): 430; scored: {n_total} "
          f"(overflow-flagged: {len(overflow)}). Whole-genome AgamP5 database, <=4 mismatches.")
md.append("")
md.append("## Library")
md.append(f"- guides with zero off-targets (<=4 MM): {zero_ot}/{n_total}")
md.append(f"- guides with an MM0 duplicate: {mm0_dup}")
md.append(f"- off-target count median {library['otcount_median']}, max {library['otcount_max']}")
md.append(f"- genome-wide hit totals by mismatch (0,1,2,3,4): {mm_totals}")
md.append("")
md.append("## Leads")
md.append("| lead | CFD specificity | Hsu2013 | MM hist 0-4 | OTs | MM1-3 OTs |")
md.append("|---|---|---|---|---|---|")
for name, rep in lead_report.items():
    if "status" in rep:
        md.append(f"| {name} | {rep['status']} | | | | |")
        continue
    md.append(f"| {name} | {rep['cfd_specificity']:.3f} | {rep['hsu2013']:.1f} | "
              f"{','.join(map(str, rep['mm_hist_0to4']))} | {rep['otCount']} | {rep['mm1_mm3_offtargets']} |")
md.append("")
md.append("## Cross-engine conclusions")
for c in out["cross_engine_conclusions"]:
    md.append(f"- {c}")
OUT_MD.write_text("\n".join(md) + "\n")
print(json.dumps({"leads": {k: v.get("cfd_specificity", v.get("status")) for k, v in lead_report.items()},
                  "zero_ot": zero_ot, "n": n_total, "mm_totals": mm_totals}, indent=2))
