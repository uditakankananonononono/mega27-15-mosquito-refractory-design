"""Chromosome-scale off-target scan for the v2 top-30 dsx gRNA candidates.

Live data (run outside CI): AgamP5 chromosome 2 (NC_064601.1, ~118 Mb),
fetched from NCBI efetch - see scripts/fetch_genome.sh. The dsx locus sits on
2R inside this chromosome, so this scan covers the entire home chromosome of
every candidate. Output: results/genome_offtargets_chr2_top30.{csv,json}.

Coordinates: `position` is 0-based on the plus strand of NC_064601.1
(accession base = position + 1).
"""
from __future__ import annotations

import csv
import json
import sys
import time
from collections import Counter

from mosqdesign.genome_scan import scan_genome_offtargets

FASTA = sys.argv[1] if len(sys.argv) > 1 else "/tmp/chr2R.fasta"
TOP_N = 30

rows = list(csv.DictReader(open("results/ranked_designs_v2.csv")))[:TOP_N]
guides = [r["protospacer"] for r in rows]
# CSV genomic_pos is the 1-based accession coordinate; convert to 0-based.
exclude = {r["protospacer"]: ("NC_064601.1", int(r["genomic_pos"]) - 1) for r in rows}

t0 = time.time()
hits = scan_genome_offtargets(guides, FASTA, block_bases=20_000_000, max_mismatches=3, exclude=exclude)
elapsed = time.time() - t0

per_guide = {g: Counter() for g in guides}
for h in hits:
    per_guide[h.guide][h.mismatches] += 1

out_rows = []
for r in rows:
    c = per_guide[r["protospacer"]]
    out_rows.append({
        "rank_v2": r["rank_v2"], "protospacer": r["protospacer"],
        "genomic_pos": r["genomic_pos"],
        "chr2_mm0": c.get(0, 0), "chr2_mm1": c.get(1, 0),
        "chr2_mm2": c.get(2, 0), "chr2_mm3": c.get(3, 0),
        "chr2_total_le3": sum(c.values()),
        "locus110kb_le3": r["offtarget_hits_le3mm_110kb"],
    })
with open("results/genome_offtargets_chr2_top30.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
    w.writeheader()
    w.writerows(out_rows)
with open("results/genome_offtargets_chr2_top30.json", "w") as fh:
    json.dump({
        "reference": "NC_064601.1 (Anopheles gambiae chromosome 2, AgamP5/idAnoGambNW_F1_1)",
        "fetch": "scripts/fetch_genome.sh (NCBI efetch, live)",
        "guides": len(guides), "max_mismatches": 3,
        "block_bases": 20_000_000, "elapsed_seconds": round(elapsed, 2),
        "totals": {str(m): sum(c.get(m, 0) for c in per_guide.values()) for m in range(4)},
        "hits": [h.__dict__ for h in hits],
    }, fh, indent=1)
print(f"scanned {len(guides)} guides x ~118 Mb chromosome 2 in {elapsed:.1f}s")
for r in out_rows[:10]:
    print(r)
