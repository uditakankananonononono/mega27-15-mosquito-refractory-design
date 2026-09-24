"""Genome-wide off-target scan for the v2 top-30 dsx gRNA candidates.

Live data (run outside CI): full AgamP5 reference (GCF_943734735.2,
idAnoGambNW_F1_1): chromosomes 2RL (NC_064601.1/OX030907.1, 118.2 Mb),
X (OX030909.1, 28.1 Mb), 3RL (OX030908.1, 99.1 Mb), MT (OX030910.2, 15.6 kb) -
fetched from NCBI efetch, see scripts/fetch_genome.sh. Output:
results/genome_offtargets_genomewide_top30.{csv,json}.

Coordinates: `position` is 0-based on the plus strand of the named accession.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from collections import Counter

from mosqdesign.genome_scan import scan_genome_offtargets

FASTA = sys.argv[1] if len(sys.argv) > 1 else "/tmp/agam_genome.fasta"
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
        "gw_mm0": c.get(0, 0), "gw_mm1": c.get(1, 0),
        "gw_mm2": c.get(2, 0), "gw_mm3": c.get(3, 0),
        "gw_total_le3": sum(c.values()),
        "locus110kb_le3": r["offtarget_hits_le3mm_110kb"],
    })
with open("results/genome_offtargets_genomewide_top30.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
    w.writeheader()
    w.writerows(out_rows)
with open("results/genome_offtargets_genomewide_top30.json", "w") as fh:
    json.dump({
        "reference": "AgamP5 GCF_943734735.2: NC_064601.1 (2RL) + OX030909.1 (X) + OX030908.1 (3RL) + OX030910.2 (MT)",
        "fetch": "scripts/fetch_genome.sh (NCBI efetch, live)",
        "guides": len(guides), "max_mismatches": 3,
        "block_bases": 20_000_000, "elapsed_seconds": round(elapsed, 2),
        "totals": {str(m): sum(c.get(m, 0) for c in per_guide.values()) for m in range(4)},
        "hits": [h.__dict__ for h in hits],
    }, fh, indent=1)
print(f"scanned {len(guides)} guides x ~246 Mb genome in {elapsed:.1f}s")
for r in out_rows[:10]:
    print(r)
