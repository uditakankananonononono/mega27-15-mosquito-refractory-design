"""V2 ranking: constraint-aware scoring (exon + splice-junction proximity bonus).

Motivation: the v1 drive-dynamics grid showed resistance generation e - controllable
through functional constraint of the target - is the binding constraint on drive
success, ahead of raw efficacy. V2 adds the a priori constraint bonus from
mosqdesign.annotation to the v1 composite.
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

from mosqdesign.annotation import parse_gene_table, exon_union, constraint_bonus
from mosqdesign.data.loaders import read_fasta, DATA_DIR
from mosqdesign.grna_scan import scan_ngg_sites, build_kmer_index, count_offtargets_fast
from mosqdesign.models.efficacy_cnn import EfficacyScorer, gc_content
from mosqdesign.grna_scan import homopolymer_run
from mosqdesign.rank import composite_score
from mosqdesign.run_analysis import KYROU_PROTOSPACER, KYROU_PAM, REGION_ACCESSION, GENE_BODY, RESULTS

REGION_OFFSET = 47600000


def run(top_n_offtarget: int = 300) -> dict:
    region = list(read_fasta(os.path.join(DATA_DIR, "dsx_locus_region.fasta")).values())[0]
    exons = exon_union(parse_gene_table(os.path.join(DATA_DIR, "dsx_gene_table.txt")))
    scorer = EfficacyScorer(os.path.join(os.path.dirname(__file__), "..", "assets",
                                         "doench2016_efficacy_cnn.pt"))

    all_sites = scan_ngg_sites(region, REGION_ACCESSION)
    eff_all = scorer.score([s.context30 for s in all_sites])

    seen: dict[tuple[str, str], int] = {}
    sites, effs = [], []
    for s, e in zip(all_sites, eff_all):
        if not (GENE_BODY[0] <= s.position <= GENE_BODY[1]):
            continue
        key = (s.protospacer, s.pam)
        if key in seen:
            if e > effs[seen[key]]:
                effs[seen[key]] = e
            continue
        seen[key] = len(sites)
        sites.append(s)
        effs.append(float(e))
    effs_arr = np.array(effs)

    # v2 composite for all sites (off-target term added for the top slice only)
    base = [composite_score(e, 0, gc_content(s.protospacer), homopolymer_run(s.protospacer))
            for s, e in zip(sites, effs_arr)]
    bonuses = [constraint_bonus(REGION_OFFSET + s.position, 23, exons) for s in sites]
    v2_pre = sorted(zip(sites, effs_arr, base, bonuses),
                    key=lambda t: t[2] + t[3][0], reverse=True)
    top = v2_pre[:top_n_offtarget]
    index = build_kmer_index(region, 5)
    ot = [count_offtargets_fast(t[0], [(REGION_ACCESSION, region, index)], max_mismatches=3)
          for t in top]
    final = sorted(zip(top, ot), key=lambda t: t[0][2] + t[0][3][0] - 0.5 * t[1], reverse=True)

    with open(os.path.join(RESULTS, "ranked_designs_v2.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank_v2", "strand", "genomic_pos", "protospacer", "pam",
                    "efficacy", "constraint_bonus", "constraint_label",
                    "offtarget_hits_le3mm_110kb", "score_v2"])
        for i, ((s, e, b, (bonus, label)), o) in enumerate(final, 1):
            w.writerow([i, s.strand, REGION_OFFSET + s.position, s.protospacer, s.pam,
                        f"{e:.4f}", f"{bonus:.2f}", label, o, f"{b + bonus - 0.5 * o:.4f}"])

    kyrou_v2 = None
    for i, ((s, e, b, (bonus, label)), o) in enumerate(final, 1):
        if s.protospacer == KYROU_PROTOSPACER and s.pam == KYROU_PAM:
            kyrou_v2 = {"rank_v2": i, "efficacy": e, "constraint_bonus": bonus,
                        "constraint_label": label, "offtargets": o,
                        "score_v2": b + bonus - 0.5 * o}
            break

    summary = {
        "n_sites": len(sites),
        "n_ranked": len(final),
        "kyrou_v2": kyrou_v2,
        "top10_v2": [{"protospacer": t[0][0].protospacer, "pam": t[0][0].pam,
                      "strand": t[0][0].strand, "genomic_pos": REGION_OFFSET + t[0][0].position,
                      "efficacy": t[0][1], "constraint": t[0][3][1],
                      "offtargets": t[1],
                      "score_v2": t[0][2] + t[0][3][0] - 0.5 * t[1]} for t in final[:10]],
    }
    with open(os.path.join(RESULTS, "analysis_v2_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary
