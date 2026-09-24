"""End-to-end item-15 analysis on the dsx locus region (NC_064601.1:47600000-47710000).

Pipeline: scan 110kb region for NGG sites -> CNN efficacy on all -> restrict to the
dsx gene body -> dedupe -> pre-rank -> fast pigeonhole off-target count (region panel)
-> final rank. Positive control: published Kyrou 2018 dsx gRNA. Drive simulations.
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

from mosqdesign.data.loaders import read_fasta, DATA_DIR
from mosqdesign.grna_scan import (scan_ngg_sites, build_kmer_index,
                                  count_offtargets_fast, reverse_complement)
from mosqdesign.models.efficacy_cnn import EfficacyScorer
from mosqdesign.rank import rank_sites
from mosqdesign.drive_sim import simulate, wright_fisher

KYROU_PROTOSPACER = "GTTTAACACAGGTCAAGCGG"  # Kyrou 2018 Nat Biotechnol nbt.4245
KYROU_PAM = "TGG"
REGION_ACCESSION = "NC_064601.1:47600000-47710000"
GENE_BODY = (10877, 100920)  # dsx gene body offsets within the region (47,610,877..47,700,920)

RESULTS = os.path.join(os.path.dirname(__file__), "..", "results")


def run(top_n_offtarget: int = 150) -> dict:
    os.makedirs(RESULTS, exist_ok=True)
    region = list(read_fasta(os.path.join(DATA_DIR, "dsx_locus_region.fasta")).values())[0]
    scorer = EfficacyScorer(os.path.join(os.path.dirname(__file__), "..", "assets",
                                         "doench2016_efficacy_cnn.pt"))

    all_sites = scan_ngg_sites(region, REGION_ACCESSION)
    eff_all = scorer.score([s.context30 for s in all_sites])

    # restrict to dsx gene body and dedupe by (protospacer, pam)
    seen: dict[tuple[str, str], int] = {}
    sites, effs = [], []
    for s, e in zip(all_sites, eff_all):
        if not (GENE_BODY[0] <= s.position <= GENE_BODY[1]):
            continue
        key = (s.protospacer, s.pam)
        if key in seen:
            if e > effs[seen[key]]:
                effs[seen[key]] = e  # keep best-context duplicate
            continue
        seen[key] = len(sites)
        sites.append(s)
        effs.append(float(e))
    effs_arr = np.array(effs)

    pre = rank_sites(sites, effs_arr, [0] * len(sites))
    top = pre[:top_n_offtarget]
    index = build_kmer_index(region, 5)
    ot = [count_offtargets_fast(r.site, [(REGION_ACCESSION, region, index)], max_mismatches=3)
          for r in top]
    final = rank_sites([r.site for r in top], np.array([r.efficacy for r in top]), ot)

    with open(os.path.join(RESULTS, "ranked_designs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "strand", "region_offset", "genomic_pos", "protospacer", "pam",
                    "efficacy", "offtarget_hits_le3mm_110kb", "gc", "homopolymer", "score"])
        for i, r in enumerate(final, 1):
            w.writerow([i, r.site.strand, r.site.position, 47600000 + r.site.position,
                        r.site.protospacer, r.site.pam, f"{r.efficacy:.4f}",
                        r.offtarget_hits, f"{r.gc:.3f}", r.homopolymer, f"{r.score:.4f}"])

    # positive control: published Kyrou 2018 dsx gRNA must be recovered and characterized
    kyrou_rank = None
    for i, r in enumerate(final, 1):
        if r.site.protospacer == KYROU_PROTOSPACER and r.site.pam == KYROU_PAM:
            kyrou_rank = {"rank_in_top150": i, "efficacy": r.efficacy,
                          "offtargets": r.offtarget_hits, "score": r.score,
                          "strand": r.site.strand, "region_offset": r.site.position}
            break
    kyrou_present = kyrou_rank is not None or any(
        s.protospacer == KYROU_PROTOSPACER for s in sites)

    grid = []
    for h in (0.80, 0.90, 0.95, 0.99):
        for e in (0.001, 0.01, 0.1, 0.5):
            traj = simulate(0.99, 0.01, 0.0, h=h, e=e, c_hom=1.0, generations=60)
            gen95 = next((g for g in range(61) if traj[g, 1] >= 0.95), None)
            grid.append({"h": h, "e": e, "peak_drive_freq": float(traj[:, 1].max()),
                         "final_resistance_freq": float(traj[-1, 2]),
                         "gen_reaches_95pct_drive": gen95})
    wf = wright_fisher(10000, 0.99, 0.01, 0.0, h=0.99, e=0.01, c_hom=1.0,
                       generations=40, seed=1)
    np.save(os.path.join(RESULTS, "wright_fisher_traj.npy"), wf)

    summary = {
        "region": REGION_ACCESSION,
        "region_len": len(region),
        "n_ngg_sites_region": len(all_sites),
        "n_ngg_sites_gene_body_dedup": len(sites),
        "efficacy_stats_gene_body": {"mean": float(effs_arr.mean()), "std": float(effs_arr.std()),
                                     "max": float(effs_arr.max()), "min": float(effs_arr.min())},
        "kyrou_positive_control": {"protospacer": KYROU_PROTOSPACER, "pam": KYROU_PAM,
                                   "present_in_gene_body_sites": bool(kyrou_present),
                                   "ranked": kyrou_rank},
        "top10_designs": [{"protospacer": r.site.protospacer, "pam": r.site.pam,
                           "strand": r.site.strand, "genomic_pos": 47600000 + r.site.position,
                           "efficacy": r.efficacy, "offtargets": r.offtarget_hits,
                           "score": r.score} for r in final[:10]],
        "drive_grid": grid,
    }
    with open(os.path.join(RESULTS, "analysis_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary
