"""Generate paper figures from analysis results."""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mosqdesign.data.loaders import read_fasta, DATA_DIR
from mosqdesign.grna_scan import scan_ngg_sites
from mosqdesign.models.efficacy_cnn import EfficacyScorer
from mosqdesign.drive_sim import simulate

RESULTS = os.path.join(os.path.dirname(__file__), "..", "results")
FIGURES = os.path.join(os.path.dirname(__file__), "..", "figures")
KYROU = "GTTTAACACAGGTCAAGCGG"


def make() -> list[str]:
    os.makedirs(FIGURES, exist_ok=True)
    summary = json.load(open(os.path.join(RESULTS, "analysis_summary.json")))
    region = list(read_fasta(os.path.join(DATA_DIR, "dsx_locus_region.fasta")).values())[0]
    scorer = EfficacyScorer(os.path.join(os.path.dirname(__file__), "..", "assets",
                                         "doench2016_efficacy_cnn.pt"))
    sites = scan_ngg_sites(region, "region")
    eff = scorer.score([s.context30 for s in sites])
    kyrou_eff = [e for s, e in zip(sites, eff) if s.protospacer == KYROU]

    # fig 1: efficacy distribution with Kyrou site marked
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.hist(eff, bins=60, color="#4472C4", edgecolor="white", linewidth=0.3)
    if kyrou_eff:
        ax.axvline(kyrou_eff[0], color="#C00000", lw=2,
                   label=f"Kyrou 2018 dsx gRNA (score {kyrou_eff[0]:.2f})")
        ax.legend(frameon=False)
    ax.set_xlabel("CNN efficacy score (Doench-2016-trained)")
    ax.set_ylabel("NGG sites in 110 kb dsx locus")
    ax.set_title("Figure 1. Predicted gRNA efficacy across the dsx locus")
    fig.tight_layout()
    p1 = os.path.join(FIGURES, "fig1_efficacy_distribution.png")
    fig.savefig(p1, dpi=200)
    plt.close(fig)

    # fig 2: top-20 ranked designs
    top = summary["top10_designs"]
    import csv
    rows = list(csv.DictReader(open(os.path.join(RESULTS, "ranked_designs.csv"))))[:20]
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [f"{r['protospacer'][:10]}... ({r['strand']}{r['genomic_pos']})" for r in rows]
    scores = [float(r["score"]) for r in rows]
    colors = ["#C00000" if r["protospacer"] == KYROU else "#4472C4" for r in rows]
    ax.barh(range(len(rows))[::-1], scores, color=colors)
    ax.set_yticks(range(len(rows))[::-1])
    ax.set_yticklabels(labels, fontsize=6, family="monospace")
    ax.set_xlabel("Composite design score")
    ax.set_title("Figure 2. Top 20 candidate drive target sites in dsx (red = published gRNA)")
    fig.tight_layout()
    p2 = os.path.join(FIGURES, "fig2_top_designs.png")
    fig.savefig(p2, dpi=200)
    plt.close(fig)

    # fig 3: drive trajectories across homing rates (e = 0.01)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for h, color in zip((0.80, 0.90, 0.95, 0.99), ("#C00000", "#ED7D31", "#70AD47", "#4472C4")):
        traj = simulate(0.99, 0.01, 0.0, h=h, e=0.01, c_hom=1.0, generations=60)
        ax.plot(traj[:, 1], color=color, label=f"h = {h}")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Drive allele frequency")
    ax.set_title("Figure 3. Drive spread vs homing rate (resistance rate e = 0.01)")
    ax.legend(frameon=False, title="homing rate")
    fig.tight_layout()
    p3 = os.path.join(FIGURES, "fig3_drive_trajectories.png")
    fig.savefig(p3, dpi=200)
    plt.close(fig)

    # fig 4: resistance accumulation vs e (h = 0.99)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for e, color in zip((0.001, 0.01, 0.1, 0.5), ("#4472C4", "#70AD47", "#ED7D31", "#C00000")):
        traj = simulate(0.99, 0.01, 0.0, h=0.99, e=e, c_hom=1.0, generations=60)
        ax.plot(traj[:, 2], color=color, label=f"e = {e}")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Resistant allele frequency")
    ax.set_title("Figure 4. End-joining resistance accumulation (homing rate h = 0.99)")
    ax.legend(frameon=False, title="resistance rate")
    fig.tight_layout()
    p4 = os.path.join(FIGURES, "fig4_resistance.png")
    fig.savefig(p4, dpi=200)
    plt.close(fig)

    # fig 5: stochastic Wright-Fisher replicates
    from mosqdesign.drive_sim import wright_fisher
    fig, ax = plt.subplots(figsize=(7, 4.2))
    det = simulate(0.99, 0.01, 0.0, h=0.99, e=0.01, c_hom=1.0, generations=40)
    ax.plot(det[:, 1], color="black", lw=2, label="deterministic")
    for seed in range(5):
        wf = wright_fisher(10000, 0.99, 0.01, 0.0, h=0.99, e=0.01, c_hom=1.0,
                           generations=40, seed=seed)
        ax.plot(wf[:, 1], color="#4472C4", alpha=0.5)
    ax.plot([], [], color="#4472C4", alpha=0.5, label="Wright-Fisher (N=10,000), 5 reps")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Drive allele frequency")
    ax.set_title("Figure 5. Deterministic vs stochastic drive spread (h = 0.99, e = 0.01)")
    ax.legend(frameon=False)
    fig.tight_layout()
    p5 = os.path.join(FIGURES, "fig5_wright_fisher.png")
    fig.savefig(p5, dpi=200)
    plt.close(fig)
    return [p1, p2, p3, p4, p5]
