"""v3 ranking: constraint-aware score fused with chromosome-wide specificity.

v2 ranks by CNN efficacy + exon/splice constraint bonus but its off-target
term came from the 110 kb locus only. The chromosome-2 scan (118 Mb) supplies
the real specificity burden per guide. v3 applies:
  - hard filter: any exact or 1-mismatch genomic off-target disqualifies
    (mm0 + mm1 == 0 required);
  - penalty: 0.50 per 2-mismatch hit, 0.05 per 3-mismatch hit
    (CFD-style decay: distal mismatches tolerate better).
The Kyrou 2018 guide is carried as a published control row.
"""
from __future__ import annotations

MM2_PENALTY = 0.50
MM3_PENALTY = 0.05


def specificity_penalty(mm2: int, mm3: int) -> float:
    return MM2_PENALTY * mm2 + MM3_PENALTY * mm3


def passes_specificity_filter(mm0: int, mm1: int) -> bool:
    return mm0 == 0 and mm1 == 0


def score_v3(efficacy: float, constraint_bonus: float, mm2: int, mm3: int) -> float:
    return efficacy + constraint_bonus - specificity_penalty(mm2, mm3)


def run(top_n: int = 30) -> dict:
    """Regenerate results/ranked_designs_v3.csv and results/v3_summary.json.

    Fuses the v2 top-`top_n` candidates (efficacy + constraint from
    results/ranked_designs_v2.csv) with exact genome-wide mismatch counts
    (results/genome_offtargets_genomewide_top30.json) under the v3 rule:
    hard filter mm0 + mm1 == 0; penalty 0.50 per mm2, 0.05 per mm3.
    The Kyrou 2018 guide is appended as a published control row, with its
    efficacy re-scored from the archived locus and its genome-wide audit
    counts taken from results/kyrou_genomewide_audit.json.
    """
    import csv
    import json
    import os

    from mosqdesign.data.loaders import read_fasta, DATA_DIR
    from mosqdesign.grna_scan import scan_ngg_sites
    from mosqdesign.models.efficacy_cnn import EfficacyScorer

    results = os.path.join(os.path.dirname(__file__), "..", "results")
    v2 = list(csv.DictReader(open(os.path.join(results, "ranked_designs_v2.csv"))))[:top_n]
    scan = json.load(open(os.path.join(results, "genome_offtargets_genomewide_top30.json")))
    mm = {p: {m: 0 for m in range(4)} for p in [r["protospacer"] for r in v2]}
    for h in scan["hits"]:
        if h["guide"] in mm and h["mismatches"] <= 3:
            mm[h["guide"]][h["mismatches"]] += 1

    rows = []
    for r in v2:
        p = r["protospacer"]
        g = mm[p]
        eff = float(r["efficacy"])
        bonus = float(r["constraint_bonus"])
        score = score_v3(eff, bonus, g[2], g[3])
        rows.append({
            "protospacer": p, "strand": r["strand"],
            "genomic_pos": r["genomic_pos"], "efficacy": eff,
            "constraint_bonus": bonus, "constraint_label": r["constraint_label"],
            "gw_mm0": g[0], "gw_mm1": g[1], "gw_mm2": g[2], "gw_mm3": g[3],
            "specificity_filter": "pass" if passes_specificity_filter(g[0], g[1]) else "fail",
            "score_v3": round(score, 4),
        })
    passing = sorted((r for r in rows if r["specificity_filter"] == "pass"),
                     key=lambda r: (-r["score_v3"], r["protospacer"]))
    for i, r in enumerate(passing, 1):
        r["rank_v3"] = i
    failing = [r for r in rows if r["specificity_filter"] == "fail"]  # v2 order
    for r in failing:
        r["specificity_filter"] = "FAIL"
        r["score_v3"] = ""
        r["rank_v3"] = ""
    rows = passing + failing

    # control: re-score the Kyrou guide from the archived locus
    region = list(read_fasta(os.path.join(DATA_DIR, "dsx_locus_region.fasta")).values())[0]
    kyrou = "GTTTAACACAGGTCAAGCGG"
    ky_sites = [s for s in scan_ngg_sites(region, "region") if s.protospacer == kyrou]
    scorer = EfficacyScorer(os.path.join(os.path.dirname(__file__), "..", "assets",
                                         "doench2016_efficacy_cnn.pt"))
    ky_eff = round(float(scorer.score([s.context30 for s in ky_sites])[0]), 3)
    audit = json.load(open(os.path.join(results, "kyrou_genomewide_audit.json")))
    kc = audit["counts"]
    control = {
        "protospacer": kyrou, "strand": "-", "genomic_pos": 47622174,
        "efficacy": ky_eff, "constraint_bonus": 0.0,
        "constraint_label": "published control (Kyrou 2018); intronic in current RefSeq model",
        "gw_mm0": int(kc["0"]), "gw_mm1": int(kc["1"]), "gw_mm2": int(kc["2"]),
        "gw_mm3": int(kc["3"]), "specificity_filter": "pass",
        "score_v3": round(score_v3(ky_eff, 0.0, int(kc["2"]), int(kc["3"])), 3),
        "rank_v3": "control",
    }

    fields = ["protospacer", "strand", "genomic_pos", "efficacy", "constraint_bonus",
              "constraint_label", "gw_mm0", "gw_mm1", "gw_mm2", "gw_mm3",
              "specificity_filter", "score_v3", "rank_v3"]
    with open(os.path.join(results, "ranked_designs_v3.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    summary = {
        "reference": "AgamP5 GCF_943734735.2 full genome (2RL+X+3RL+MT, ~246 Mb)",
        "control_kyrou": control,
        "n_pass": len(passing),
        "n_fail": len(rows) - len(passing),
        "scoring": "v3 = efficacy + constraint_bonus - 0.50*mm2 - 0.05*mm3; hard filter mm0+mm1==0; exact genome-wide enumeration",
    }
    with open(os.path.join(results, "v3_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    return summary
