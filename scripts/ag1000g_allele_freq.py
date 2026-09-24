"""Per-site allele-frequency analysis of Ag1000G on-target variation at the
dsx guide targets (AgamP4 coordinates, minus-strand guides).

Reads the per-sample JSONL produced by scripts/ag1000g_scan_all.py /
ag1000g_on_target.py (one row per sample; every non-reference genotype call
inside the fetch window 2R:48,711,450-48,714,700) and computes, per target
window site:

  * allele counts (AC) and allele frequencies (AF) per ALT allele, with an
    all-sites-VCF denominator: samples with no record at a site are
    homozygous-reference (window rows are complete: n_rows == 3250);
  * guide-allele classification: at AgamP4/guide mismatch positions the ALT
    matching the guide is the *intact* (cleavable) allele, so the AgamP4
    reference allele itself scores as a resistance (mismatch) allele;
  * functional role: PAM N position (tolerated), PAM G2/G3 (disruptive),
    protospacer seed (<=7 nt from PAM) vs distal protospacer;
  * per-sample, per-guide compromise status (intact / het / hom) defined as:
    an allele is compromised if it disrupts the PAM GG dinucleotide or
    creates any mismatch vs the guide within the seed.

Outputs:
  results/ag1000g/allele_freq.json  - full per-site and per-sample records
  results/ag1000g/allele_freq.md    - paper-ready tables (section 4.7)

Verified reference bases (data/agamp4/AgamP4_2R_dsx_region.fasta, re-checked
against results/agamp4_reannotation.json aligned_sequence fields):
  dsx-v3-1 proto plus ACCCTACCGCTTACTACCCA, PAM plus CCC (minus GGG)
  dsx-v3-2 proto plus ACGCTTCGTAGGTCTTAATG, PAM plus CCA (minus TGG)
  kyrou    proto plus CCGCTTGACCTGTGTTAAAC, PAM plus CCA (minus TGG)
"""
import argparse
import collections
import json
import os

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


# 1-based inclusive AgamP4 2R windows (verified against the reannotation).
GUIDES = {
    "dsx-v3-1": {
        "guide": "TGGGCAGTATGCGTTAGGGT",      # minus strand, 5'->3'
        "ref_minus": "TGGGTAGTAAGCGGTAGGGT",  # AgamP4 aligned sequence
        "proto": (48711501, 48711520),
        "pam": (48711498, 48711500),          # minus-strand 5'-GGG-3'
    },
    "dsx-v3-2": {
        "guide": "CATTAAGACCTACGAAGCGC",
        "ref_minus": "CATTAAGACCTACGAAGCGT",
        "proto": (48712765, 48712784),
        "pam": (48712762, 48712764),          # minus-strand 5'-TGG-3'
    },
    "kyrou": {
        "guide": "GTTTAACACAGGTCAAGCGG",
        "ref_minus": "GTTTAACACAGGTCAAGCGG",  # exact match to AgamP4
        "proto": (48714640, 48714659),
        "pam": (48714637, 48714639),          # minus-strand 5'-TGG-3'
    },
}

SEED_MAX_DIST = 7  # first 8 PAM-proximal protospacer bases


def position_map(gname):
    """Per-genomic-position annotation for one guide.

    Returns {pos: dict(role, ref_plus, guide_plus, dist_from_pam)}.
    guide_plus is None at PAM positions. For protospacer positions the
    PAM-proximal end is the low-coordinate end (minus-strand guide, PAM
    upstream); dist_from_pam = pos - proto_start.
    """
    g = GUIDES[gname]
    pa, pb = g["proto"]
    out = {}
    for pos in range(pa, pb + 1):
        i = pb - pos  # guide index: 5' guide base at highest coordinate
        out[pos] = {
            "role": "seed" if pos - pa <= SEED_MAX_DIST else "proto",
            "ref_plus": g["ref_minus"][i].translate(COMP),
            "guide_plus": g["guide"][i].translate(COMP),
            "dist_from_pam": pos - pa,
        }
    sa, sb = g["pam"]
    # minus-strand PAM 5'-N G G-3': N at highest coord, G2 next, G3 lowest
    for pos in range(sa, sb + 1):
        k = sb - pos  # 0 -> N, 1 -> G2, 2 -> G3
        role = "PAM-N" if k == 0 else "PAM-G"
        out[pos] = {
            "role": role,
            "ref_plus": "C",  # N may differ per guide; G is always C on plus
            "guide_plus": None,
            "dist_from_pam": None,
        }
    # fix the N-position ref base from the actual PAM sequence is unnecessary
    # for classification (any ALT at N is tolerated); keep recorded ref.
    return out


def parse_gt(gt, alt_field):
    """Decode a diploid VCF genotype into two allele bases (0 = ref)."""
    alts = alt_field.split(",")
    sep = "|" if "|" in gt else "/"
    idx = gt.split(sep)
    if len(idx) != 2:
        return None
    out = []
    for x in idx:
        if x == ".":
            return None
        out.append(int(x))
    return out  # indices into [ref] + alts


def load_records(path):
    recs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                recs.append(json.loads(line))
    return recs


def compute(records, min_gq=20, min_dp=5):
    """Aggregate per-site allele counts and per-sample compromise status.

    Compromise status is unphased: a sample het-disruptive at two sites
    in the same guide window is called het (max over positions), a lower
    bound on the true compromised-allele count.
    """
    n_samples = len(records)
    # site -> stats
    sites = collections.defaultdict(lambda: {
        "alt_ac": collections.Counter(), "ref_from_calls": 0,
        "n_called": 0, "n_lowqual": 0, "ref": None,
    })
    # sample -> guide -> compromised allele count (0,1,2)
    sample_status = {r["sample"]: {g: 0 for g in GUIDES} for r in records}
    pos2guide = {}
    for gname in GUIDES:
        for pos, meta in position_map(gname).items():
            pos2guide[pos] = (gname, meta)

    for rec in records:
        sid = rec["sample"]
        for hit in rec.get("target_hits", []):
            pos = hit["pos"]
            if pos not in pos2guide:
                continue
            st = sites[pos]
            gq, dp = hit.get("gq", 0), hit.get("dp", 0)
            gt = parse_gt(hit["gt"], hit["alt"])
            if gt is None or gq < min_gq or dp < min_dp:
                st["n_lowqual"] += 1
                continue
            alts = hit["alt"].split(",")
            st["n_called"] += 1
            st["ref"] = hit["ref"]
            gname, meta = pos2guide[pos]
            for a in gt:
                if a == 0:
                    st["ref_from_calls"] += 1
                else:
                    st["alt_ac"][alts[a - 1]] += 1
            # compromise status per allele
            comp = 0
            for a in gt:
                base = hit["ref"] if a == 0 else alts[a - 1]
                if allele_disrupts(base, meta):
                    comp += 1
            sample_status[sid][gname] = max(sample_status[sid][gname], comp)

    # finalize per-site AF
    out_sites = {}
    for pos, st in sorted(sites.items()):
        gname, meta = pos2guide[pos]
        an = 2 * (n_samples - st["n_lowqual"])
        n_unobserved = n_samples - st["n_called"] - st["n_lowqual"]
        ref_ac = 2 * n_unobserved + st["ref_from_calls"]
        alt_af = {a: c / an for a, c in sorted(st["alt_ac"].items())}
        entry = {
            "pos": pos, "guide": gname, "role": meta["role"],
            "ref": st["ref"], "alt_ac": dict(st["alt_ac"]), "alt_af": alt_af,
            "ref_ac": ref_ac, "an": an,
            "n_called": st["n_called"], "n_lowqual": st["n_lowqual"],
            "dist_from_pam": meta["dist_from_pam"],
            "guide_plus": meta["guide_plus"],
        }
        out_sites[pos] = entry
    return {"n_samples": n_samples, "sites": out_sites,
            "sample_status": sample_status}


def allele_disrupts(base, meta):
    """Does this allele disrupt cleavage at this position?"""
    if meta["role"] == "PAM-N":
        return False  # any base tolerated at N
    if meta["role"] == "PAM-G":
        return base != "C"  # G on minus strand = C on plus
    return base != meta["guide_plus"]  # any mismatch vs guide


def guide_summary(res):
    """Per-guide headline numbers: carrier/compromise frequencies."""
    n = res["n_samples"]
    out = {}
    for g in GUIDES:
        counts = collections.Counter(res["sample_status"][s][g]
                                     for s in res["sample_status"])
        compromised_alleles = sum(res["sample_status"][s][g]
                                  for s in res["sample_status"])
        out[g] = {
            "hom_compromised": counts[2], "het": counts[1],
            "intact": counts[0],
            "carrier_fraction": (counts[1] + counts[2]) / n,
            "compromised_allele_fraction": compromised_alleles / (2 * n),
        }
    return out


def render_md(res, summary):
    n = res["n_samples"]
    lines = [
        f"## Ag1000G on-target polymorphism at dsx guide targets "
        f"(n = {n} samples)",
        "",
        "Unphased genotypes: het-at-two-sites is called het "
        "(lower bound on compromised alleles).",
        "",
        "Per-guide compromise status (allele compromised = PAM GG disrupted "
        "or any mismatch vs guide within the 8-nt PAM-proximal seed):",
        "",
        "| guide | intact | het | hom-compromised | carrier fraction | "
        "compromised AF |",
        "|---|---|---|---|---|---|",
    ]
    for g, s in summary.items():
        lines.append(
            f"| {g} | {s['intact']} | {s['het']} | {s['hom_compromised']} "
            f"| {s['carrier_fraction']:.4f} | "
            f"{s['compromised_allele_fraction']:.4f} |")
    lines += ["", "### Variant sites inside target windows", "",
              "| guide | pos | role | dist from PAM | ref->alt (AC, AF) |",
              "|---|---|---|---|---|"]
    for pos, e in res["sites"].items():
        alts = "; ".join(f"{a} (AC {c}, AF {e['alt_af'][a]:.4f})"
                         for a, c in sorted(e["alt_ac"].items()))
        d = e["dist_from_pam"]
        lines.append(f"| {e['guide']} | {pos} | {e['role']} | "
                     f"{d if d is not None else '-'} | {alts} |")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="results/ag1000g/on_target.jsonl")
    ap.add_argument("--outdir", default="results/ag1000g")
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--min-dp", type=int, default=5)
    args = ap.parse_args()
    records = load_records(args.input)
    res = compute(records, args.min_gq, args.min_dp)
    summary = guide_summary(res)
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "allele_freq.json"), "w") as fh:
        json.dump({"n_samples": res["n_samples"],
                   "sites": {str(k): v for k, v in res["sites"].items()},
                   "guide_summary": summary,
                   "sample_status": res["sample_status"]}, fh, indent=1)
    with open(os.path.join(args.outdir, "allele_freq.md"), "w") as fh:
        fh.write(render_md(res, summary))
    print(f"n_samples={res['n_samples']} variant_sites={len(res['sites'])}")
    for g, s in summary.items():
        print(f"{g}: intact={s['intact']} het={s['het']} "
              f"hom={s['hom_compromised']} "
              f"compromised_AF={s['compromised_allele_fraction']:.4f}")


if __name__ == "__main__":
    main()
