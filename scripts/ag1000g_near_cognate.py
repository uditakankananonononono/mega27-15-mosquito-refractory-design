"""Near-cognate target enumeration at the dsx guide windows from Ag1000G
on-target calls (AgamP4 coordinates, minus-strand guides).

Beyond the binary intact/compromised classification of
scripts/ag1000g_allele_freq.py, this script enumerates the distinct
near-cognate protospacer/PAM alleles actually standing in the population:
every observed ALT allele at a target-window site defines a single-variant
near-cognate allele whose cleavage phenotype (PAM disruption, seed
mismatch, distal mismatch, tolerated PAM-N change) determines whether it
is pre-existing resistance substrate.

Phasing discipline (unphased VCF):
  * A single-variant allele is confirmed present on at least one
    chromosome whenever its ALT is called (het or hom); its chromosome
    count is the VCF allele count (AC).
  * Multi-variant haplotypes are enumerated ONLY when phase is certain:
    a sample homozygous-ALT at every variant position it carries inside
    the window contributes exactly one multi-variant allele (count 2).
    Het-at-multiple-sites samples are NOT phased; their constituent
    single-variant alleles are already counted and the residual ambiguity
    is reported as n_ambiguous_multisite_samples.

Per enumerated allele we report, in guide (minus-strand) orientation:
  * protospacer allele sequence (20 nt) and PAM allele (3 nt),
  * mismatch count vs guide, mismatched guide indices (1-20, 1 = 5' end),
    seed mismatches (8 PAM-proximal bases = guide indices 13-20),
  * PAM class: intact / N-variant (tolerated) / GG-disrupted,
  * predicted class: pam_disrupted | seed_mismatch | distal_mismatch |
    pam_n_variant | wild_type_equivalent,
  * chromosome count (AC) and AF denominator note.

Outputs:
  results/ag1000g/near_cognate.json - full enumeration records
  results/ag1000g/near_cognate.md   - paper-ready table
"""
import argparse
import collections
import json

from ag1000g_allele_freq import GUIDES, position_map, parse_gt, load_records

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def minus_allele(base_plus):
    """Plus-strand base to minus-strand complement."""
    return base_plus.translate(COMP)


def classify_variant(pos, alt_base_plus, meta):
    """Functional class of one single-base variant at a window position."""
    role = meta["role"]
    if role == "PAM-N":
        return "pam_n_variant"
    if role == "PAM-G":
        return "pam_disrupted" if alt_base_plus != "C" else "wild_type"
    # protospacer
    if alt_base_plus == meta["guide_plus"]:
        return "wild_type"  # ALT restores guide match (assembly-mismatch pos)
    return "seed_mismatch" if role == "seed" else "distal_mismatch"


def guide_index(pos, gname):
    """Guide base index 1-20 (1 = 5' end of guide) for a proto position."""
    pa, pb = GUIDES[gname]["proto"]
    return pb - pos + 1


def enumerate_alleles(records, min_gq=20, min_dp=5):
    """Enumerate confirmed near-cognate alleles per guide window."""
    pos2guide = {}
    for gname in GUIDES:
        for pos, meta in position_map(gname).items():
            pos2guide[pos] = (gname, meta)

    # guide -> pos -> alt_base -> AC
    alt_counts = {g: collections.defaultdict(collections.Counter)
                  for g in GUIDES}
    # guide -> sample -> list of (pos, gt_indices, alts, gq, dp)
    per_sample = {g: collections.defaultdict(list) for g in GUIDES}
    n_lowqual = 0

    for rec in records:
        sid = rec["sample"]
        for hit in rec.get("target_hits", []):
            pos = hit["pos"]
            if pos not in pos2guide:
                continue
            gname, meta = pos2guide[pos]
            gq, dp = hit.get("gq", 0), hit.get("dp", 0)
            gt = parse_gt(hit["gt"], hit["alt"])
            if gt is None or gq < min_gq or dp < min_dp:
                n_lowqual += 1
                continue
            alts = hit["alt"].split(",")
            for a in gt:
                if a == 0:
                    continue
                alt_counts[gname][pos][alts[a - 1]] += 1
            per_sample[gname][sid].append(
                (pos, gt, alts, hit["ref"]))

    enumerated = {}
    for gname in GUIDES:
        g = GUIDES[gname]
        pa, pb = g["proto"]
        sa, sb = g["pam"]
        alleles = []
        # single-variant alleles (confirmed; AC = chromosome count)
        for pos, counter in sorted(alt_counts[gname].items()):
            meta = position_map(gname)[pos]
            for alt_base, ac in sorted(counter.items()):
                cls = classify_variant(pos, alt_base, meta)
                # build allele sequences in guide orientation
                if meta["role"].startswith("PAM"):
                    proto_alt = g["guide"]  # proto unchanged
                    k = sb - pos  # 0 N, 1 G2, 2 G3
                    pam_minus = list(rc(_pam_plus_ref(gname)))
                    pam_minus[k] = minus_allele(alt_base)
                    pam_alt = "".join(pam_minus)
                    mm_idx, seed_mm = [], 0
                else:
                    i = pb - pos
                    proto_alt = (g["ref_minus"][:i] +
                                 minus_allele(alt_base) +
                                 g["ref_minus"][i + 1:])
                    pam_alt = rc(_pam_plus_ref(gname))
                    idx = guide_index(pos, gname)
                    mm_idx = [idx] if cls in ("seed_mismatch",
                                              "distal_mismatch") else []
                    seed_mm = 1 if cls == "seed_mismatch" else 0
                n_mm = (len(mm_idx) +
                        (1 if cls == "pam_disrupted" else 0))
                alleles.append({
                    "kind": "single_variant",
                    "pos": pos, "ref_plus": meta["ref_plus"],
                    "alt_plus": alt_base, "ac": ac,
                    "proto_allele_minus": proto_alt,
                    "pam_allele_minus": pam_alt,
                    "n_mismatch_vs_guide": n_mm,
                    "mismatch_guide_indices": mm_idx,
                    "seed_mismatches": seed_mm,
                    "class": cls,
                })
        # phase-certain multi-variant alleles: hom-ALT at every variant
        # position the sample carries in this window
        n_ambiguous = 0
        multi = []
        for sid, hits in per_sample[gname].items():
            if len(hits) < 2:
                continue
            all_hom = all(gt[0] == gt[1] and gt[0] != 0
                          for _, gt, _, _ in hits)
            if not all_hom:
                n_ambiguous += 1
                continue
            proto_alt = list(g["ref_minus"])
            pam_plus = list(_pam_plus_ref(gname))
            mm_idx, seed_mm, pam_disrupted = [], 0, False
            classes = []
            for pos, gt, alts, ref in hits:
                meta = position_map(gname)[pos]
                alt_base = alts[gt[0] - 1]
                cls = classify_variant(pos, alt_base, meta)
                classes.append(cls)
                if meta["role"].startswith("PAM"):
                    k_plus = pos - sa  # plus-strand PAM index
                    pam_plus[k_plus] = alt_base
                    if cls == "pam_disrupted":
                        pam_disrupted = True
                else:
                    i = pb - pos
                    proto_alt[i] = minus_allele(alt_base)
                    if cls in ("seed_mismatch", "distal_mismatch"):
                        mm_idx.append(guide_index(pos, gname))
                        if cls == "seed_mismatch":
                            seed_mm += 1
            cls = ("pam_disrupted" if pam_disrupted else
                   "seed_mismatch" if seed_mm else
                   "distal_mismatch" if mm_idx else
                   "pam_n_variant" if "pam_n_variant" in classes else
                   "wild_type_equivalent")
            multi.append({
                "kind": "multi_variant_phase_certain",
                "sample": sid, "ac": 2,
                "variant_positions": [h[0] for h in hits],
                "proto_allele_minus": "".join(proto_alt),
                "pam_allele_minus": rc("".join(pam_plus)),
                "n_mismatch_vs_guide": len(mm_idx) + (1 if pam_disrupted
                                                      else 0),
                "mismatch_guide_indices": sorted(mm_idx),
                "seed_mismatches": seed_mm,
                "class": cls,
            })
        enumerated[gname] = {
            "single_variant_alleles": alleles,
            "multi_variant_phase_certain": multi,
            "n_ambiguous_multisite_samples": n_ambiguous,
        }
    return {"n_samples": len(records), "n_lowqual_calls": n_lowqual,
            "guides": enumerated}


# plus-strand PAM reference sequences (verified window constants)
_PAM_PLUS = {
    "dsx-v3-1": "CCC",   # 48711498-48711500, minus 5'-GGG-3'
    "dsx-v3-2": "CCA",   # 48712762-48712764, minus 5'-TGG-3'
    "kyrou":    "CCA",   # 48714637-48714639, minus 5'-TGG-3'
}


def _pam_plus_ref(gname):
    return _PAM_PLUS[gname]


def summarize(res):
    out = {}
    for gname, g in res["guides"].items():
        sv = g["single_variant_alleles"]
        cls_counts = collections.Counter(a["class"] for a in sv)
        cls_ac = collections.Counter()
        for a in sv:
            cls_ac[a["class"]] += a["ac"]
        out[gname] = {
            "n_single_variant_alleles": len(sv),
            "by_class_counts": dict(cls_counts),
            "by_class_chromosomes": dict(cls_ac),
            "n_multi_variant_phase_certain":
                len(g["multi_variant_phase_certain"]),
            "n_ambiguous_multisite_samples":
                g["n_ambiguous_multisite_samples"],
        }
    return out


def render_md(res, summary):
    n = res["n_samples"]
    lines = [
        f"## Near-cognate allele enumeration at dsx guide windows "
        f"(n = {n} samples)", "",
        "Single-variant alleles confirmed by VCF allele counts; "
        "multi-variant haplotypes enumerated only when phase is certain "
        "(homozygous-ALT at every carried variant). Het-multisite samples "
        "are reported as ambiguous, never phased by assumption.", "",
        "| guide | class | distinct alleles | chromosomes (AC) |",
        "|---|---|---|---|",
    ]
    order = ["pam_disrupted", "seed_mismatch", "distal_mismatch",
             "pam_n_variant", "wild_type"]
    for gname, s in summary.items():
        for cls in order:
            c = s["by_class_counts"].get(cls, 0)
            ac = s["by_class_chromosomes"].get(cls, 0)
            if c:
                lines.append(f"| {gname} | {cls} | {c} | {ac} |")
        lines.append(
            f"| {gname} | multi-variant (phase-certain) | "
            f"{s['n_multi_variant_phase_certain']} | - |")
        lines.append(
            f"| {gname} | ambiguous multisite samples (unphased) | "
            f"{s['n_ambiguous_multisite_samples']} | - |")
    lines += ["", "### Enumerated single-variant alleles", "",
              "| guide | pos | ref>alt | class | guide idx | AC |",
              "|---|---|---|---|---|---|"]
    for gname, g in res["guides"].items():
        for a in g["single_variant_alleles"]:
            idx = ",".join(str(i) for i in a["mismatch_guide_indices"]) or "-"
            lines.append(
                f"| {gname} | {a['pos']} | {a['ref_plus']}>{a['alt_plus']} "
                f"| {a['class']} | {idx} | {a['ac']} |")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="results/ag1000g/on_target.jsonl")
    ap.add_argument("--out-json", default="results/ag1000g/near_cognate.json")
    ap.add_argument("--out-md", default="results/ag1000g/near_cognate.md")
    args = ap.parse_args()
    records = load_records(args.input)
    res = enumerate_alleles(records)
    summary = summarize(res)
    with open(args.out_json, "w") as fh:
        json.dump({"summary": summary, "detail": res}, fh, indent=2)
    with open(args.out_md, "w") as fh:
        fh.write(render_md(res, summary))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
