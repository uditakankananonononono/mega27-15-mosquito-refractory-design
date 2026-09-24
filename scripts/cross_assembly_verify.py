"""Cross-assembly verification of dsx guide targets: AgamP4 vs AgamP5.

The honest replication path for the on-target audit: every specificity and
window claim made on AgamP4 (GCA_000005575.1) is re-checked on the AgamP5
reference (NC_064601.1 2R; fetched by scripts/fetch_genome.sh from NCBI
RefSeq). For each guide we search AgamP5 2R (both strands) for the 23-nt
protospacer+PAM unit in plus-strand orientation and report:

  * exact matches (count, positions),
  * near matches at <=1, <=2, <=3 mismatches (pigeonhole search: the
    23-mer is split into 4 blocks; any <=3-mismatch hit has >=1 intact
    block, so exact block hits seed a full Hamming check),
  * whether the AgamP4 target allele is reproduced identically in AgamP5
    (assembly-consistency of the on-target sequence itself).

No alignment tool is invoked; the method is exact enumeration with a
combinatorial filter, so results are reproducible and auditable.

Usage: python3 scripts/cross_assembly_verify.py --agamp5 /path/NC_064601.fa
Output: results/cross_assembly_verify.{json,md}
"""
import argparse
import json

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


# plus-strand 23-mers (proto+PAM) for each minus-strand guide, from the
# twice-verified window constants (scripts/ag1000g_allele_freq.py GUIDES):
# plus proto = rc(ref_minus), plus PAM = rc(minus PAM).
def plus_units():
    # Minus-strand targets: minus 5'->3' over the 23-nt window is
    # ref_minus(proto) + pam_minus (high to low plus coords). The
    # plus-strand 5'->3' unit is therefore rc(ref_minus + pam_minus),
    # i.e. CCN-PAM first, then protospacer.
    from ag1000g_allele_freq import GUIDES
    pam_minus = {"dsx-v3-1": "GGG", "dsx-v3-2": "TGG", "kyrou": "TGG"}
    return {g: rc(d["ref_minus"] + pam_minus[g])
            for g, d in GUIDES.items()}


def guide_units():
    """Guide-derived 23-mers: what a perfect target allele looks like.
    Same orientation logic as plus_units but with the guide sequence."""
    from ag1000g_allele_freq import GUIDES
    pam_minus = {"dsx-v3-1": "GGG", "dsx-v3-2": "TGG", "kyrou": "TGG"}
    return {g: rc(d["guide"] + pam_minus[g]) for g, d in GUIDES.items()}


def blocks(seq23):
    """4 near-equal blocks covering the 23-mer."""
    return [(0, 6), (6, 12), (12, 18), (18, 23)]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def search_unit(genome, unit, max_mm=3):
    """Exact + near-match search of a 23-mer over one strand (string)."""
    n = len(genome)
    L = len(unit)
    exact = find_all(genome, unit)
    cand = set()
    for a, b in blocks(unit):
        for pos in find_all(genome, unit[a:b]):
            start = pos - a
            if 0 <= start <= n - L:
                cand.add(start)
    exact_set = set(exact)
    near = {1: [], 2: [], 3: []}
    for s in cand:
        if s in exact_set:
            continue
        mm = hamming(genome[s:s + L], unit)
        if 1 <= mm <= max_mm:
            near[mm].append(s)
    return {"exact": exact, "near": near}


def load_fasta(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    return seq.upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--agamp5", required=True,
                    help="AgamP5 2R FASTA (NC_064601.1)")
    ap.add_argument("--out-json",
                    default="results/cross_assembly_verify.json")
    ap.add_argument("--out-md", default="results/cross_assembly_verify.md")
    args = ap.parse_args()
    genome = load_fasta(args.agamp5)
    genome_rc = rc(genome)
    n = len(genome)
    units = plus_units()
    gunits = guide_units()
    report = {}
    for g, unit in units.items():
        gunit = gunits[g]
        gfwd = search_unit(genome, gunit)
        grev = search_unit(genome_rc, gunit)
        fwd = search_unit(genome, unit)
        rev = search_unit(genome_rc, unit)
        # rc-strand hit at index s corresponds to plus coord n-s-23
        rev_plus = {mm: [n - s - 23 for s in v] for mm, v in rev["near"].items()}
        report[g] = {
            "unit_plus": unit,
            "guide_unit_plus": gunit,
            "guide_exact_plus": gfwd["exact"],
            "guide_exact_minus": [n - s - 23 for s in grev["exact"]],
            "guide_near": {
                "plus": {str(k): v for k, v in gfwd["near"].items()},
                "minus": {str(k): [n - s - 23 for s in grev["near"][k]]
                          for k in (1, 2, 3)},
            },
            "agamp5_len": n,
            "exact_plus": fwd["exact"],
            "exact_minus": [n - s - 23 for s in rev["exact"]],
            "near_plus": {str(k): v for k, v in fwd["near"].items()},
            "near_minus": {str(k): rev_plus[k] for k in (1, 2, 3)},
        }
    with open(args.out_json, "w") as fh:
        json.dump({"method": "exact + pigeonhole near-match (4 blocks, "
                             "<=3 mismatches), both strands, no aligner",
                   "agamp5_query": args.agamp5, "guides": report},
                  fh, indent=1)
    lines = ["## Cross-assembly verification: guide targets on AgamP5 "
             "(NC_064601.1 2R)", "",
             "Exact and <=3-mismatch matches of each 23-nt protospacer+PAM "
             "unit (plus-strand orientation) on both strands of AgamP5 2R. "
             "Method: exact substring + pigeonhole near-match (4 blocks), "
             "no aligner.", "",
             "AgamP4-reference units and guide-derived units searched "
             "separately: AgamP4 can carry a minor allele (assembly "
             "mismatch); the guide unit is the intended target allele.", "",
             "| guide | AgamP4 unit exact +/- | guide unit exact +/- | "
             "guide mm1 | guide mm2 | guide mm3 |",
             "|---|---|---|---|---|---|---|"]
    for g, r in report.items():
        gn = r["guide_near"]
        lines.append(
            f"| {g} | {len(r['exact_plus'])}/{len(r['exact_minus'])} "
            f"| {len(r['guide_exact_plus'])}/{len(r['guide_exact_minus'])} "
            f"| {len(gn['plus']['1']) + len(gn['minus']['1'])} | "
            f"{len(gn['plus']['2']) + len(gn['minus']['2'])} | "
            f"{len(gn['plus']['3']) + len(gn['minus']['3'])} |")
    with open(args.out_md, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(json.dumps({g: {"ref_exact": len(r["exact_plus"]) +
                                len(r["exact_minus"]),
                          "guide_exact": len(r["guide_exact_plus"]) +
                                len(r["guide_exact_minus"]),
                          "guide_mm1": len(r["guide_near"]["plus"]["1"]) +
                                len(r["guide_near"]["minus"]["1"]),
                          "guide_mm2": len(r["guide_near"]["plus"]["2"]) +
                                len(r["guide_near"]["minus"]["2"]),
                          "guide_mm3": len(r["guide_near"]["plus"]["3"]) +
                                len(r["guide_near"]["minus"]["3"])}
                      for g, r in report.items()}, indent=1))


if __name__ == "__main__":
    main()
