"""Cross-assembly verification of dsx guide targets: AgamP4 vs AgamP5.

The honest replication path for the on-target audit: every specificity and
window claim made on AgamP4 (GCA_000005575.1) is re-checked on the full
AgamP5 reference (NC_064601.1 chr2, OX030908.1 chr3, OX030909.1 chrX,
OX030910.2 MT; fetched by scripts/fetch_genome.sh from NCBI RefSeq). For
each guide we search every molecule (both strands) for the 23-nt
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


def _search_one_molecule(genome, unit, gunit):
    """Search one molecule (both strands) with ref and guide units."""
    genome_rc = rc(genome)
    n = len(genome)
    gfwd = search_unit(genome, gunit)
    grev = search_unit(genome_rc, gunit)
    fwd = search_unit(genome, unit)
    rev = search_unit(genome_rc, unit)
    return {
        "length": n,
        "exact_plus": fwd["exact"],
        "exact_minus": [n - s - 23 for s in rev["exact"]],
        "near_plus": {str(k): v for k, v in fwd["near"].items()},
        "near_minus": {str(k): [n - s - 23 for s in rev["near"][k]]
                       for k in (1, 2, 3)},
        "guide_exact_plus": gfwd["exact"],
        "guide_exact_minus": [n - s - 23 for s in grev["exact"]],
        "guide_near": {
            "plus": {str(k): v for k, v in gfwd["near"].items()},
            "minus": {str(k): [n - s - 23 for s in grev["near"][k]]
                      for k in (1, 2, 3)},
        },
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--agamp5", required=True, nargs="+",
                    help="AgamP5 molecule FASTA(s) (e.g. NC_064601.1)")
    ap.add_argument("--out-json",
                    default="results/cross_assembly_verify.json")
    ap.add_argument("--out-md", default="results/cross_assembly_verify.md")
    args = ap.parse_args()
    units = plus_units()
    gunits = guide_units()
    molecules = {}
    for path in args.agamp5:
        with open(path) as fh:
            header = fh.readline().strip().lstrip(">")
        molecules[header.split()[0]] = (path, load_fasta(path))
    report = {}
    for g, unit in units.items():
        per_mol = {}
        for mol, (path, genome) in molecules.items():
            per_mol[mol] = _search_one_molecule(genome, unit, gunits[g])
        report[g] = {
            "unit_plus": unit,
            "guide_unit_plus": gunits[g],
            "molecules": per_mol,
        }
    _write_outputs(args, molecules, report)


def _agg(mol_dict):
    """Genome-wide totals across molecules for one guide."""
    t = {"ref_exact": 0, "guide_exact": 0,
         "guide_mm1": 0, "guide_mm2": 0, "guide_mm3": 0}
    for m in mol_dict.values():
        t["ref_exact"] += len(m["exact_plus"]) + len(m["exact_minus"])
        t["guide_exact"] += (len(m["guide_exact_plus"]) +
                             len(m["guide_exact_minus"]))
        gn = m["guide_near"]
        for k in (1, 2, 3):
            t[f"guide_mm{k}"] += (len(gn["plus"][str(k)]) +
                                   len(gn["minus"][str(k)]))
    return t


def _write_outputs(args, molecules, report):
    with open(args.out_json, "w") as fh:
        json.dump({"method": "exact + pigeonhole near-match (4 blocks, "
                             "<=3 mismatches), both strands, no aligner",
                   "molecules": {m: v[1] and len(v[1])
                                 for m, v in molecules.items()},
                   "agamp5_query": args.agamp5, "guides": report},
                  fh, indent=1)
    lines = ["## Cross-assembly verification: guide targets on AgamP5", "",
             "Exact and <=3-mismatch matches of each 23-nt protospacer+PAM "
             "unit (plus-strand orientation) on both strands of every "
             "AgamP5 molecule (NC_064601.1 chr2, OX030908.1 chr3, "
             "OX030909.1 chrX, OX030910.2 MT). Method: exact substring + "
             "pigeonhole near-match (4 blocks), no aligner.", "",
             "AgamP4-reference units and guide-derived units searched "
             "separately: AgamP4 can carry a minor allele (assembly "
             "mismatch); the guide unit is the intended target allele.", "",
             "| guide | AgamP4 unit exact | guide unit exact | "
             "guide mm1 | guide mm2 | guide mm3 |",
             "|---|---|---|---|---|---|---|"]
    for g, r in report.items():
        t = _agg(r["molecules"])
        lines.append(
            f"| {g} | {t['ref_exact']} | {t['guide_exact']} "
            f"| {t['guide_mm1']} | {t['guide_mm2']} | {t['guide_mm3']} |")
    lines += ["", "### Per-molecule exact hits (plus-strand 1-based = "
              "0-based index + 1)", "",
              "| guide | molecule | ref exact + | guide exact + |", 
              "|---|---|---|---|"]
    for g, r in report.items():
        for mol, m in r["molecules"].items():
            rp = [x + 1 for x in m["exact_plus"]]
            gp = [x + 1 for x in m["guide_exact_plus"]]
            if rp or gp:
                lines.append(f"| {g} | {mol} | {rp} | {gp} |")
    with open(args.out_md, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(json.dumps({g: _agg(r["molecules"])
                      for g, r in report.items()}, indent=1))


if __name__ == "__main__":
    main()
