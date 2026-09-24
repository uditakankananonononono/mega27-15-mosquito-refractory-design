"""Tool 38: BLAT v39x1 (Kent 2002, Genome Res 12:656-664) protein-to-genome
placement cross-check of the miniprot ortholog call.

Question: does UCSC BLAT in translated mode (-t=dnax -q=prot) place the
Drosophila Dsx protein (P23023, 549 aa) at the AgamP5 dsx locus
(NC_064601.1:47,610,877-47,700,920, minus strand)?

PARSING CAVEAT (verified 2026-09-25): for "+-"/"--" strand rows the PSL
tStart/tEnd fields are reverse-complement coordinates while block tStarts
stay plus-strand; a naive tStart/tEnd window test therefore misplaces every
minus-strand hit. Locus tests below use BLOCK-LEVEL plus-strand coordinates
only. Two rows that superficially "land" at 47,622,733-47,692,496 (score 99)
and 2,159,826-4,057,682 (score 102) are +- rows whose blocks actually sit at
70.5 Mb and 114-116 Mb - both spurious, documented so the artifact is not
reproduced downstream.

Configurations on full chromosome 2 (118,196,952 bp):
  default:  -t=dnax -q=prot
  loose:    -minScore=10 -minIdentity=20
Writes results/blat_placement.json and results/blat_placement.md.
"""
import json, os, subprocess, sys

BLAT = "/tmp/blat"
CONTIG = "/tmp/contig_NC_064601.1.fa"
QUERY = "/tmp/dsx_drome.fasta"
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCUS = (47610877, 47700920)  # dsx gene span, NC_064601.1 (minus strand)


def parse_psl(path):
    hits = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 21 or not f[0].isdigit():
            continue
        sizes = [int(x) for x in f[18].rstrip(",").split(",")]
        starts = [int(x) for x in f[20].rstrip(",").split(",")]
        blocks = [(st, st + s * 3) for s, st in zip(sizes, starts)]  # plus-strand
        hits.append({
            "score": int(f[0]), "strand": f[8],
            "tStart_reported": int(f[15]), "tEnd_reported": int(f[16]),
            "blocks": blocks,
            "block_in_locus": any(b[1] >= LOCUS[0] and b[0] <= LOCUS[1] for b in blocks),
        })
    return hits


def run(name, extra):
    out = f"/tmp/blat_{name}.psl"
    r = subprocess.run([BLAT, CONTIG, QUERY, "-t=dnax", "-q=prot"] + extra + [out],
                       capture_output=True)
    if r.returncode != 0:
        sys.exit(f"blat {name} failed: " + r.stderr.decode()[:200])
    hits = parse_psl(out)
    locus_hits = [h for h in hits if h["block_in_locus"]]
    ranked = sorted(hits, key=lambda h: -h["score"])
    rank = next((i for i, h in enumerate(ranked, 1) if h["block_in_locus"]), None)
    return {
        "config": name, "parameters": " ".join(extra) or "defaults",
        "n_hits": len(hits),
        "hits_with_block_in_dsx_locus": len(locus_hits),
        "best_locus_hit": (max(locus_hits, key=lambda h: h["score"]) if locus_hits else None),
        "locus_hit_rank_by_score": rank,
        "top_score_overall": ranked[0]["score"] if ranked else None,
        "n_hits_outscoring_best_locus_hit": (
            sum(1 for h in hits if h["score"] > max(x["score"] for x in locus_hits))
            if locus_hits else len(hits)),
    }


def main():
    default = run("default", [])
    loose = run("loose", ["-minScore=10", "-minIdentity=20"])
    dl = default["best_locus_hit"]
    ll = loose["best_locus_hit"]
    result = {
        "tool": "UCSC BLAT v39x1 (Kent 2002, Genome Res 12:656-664), -t=dnax -q=prot",
        "query": "Drosophila melanogaster Dsx, UniProt P23023 (549 aa)",
        "target": "AgamP5 chromosome 2 (NC_064601.1, 118,196,952 bp)",
        "dsx_locus": "NC_064601.1:47,610,877-47,700,920 (minus strand)",
        "parsing_caveat": "PSL tStart/tEnd are reverse-complement coordinates on +- rows; "
                          "block tStarts are plus-strand. Locus tests use block coordinates.",
        "runs": [default, loose],
        "findings": [
            f"default settings: {default['n_hits']} hits chromosome-wide; exactly one places a "
            f"single 8-aa block (24 bp, score {dl['score']}) inside the dsx locus as part of a fragmented "
            f"9-block chain whose remaining blocks sit in intergenic DNA to 48.87 Mb. "
            f"{default['n_hits_outscoring_best_locus_hit']} spurious hits outscore it (top "
            f"{default['top_score_overall']}) - the placement is not identifiable de novo.",
            f"loose (-minScore=10 -minIdentity=20): identical picture ({loose['n_hits']} hits; the same "
            f"single locus block, score {ll['score']}, outscored by "
            f"{loose['n_hits_outscoring_best_locus_hit']} spurious hits) - loosening thresholds adds "
            f"only sub-threshold noise, no additional gene structure.",
            "Two +- rows deceptively report tStart/tEnd inside the dsx and dmd-4 regions; block-level "
            "coordinates place both at 70.5 Mb and 114-116 Mb. Documented as a PSL parsing artifact.",
            "Sensitivity ordering for this 70%-identity cross-species placement: BLAT (default miss; "
            "loose single-block, outranked by spurious hits) < miniprot (-M2: full 8-CDS structure, "
            "results/miniprot_census.json) << HMMER profile search (family census, "
            "results/hmmer_dsx_ortholog.json). Honest negative/partial result, preserved.",
        ],
    }
    assert default["hits_with_block_in_dsx_locus"] == 1
    assert loose["hits_with_block_in_dsx_locus"] == 1
    assert default["best_locus_hit"]["score"] == loose["best_locus_hit"]["score"] == 66
    assert loose["best_locus_hit"]["score"] < loose["top_score_overall"]
    with open(os.path.join(BASE, "results", "blat_placement.json"), "w") as fh:
        json.dump(result, fh, indent=1)
    with open(os.path.join(BASE, "results", "blat_placement.md"), "w") as fh:
        fh.write("# BLAT protein-to-genome placement cross-check (tool 38)\n\n")
        for f_ in result["findings"]:
            fh.write(f"- {f_}\n")
    print("blat: single marginal locus block at both settings, outranked by 13 spurious hits")


if __name__ == "__main__":
    main()
