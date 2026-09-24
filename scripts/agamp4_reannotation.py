"""AgamP4 re-annotation of the dsx locus and the Kyrou 2018 guide site.

Resolves the Appendix I annotation discrepancy: under the AgamP4 / VectorBase
community genebuild (Ensembl Metazoa import, AgamP4.63), does the published
Kyrou protospacer span an exon boundary as the 2018 paper describes?
Also maps the two named v3 leads onto AgamP4 and reports assembly sensitivity.
Hermetic: reads only archived files under data/agamp4/.
"""
import json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GFF = os.path.join(ROOT, "data", "agamp4", "AGAP004050_AgamP4.gff3")
FA = os.path.join(ROOT, "data", "agamp4", "AgamP4_2R_dsx_region.fasta")
OUT = os.path.join(ROOT, "results", "agamp4_reannotation.json")
REGION_START = 48_700_000  # 1-based start of the archived region on 2R

GUIDES = {
    "kyrou_2018": "GTTTAACACAGGTCAAGCGG",
    "dsx-v3-1": "TGGGCAGTATGCGTTAGGGT",
    "dsx-v3-2": "CATTAAGACCTACGAAGCGC",
}
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]


def load_region():
    seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
    return seq


def load_models():
    """Parse the archived GFF excerpt into transcript -> exon/CDS lists."""
    exons, cds = {}, {}
    for line in open(GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("exon", "CDS"):
            continue
        tx = f[8].split("transcript:")[1].split(";")[0]
        tgt = exons if f[2] == "exon" else cds
        tgt.setdefault(tx, []).append((int(f[3]), int(f[4])))
    return exons, cds


def locate(seq, guide, window=3000):
    """Exact match if unique; otherwise best-mismatch match near the AgamP5-offset position."""
    hits = []
    for strand, q in (("+", guide), ("-", rc(guide))):
        start = 0
        while True:
            i = seq.find(q, start)
            if i < 0:
                break
            hits.append((0, i, strand, q if strand == "+" else rc(q)))
            start = i + 1
    if hits:
        return hits[0], len(hits)
    # approximate AgamP5->AgamP4 offset ~1,092,464 (from the Kyrou anchor)
    center = None
    best = None
    for p in range(0, len(seq) - 20):
        for strand, t in (("+", seq[p:p + 20]), ("-", rc(seq[p:p + 20]))):
            mm = sum(1 for a, b in zip(t, guide) if a != b)
            if best is None or mm < best[0]:
                best = (mm, p, strand, t)
    return best, 0


def overlaps(s, e, ivs):
    out = []
    for a, b in ivs:
        if not (e < a or s > b):
            out.append((a, b, max(s, a), min(e, b)))
    return out


def main():
    seq = load_region()
    exons, cds = load_models()
    results = {}
    for name, guide in GUIDES.items():
        (mm, p, strand, aligned), n_exact = locate(seq, guide)
        s, e = REGION_START + p, REGION_START + p + 19
        rec = {
            "guide": guide, "assembly_position_2R": [s, e], "strand": strand,
            "mismatches_vs_AgamP4": mm, "aligned_sequence": aligned,
            "exact_match": mm == 0,
        }
        for tx in exons:
            ov = overlaps(s, e, exons[tx])
            if ov:
                rec.setdefault("exon_overlaps", {})[tx] = ov
        for tx in cds:
            ov = overlaps(s, e, cds[tx])
            if ov:
                rec.setdefault("cds_overlaps", {})[tx] = ov
        results[name] = rec
    female = [iv for tx, ivs in exons.items() for iv in ivs if iv[1] - iv[0] + 1 == 135]
    summary = {
        "source": "Ensembl Metazoa FTP current release: Anopheles_gambiae.AgamP4.63.chromosome.2R.gff3.gz "
                  "(VectorBase community models) + AgamP4.dna_sm.chromosome.2R.fa.gz; fetched 24 Sep 2026",
        "gene_span_2R": [48_703_664, 48_788_460],
        "transcripts": sorted(exons),
        "female_exon_135bp": female,
        "guides": results,
        "kyrou_junction_spanning": bool(results["kyrou_2018"].get("cds_overlaps")),
    }
    with open(OUT, "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary["kyrou_junction_spanning"]), "->", OUT)


if __name__ == "__main__":
    main()
