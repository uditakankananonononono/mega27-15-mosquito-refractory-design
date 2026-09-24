"""Genome-wide miniprot DM-domain census (tool 33 extension) with memory engineering.

Whole-genome and default-sensitivity per-contig miniprot runs OOM-kill this 2 GB
sandbox. Per-contig -M2 (1/4 modimiser sampling) fits (peaks 0.034-0.943 GB) and
recovers the full dsx gene structure on chr2; -M3 (1/8) loses even the known dsx
locus. Both HMMER-detected DM-family paralogs (dmd-4, DmrtA2) are missed by every
miniprot configuration tried, including full-sensitivity locus-level runs with
loose chaining; native-protein controls map both loci at 100% identity, so the
miss is a documented seed-chaining sensitivity limit, not absence.

Commands (executed; regenerable):
  seqkit/awk split of /tmp/gcf_genomic.fna -> /tmp/contig_*.fa, /tmp/gcf_scaffolds.fa
  miniprot -t1 -M2 --gff /tmp/contig_<acc>.fa /tmp/dsx_drome.fasta   # per contig
  miniprot -t1 -M2 --gff /tmp/gcf_scaffolds.fa /tmp/dsx_drome.fasta  # 187 scaffolds
  miniprot -t1 --gff /tmp/locus_<gene>.fa /tmp/dsx_drome.fasta       # locus, default
  miniprot -t1 -n1 -m20 /tmp/locus_<gene>.fa /tmp/dsx_drome.fasta    # locus, loose
  miniprot -t1 --gff /tmp/locus_<gene>.fa /tmp/q_<gene>.faa          # native control

Outputs: results/miniprot_census.json (+ .md companion, hand-written)
"""
import json, subprocess, os, re

MP = "/tmp/miniprot-src/miniprot"
DSX = "/tmp/dsx_drome.fasta"
CONTIGS = {"NC_064601.1": "chr2", "NC_064602.1": "chr3",
           "NC_064600.1": "chrX", "NC_083487.1": "mito"}
# contig fastas for the new accessions come from the complete GCF file;
# old-accession files (OX*) are byte-identical and reused for X/3/mito.
CONTIG_FA = {"NC_064601.1": "/tmp/contig_NC_064601.1.fa",
             "NC_064602.1": "/tmp/contig_OX030908.1.fa",
             "NC_064600.1": "/tmp/contig_OX030909.1.fa",
             "NC_083487.1": "/tmp/contig_OX030910.2.fa"}

def run(args):
    r = subprocess.run(args, capture_output=True, text=True)
    peak = None
    m = re.search(r"Peak RSS: ([0-9.]+) GB", r.stderr)
    if m: peak = float(m.group(1))
    return r.stdout, peak

def gff_mrnas(gff):
    out = []
    for line in gff.splitlines():
        if line.startswith("#"): continue
        f = line.split("\t")
        if len(f) > 8 and f[2] == "mRNA":
            out.append({"seqid": f[0], "start": int(f[3]), "end": int(f[4]),
                        "score": f[5], "strand": f[6], "attrs": f[8]})
    return out

def main():
    scans = {}
    for acc, label in CONTIGS.items():
        gff, peak = run([MP, "-t1", "-M2", "--gff", CONTIG_FA[acc], DSX])
        scans[acc] = {"chromosome": label, "mode": "-M2 (1/4 modimiser sampling)",
                      "peak_rss_gb": peak, "mrna_features": gff_mrnas(gff)}
    gff_s, peak_s = run([MP, "-t1", "-M2", "--gff", "/tmp/gcf_scaffolds.fa", DSX])
    scans["NW_scaffolds_187"] = {"chromosome": "unplaced (187 scaffolds, 19.0 Mb)",
                                 "mode": "-M2", "peak_rss_gb": peak_s,
                                 "mrna_features": gff_mrnas(gff_s)}

    loci = {}
    for gene, acc_region in [("dmd4", "NC_064601.1:4000000-4200000"),
                             ("dmrta2", "NC_064600.1:7600000-7800000")]:
        gff_d, _ = run([MP, "-t1", "--gff", f"/tmp/locus_{gene}.fa", DSX])
        paf_l, _ = run([MP, "-t1", "-n1", "-m20", f"/tmp/locus_{gene}.fa", DSX])
        gff_n, _ = run([MP, "-t1", "--gff", f"/tmp/locus_{gene}.fa", f"/tmp/q_{gene}.faa"])
        loci[gene] = {
            "locus_window": acc_region,
            "dsx_query_default_sensitivity_mrnas": gff_mrnas(gff_d),
            "dsx_query_loose_chaining_paf_lines": len([l for l in paf_l.splitlines() if l.strip()]),
            "native_protein_control_mrnas": gff_mrnas(gff_n),
        }

    out = {
        "tool": "miniprot 0.18-r281 (lh3/miniprot 81f9b93; Li 2023 Bioinformatics 39:btad014)",
        "query": "D. melanogaster Dsx (UniProt P23023)",
        "reference": "complete GCF_943734735.2 idAnoGambNW_F1_1, per-contig",
        "memory_engineering": {
            "whole_genome": "OOM-killed twice at default settings (2 GB sandbox)",
            "chr2_default_M1": "OOM-killed during index build (923,414 blocks)",
            "M2": "fits: peaks 0.034-0.943 GB per contig",
            "M3": "fits trivially (0.613 GB chr2) but loses ALL hits incl. known dsx locus - sampling sensitivity floor",
        },
        "per_contig_scans": scans,
        "paralog_locus_checks": {
            "dmd4_annotation": "LOC1281789, NC_064601.1:4056166-4057948 (chr2)",
            "dmrta2_annotation": "LOC1271814, NC_064600.1:7711435-7726245 (chrX)",
            **loci,
        },
        "verdict": "miniprot genome census with cross-species Dsx query detects 1/3 DM-family genes (dsx only, full 8-CDS structure on chr2); dmd-4 and DmrtA2 missed at every configuration incl. full-sensitivity locus runs - documented seed-chaining sensitivity limit, controls map both loci at 100% identity. HMMER (profile HMM) remains the 3/3 family census.",
    }
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/miniprot_census.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k in ("per_contig_scans", "verdict")}, indent=1)[:1500])

if __name__ == "__main__":
    main()
