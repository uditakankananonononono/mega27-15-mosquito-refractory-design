"""Functional context of bowtie2 near-cognate sites (bedtools 2.31.1, tool 30/40).

Question: bowtie2 enumerated 24,278 ungapped <=3-mismatch alignments of all
421 dsx-window guides against AgamP5 (PAM-agnostic). How many of those
near-cognate sites overlap annotated protein-coding sequence - i.e. carry
functional off-target potential if a PAM happened to flank them - and does the
answer change the safety picture for the three leads and the Kyrou control?

Inputs (regenerable):
  /tmp/bt2_hits2.sam     - of-record bowtie2 alignments (scripts/bowtie2_run.py)
  /tmp/agamp5.gff.gz     - NCBI RefSeq annotation of AgamP5 = idAnoGambNW_F1_1
                           (GCF_943734735.2, release RS_2023_12), fetched from
                           ftp.ncbi.nlm.nih.gov/genomes/all/GCF/943/734/735/

Method: SAM -> BED (CIGAR M/D/N span, NM tag as mismatch count); GFF CDS
features -> BED (1-based inclusive -> 0-based half-open); bedtools intersect
-sorted -wa -wb. A hit is "in CDS" if it overlaps >=1 CDS base. self-hits
(NM=0 at the guide's own locus) are reported separately from off-targets.

Outputs: results/offtarget_context_bedtools.json + .md
"""
import json, re, subprocess, collections, os

BEDTOOLS = "/home/sandbox/tools/bedtools2/bin/bedtools"
SAM = "/tmp/bt2_hits2.sam"
GFF = "/tmp/agamp5.gff.gz"
HITS_BED = "/tmp/bt2_hits.sorted.bed"
CDS_BED = "/tmp/agamp5_cds.bed"

LEADS = {  # guide id prefix in bt2_guides.fa -> display name
    "dsx_W1_340_363": "dsx-v3-1", "dsx_W1_1604_1627": "dsx-v3-2",
    "dsx_W2_796_819": "dsx-v3-3", "dsx_W2_174_197": "kyrou",
}

def build_beds():
    import gzip
    n_cds = 0
    with gzip.open(GFF, "rt") as fh, open(CDS_BED, "w") as out:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            m = re.search(r"gene=([^;]+)", f[8])
            out.write(f"{f[0]}\t{int(f[3])-1}\t{f[4]}\t{m.group(1) if m else '.'}\n")
            n_cds += 1
    n = 0
    with open(SAM) as fh, open("/tmp/bt2_hits.bed", "w") as out:
        for line in fh:
            if line.startswith("@"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or f[2] == "*": continue
            nm = [t for t in f[11:] if t.startswith("NM:i:")]
            nm = int(nm[0][5:]) if nm else -1
            span = sum(int(x) for x, l in re.findall(r"(\d+)([MDN])", f[5]))
            s0 = int(f[3]) - 1
            out.write(f"{f[2]}\t{s0}\t{s0+span}\t{f[0]}\t{nm}\n"); n += 1
    subprocess.run(f"sort -k1,1 -k2,2n /tmp/bt2_hits.bed > {HITS_BED}",
                   shell=True, check=True)
    return n, n_cds

def main():
    n_aln, n_cds = build_beds()
    # CDS BED must be sorted for -sorted
    subprocess.run(f"sort -k1,1 -k2,2n {CDS_BED} > /tmp/agamp5_cds.sorted.bed",
                   shell=True, check=True)
    r = subprocess.run([BEDTOOLS, "intersect", "-sorted", "-wa", "-wb",
                        "-a", HITS_BED, "-b", "/tmp/agamp5_cds.sorted.bed"],
                       capture_output=True, text=True, check=True)
    in_cds = {}   # (contig,start) -> set of genes
    for line in r.stdout.splitlines():
        f = line.split("\t")
        in_cds.setdefault((f[0], int(f[1])), set()).add(f[8])

    # walk the hits BED for full accounting
    per_tier = collections.defaultdict(lambda: [0, 0])  # nm -> [total, in_cds]
    per_guide = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    for line in open(HITS_BED):
        c, s, e, qid, nm = line.split("\t")
        nm = int(nm)
        key = (c, int(s))
        hit_cds = key in in_cds
        genes = sorted(in_cds.get(key, ()))
        per_tier[nm][0] += 1
        per_tier[nm][1] += hit_cds
        per_guide[qid][nm][0] += 1
        per_guide[qid][nm][1] += hit_cds
        if hit_cds and nm >= 1:
            per_guide[qid].setdefault("cds_genes", set()).update(g for g in genes if g != ".")

    leads = {}
    for prefix, name in LEADS.items():
        rows = {}
        for qid, tiers in per_guide.items():
            if qid.startswith(prefix):
                for nm, pair in tiers.items():
                    if nm == "cds_genes": continue
                    tot, cds = pair
                    rows.setdefault(nm, [0, 0])
                    rows[nm][0] += tot; rows[nm][1] += cds
                genes = tiers.get("cds_genes", set())
                leads[name] = {"mm_tiers_total_inCDS": {str(k): v for k, v in sorted(rows.items())},
                               "cds_overlap_genes": sorted(genes) if isinstance(genes, set) else []}

    guides_with_cds_ot = sum(1 for q, t in per_guide.items()
                             if any(isinstance(nm, int) and nm >= 1 and c[1] > 0
                                    for nm, c in t.items()))
    # PAM-adjacency resolution for lead CDS-overlapping near-cognates:
    # both are minus-strand alignments; an SpCas9 PAM would be NGG on the
    # guide strand. Checked both flanks in AgamP5 (data via Biopython):
    #   dsx-v3-3 MM3 chr2:38,008,162-38,008,181 in LOC1270942 (VPS13B):
    #     flanks GCG|GTG - no NGG/CCN -> Cas9-inert
    #   kyrou    MM3 chr2:41,285,163-41,285,182 in LOC1274366 (uncharacterized):
    #     flanks GTT|TAC - no NGG/CCN -> Cas9-inert
    lead_cds_pam_resolution = {
        "dsx-v3-3": {"site": "NC_064601.1:38008162-38008181", "mm": 3,
                     "gene": "LOC1270942 (VPS13B, intermembrane lipid transfer protein)",
                     "plus_flanks": "GCG|GTG", "pam": None,
                     "verdict": "no adjacent NGG/CCN - Cas9-inert"},
        "kyrou": {"site": "NC_064601.1:41285163-41285182", "mm": 3,
                  "gene": "LOC1274366 (uncharacterized)",
                  "plus_flanks": "GTT|TAC", "pam": None,
                  "verdict": "no adjacent NGG/CCN - Cas9-inert"},
    }
    out = {
        "tool": "bedtools 2.31.1 intersect (Quinlan & Hall 2010 Bioinformatics 26:841)",
        "lead_cds_hits_pam_resolution": lead_cds_pam_resolution,
        "annotation": "NCBI RefSeq GCF_943734735.2 idAnoGambNW_F1_1 (AgamP5), RS_2023_12",
        "cds_features": n_cds, "alignments": n_aln,
        "per_mm_tier_total_inCDS": {str(k): v for k, v in sorted(per_tier.items())},
        "library_guides_with_at_least_one_CDS_overlapping_near_cognate": guides_with_cds_ot,
        "leads": leads,
    }
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/offtarget_context_bedtools.json", "w"), indent=1)
    print(json.dumps(out, indent=1)[:2000])

if __name__ == "__main__":
    main()
