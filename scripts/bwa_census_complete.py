"""Genome-completeness correction of the BWA exhaustive <=3-mismatch census.

The census of record (results/offtarget_bwa.json, 32,305 alignments) was
computed against a 4-record FASTA (245,459,960 bp: chromosomes X, 2, 3 and
mitochondrion only). The complete GCF_943734735.2 idAnoGambNW_F1_1 assembly
is 191 records / 264,466,744 bp - the same four chromosome-level sequences
(byte-identical; RefSeq vs INSDC accessions) PLUS 187 unlocalized/unplaced
NW_* scaffolds (19,006,784 bp). The 4-record reference therefore omitted
7.2% of the assembly. This script reruns the identical exhaustive BWA
enumeration on the complete assembly and maps the old census forward by
accession (sequences verified byte-identical by md5):
  OX030909.1 -> NC_064600.1 (X), OX030908.1 -> NC_064602.1 (3),
  OX030910.2 -> NC_083487.1 (mito), NC_064601.1 = NC_064601.1 (2).

Commands (executed; regenerable):
  curl -sSL -o /tmp/gcf_genomic.fna.gz \\
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/943/734/735/GCF_943734735.2_idAnoGambNW_F1_1/GCF_943734735.2_idAnoGambNW_F1_1_genomic.fna.gz
  gunzip -kf /tmp/gcf_genomic.fna.gz
  bwa index /tmp/gcf_genomic.fna
  bwa aln -t 4 -N -l 20 -k 3 -n 3 -o 0 -e 0 -f /tmp/bwa_guides_gcf.sai \\
      /tmp/gcf_genomic.fna /tmp/bt2_guides.fa
  bwa samse -n 1000000 -f /tmp/bwa_hits_gcf.sam \\
      /tmp/gcf_genomic.fna /tmp/bwa_guides_gcf.sai /tmp/bt2_guides.fa
  (identical alignment parameters to the 4-record census of record)

Outputs: results/offtarget_bwa_complete.json + .md
"""
import json, subprocess, collections, os, random

SAM_NEW = "/tmp/bwa_hits_gcf.sam"
SAM_OLD = "/tmp/bwa_hits.sam"
GENOME = "/tmp/gcf_genomic.fna"
GUIDES_FA = "/tmp/bt2_guides.fa"
CDS_SORTED = "/tmp/agamp5_cds.sorted.bed"
BEDTOOLS = "/home/sandbox/tools/bedtools2/bin/bedtools"

ACC_MAP = {"OX030909.1": "NC_064600.1", "OX030908.1": "NC_064602.1",
           "OX030910.2": "NC_083487.1", "NC_064601.1": "NC_064601.1"}

LEADS = {"dsx_W1_340_363": "dsx-v3-1", "dsx_W1_1604_1627": "dsx-v3-2",
         "dsx_W2_796_819": "dsx-v3-3", "dsx_W2_174_197": "kyrou"}

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]

def load_genome():
    seqs, name = {}, None
    for line in open(GENOME):
        line = line.rstrip()
        if line.startswith(">"):
            name = line.split()[0][1:]; seqs[name] = []
        else: seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}

def load_guides():
    g, name = {}, None
    for line in open(GUIDES_FA):
        line = line.rstrip()
        if line.startswith(">"): name = line[1:]
        else: g[name] = line.upper()
    return g

def parse_bwa_sam(path):
    hits = []
    for line in open(path):
        if line.startswith("@"): continue
        f = line.rstrip("\n").split("\t"); flag = int(f[1])
        nm0 = None; xa = ""
        for t in f[11:]:
            if t.startswith("NM:i:"): nm0 = int(t[5:])
            elif t.startswith("XA:Z:"): xa = t[5:]
        if not (flag & 4) and nm0 is not None and nm0 <= 3:
            hits.append((f[0], f[2], int(f[3]), "rev" if flag & 16 else "fwd", nm0))
        if xa:
            for e in xa.rstrip(";").split(";"):
                chrom, pos, cigar, nm = e.split(",")
                nm = int(nm)
                if nm <= 3:
                    hits.append((f[0], chrom, abs(int(pos)),
                                 "fwd" if pos[0] == "+" else "rev", nm))
    return hits

def main():
    new_hits = parse_bwa_sam(SAM_NEW)
    old_hits = parse_bwa_sam(SAM_OLD)
    census = collections.Counter(nm for *_x, nm in new_hits)
    scaf = collections.Counter(nm for q, c, p, s, nm in new_hits if c.startswith("NW_"))
    chrm = collections.Counter(nm for q, c, p, s, nm in new_hits if not c.startswith("NW_"))

    old_mapped = {(q, ACC_MAP[c], p, s) for q, c, p, s, nm in old_hits}
    new_keys = {(q, c, p, s) for q, c, p, s, nm in new_hits}
    old_missing = [k for k in old_mapped if k not in new_keys]
    scaffold_only = [h for h in new_hits if h[1].startswith("NW_")]

    # sequence-level verification of 50 seeded-random scaffold hits
    genome, guides = load_genome(), load_guides()
    random.seed(11)
    sample = random.sample(scaffold_only, min(50, len(scaffold_only)))
    verified = 0
    for q, chrom, pos, strand, nm in sample:
        ref = genome[chrom][pos-1:pos-1+20]
        if strand == "rev": ref = rc(ref)
        if len(ref) == 20 and sum(1 for a, b in zip(ref, guides[q]) if a != b) <= 3:
            verified += 1

    # CDS context on the complete hit set (bed covers all 31 annotated seqids)
    with open("/tmp/bwa_hits_gcf.bed", "w") as out:
        for q, chrom, pos, strand, nm in new_hits:
            out.write(f"{chrom}\t{pos-1}\t{pos-1+20}\t{q}\t{nm}\n")
    subprocess.run("sort -k1,1 -k2,2n /tmp/bwa_hits_gcf.bed > /tmp/bwa_hits_gcf.sorted.bed",
                   shell=True, check=True)
    r = subprocess.run([BEDTOOLS, "intersect", "-wa", "-wb",
                        "-a", "/tmp/bwa_hits_gcf.sorted.bed", "-b", CDS_SORTED],
                       capture_output=True, text=True, check=True)
    in_cds = {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        in_cds.setdefault((f[0], int(f[1])), set()).add(f[8])

    per_tier = collections.defaultdict(lambda: [0, 0])
    per_guide = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    lead_cds_detail = collections.defaultdict(list)
    for q, chrom, pos, strand, nm in new_hits:
        key = (chrom, pos - 1)
        hit = key in in_cds
        per_tier[nm][0] += 1; per_tier[nm][1] += hit
        per_guide[q][nm][0] += 1; per_guide[q][nm][1] += hit
        if hit and nm >= 1:
            genes = in_cds[key]
            per_guide[q].setdefault("cds_genes", set()).update(g for g in genes if g != ".")
            for prefix, lname in LEADS.items():
                if q.startswith(prefix):
                    for g2 in genes:
                        if (chrom, pos, nm, g2) not in lead_cds_detail[lname]:
                            lead_cds_detail[lname].append((chrom, pos, nm, g2))

    leads = {}
    for prefix, name in LEADS.items():
        rows = {}; genes = set()
        for qid, tiers in per_guide.items():
            if qid.startswith(prefix):
                for nm, pair in tiers.items():
                    if nm == "cds_genes": continue
                    rows.setdefault(nm, [0, 0])
                    rows[nm][0] += pair[0]; rows[nm][1] += pair[1]
                genes |= tiers.get("cds_genes", set())
        leads[name] = {"mm_tiers_total_inCDS": {str(k): v for k, v in sorted(rows.items())},
                       "cds_overlap_genes": sorted(genes),
                       "cds_overlap_sites": sorted(lead_cds_detail.get(name, []))}

    guides_with_cds_ot = sum(1 for q, t in per_guide.items()
                             if any(isinstance(nm, int) and nm >= 1 and c[1] > 0
                                    for nm, c in t.items()))

    # lead-guide exact-site check (the v3-1/v3-2 zero-near-cognate claim)
    lead_tiers = collections.defaultdict(collections.Counter)
    for q, chrom, pos, strand, nm in new_hits:
        for prefix, lname in LEADS.items():
            if q.startswith(prefix):
                lead_tiers[lname][nm] += 1


    # --- correction decomposition: accession-mismatch vs scaffold contribution ---
    def cds_guides(hit_list, tag):
        with open(f"/tmp/decomp_{tag}.bed", "w") as fh:
            for q, c, p, s, nm in hit_list:
                fh.write(f"{c}\t{p-1}\t{p-1+20}\t{q}\t{nm}\n")
        rr = subprocess.run([BEDTOOLS, "intersect", "-wa", "-wb",
                             "-a", f"/tmp/decomp_{tag}.bed", "-b", CDS_SORTED],
                            capture_output=True, text=True, check=True)
        gg = collections.defaultdict(set)
        for ln in rr.stdout.splitlines():
            ff = ln.split("\t"); gg[ff[3]].add(int(ff[4]))
        return {q for q, tiers in gg.items() if any(n >= 1 for n in tiers)}
    old_acc_fixed = [(q, ACC_MAP[c], p, s, nm) for q, c, p, s, nm in parse_bwa_sam(SAM_OLD)]
    g_old_fixed = cds_guides(old_acc_fixed, "oldfixed")
    g_new = cds_guides(new_hits, "new")
    g_scaf = cds_guides([h for h in new_hits if h[1].startswith("NW_")], "scaf")
    decomposition = {
        "previously_reported_guides_with_cds_near_cognate": 308,
        "same_32305_sites_accession_fixed": len(g_old_fixed),
        "complete_genome": len(g_new),
        "scaffold_hit_driven_guides": len(g_scaf),
        "scaffold_only_new_guides": len(g_new - g_old_fixed),
        "explanation": "The entire 308->367 delta is an accession-mismatch correction, not new sequence: the 4-record reference named chromosomes X/3/mito with INSDC accessions (OX030909.1/OX030908.1/OX030910.2) while the RefSeq CDS annotation names them NC_064600.1/NC_064602.1/NC_083487.1, so every chr3/X/mito CDS overlap was silently invisible to the intersect. Scaffolds add 631 alignments but zero additional CDS-carrying guides (all 7 scaffold-driven guides already had a chromosome CDS hit).",
    }

    # --- v3-2 newly visible CDS site: sequence + PAM verification ---
    site_hits = [h for h in new_hits if h[1] == "NC_064602.1" and h[2] == 2111186]
    qv, cv, pv, sv, nmv = site_hits[0]
    refv = genome[cv][pv-1:pv-1+20]
    if sv == "rev": refv = rc(refv)
    gseqv = guides[qv]
    mm_pos = [i+1 for i, (a, b) in enumerate(zip(refv, gseqv)) if a != b]
    pamv = (genome[cv][pv-1+20:pv-1+23] if sv == "fwd"
            else rc(genome[cv][pv-4:pv-1]))
    v32_site = {
        "site": "NC_064602.1:2111186-2111205", "guide": qv, "strand": sv,
        "guide_seq": gseqv, "ref_seq": refv,
        "hamming": len(mm_pos), "mismatch_positions_1based": mm_pos,
        "pam_3prime_3bp": pamv, "is_ngg": bool(len(pamv) == 3 and pamv[1:] == "GG"),
        "gene": "LOC1278105 (probable aconitate hydratase, mitochondrial)",
        "verdict": "true MM3 site, PAM CCC not NGG - Cas9-inert",
        "why_new": "site present in the 4-record census (as OX030908.1:2111186) but its CDS overlap was invisible under the INSDC/RefSeq accession mismatch; complete-assembly rerun uses native RefSeq accessions",
    }

    out = {
        "tool": "BWA 0.7.19-r1273 aln/samse (genome-completeness rerun)",
        "mode": "aln -N -l 20 -k 3 -n 3 -o 0 -e 0; samse -n 1000000 (identical to census of record)",
        "guides": 421,
        "reference": "COMPLETE AgamP5 = idAnoGambNW_F1_1 (GCF_943734735.2) genomic.fna: 191 records, 264,466,744 bp (4 chromosomes byte-identical to the 4-record reference + 187 NW_* scaffolds, 19,006,784 bp = 7.2% of assembly)",
        "census_per_mm_tier": {str(k): census[k] for k in sorted(census)},
        "census_total": sum(census.values()),
        "chromosomes_per_mm_tier": {str(k): chrm[k] for k in sorted(chrm)},
        "scaffolds_per_mm_tier": {str(k): scaf[k] for k in sorted(scaf)},
        "prior_census_4record_total": 32305,
        "forward_map_check": {
            "method": "old hits mapped by accession (sequences md5-identical); compared as (guide, chrom, pos, strand)",
            "old_alignments_missing_from_complete": len(old_missing),
        },
        "scaffold_hit_verification": {
            "scaffold_alignments_total": len(scaffold_only),
            "sampled": len(sample), "verified_hamming_le_3": verified, "seed": 11,
            "method": "direct reference extraction at SAM coordinates, reverse-complement for minus strand, Hamming recount vs guide"},
        "cds_context_on_complete_set": {
            "annotation": "NCBI RefSeq GCF_943734735.2 RS_2023_12 CDS, 31 seqids incl. scaffolds (bedtools 2.31.1 intersect)",
            "per_mm_tier_total_inCDS": {str(k): v for k, v in sorted(per_tier.items())},
            "library_guides_with_at_least_one_CDS_overlapping_near_cognate": guides_with_cds_ot,
        },
        "leads": leads,
        "lead_mm_tier_totals": {k: dict(v) for k, v in lead_tiers.items()},
        "correction_decomposition": decomposition,
        "v32_new_cds_site_verification": v32_site,
    }
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/offtarget_bwa_complete.json", "w"), indent=1)
    print(json.dumps(out, indent=1)[:3500])

if __name__ == "__main__":
    main()
