"""Exhaustive ungapped <=3-mismatch census of 421 dsx guides vs AgamP5 with
BWA 0.7.19-r1273 (tool 31/40) + CDS context on the corrected hit set.

Motivation: CRISPOR's engine is BWA (Li & Durbin 2009 Bioinformatics
25:1754-60); running it first-party on our own data both replicates that
engine and cross-validates the bowtie2 census (tool 29/40). Result: BWA is
a strict superset - 0 of 24,278 bowtie2 alignments missing, 8,027 bwa-only
alignments (+142 MM2, +7,885 MM3). 50 seeded-random bwa-only hits were
verified by direct reference extraction: 50/50 are true Hamming <=3 sites.
The bowtie2 of-record census (24,278) therefore under-reports by 25%
(32,305 true); bowtie2's -a mode still score-trims high-mismatch ungapped
hits even after seed tuning. The corrected census is adopted; the bowtie2
table is retained as a documented under-count (negative result preserved).

Commands (executed; regenerable):
  bwa index /tmp/agam_genome.fasta            # 278.7 s
  bwa aln -t 4 -N -l 20 -k 3 -n 3 -o 0 -e 0 -f /tmp/bwa_guides.sai \
      /tmp/agam_genome.fasta /tmp/bt2_guides.fa
  bwa samse -n 1000000 -f /tmp/bwa_hits.sam \
      /tmp/agam_genome.fasta /tmp/bwa_guides.sai /tmp/bt2_guides.fa
  (-N = non-iterative: search for ALL n-difference hits; -l 20 full-length
   seed; -k 3 seed differences; -n 3 max differences; -o 0 -e 0 = no gaps,
   i.e. pure-Hamming enumeration; samse -n keeps every hit in XA.)

Outputs: results/offtarget_bwa.json + .md
"""
import json, re, subprocess, collections, os, random

SAM = "/tmp/bwa_hits.sam"
BT2_SAM = "/tmp/bt2_hits2.sam"
GENOME = "/tmp/agam_genome.fasta"
GUIDES_FA = "/tmp/bt2_guides.fa"
CDS_SORTED = "/tmp/agamp5_cds.sorted.bed"
BEDTOOLS = "/home/sandbox/tools/bedtools2/bin/bedtools"

LEADS = {"dsx_W1_340_363": "dsx-v3-1", "dsx_W1_1604_1627": "dsx-v3-2",
         "dsx_W2_796_819": "dsx-v3-3", "dsx_W2_174_197": "kyrou"}

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

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]

def parse_bwa_sam():
    """Return list of (qid, chrom, pos1, strand, nm) for all primary+XA hits."""
    hits = []
    for line in open(SAM):
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

def parse_bt2_set():
    s = set()
    for line in open(BT2_SAM):
        if line.startswith("@"): continue
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & 4: continue
        nm = [t for t in f[11:] if t.startswith("NM:i:")]
        if nm and int(nm[0][5:]) <= 3:
            s.add((f[0], f[2], int(f[3]), "rev" if int(f[1]) & 16 else "fwd"))
    return s

def main():
    hits = parse_bwa_sam()
    census = collections.Counter(nm for *_x, nm in hits)

    bt2 = parse_bt2_set()
    bwa_keys = {(q, c, p, s) for q, c, p, s, nm in hits}
    bt2_missing = [k for k in bt2 if k not in bwa_keys]
    bwa_only = [h for h in hits if (h[0], h[1], h[2], h[3]) not in bt2]
    delta_by_tier = collections.Counter(h[4] for h in bwa_only)

    # sequence-level verification of 50 seeded-random bwa-only hits
    genome, guides = load_genome(), load_guides()
    random.seed(7)
    sample = random.sample(bwa_only, 50)
    verified = 0
    for q, chrom, pos, strand, nm in sample:
        ref = genome[chrom][pos-1:pos-1+20]
        if strand == "rev": ref = rc(ref)
        if len(ref) == 20 and sum(1 for a, b in zip(ref, guides[q]) if a != b) <= 3:
            verified += 1

    # CDS context on the corrected hit set
    with open("/tmp/bwa_hits.bed", "w") as out:
        for q, chrom, pos, strand, nm in hits:
            out.write(f"{chrom}\t{pos-1}\t{pos-1+20}\t{q}\t{nm}\n")
    subprocess.run("sort -k1,1 -k2,2n /tmp/bwa_hits.bed > /tmp/bwa_hits.sorted.bed",
                   shell=True, check=True)
    r = subprocess.run([BEDTOOLS, "intersect", "-sorted", "-wa", "-wb",
                        "-a", "/tmp/bwa_hits.sorted.bed", "-b", CDS_SORTED],
                       capture_output=True, text=True, check=True)
    in_cds = {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        in_cds.setdefault((f[0], int(f[1])), set()).add(f[8])

    per_tier = collections.defaultdict(lambda: [0, 0])
    per_guide = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    lead_cds_detail = collections.defaultdict(list)
    for q, chrom, pos, strand, nm in hits:
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

    out = {
        "tool": "BWA 0.7.19-r1273 aln/samse (lh3/bwa d82444c, built from source; Li & Durbin 2009 Bioinformatics 25:1754-60)",
        "mode": "aln -N -l 20 -k 3 -n 3 -o 0 -e 0: exhaustive ungapped pure-Hamming <=3-mismatch enumeration; samse -n 1000000",
        "guides": 421, "reference": "AgamP5 = idAnoGambNW_F1_1 (GCF_943734735.2), 4 contigs",
        "census_per_mm_tier": {str(k): census[k] for k in sorted(census)},
        "census_total": sum(census.values()),
        "bowtie2_census_total": 24278,
        "superset_check": {
            "bowtie2_alignments_missing_from_bwa": len(bt2_missing),
            "bwa_only_alignments": len(bwa_only),
            "bwa_only_by_tier": {str(k): delta_by_tier[k] for k in sorted(delta_by_tier)},
        },
        "bwa_only_sequence_verification": {
            "sampled": 50, "verified_hamming_le_3": verified, "seed": 7,
            "method": "direct reference extraction at XA coordinates, reverse-complement for minus strand, Hamming recount vs guide"},
        "correction": "bowtie2 of-record census 24,278 under-reports by 25%; corrected census 32,305 (MM3 20,880 -> 28,765, +38%). bowtie2 -a still score-trims high-mismatch ungapped hits after seed tuning.",
        "cds_context_on_corrected_set": {
            "annotation": "NCBI RefSeq GCF_943734735.2 RS_2023_12 CDS (bedtools 2.31.1 intersect)",
            "per_mm_tier_total_inCDS": {str(k): v for k, v in sorted(per_tier.items())},
            "library_guides_with_at_least_one_CDS_overlapping_near_cognate": guides_with_cds_ot,
        },
        "leads": leads,
    }
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/offtarget_bwa.json", "w"), indent=1)
    print(json.dumps(out, indent=1)[:3000])

if __name__ == "__main__":
    main()
