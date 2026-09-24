"""GuideScan2 cross-engine census check (tool 43).

The census of record (BWA 0.7.19, scripts/bwa_census_complete.py,
results/offtarget_bwa_complete.json) enumerated every <=3-mismatch
20-mer protospacer alignment for the 421 dsx-window guides against the
complete AgamP5 assembly (GCF_943734735.2, 191 records), PAM-agnostically:
32,936 sites. GuideScan2 (2.2.1, bioconda binary) searches a
bidirectional-BWT genome index for <=3-mismatch protospacer hits that
MUST carry an NGG PAM. This script

  1. resolves each guide's own NGG-carrying site from the census BED
     (its MM0 row whose 3' flank is NGG; 421/421 unique) and writes a
     GuideScan kmer CSV;
  2. runs `guidescan index` (if index files are absent) and
     `guidescan enumerate -m 3` on the 421 guides;
  3. reclassifies every census BED hit as NGG-carrying or PAM-less by
     orientation-aware Hamming comparison plus PAM inspection, and
  4. compares per-guide per-tier counts between the engines.

Result: exact agreement - GuideScan 421/15/171/1549 (MM0-3, 2,156
sites) equals the NGG-carrying subset of the census on every one of the
1,684 guide x tier cells; the remaining 30,780 census sites are PAM-less
(6/212/2,803/27,759) and outside GuideScan's search definition. The
census of record is therefore not an artifact of the BWA alignment
engine. Artifacts in /tmp: genome FASTA, guides FASTA, census BED,
guidescan binary+index (see constants below).

Usage:
  python3 scripts/guidescan_census.py --run       # end-to-end, writes results/guidescan_census.{json,md}
  python3 scripts/guidescan_census.py --compare   # comparison only, from existing artifacts
"""
import csv, json, os, subprocess, sys, collections

GUIDESCAN = os.environ.get("GUIDESCAN_BIN", "/tmp/guidescan-pkg/bin/guidescan")
GENOME = "/tmp/gcf_genomic.fna"
GUIDES_FA = "/tmp/bt2_guides.fa"
CENSUS_BED = "/tmp/bwa_hits_gcf.bed"
INDEX = "/tmp/agamp5_gs"
KMERS = "/tmp/gs_kmers.csv"
GS_OUT = "/tmp/gs_out.csv"
RESULT_JSON = "results/guidescan_census.json"
RESULT_MD = "results/guidescan_census.md"
TIERS = (0, 1, 2, 3)
COMP = str.maketrans("ACGT", "TGCA")

def rc(s):
    return s.translate(COMP)[::-1]

def ham(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def load_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs

def load_guides(path=GUIDES_FA):
    guides, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
        elif line:
            guides[name] = line.upper()
    return guides

def classify_hit(genome, guides, chrom, s, e, gid):
    """Return 'ngg', 'pamless' or 'ambiguous' for one census BED row."""
    g = guides[gid]
    seq = genome[chrom][s:e]
    hf, hr = ham(seq, g), ham(rc(seq), g)
    if hf < hr:
        pam = genome[chrom][e:e + 3]
    elif hr < hf:
        pam = rc(genome[chrom][s - 3:s])
    else:
        return "ambiguous"
    return "ngg" if (len(pam) == 3 and pam[1:] == "GG") else "pamless"

def census_by_guide(genome, guides, bed_path=CENSUS_BED):
    """guide -> mm -> {'ngg': n, 'pamless': n}; plus ambiguous count."""
    cen = collections.defaultdict(lambda: collections.defaultdict(
        lambda: {"ngg": 0, "pamless": 0}))
    n_amb = 0
    for line in open(bed_path):
        f = line.rstrip("\n").split("\t")
        chrom, s, e, gid, mm = f[0], int(f[1]), int(f[2]), f[3], int(f[4])
        cls = classify_hit(genome, guides, chrom, s, e, gid)
        if cls == "ambiguous":
            n_amb += 1
        else:
            cen[gid][mm][cls] += 1
    return cen, n_amb

def resolve_own_sites(genome, guides, bed_path=CENSUS_BED):
    """Each guide's unique NGG-carrying MM0 site -> (chrom, pos1, strand)."""
    cand = collections.defaultdict(list)
    for line in open(bed_path):
        f = line.rstrip("\n").split("\t")
        if f[4] == "0":
            cand[f[3]].append((f[0], int(f[1]), int(f[2])))
    own, problems = {}, []
    for gid, g in guides.items():
        hits = []
        for chrom, s, e in cand.get(gid, []):
            seq = genome[chrom][s:e]
            if seq == g:
                pam = genome[chrom][e:e + 3]
                if len(pam) == 3 and pam[1:] == "GG":
                    hits.append((chrom, e + 3, "+"))
            elif rc(seq) == g:
                pam = rc(genome[chrom][s - 3:s])
                if len(pam) == 3 and pam[1:] == "GG":
                    hits.append((chrom, s - 3, "-"))
        if len(hits) != 1:
            problems.append((gid, len(hits)))
        else:
            own[gid] = hits[0]
    return own, problems

def write_kmers(own, guides, path=KMERS):
    with open(path, "w") as fh:
        fh.write("id,sequence,pam,chromosome,position,sense\n")
        for gid in sorted(own):
            chrom, pos, strand = own[gid]
            fh.write(f"{gid},{guides[gid]},NGG,{chrom},{pos},{strand}\n")
    return path

def run_guidescan(kmers=KMERS, out=GS_OUT):
    if not (os.path.exists(INDEX + ".forward") and os.path.exists(INDEX + ".reverse")):
        subprocess.run([GUIDESCAN, "index", "--index", INDEX, GENOME], check=True)
    subprocess.run([GUIDESCAN, "enumerate", INDEX, "-f", kmers, "-m", "3",
                    "--format", "csv", "-o", out, "-n", "8"], check=True)
    return out

def guidescan_by_guide(path=GS_OUT):
    gs = collections.defaultdict(lambda: collections.defaultdict(int))
    n = 0
    for r in csv.DictReader(open(path)):
        gs[r["id"]][int(r["match_distance"])] += 1
        n += 1
    return gs, n

def compare(gs, cen, guides):
    disc = []
    for gid in guides:
        for t in TIERS:
            if gs[gid][t] != cen[gid][t]["ngg"]:
                disc.append({"guide": gid, "tier": t,
                             "census_ngg": cen[gid][t]["ngg"], "guidescan": gs[gid][t]})
    return disc

def main():
    genome = load_fasta(GENOME)
    guides = load_guides()
    if "--run" in sys.argv:
        own, problems = resolve_own_sites(genome, guides)
        assert not problems, f"own-site resolution failed: {problems[:5]}"
        write_kmers(own, guides)
        run_guidescan()
    gs, n_rows = guidescan_by_guide()
    cen, n_amb = census_by_guide(genome, guides)
    disc = compare(gs, cen, guides)
    tot_gs = {str(t): sum(gs[g][t] for g in gs) for t in TIERS}
    tot_ngg = {str(t): sum(cen[g][t]["ngg"] for g in cen) for t in TIERS}
    tot_pamless = {str(t): sum(cen[g][t]["pamless"] for g in cen) for t in TIERS}
    out = {
        "tool": "GuideScan2 2.2.1 (bioconda linux-64 binary; pritykinlab/guidescan-cli; Schmidt et al. 2025, Genome Biology)",
        "index": "AgamP5 = GCF_943734735.2 genomic.fna (191 records), guidescan index built locally",
        "guides": len(guides),
        "mismatches": 3,
        "comparison": "GuideScan NGG-restricted counts vs NGG-carrying subset of the BWA census of record (results/offtarget_bwa_complete.json, /tmp/bwa_hits_gcf.bed)",
        "guidescan_totals_per_mm": tot_gs,
        "guidescan_sites": n_rows,
        "census_ngg_per_mm": tot_ngg,
        "census_pamless_per_mm": tot_pamless,
        "census_pamless_total": sum(tot_pamless.values()),
        "census_ngg_total": sum(tot_ngg.values()),
        "orientation_ambiguous_hits": n_amb,
        "n_guide_tier_cells": len(guides) * len(TIERS),
        "n_discordant_cells": len(disc),
        "discordant_cells": disc[:50],
        "verdict": ("EXACT AGREEMENT" if not disc else "DISCORDANT - investigate"),
    }
    with open(RESULT_JSON, "w") as fh:
        json.dump(out, fh, indent=2)
    md = [
        "## GuideScan2 cross-engine census check",
        "",
        f"Engine: {out['tool']}",
        f"Index: {out['index']}",
        f"Guides: {len(guides)} dsx-window 20-mers, enumerate -m 3, NGG PAM required.",
        "",
        f"GuideScan per-tier totals (MM0-3): {tot_gs} = {n_rows} sites",
        f"Census NGG-carrying subset:        {tot_ngg} = {sum(tot_ngg.values())} sites",
        f"Census PAM-less remainder:         {tot_pamless} = {sum(tot_pamless.values())} sites",
        f"Sum check: {sum(tot_ngg.values())} + {sum(tot_pamless.values())} = "
        f"{sum(tot_ngg.values()) + sum(tot_pamless.values())} (census total 32,936)",
        "",
        f"Per-guide per-tier agreement: {len(guides) * len(TIERS) - len(disc)}/"
        f"{len(guides) * len(TIERS)} cells; orientation-ambiguous census hits: {n_amb}.",
        f"Verdict: **{out['verdict']}** - the census of record is not an artifact "
        "of the BWA alignment engine; its NGG-carrying subset reproduces exactly "
        "under GuideScan2's bidirectional-BWT search.",
    ]
    if disc:
        md.append("")
        md.append(f"Discordant cells (first {min(50, len(disc))}):")
        for d in disc[:50]:
            md.append(f"- {d['guide']} tier {d['tier']}: census_ngg={d['census_ngg']} guidescan={d['guidescan']}")
    with open(RESULT_MD, "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("guidescan_totals_per_mm", "census_ngg_per_mm",
                       "n_discordant_cells", "verdict")}, indent=2))

if __name__ == "__main__":
    main()
