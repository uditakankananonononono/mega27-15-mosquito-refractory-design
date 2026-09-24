"""SnpEff functional annotation of Ag1000G dsx-window variant sites (tool 44).

The Ag1000G scan (results/ag1000g/on_target.jsonl, 4,693 per-sample VCFs)
records every non-reference call in the 3.25 kb dsx window
(2R:48,711,450-48,714,700, AgamP4 = NC_064601.1 in AgamP5 naming). This
script collapses those calls to distinct (pos, ref, alt) sites, maps each
sample's GT allele index to its ALT nucleotide, writes a VCF, and annotates
it with SnpEff 4.3t against a purpose-built AgamP5 chromosome-2 database
(RefSeq GCF_943734735.2 GFF + FASTA; DB scope = NC_064601.1 only, a
memory constraint of the build box, stated in the paper).

Question answered: which wild-population variant sites in and around the
three guide windows are protein-altering (missense / stop / splice), which
are synonymous or intronic, and how the functional classes distribute
between guide protospacer/PAM positions and the flanking window - plus the
nested aminopeptidase N-like gene (LOC11175624) the Prodigal analysis
surfaced inside a dsx intron.

Usage:
  python3 scripts/snpeff_window_ann.py --vcf      # build VCF of distinct sites
  python3 scripts/snpeff_window_ann.py --run      # run SnpEff (needs the local DB)
  python3 scripts/snpeff_window_ann.py --summarize  # results/snpeff_window_ann.{json,md}
"""
import json, os, subprocess, sys, collections

SCAN_IN = "results/ag1000g/on_target.jsonl"
VCF = "/tmp/ag1000g_dsx_window_sites.vcf"
ANN_VCF = "/tmp/ag1000g_dsx_window_sites.ann.vcf"
SNPEFF = "/tmp/snpEff/snpEff.jar"
DB = "agamp4_2r"
RESULT_JSON = "results/snpeff_window_ann.json"
RESULT_MD = "results/snpeff_window_ann.md"
TARGETS = {
    "dsx-v3-1_proto": (48711501, 48711520),
    "dsx-v3-1_PAM":   (48711498, 48711500),
    "dsx-v3-2_proto": (48712765, 48712784),
    "dsx-v3-2_PAM":   (48712762, 48712764),
    "kyrou_proto":    (48714640, 48714659),
    "kyrou_PAM":      (48714637, 48714639),
}
NESTED_GENE = ("LOC11175624", 47636160, 47639924)  # Prodigal/miniprot finding

def in_targets(pos):
    return [n for n, (a, b) in TARGETS.items() if a <= pos <= b]

def distinct_sites(path=SCAN_IN):
    """distinct (pos, ref, alt) -> carrier-chromosome tally (GT allele copies)."""
    sites = collections.Counter()
    for line in open(path):
        r = json.loads(line)
        if "error" in r:
            continue
        hits = r.get("target_hits", r.get("hits", []))
        for h in hits:
            alts = h["alt"].split(",")
            for a in h["gt"].replace("|", "/").split("/"):
                if a in (".", "0"):
                    continue
                try:
                    alt = alts[int(a) - 1]
                except (ValueError, IndexError):
                    continue
                sites[(h["pos"], h["ref"], alt)] += 1
    return sites

def write_vcf(sites, path=VCF):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=2R>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for (pos, ref, alt), ac in sorted(sites.items()):
            fh.write(f"2R\t{pos}\t.\t{ref}\t{alt}\t.\t.\tAC={ac}\n")
    return path

def run_snpeff(vcf=VCF, out=ANN_VCF):
    with open(out, "w") as fh:
        subprocess.run(["java", "-Xmx900m", "-jar", SNPEFF, "ann",
                        "-noStats", DB, vcf],
                       stdout=fh, stderr=subprocess.PIPE, check=True)
    return out

def parse_ann(path=ANN_VCF):
    """Parse ANN fields -> per-site records."""
    recs = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        pos, ref, alt = int(f[1]), f[3], f[4]
        ac = 0
        anns = []
        for kv in f[7].split(";"):
            if kv.startswith("AC="):
                ac = int(kv[3:])
            elif kv.startswith("ANN="):
                anns = kv[4:].split(",")
        effects = []
        for a in anns:
            p = a.split("|")
            if len(p) > 3:
                effects.append({"effect": p[1], "impact": p[2], "gene": p[3],
                                "transcript": p[6] if len(p) > 6 else ""})
        recs.append({"pos": pos, "ref": ref, "alt": alt, "ac": ac,
                     "targets": in_targets(pos), "effects": effects})
    return recs

def summarize(recs):
    worst_rank = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
    def worst(rec):
        if not rec["effects"]:
            return ("intergenic_or_unannotated", "NONE", "")
        w = max(rec["effects"], key=lambda e: worst_rank.get(e["impact"], -1))
        return (w["effect"], w["impact"], w["gene"])
    by_class = collections.Counter()
    tgt_sites, nested_sites = [], []
    for r in recs:
        eff, impact, gene = worst(r)
        by_class[eff] += 1
        if r["targets"]:
            tgt_sites.append({**r, "worst_effect": eff, "worst_impact": impact, "gene": gene})
    protein_altering = [t for t in tgt_sites if t["worst_impact"] in ("HIGH", "MODERATE")]
    per_guide = collections.Counter()
    for t in tgt_sites:
        for name in t["targets"]:
            guide = name.rsplit("_", 1)[0]
            per_guide[(guide, t["worst_impact"])] += 1
    return {
        "n_sites": len(recs),
        "effect_classes": dict(by_class.most_common()),
        "n_target_window_sites": len(tgt_sites),
        "target_window_sites": tgt_sites,
        "n_protein_altering_target_sites": len(protein_altering),
        "protein_altering_target_sites": protein_altering,
        "per_guide_impact_of_variant_sites": {
            f"{g}|{i}": n for (g, i), n in sorted(per_guide.items())},
        "nested_gene_note": ("the nested aminopeptidase N-like gene LOC11175624 "
                             "(47,636,160-47,639,924) lies ~97 kb outside the scanned "
                             "window, so this annotation makes no statement about it; "
                             "it remains a Prodigal/miniprot design caveat only"),
    }

def main():
    if "--vcf" in sys.argv or "--run" in sys.argv:
        sites = distinct_sites()
        write_vcf(sites)
        print(f"distinct sites: {len(sites)}")
    if "--run" in sys.argv:
        run_snpeff()
        print("SnpEff annotation done")
    if "--summarize" in sys.argv or "--run" in sys.argv:
        recs = parse_ann()
        out = {
            "tool": "SnpEff 4.3t (pcingola/SnpEff; Cingolani et al. 2012, Fly)",
            "database": "custom AgamP4 2R DB (Ensembl Metazoa AgamP4.63 GFF3 + dna_sm FASTA, chromosome 2R only - build-box memory constraint, documented); annotation in the Ag1000G-native AgamP4 coordinate system after an AgamP5-coordinate first pass was rejected as a namespace mismatch",
            "source": "distinct (pos, ref, alt) sites from the 4,693-sample Ag1000G dsx-window scan (results/ag1000g/on_target.jsonl), GT allele-index resolved",
            **summarize(recs),
        }
        with open(RESULT_JSON, "w") as fh:
            json.dump(out, fh, indent=2)
        md = [
            "## SnpEff functional annotation of the Ag1000G dsx-window variant sites",
            "",
            f"Tool: {out['tool']}",
            f"Database: {out['database']}",
            "",
            f"Distinct variant sites annotated: {out['n_sites']}",
            f"Effect classes (worst effect per site): {out['effect_classes']}",
            "",
            f"Sites inside guide protospacer/PAM targets: {out['n_target_window_sites']}",
            f"Protein-altering (HIGH/MODERATE) target sites: {out['n_protein_altering_target_sites']}",
            f"Per-guide worst-impact distribution of variant sites: "
            f"{out['per_guide_impact_of_variant_sites']}",
            f"Note: {out['nested_gene_note']}",
        ]
        with open(RESULT_MD, "w") as fh:
            fh.write("\n".join(md) + "\n")
        print(json.dumps({"n_sites": out["n_sites"],
                          "effect_classes": out["effect_classes"],
                          "n_target_window_sites": out["n_target_window_sites"],
                          "n_protein_altering_target_sites": out["n_protein_altering_target_sites"]},
                         indent=2))

if __name__ == "__main__":
    main()
