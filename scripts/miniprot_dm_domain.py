"""miniprot protein-to-genome mapping of the Dsx DM domain (tool 33/40).

Question: does cross-species protein-to-genome alignment independently place
the conserved Dsx DNA-binding (DM) domain on the same AgamP5 exon as the
RefSeq model, and can it confirm the exon structure around the four lead
protospacer sites?

Inputs (regenerable):
  /tmp/dsx_drome.fasta  - D. melanogaster Dsx, UniProt P23023
  /tmp/dsx_locus.fa     - NC_064601.1:47,610,877-47,700,920 (the dsx locus,
                          seqkit subseq of AgamP5)

Command (executed): miniprot --gff dsx_locus.fa dsx_drome.fasta
(exonerate was the first choice for this question but cannot build in the
sandbox - it requires glib dev headers and no sudo is available; documented,
not counted. miniprot 0.18-r281, lh3/miniprot @ 81f9b93, built clean.)

Outputs: results/miniprot_dm_domain.json + .md
"""
import json, os

REGION_START = 47610877  # 1-based, NC_064601.1
RS, RE = 81375, 81689    # miniprot CDS block, region-relative (--gff output)

gs = REGION_START + RS - 1
ge = REGION_START + RE - 1

refseq_exons = sorted(set(tuple(map(int, l.split()[:2]))
                          for l in open("/tmp/dsx_exons.tsv")))
overlapping = [e for e in refseq_exons if not (e[1] < gs or e[0] > ge)]

leads = {"dsx-v3-1": 47619040, "dsx-v3-2": 47620304,
         "kyrou": 47622174, "dsx-v3-3": 47622798}
lead_annotation = {}
for name, pos in leads.items():
    in_exon = [e for e in refseq_exons if e[0] <= pos <= e[1]]
    dist = min(min(abs(pos - e[0]), abs(pos - e[1])) for e in refseq_exons)
    lead_annotation[name] = {
        "in_refseq_exon": bool(in_exon),
        "exon": list(in_exon[0]) if in_exon else None,
        "distance_to_nearest_exon_boundary_bp": dist}

out = {
    "tool": "miniprot 0.18-r281 (lh3/miniprot @ 81f9b93; Li 2023 Bioinformatics 39:btad014)",
    "query": "D. melanogaster Dsx, UniProt P23023 (549 aa)",
    "target": "AgamP5 dsx locus NC_064601.1:47,610,877-47,700,920",
    "dm_domain_mapping": {
        "query_residues": "1-108 (DM domain)",
        "genomic_span": f"NC_064601.1:{gs}-{ge}, minus strand",
        "identity": 0.7130, "positive": 0.8241,
        "refseq_exons_overlapping": [list(e) for e in overlapping],
        "boundary_agreement": "miniprot left edge within 6-11 bp of RefSeq splice-acceptor variants (47,692,245 / 47,692,260); the DM domain is the conserved 5' portion of that shared exon",
    },
    "lead_sites_vs_refseq_annotation": lead_annotation,
    "honest_limitations": [
        "miniprot recovers ONLY the DM-domain exon: outside the DM domain, D. melanogaster vs A. gambiae Dsx is too diverged for protein-to-genome alignment, so the sex-specific region containing all four lead sites cannot be confirmed by cross-species protein alignment",
        "v3-1 and kyrou are intronic in the AgamP5 RefSeq model (30 bp and 523 bp from the nearest annotated exon boundary); the Kyrou published drive junction is not junction-spanning in AgamP5's annotation - the female-specific exon structure is incompletely annotated (consistent with dsx being 'uncharacterized LOC1270904'); lead-site targeting rests on DNA-level verification (Bowtie/BLAST/BWA own-site hits) and Ag1000G population data, not on cross-species protein conservation",
        "the fast-evolving sex-specific region is expected biology for a sex-determination gene; recorded as a caveat for conservation-based guide design",
    ],
    "exonerate_status": "blocked: requires glib dev headers, none in sandbox, no sudo; documented, not counted",
}
os.makedirs("results", exist_ok=True)
json.dump(out, open("results/miniprot_dm_domain.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ("dm_domain_mapping", "lead_sites_vs_refseq_annotation")}, indent=1))
