# Genome-wide miniprot DM-domain census (tool 33 extension)

**Tool:** miniprot 0.18-r281 (protein-to-genome aligner; Li 2023)
**Query:** D. melanogaster Dsx (UniProt P23023)
**Reference:** complete GCF_943734735.2 (191 records, 264.5 Mb), per-contig

## Memory engineering

Whole-genome indexing OOM-kills this 2 GB sandbox twice (documented in
`miniprot_dm_domain.json`). Per-contig runs still OOM at default `-M1` on chr2 (118 Mb).
Reducing modimiser sampling to `-M2` (1/4 k-mer sampling) fits: chr2 peak RSS 0.943 GB,
chr3 0.791 GB, chrX 0.251 GB, mito 0.034 GB, all-187-scaffolds 0.104 GB. `-M3` (1/8) fits
trivially but loses ALL hits including the known dsx locus — sampling sensitivity floor.

## Result

- **chr2 (NC_064601.1): dsx recovered, full gene structure** — mRNA 47,545,556-47,692,565
  (minus strand), 8 CDS blocks including the DM-domain exon at 47,692,245-47,692,565
  (70.0% identity). Better than the 90 kb locus run, which recovered only the DM-domain exon.
- **chrX, chr3, mito, 187 scaffolds: zero hits.**
- **Both DM-family paralogs missed.** dmd-4 (LOC1281789, chr2:4,056,166-4,057,948) and
  DmrtA2 (LOC1271814, chrX:7,711,435-7,726,245) are NOT detected — not by the -M2 contig
  scans, and not even by locus-level runs (±100 kb windows) at full default sensitivity
  with loose chaining (-n1 -m20). Native-protein controls map both loci at 100% identity,
  so the loci and extraction are correct: miniprot's exact 6-mer seed chaining simply has
  no signal for homologs at dmd-4/DmrtA2 divergence from Dsx (HMMER E = 5.9e-18 / 4.1e-21
  vs dsx 3.6e-51).

## Honest verdict

Cross-tool family census: HMMER (profile HMM) finds 3/3 DM-domain genes; miniprot
(seed-chained protein-to-genome) finds 1/3. The miss is a documented sensitivity limit of
seed-based alignment for diverged homologs, not evidence of absence — the paralogs are real
(100%-identity native mappings) and HMMER-detected. For single-ortholog confirmation
(miniprot's strength, splice-aware placement) the tools agree: dsx is the only DM-domain
gene miniprot-detectable by cross-species Dsx query, consistent with HMMER's 30-orders-of-
magnitude score gap to the best paralog.

*Regenerate:* `python3 scripts/miniprot_census.py`
