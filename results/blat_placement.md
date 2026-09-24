# BLAT protein-to-genome placement cross-check (tool 38)

- default settings: 31 hits chromosome-wide; exactly one places a single 8-aa block (24 bp, score 66) inside the dsx locus as part of a fragmented 9-block chain whose remaining blocks sit in intergenic DNA to 48.87 Mb. 13 spurious hits outscore it (top 102) - the placement is not identifiable de novo.
- loose (-minScore=10 -minIdentity=20): identical picture (32 hits; the same single locus block, score 66, outscored by 13 spurious hits) - loosening thresholds adds only sub-threshold noise, no additional gene structure.
- Two +- rows deceptively report tStart/tEnd inside the dsx and dmd-4 regions; block-level coordinates place both at 70.5 Mb and 114-116 Mb. Documented as a PSL parsing artifact.
- Sensitivity ordering for this 70%-identity cross-species placement: BLAT (default miss; loose single-block, outranked by spurious hits) < miniprot (-M2: full 8-CDS structure, results/miniprot_census.json) << HMMER profile search (family census, results/hmmer_dsx_ortholog.json). Honest negative/partial result, preserved.
