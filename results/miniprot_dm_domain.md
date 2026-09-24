# miniprot protein-to-genome mapping of the Dsx DM domain (tool 33/40)

**Tool**: miniprot 0.18-r281 (Li 2023, *Bioinformatics* 39:btad014), built from
source. **Query**: D. melanogaster Dsx (UniProt P23023, 549 aa). **Target**:
the AgamP5 dsx locus, NC_064601.1:47,610,877-47,700,920.

(exonerate was the first choice for splice-aware alignment but cannot build in
this sandbox - it requires glib dev headers and no sudo is available.
Documented, not counted; miniprot answers the same protein-to-genome question.)

## Results

- **DM domain placed independently**: miniprot maps P23023 residues 1-108 (the
  DM DNA-binding domain) to NC_064601.1:47,692,251-47,692,565 (minus strand),
  71.3% identity / 82.4% positive. The block overlaps three RefSeq exon
  annotations sharing the 47,693,483 donor edge; miniprot's left edge agrees
  with the RefSeq splice-acceptor variants (47,692,245 / 47,692,260) within
  6-11 bp. Protein-level and annotation-level exon placement agree.
- **Honest limitation**: miniprot recovers ONLY the DM-domain exon. Outside it,
  D. melanogaster vs A. gambiae Dsx is too diverged for protein-to-genome
  alignment - so the sex-specific region containing all four lead sites cannot
  be confirmed by cross-species protein alignment.
- **Lead sites vs the AgamP5 RefSeq model** (GFF exon audit, this analysis):
  dsx-v3-2 is exonic (47,619,070-47,620,336; 32 bp from the boundary) and
  dsx-v3-3 is exonic (47,622,697-47,622,831; 33 bp); **dsx-v3-1 and kyrou are
  intronic in the AgamP5 model** - 30 bp and 523 bp from the nearest annotated
  boundary. The published Kyrou drive junction is therefore not
  junction-spanning in AgamP5's annotation: the female-specific exon structure
  is incompletely annotated, consistent with dsx being 'uncharacterized
  LOC1270904' (tool 32 finding). Lead targeting rests on DNA-level
  verification (Bowtie/BLAST/BWA own-site hits) and Ag1000G population data,
  not on cross-species protein conservation - expected for the fast-evolving
  sex-specific region of a sex-determination gene, and recorded as a caveat.
