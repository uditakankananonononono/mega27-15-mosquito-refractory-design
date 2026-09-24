# HMMER profile-homology verification of the dsx locus (tool 32/40)

**Tool**: HMMER 3.4 `phmmer` / `jackhmmer` (Eddy 2011, *PLoS Comput Biol*
7:e1002195), built from source. **Query**: D. melanogaster Dsx (UniProt
P23023). **Database**: AgamP5 RefSeq proteome (GCF_943734735.2, 30,505 proteins).

## Why

The AgamP5 RefSeq annotation labels the dsx locus **"uncharacterized protein
LOC1270904"**. Every targeting claim in this project therefore rests on our own
alignment-based checks (Bowtie/BLAST own-site hits). Profile-level homology to
the Drosophila Dsx reference provides an independent, annotation-free
confirmation - and a census of the DM-domain paralogs a drive could
inadvertently perturb.

## Results

- **dsx confirmed**: 9 isoforms of LOC1270904 top the phmmer search,
  E = 3.6e-51 to 6.6e-48 (best 175.3 bits). The gene spans
  NC_064601.1:47,610,877-47,700,920 (minus strand) and **contains all four
  lead protospacer loci** (47,619,040 / 47,620,304 / 47,622,174 / 47,622,798).
- **Unique ortholog**: the DM-domain family in AgamP5 numbers exactly 3 genes -
  dsx (LOC1270904), DmrtA2 (2 isoforms, best E = 4.1e-21) and dmd-4
  (E = 5.9e-18). The best non-dsx paralog trails by **30.1 orders of
  magnitude**. There is no second dsx-like locus for a drive to hit.
- **jackhmmer** (2 iterations): dsx at E = 4.1e-104.
- **Iterative-profile drift caveat (honest)**: jackhmmer iteration 2 boosts a
  spurious low-complexity ser/thr kinase hit (DDB_G0282963, phmmer
  E = 5.7e-08) to E = 1.1e-62 - above the true Dmrt paralogs. The phmmer
  ranking is the reliable one; iterative profiles can be contaminated by
  repetitive sequence. Recorded as a caveat, not used for any conclusion.
- **Annotation gap documented**: AgamP5 RefSeq does not name dsx; the project's
  target naming is correct but not annotation-derived.
