# bedtools 2.31.1 functional-context audit of bowtie2 near-cognates (tool 30/40)

Annotation: NCBI RefSeq GCF_943734735.2 idAnoGambNW_F1_1 (AgamP5), release RS_2023_12; 206,438 CDS features.
Input: all 24,278 ungapped <=3-mismatch bowtie2 alignments of the 421 dsx-window guides (results/offtarget_bowtie2.json; PAM-agnostic).

## Library-wide CDS overlap by mismatch tier (total, in-CDS)

| MM | total alignments | overlapping CDS | fraction |
|---|---|---|---|
| 0 | 427 | 244 | 57.1% |
| 1 | 227 | 29 | 12.8% |
| 2 | 2,744 | 390 | 14.2% |
| 3 | 20,880 | 2,624 | 12.6% |

- MM0 fraction is high because guides were designed inside the dsx locus (self-hits, 421 of 427).
- 239 of 421 library guides (56.8%) have at least one CDS-overlapping near-cognate (MM 1-3): a library-scale reminder that PAM-agnostic near-cognates are common, and the NGG filter plus mismatch intolerance does the protective work.

## Leads (MM tier: [total, in-CDS])

| lead | MM0 | MM1 | MM2 | MM3 | CDS-overlapping near-cognates |
|---|---|---|---|---|---|
| dsx-v3-1 | [1, 0] | 0 | 0 | [3, 0] | none |
| dsx-v3-2 | [1, 1] | 0 | 0 | [3, 0] | none (self-hit in dsx CDS, expected) |
| dsx-v3-3 | [1, 1] | 0 | 0 | [4, 1] | one: chr2:38,008,162-38,008,181 in LOC1270942 (VPS13B) |
| kyrou | [1, 0] | 0 | 0 | [2, 1] | one: chr2:41,285,163-41,285,182 in LOC1274366 (uncharacterized) |

## PAM resolution of the two lead CDS overlaps

Both sites are minus-strand alignments. Neither has an adjacent NGG (or CCN in the
reverse orientation) in AgamP5:

- dsx-v3-3 MM3 in VPS13B: flanks GCG | GTG - Cas9-inert.
- kyrou MM3 in LOC1274366: flanks GTT | TAC - Cas9-inert.

Verdict: the only CDS-overlapping near-cognates of any lead are PAM-less and
therefore not cleavable by SpCas9 - the functional-context audit adds gene-level
annotation without changing the leads' safety picture; it also quantifies the
background the PAM filter protects against (one in eight MM1-3 near-cognates
genome-wide overlaps coding sequence).
