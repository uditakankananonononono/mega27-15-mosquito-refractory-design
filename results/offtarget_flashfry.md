# FlashFry 2.0 library-scale off-target audit (AgamP5)

Guides discovered in dsx_W1+dsx_W2 (3,373 bp): 430; scored: 421 (overflow-flagged: 0). Whole-genome AgamP5 database, <=4 mismatches.

## Library
- guides with zero off-targets (<=4 MM): 0/421
- guides with an MM0 duplicate: 0
- off-target count median 14, max 1570
- genome-wide hit totals by mismatch (0,1,2,3,4): [421, 15, 171, 1491, 11276]

## Leads
| lead | CFD specificity | Hsu2013 | MM hist 0-4 | OTs | MM1-3 OTs |
|---|---|---|---|---|---|
| dsx-v3-1 | 0.702 | 98.8 | 1,0,0,1,5 | 7 | 1 |
| dsx-v3-2 | 0.782 | 99.0 | 1,0,0,0,3 | 4 | 0 |
| dsx-v3-3 | 0.567 | 97.4 | 1,0,0,2,7 | 10 | 2 |
| kyrou-2018 | 0.706 | 98.9 | 1,0,0,1,5 | 7 | 1 |

## Cross-engine conclusions
- v3-2 shows ZERO 1-3-mismatch off-targets in AgamP5 (mm_hist 1,0,0,0,3): the one-mismatch site CRISPOR found in AgamP4 and CRISPRoff flagged at 77% of on-target energy is an AgamP4-reference-allele artifact, absent from the AgamP5 assembly; the flag tracks the assembly, not the guide.
- No lead carries an MM0 duplicate in AgamP5 (all mm0 counts = 1, the on-target): CHOPCHOP's MM0=2 flag on a lead is not reproduced by FlashFry, consistent with the CRISPOR cross-engine resolution.
- v3-3 is the least specific lead on this third independent axis too (CFD specificity lowest, 10 off-targets <=4 MM, 2 at 3 MM), adding to the ViennaRNA structural outlier flag.
- v3-2 is the most sequence-specific lead in AgamP5 (CFD specificity highest, only 3 MM4 hits): its withdrawal rests on population-diversity evidence, not reference-genome off-target load.
