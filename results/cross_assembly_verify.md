## Cross-assembly verification: guide targets on AgamP5

Exact and <=3-mismatch matches of each 23-nt protospacer+PAM unit (plus-strand orientation) on both strands of every AgamP5 molecule (NC_064601.1 chr2, OX030908.1 chr3, OX030909.1 chrX, OX030910.2 MT). Method: exact substring + pigeonhole near-match (4 blocks), no aligner.

AgamP4-reference units and guide-derived units searched separately: AgamP4 can carry a minor allele (assembly mismatch); the guide unit is the intended target allele.

| guide | AgamP4 unit exact | guide unit exact | guide mm1 | guide mm2 | guide mm3 |
|---|---|---|---|---|---|---|
| dsx-v3-1 | 0 | 0 | 1 | 0 | 0 |
| dsx-v3-2 | 0 | 1 | 0 | 0 | 0 |
| kyrou | 1 | 1 | 0 | 0 | 0 |

### Per-molecule exact hits (plus-strand 1-based = 0-based index + 1)

| guide | molecule | ref exact + | guide exact + |
|---|---|---|---|
| dsx-v3-2 | NC_064601.1 | [] | [47620304] |
| kyrou | NC_064601.1 | [47622174] | [47622174] |
