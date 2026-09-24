# CRISPOR independent specificity cross-check (tool 26/40)

**Tool:** CRISPOR (Haeussler et al. 2016, Genome Biol 17:148), crispor.gi.ucsc.edu, accessed 2026-09-25.
**Input:** two AgamP5 dsx windows (CRISPOR caps input at 2,300 bp): W1 NC_064601.1:47,618,700-47,620,999 (2,300 bp; dsx-v3-1, dsx-v3-2), W2 47,622,000-47,623,000 (1,001 bp; Kyrou control, dsx-v3-3). All four spacers verified present in the submitted sequence by exact string match pre-submission. PAM = NGG.
**Genomes:** (a) ensAnoGam = A. gambiae EnsemblMetazoa 76 (AgamP4, the assembly CHOPCHOP used); (b) GCF_943734685.1 = A. coluzzii AcolN3 (cross-species check). 345 (W1) and 85 (W2) window guides ranked per genome.
**Job URLs:** AgamP4 W1 batchId=W3xteqHaUZBayGdaDFz5, W2 batchId=5ZlWiVBOsR2HZBCJ9mzO; AcolN3 W1 batchId=FuBYOfv71qhlKQGVvBop, W2 batchId=YR58KOD5Z3WnrFqhmRyu (https://crispor.gi.ucsc.edu/crispor.py?batchId=...).

## Lead results (MM = off-target mismatch bins, genome-wide)

| lead | genome | MIT spec | CFD spec | MM0 | MM1 | MM2 | MM3 | MM4 |
|---|---|---|---|---|---|---|---|---|
| dsx-v3-1 | AgamP4 | 99 | 100 | 0 | 0 | 0 | 2 | 4 |
| dsx-v3-1 | AcolN3 | 98 | 99 | 0 | 0 | 0 | 2 | 7 |
| dsx-v3-2 | AgamP4 | 70 | 99 | 0 | 1 | 0 | 0 | 5 |
| dsx-v3-2 | AcolN3 | 100 | 100 | 0 | 0 | 0 | 0 | 2 |
| dsx-v3-3 | AgamP4 | 97 | 100 | 0 | 0 | 0 | 2 | 6 |
| dsx-v3-3 | AcolN3 | 97 | 99 | 0 | 0 | 0 | 2 | 6 |
| kyrou | AgamP4 | 98 | 100 | 0 | 0 | 0 | 0 | 7 |
| kyrou | AcolN3 | 98 | 99 | 0 | 0 | 0 | 0 | 7 |

## Findings

1. **v3-2's only sub-MM4 hit is its own locus.** The single MM1 in AgamP4 sits at 2R:48,712,762-48,712,784(-), exon AGAP004050 (dsx), exactly the guide's target site: the AgamP4 reference carries allele ...CGT where the guide is ...CGC (PAM-distal position 20 C>T; mit 41.70, CFD 0.300). Our own earlier reannotation (results/agamp4_reannotation.json) recorded this same AgamP4 allele at the same coordinates. It is a target-site allele difference between reference assemblies, not a dispersed off-target. In A. coluzzii the reference matches the guide's C allele and the lead has zero hits through 3 mismatches (MIT spec 100, CFD spec 100).
2. **The v3-3 duplicate flag is not reproduced.** CHOPCHOP's AgamP4 index returned MM0=2 for dsx-v3-3 (own site + one exact duplicate). CRISPOR's independent AgamP4 index returns MM0=0, and A. coluzzii also MM0=0. The flag is downgraded from "confirmed assembly-dependent duplicate" to "engine-specific, unresolved"; it remains a caveat on the third lead, which the population screen already withdrew.
3. **v3-1 and Kyrou agree across engines** (v3-1: MM0-2=0, 2 MM3, both engines; Kyrou: CHOPCHOP's MM0=1 is its conserved own site, which CRISPOR excludes by convention - nothing nearer than MM4 in either engine).
4. **Cross-species conservation.** All four on-target sites are recovered exactly in A. coluzzii at the dsx ortholog (XM_049607256.1 etc.): the drive target sites are conserved across the two sympatric, hybridizing members of the gambiae complex. dsx-v3-2 and the Kyrou guide are clean through MM3 in coluzzii; dsx-v3-1 and dsx-v3-3 clean through MM2.
