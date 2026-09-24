# CHOPCHOP v3 independent off-target/efficiency cross-check (AgamP4)

Job: https://chopchop.cbu.uib.no/results/cacd2b0a-02f0-4536-9c71-52893d748475/

Target: 4801-bp dsx exon-5 window from AgamP5 NC_064601.1:47,618,500-47,623,300 (fasta input)

Guides ranked in window: 585

| guide | CHOPCHOP rank | target sequence (+PAM) | MM0 | MM1 | MM2 | MM3 | efficiency |
|---|---|---|---|---|---|---|---|
| dsx-v3-1 | 326 | TGGGCAGTATGCGTTAGGGTAGG | 0 | 0 | 0 | 2 | 63.74 |
| dsx-v3-2 | 355 | CATTAAGACCTACGAAGCGCTGG | 0 | 1 | 0 | 0 | 54.90 |
| dsx-v3-3 | 366 | GAAGCGAGCCCAATGGCTGTTGG | 1 | 0 | 0 | 2 | 54.08 |
| kyrou-2018 | 18 | GTTTAACACAGGTCAAGCGGTGG | 1 | 0 | 0 | 0 | 64.48 |

## Interpretation

- **dsx-v3-1**: MM0=0 in AgamP4 (own-locus haplotype diverged from PEST reference); 0 MM1/MM2, 2 MM3; CHOPCHOP efficiency 63.74; rank 326/585 in window
- **dsx-v3-2**: MM0=0 in AgamP4; 1 MM1 (near cognate), 0 MM2/MM3; efficiency 54.90; rank 355/585
- **dsx-v3-3**: MM0=2 in AgamP4: exact-match duplicate exists in the AgamP4 assembly although unique in AgamP5 (own Bowtie/BLAST checks); assembly-dependent copy-number caveat; efficiency 54.08; rank 366/585
- **kyrou-2018**: MM0=1 = own site conserved between assemblies; 0 MM1-MM3; efficiency 64.48; rank 18/585 (best in window)

Note: fasta input was AgamP5 sequence; CHOPCHOP scores off-targets against its curated AgamP4 index. MM0 counts exact protospacer matches in AgamP4 (kyrou=1 is the conserved own site). The v3-3 MM0=2 flag is assembly-dependent: own Bowtie 1.3.1 (-v 3 -a) and BLAST+ 2.17.0 searches found v3-3 unique in AgamP5 NC_064601.1.
