# Cross-assembly check: minimap2 (orthogonal aligner)

minimap2 2.28-r1209 (preset asm5, -c), static binary from GitHub release

query: data/agamp4/AgamP4_2R_dsx_region.fasta (AgamP4 2R:48,700,000-48,790,000, 90,001 bp)
target: AgamP5 NC_064601.1 (2RL, 119,885,574 bp)

| guide | AgamP4 PAM start | minimap2 AgamP5 | aligner-free AgamP5 | delta |
|---|---|---|---|---|
| dsx-v3-1 | 48,711,498 | 47,619,047 | 47,619,040 | +7 bp |
| dsx-v3-2 | 48,712,762 | 47,620,311 | 47,620,304 | +7 bp |
| kyrou | 48,714,637 | 47,622,186 | 47,622,174 | +12 bp |

An aligner with gap placement (minimap2 asm5 chain) independently places all three guide windows on the same AgamP5 chromosome-2 strand within 7-12 bp of the aligner-free exact-substring coordinates. The uniform positive deltas are gap-placement effects of the chained alignment (largest at the kyrou window, consistent with the indel-rich RB-exon region of Appendix I). Cross-assembly concordance is confirmed at block resolution; the exact-substring coordinates remain the base-precise claim.
