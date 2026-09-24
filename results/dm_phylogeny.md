# DM-domain family phylogeny (tools 39-40: MUSCLE + FastTree)

Tree: `(dsx_AGAM:0.53798,Dsx_DROME:0.45880,(DmrtA2_AGAM:0.74733,dmd4_AGAM:0.48946)1.000:0.63690);`

- FastTree groups AgamP5 dsx with Drosophila Dsx (local support 1.0), excluding both paralogs - the tree-level engine agrees with HMMER (results/hmmer_dsx_ortholog.json) and DIAMOND (results/diamond_census.json).
- Patristic distance dsx_AGAM-Dsx_DROME = 0.9968; nearest paralog (DmrtA2) = 1.6643 (1.67x further). The ortholog call is quantitatively separated at the tree level, not only by search E-values.
- Alignment caveat: full-length proteins; only the DM domain is conserved across all four, so the conserved domain drives the signal. Documented, not smoothed over.
