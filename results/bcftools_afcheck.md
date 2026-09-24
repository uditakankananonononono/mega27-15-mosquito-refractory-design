## bcftools independent recheck of the Ag1000G dsx-window scan

Engine: bcftools 1.21 (htslib full VCF parse) vs pysam.TabixFile (raw text slice)
Subsample: seed 42, 400 carriers + 300 zero-hit samples of 4,693.
Compared: 700 samples (0 fetch errors).
Target-window row concordance (pos/ref/alt/GT): **677/700 = 0.967143**
Whole-window row-count concordance (n_rows-carrying schemas only, 453 samples): 0/453 = 0.0
Off-by-one row-count delta (bcftools = pysam + 1): 453/453 - coordinate convention: pysam fetch() is 0-based half-open, bcftools -r is 1-based inclusive, so bcftools also returns the row at 2R:48,711,450 (outside all six targets).

Shared transport (both link htslib for bgzf/tabix/HTTPS); independent parse layer
(raw text slicing vs full VCF parse + re-emit).

Scan-vintage split: 453 scan-of-record rows (unfiltered) vs 247 legacy rows (ag1000g_on_target.py, GQ>=20/DP>=5 filter).
Target concordance on scan-of-record rows: **453/453 = 1.0**
Legacy-row discordance: 23 samples, of which 23 are fully explained by the legacy quality filter (bcftools recovers exactly the calls that filter drops).

Discordant samples (first 50): 23 total
- AB0107-C: only_pysam=[] only_bcftools=[(48712765, 'A', 'C,T,G', '3/3')]
- AB0349-C: only_pysam=[] only_bcftools=[(48712765, 'A', 'C,T,G', '3/3')]
- AB0262-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '2/2'), (48711507, 'C', 'A,T,G', '2/2'), (48711510, 'C', 'A,T,G', '2/2'), (48711511, 'T', 'A,C,G', '1/1')]
- AB0303-C: only_pysam=[] only_bcftools=[(48711516, 'A', 'C,T,G', '0/3')]
- AB0154-Cx: only_pysam=[] only_bcftools=[(48712765, 'A', 'C,T,G', '3/3')]
- AB0216-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2')]
- AB0156-Cx: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '2/2')]
- AB0108-C: only_pysam=[] only_bcftools=[(48711501, 'A', 'C,T,G', '0/1'), (48711516, 'A', 'C,T,G', '0/3')]
- AB0304-CW: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2'), (48711507, 'C', 'A,T,G', '0/1'), (48711511, 'T', 'A,C,G', '0/1'), (48711516, 'A', 'C,T,G', '0/3'), (48712765, 'A', 'C,T,G', '3/3')]
- AB0213-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2'), (48711510, 'C', 'A,T,G', '2/2'), (48711511, 'T', 'A,C,G', '1/1')]
- AB0216-CW: only_pysam=[] only_bcftools=[(48711511, 'T', 'A,C,G', '1/1')]
- AB0107-Cx: only_pysam=[] only_bcftools=[(48712765, 'A', 'C,T,G', '3/3')]
- AB0310-CW: only_pysam=[] only_bcftools=[(48711498, 'C', 'A,T,G', '0/1'), (48711500, 'C', 'A,T,G', '2/2'), (48711511, 'T', 'A,C,G', '1/1'), (48711516, 'A', 'C,T,G', '3/3'), (48712765, 'A', 'C,T,G', '3/3')]
- AB0305-CW: only_pysam=[] only_bcftools=[(48712765, 'A', 'C,T,G', '3/3')]
- AB0182-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2'), (48711507, 'C', 'A,T,G', '0/2')]
- AA0109-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2'), (48711507, 'C', 'A,T,G', '0/2'), (48711510, 'C', 'A,T,G', '0/2')]
- AB0141-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '2/2')]
- AB0132-CW: only_pysam=[] only_bcftools=[(48711507, 'C', 'A,T,G', '0/2'), (48711510, 'C', 'A,T,G', '0/2')]
- AB0247-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '2/2'), (48711507, 'C', 'A,T,G', '2/2'), (48711510, 'C', 'A,T,G', '2/2'), (48711511, 'T', 'A,C,G', '1/1')]
- AB0117-Cx: only_pysam=[] only_bcftools=[(48711499, 'C', 'A,T,G', '0/1')]
- AB0114-Cx: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2'), (48711507, 'C', 'A,T,G', '0/2'), (48711509, 'G', 'A,C,T', '0/2')]
- AB0184-C: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '0/2'), (48711510, 'C', 'A,T,G', '0/2')]
- AB0121-Cx: only_pysam=[] only_bcftools=[(48711500, 'C', 'A,T,G', '2/2'), (48711507, 'C', 'A,T,G', '2/2'), (48711510, 'C', 'A,T,G', '2/2'), (48711511, 'T', 'A,C,G', '1/1')]
