# Population allele redundancy audit (MMseqs2)

MMseqs2 (release 057db43, static AVX2 build) cluster, --min-seq-id 0.87 -c 0.9 --cov-mode 1 -k 6

43 allele 23-mers -> 27 clusters; 5 with >1 member.

- No allele sequence duplicates across guides: every multi-member cluster is intra-guide, so the three inventories are sequence-independent.
- Multi-member clusters are exactly single variants at adjacent positions of the same 23-mer unit (>=20/23 identity by construction), quantifying how concentrated population variation is around the PAM/seed region.
- dsx-v3-2's guide-matching restoration allele (AC 7,640) sits one cluster away from compromised neighbors (PAM-G disruption AC 11, seed AC 9): the functional and resistant haplotypes differ by 1-2 SNVs.
