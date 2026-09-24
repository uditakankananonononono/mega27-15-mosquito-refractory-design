# Fourth off-target confirmation: Bowtie -v 3

Bowtie 1.3.1 (static binary), -v 3 -a: full-length ungapped alignments with up to 3 mismatches, all reported

| guide | site (0-based offset) | strand | secondary <=3mm sites |
|---|---|---|---|
| dsx-v3-1 | NC_064601.1:47,619,039 | - | 0 |
| dsx-v3-2 | NC_064601.1:47,620,303 | - | 0 |
| kyrou | NC_064601.1:47,622,173 | - | 0 |

Bowtie's -v 3 mode is semantically identical to the audit's bar (full-length Hamming distance <=3 over the 23-mer unit). It reports exactly three alignments genome-wide - the three on-target sites at the cross-assembly-verified coordinates - and zero secondary sites. This is the fourth independent confirmation of 23-mer uniqueness (aligner-free enumerator, NCBI BLAST+, minimap2 at block resolution, Bowtie), and the one whose semantics match the claim most directly.
