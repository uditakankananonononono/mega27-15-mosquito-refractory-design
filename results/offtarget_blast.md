# Orthogonal off-target confirmation: NCBI BLAST+

NCBI BLAST+ 2.17.0 (blastn -task blastn-short, word_size 7, dust off, evalue 1000), local DB from the four archived AgamP5 molecules (makeblastdb)

| guide | exact full-length hits | near full-length hits (1-6+ mm) | site |
|---|---|---|---|
| dsx-v3-1 | 1 | 0 | NC_064601.1:47,619,040-47,619,062 (minus) |
| dsx-v3-2 | 1 | 0 | NC_064601.1:47,620,304-47,620,326 (minus) |
| kyrou | 1 | 0 | NC_064601.1:47,622,174-47,622,196 (minus) |

An alignment-based search engine with independent indexing and scoring confirms the aligner-free enumerator's central specificity claim: over the full 23-nt cleavage-relevant unit, every guide has exactly one ungapped full-length genomic hit - its own site, at the cross-assembly-verified coordinate - and zero full-length hits at any non-zero mismatch count. The 2,171 shorter/gapped HSPs are sub-23-bp partial matches, below the cleavage-relevant unit, and do not contradict the <=3-mismatch 23-mer uniqueness result of Section 4.10.
