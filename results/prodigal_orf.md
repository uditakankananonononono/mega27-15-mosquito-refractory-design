# Prodigal annotation-independent coding-potency probe (tool 41)

- ORF COUNT does not discriminate: 86 ORFs in the dsx locus vs 79 in the gene-free control - expected for a prokaryotic tool on eukaryotic sequence.
- ORF SCORE does: the two strongest de novo ORFs in the locus (51.96 and 36.35) are the two real genes, both above 2x the control maximum (15.82).
- Strongest ORF (score 51.96, conf 100.0, + strand, locus 26981-28792 = genomic 47,637,857-47,639,668) is LOC11175624, an aminopeptidase N-like gene NESTED INSIDE A DSX INTRON - recovered with zero annotation input. Design note: any drive-cargo placement in dsx introns downstream of the lead cluster must account for this nested gene (the four leads sit ~13 kb upstream of it).
- Second-strongest (score 36.35, complement(8265..9554)) covers the RefSeq constitutive exon 8265-9460 with an EXACT shared stop-boundary at 8265.
- The lead-target exon map agrees with the miniprot exon audit: dsx-v3-2 (locus 9428) and dsx-v3-3 (11,922) fall inside RefSeq exons; dsx-v3-1 (8,164) and Kyrou (11,298) are intronic - independently corroborated here without using the annotation at prediction time.
