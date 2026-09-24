# Asset provenance

`doench2016_efficacy_cnn.pt` - GuideEfficacyCNN trained on Doench 2016 FC+RES
(5,310 experimentally measured 30-mer guides; MicrosoftResearch/azimuth mirror).
Held-out Spearman 0.556; cross-dataset on Doench V1 0.609 (see metrics JSON).
Architecture and training code: MEGA27-06 crisprgap package (same program).

## dsx locus identity verification
NCBI Gene 1270904 is labeled "uncharacterized LOC1270904" in current RefSeq, but it is
the record NCBI Gene returns for the alias AGAP004050 (VectorBase doublesex id);
RefSeq naming lags VectorBase. Identity independently confirmed in silico: the
published Kyrou 2018 dsx drive gRNA (5'-GTTTAACACAGGTCAAGCGG TGG-3', from the
paper's full text, PMC6871539) matches exactly one site inside the fetched locus
(NC_064601.1:47600000-47710000, minus strand, consistent with the gene's annotated
strand and position).
