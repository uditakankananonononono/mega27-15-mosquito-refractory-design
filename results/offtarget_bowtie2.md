# bowtie2 2.5.4 whole-genome off-target audit (tool 29/40)

- guides aligned: 421 (unique 20-mers), self-hit found for 421
- library off-targets <=3 MM (PAM-agnostic): 23,851; MM hist {1: 227, 2: 2744, 3: 20880}
- median 9, max 6598 off-targets per guide

## Leads

- v3-1: 3 OTs <=3MM (hist {'3': 3}), self-hits 1; FlashFry PAM-adjacent otCount 7, FF MM0-4 1,0,0,1,5
- v3-2: 3 OTs <=3MM (hist {'3': 3}), self-hits 1; FlashFry PAM-adjacent otCount 4, FF MM0-4 1,0,0,0,3
- v3-3: 4 OTs <=3MM (hist {'3': 4}), self-hits 1; FlashFry PAM-adjacent otCount 10, FF MM0-4 1,0,0,2,7
- kyrou: 2 OTs <=3MM (hist {'3': 2}), self-hits 1; FlashFry PAM-adjacent otCount 7, FF MM0-4 1,0,0,1,5

## MM0 duplicate resolution

427 MM0 alignments / 421 guides. The 6 extra copies (3 guides: dsx_W1_876_899_FWD chr2:21,601,508; dsx_W1_578_601_RVS chr2:107,404,399(-), chr2:41,326,168, OX030908.1:7,545,406; dsx_W1_1445_1468_FWD OX030909.1:2,997,367, OX030909.1:3,004,474(-)) all lack an adjacent NGG PAM -> Cas9-inert. FlashFry's zero-MM0-duplicate result for NGG-adjacent Cas9 sites confirmed; 3/421 guides (0.7%) carry PAM-less perfect copies visible only to PAM-agnostic aligners.

## Sensitivity caveat (honest negative result)

A first pass with bowtie2's default --sensitive preset (-D 15 -R 2 -L 22) under-reported the library: 660 alignments, zero 2-3-mismatch sites - its 22 bp seed cannot tolerate 2-3 mismatches in a 20-mer read. The of-record run uses --very-sensitive -D 100 -R 5 -L 14 -N 1 -i S,1,0.50 (24,278 alignments, MM hist 427/227/2,744/20,880). PAM-agnostic counts exceed FlashFry's PAM-adjacent counts ~14x uniformly across MM1-3, quantifying the NGG constraint's filtering power on this genome.
