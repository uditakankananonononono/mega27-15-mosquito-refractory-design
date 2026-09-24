# Guide-spacer secondary-structure audit (ViennaRNA 2.7.2)

Null: 500 random 2R 20-mers, mean MFE -1.104 kcal/mol, sd 1.703.

| guide | spacer | MFE (kcal/mol) | null z | structured (z<-2) |
|---|---|---|---|---|
| v3-1 | TGGGCAGTATGCGTTAGGGT | -1.20 | -0.056 | False |
| v3-2 | CATTAAGACCTACGAAGCGC | 0.00 | 0.648 | False |
| v3-3 | GAAGCGAGCCCAATGGCTGT | -4.70 | -2.111 | True |
| kyrou2018-control | GTTTAACACAGGTCAAGCGG | -1.80 | -0.409 | False |

Verdict: at least one lead spacer is unusually structured.
Method: RNA.fold MFE on the 20-nt spacer; spacers treated as RNA (the species loaded into Cas9).
