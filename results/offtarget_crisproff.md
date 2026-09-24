# CRISPRoff/CRISPRspec energy-based specificity audit (tool 27/40)

- Tool: CRISPRoff/CRISPRspec pipeline v1.1.2 (github.com/RTH-tools/crisproff; Alkan et al. 2018, Genome Biol 19:177)
- Validation before use: repo test.sh reproduces test_data.example.out diff-clean; RNAfold requirement met by a shim over the ViennaRNA 2.7.2 python bindings (same library as results/guide_structure_vienna.json).
- Inputs: CRISPOR AgamP4 off-target sets per lead (lineage results/offtarget_crispor.json); each on-target 23-mer spacer-verified before the run.

## CRISPRspec specificity scores (energy-based; higher = more specific given the input off-target landscape)

| rank | lead | CRISPRspec | input off-targets (MM hist) |
|---|---|---|---|
| 1 | kyrou | 14.572 | 7 ({'4': 7}) |
| 2 | dsx-v3-1 | 11.947 | 6 ({'3': 2, '4': 4}) |
| 3 | dsx-v3-3 | 11.661 | 8 ({'3': 2, '4': 6}) |
| 4 | dsx-v3-2 | 3.716 | 6 ({'1': 1, '4': 5}) |

## Key corroborations

- dsx-v3-2 scores LOWEST (3.716) - the only lead with a <=1-mismatch near-cognate in its input set (the own-locus AgamP4 pos-20 C>T allele site resolved in results/offtarget_crispor.json / results/agamp4_reannotation.json). The energy model independently penalizes exactly that site.
- That MM1 site (CATTAAGACCTACGAAGCGTTGG) carries per-site CRISPRoff score 18.09 vs 23.37 for the on-target (77% of on-target binding energy) - quantifying how strongly the variant allele would bind if present in a target mosquito; reinforces the population-screen withdrawal of v3-2 as lead.
- kyrou ranks first (14.572), dsx-v3-1 second (11.947), dsx-v3-3 third (11.661); all three are separated from their nearest near-cognates by >=3 mismatches, consistent with the CHOPCHOP/CRISPOR cross-engine agreement.
- Independent method family: nucleic-acid duplex energy parameters (RNA/DNA hybrid thermodynamics), not alignment counts or ML scores - a third, physics-based axis agreeing with the alignment engines.

Method note: input sets are the CRISPOR AgamP4 near-cognate lists (up to 4 mismatches); CRISPRoff scores are therefore conditional on that input landscape, and CRISPRspec values are comparable across leads here because inputs come from the same engine, genome and thresholds.
