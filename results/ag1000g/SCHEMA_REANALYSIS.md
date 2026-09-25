# Ag1000G mixed-schema correction (2026-09-25)

Primary analysis uses 4,106 complete, unfiltered, 3,250-row-per-window VCF slices. There are 575 older `hits` rows already prefiltered by GQ>=20 and DP>=5, analyzed separately in `allele_freq_all_vintages.json`. Twelve fetch-error rows remain in the manifest audit but supply no genotype.

The original aggregator read only `target_hits`, ignored 575 legacy `hits` records, included twelve fetch failures in denominators, and initialized all missing variant records as guide-matching instead of using their AgamP4 reference bases. The last defect is severe at assembly-discordant guide sites. The old bcftools selector also ignored `hits`: its 300 supposed zero-hit controls included legacy carrier rows. The corrected selector finds only 110 genuinely zero-target-call rows among all successful records at the same selection threshold; no new 300-control validation is claimed. The 453 modern-row parse comparison still checks positions and genotype calls, not this reference-background interpretation.

Unphased multisite het calls are a lower bound on allele compromise; low-quality sites make guide status unknown unless a different site already establishes two compromised alleles. The primary and sensitivity rates differ, especially for v3-2 and Kyrou, so neither is independent population replication.

| Guide | Primary callable n | Primary compromised AF | Mixed-vintage callable n | Mixed compromised AF |
|---|---:|---:|---:|---:|
| dsx-v3-1 | 4070 | 0.923710 | 4645 | 0.922713 |
| dsx-v3-2 | 3939 | 0.065626 | 4514 | 0.061918 |
| kyrou | 4071 | 0.035863 | 4646 | 0.031640 |

The prior 40.1% / 3.0% / 3.1% headline and 11-12-generation sieve law are withdrawn. Corrected primary-input simulation starts v3-1 above 50% at generation zero; the two other target windows cross at generation twelve under the tested model endpoints. This is simulation, not observed drive behavior.
