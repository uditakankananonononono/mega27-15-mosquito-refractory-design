# Genome-completeness correction of the BWA census of record

**Tool:** BWA 0.7.19-r1273 aln/samse (rerun; identical parameters to the census of record)
**Reference:** COMPLETE AgamP5 = idAnoGambNW_F1_1 (GCF_943734735.2) `genomic.fna` — 191 records, 264,466,744 bp

## Why

The census of record (`offtarget_bwa.json`, 32,305 alignments) was computed against a
4-record FASTA (245,459,960 bp): chromosomes X, 2, 3 and mitochondrion only. The
complete assembly adds 187 unlocalized/unplaced `NW_*` scaffolds (19,006,784 bp = 7.2%
of the assembly). The four chromosome-level sequences are **byte-identical** (md5-verified)
between the two files; only the accessions differ (INSDC `OX030909.1/OX030908.1/OX030910.2`
= RefSeq `NC_064600.1/NC_064602.1/NC_083487.1`; chr2 is `NC_064601.1` in both).

## Census on the complete assembly

| tier | MM0 | MM1 | MM2 | MM3 | total |
|---|---|---|---|---|---|
| chromosomes | 427 | 227 | 2,886 | 28,765 | 32,305 |
| scaffolds | 0 | 0 | 88 | 543 | 631 |
| **complete** | **427** | **227** | **2,974** | **29,308** | **32,936** |

- Forward-map check: 0 of the 32,305 accession-mapped old alignments are missing from the complete rerun.
- Scaffolds contribute **no MM0 and no MM1 site anywhere** — no new exact or one-mismatch near-cognate.
- 50 seeded-random (seed 11) scaffold hits sequence-verified: 50/50 true Hamming ≤3 sites.
- **Corrected census of record: 32,936.**

## Correction decomposition (308 → 367 CDS-carrying guides)

The guide-level CDS-context count change is **entirely an accession-mismatch fix**, not new sequence:
the 4-record reference named chromosomes X/3/mito with INSDC accessions while the RefSeq CDS
annotation names them `NC_*`, so every chr3/X/mito CDS overlap was silently invisible to the
prior intersects (affected: the bedtools tenth check and the BWA eleventh check).

- previously reported: 308/421 guides with ≥1 CDS-overlapping near-cognate
- same 32,305 sites, accession-fixed: **367/421 (87.2%)**
- complete genome: 367/421 — scaffolds add **zero** additional CDS-carrying guides
  (7 guides have scaffold CDS hits, all already counted via chromosome hits)
- per-tier CDS overlap (complete set): MM0 247/427, MM1 60/227, MM2 816/2,974, MM3 7,136/29,308

## Lead-level impact

| lead | MM tiers (complete) | CDS-overlapping near-cognates |
|---|---|---|
| dsx-v3-1 | MM0=1 (own site), MM3=4 | none |
| dsx-v3-2 | MM0=1, MM3=9 | **1 — NEW (see below)** |
| dsx-v3-3 | MM0=1, MM3=10 | 1 (LOC1270942/VPS13B, chr2:38,008,162, PAM-less) — unchanged |
| kyrou | MM0=1, MM3=6 | 1 (LOC1274366, chr2:41,285,163, PAM-less) — unchanged |

**dsx-v3-2 newly visible CDS site** (present in the old census as `OX030908.1:2,111,186`, invisible
under the accession mismatch): guide `dsx_W1_1604_1627_RVS`, chr3 (`NC_064602.1`):2,111,186-2,111,205,
minus strand, Hamming = 3 (positions 4, 5, 7; guide `CATTAAGACCTACGAAGCGC` → ref `CATGGACACCTACGAAGCGC`),
inside **LOC1278105** (probable aconitate hydratase, mitochondrial). 3′ flank `CCC` — **no NGG PAM,
Cas9-inert**. All three lead CDS-overlapping near-cognates are now MM3 + PAM-less, a uniform picture:
no lead carries a CDS-overlapping near-cognate that Cas9 could actually cut.

*Regenerate:* `python3 scripts/bwa_census_complete.py` (commands in the script docstring).
