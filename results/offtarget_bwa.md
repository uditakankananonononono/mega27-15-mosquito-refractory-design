# BWA exhaustive ungapped <=3-mismatch census (tool 31/40) - and a correction

**Tool**: BWA 0.7.19-r1273 `aln`/`samse` (lh3/bwa @ d82444c, built from source;
Li & Durbin 2009, *Bioinformatics* 25:1754-60). CRISPOR's engine is BWA; this
runs that engine first-party on our own data and cross-validates the bowtie2
census (tool 29/40).

**Mode**: `bwa aln -N -l 20 -k 3 -n 3 -o 0 -e 0` = non-iterative search for
ALL n-difference hits, full-length 20 nt seed, <=3 differences, gaps disabled
(pure-Hamming enumeration); `bwa samse -n 1000000` keeps every hit.

## Corrected census (421 guides vs AgamP5 = GCF_943734735.2)

| tier | bowtie2 (of record, tool 29) | BWA (this work) | delta |
|---|---|---|---|
| MM0 | 427 | 427 | 0 |
| MM1 | 227 | 227 | 0 |
| MM2 | 2,744 | 2,886 | +142 |
| MM3 | 20,880 | 28,765 | +7,885 |
| **total** | **24,278** | **32,305** | **+8,027 (+33%)** |

## Superset verification

- 0 of 24,278 bowtie2 alignments are missing from the BWA set (strict superset).
- 8,027 bwa-only alignments, all MM2/MM3 (table above).
- 50 seeded-random (seed 7) bwa-only hits re-verified by direct reference
  extraction at XA coordinates (reverse-complement on minus strand) and
  Hamming recount against the guide: **50/50 are true Hamming <=3 sites**.

## Verdict (honest, and it stings our own prior table)

The bowtie2 census under-reports true <=3-mismatch ungapped near-cognates by
25% (MM3 by 38%): `bowtie2 -a` still score-trims high-mismatch hits even after
seed tuning (`-L 20 -N 1`). MM0/MM1 agree exactly, so the failure is confined
to the highest-mismatch tiers. The corrected census of record is **32,305**
(BWA); the bowtie2 table is retained in Appendix/below as a documented
under-count (negative result preserved, per program rules). This is the
second bowtie2 sensitivity failure this project has surfaced (the first:
default presets miss 2-3 MM sites in 20-mers entirely) - lesson hardened:
**for guide enumeration, verify aligner recall against a second engine, not
just its own settings sweep.**

## CDS context on the corrected set (bedtools 2.31.1, AgamP5 RefSeq RS_2023_12)

| tier | total | in CDS | frac |
|---|---|---|---|
| MM0 | 427 | 244 | 57% (self-hits in dsx) |
| MM1 | 227 | 29 | 12.8% |
| MM2 | 2,886 | 412 | 14.3% |
| MM3 | 28,765 | 3,452 | 12.0% |

Library: **308/421 guides** carry >=1 CDS-overlapping near-cognate (was 239/421
on the under-counted set). Superset adds no new MM1 sites anywhere.

## Leads (unchanged safety picture)

- **dsx-v3-1**: MM3 4 total, 0 in CDS; no MM1/MM2 anywhere.
- **dsx-v3-2**: MM3 9 total, 0 in CDS; no MM1/MM2 anywhere. (MM0 self-hit in dsx.)
- **dsx-v3-3**: MM3 10 total; the only CDS overlap remains the known site
  NC_064601.1:38,008,162-38,008,181 in LOC1270942 (VPS13B) - flanks GCG|GTG,
  no adjacent NGG/CCN in either orientation -> **Cas9-inert**.
- **kyrou**: MM3 6 total; the only CDS overlap remains
  NC_064601.1:41,285,163-41,285,182 in LOC1274366 (uncharacterized) - flanks
  GTT|TAC, no adjacent NGG/CCN -> **Cas9-inert**.

The superset's extra 8,027 sites add no new CDS-overlapping near-cognate for
any lead; the two known lead CDS overlaps are unchanged and already resolved.
