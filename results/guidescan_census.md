## GuideScan2 cross-engine census check

Engine: GuideScan2 2.2.1 (bioconda linux-64 binary; pritykinlab/guidescan-cli; Schmidt et al. 2025, Genome Biology)
Index: AgamP5 = GCF_943734735.2 genomic.fna (191 records), guidescan index built locally
Guides: 421 dsx-window 20-mers, enumerate -m 3, NGG PAM required.

GuideScan per-tier totals (MM0-3): {'0': 421, '1': 15, '2': 171, '3': 1549} = 2156 sites
Census NGG-carrying subset:        {'0': 421, '1': 15, '2': 171, '3': 1549} = 2156 sites
Census PAM-less remainder:         {'0': 6, '1': 212, '2': 2803, '3': 27759} = 30780 sites
Sum check: 2156 + 30780 = 32936 (census total 32,936)

Per-guide per-tier agreement: 1684/1684 cells; orientation-ambiguous census hits: 0.
Verdict: **EXACT AGREEMENT** - the census of record is not an artifact of the BWA alignment engine; its NGG-carrying subset reproduces exactly under GuideScan2's bidirectional-BWT search.
