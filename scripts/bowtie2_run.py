#!/usr/bin/env python3
"""bowtie2 2.5.4 whole-genome off-target audit of the 421-guide dsx library (tool 29/40).

Independent FM-index aligner cross-check of FlashFry 2.0's AgamP5 off-target
histograms (tool 28/40) and the earlier bowtie1 sweep (offtarget_bowtie.json).

Pipeline (exact commands; index build ~9 min on the 2-core sandbox):
  bowtie2-build /tmp/agam_genome.fasta /tmp/bowtie2_idx/agamp5 --threads 2
  bowtie2 -x /tmp/bowtie2_idx/agamp5 -f -U /tmp/bt2_guides.fa --end-to-end \
      --score-min L,0,-0.9 --rdg 1000,1000 --rfg 1000,1000\n      --very-sensitive -D 100 -R 5 -L 14 -N 1 -i S,1,0.50 -a -p 2 \
      --no-unal -S /tmp/bt2_hits.sam
  bowtie2 has no -v mode; --score-min L,0,-0.9 with default --mp 6 caps a
  20-mer at exactly 3 mismatches (-18), and prohibitive gap penalties
  (--rdg/--rfg 1000,1000) keep the search ungapped -- the same search space
  as the bowtie1 -v 3 sweep and the MM0-3 slice of FlashFry's histogram.
    -a        report ALL alignments (like bowtie1 -a)
Input: results/flashfry_scored_guides.tsv -> 421 unique 20-mer protospacers
(PAM stripped; window-relative ids contig_start_stop_orientation).
Search space is PAM-agnostic (like bowtie1), so per-guide counts are expected
to be >= FlashFry's PAM-adjacent counts; the comparison of record is bowtie1.
"""
import csv, json, collections, os

SAM = '/tmp/bt2_hits2.sam'
TSV = os.path.join(os.path.dirname(__file__), '..', 'results', 'flashfry_scored_guides.tsv')
OUTJ = os.path.join(os.path.dirname(__file__), '..', 'results', 'offtarget_bowtie2.json')
OUTM = os.path.join(os.path.dirname(__file__), '..', 'results', 'offtarget_bowtie2.md')

# on-target lookup: window-relative id -> (contig, start) from the TSV
trows = list(csv.DictReader(open(TSV), delimiter='\t'))
by_id = {f"{r['contig']}_{r['start']}_{r['stop']}_{r['orientation']}": r for r in trows}

hits = collections.defaultdict(list)   # gid -> list of (nm, rname, pos)
with open(SAM) as fh:
    for line in fh:
        if line.startswith('@'):
            continue
        f = line.rstrip('\n').split('\t')
        gid, flag, rname, pos = f[0], int(f[1]), f[2], int(f[3])
        nm = next((int(t.split(':')[-1]) for t in f[11:] if t.startswith('NM:i:')), None)
        hits[gid].append((nm, rname, pos))

per_guide = {}
tot_mm = collections.Counter()
no_self = []
for gid, r in by_id.items():
    hs = hits.get(gid, [])
    mm_hist = collections.Counter(nm for nm, _, _ in hs)
    self_hits = [h for h in hs if h[0] == 0]
    if not self_hits:
        no_self.append(gid)
    ot = len(hs) - len(self_hits)          # PAM-agnostic off-targets, <=3 MM
    ot_mm = collections.Counter(nm for nm, _, _ in hs if not (nm == 0))
    for nm, c in ot_mm.items():
        tot_mm[nm] += c
    per_guide[gid] = {
        'contig': r['contig'], 'start': int(r['start']), 'orientation': r['orientation'],
        'self_hits_mm0': len(self_hits), 'offtargets_le3mm_pam_agnostic': ot,
        'ot_mm_hist': {str(k): ot_mm[k] for k in sorted(ot_mm)},
        'flashfry_otcount': int(r['otCount']),
        'flashfry_mm_0_4': r['0-1-2-3-4_mismatch'],
    }

# lead lookup by spacer
sp2gid = {r['target'][:20]: f"{r['contig']}_{r['start']}_{r['stop']}_{r['orientation']}" for r in trows}
lead_spacers = {'v3-1': 'TGGGCAGTATGCGTTAGGGT', 'v3-2': 'CATTAAGACCTACGAAGCGC',
                'v3-3': 'GAAGCGAGCCCAATGGCTGT', 'kyrou': 'GTTTAACACAGGTCAAGCGG'}
lead_rows = {}
for name, sp in lead_spacers.items():
    gid = sp2gid.get(sp)
    lead_rows[name] = {'guide_id': gid, **(per_guide.get(gid) or {})}

lib_ot = sum(p['offtargets_le3mm_pam_agnostic'] for p in per_guide.values())
summary = {
    'tool': 'bowtie2 2.5.4 (BenLangmead bowtie2, FM-index; Langmead & Salzberg 2012 Nat Methods 9:357)',
    'index': 'AgamP5 (4 molecules incl. NC_064601.1 whole chr 2), bowtie2-build --threads 2',
    'align': '-f -U --end-to-end --score-min L,0,-0.9 --rdg/--rfg 1000,1000 --very-sensitive -D 100 -R 5 -L 14 -N 1 -i S,1,0.50 -a -p 2 --no-unal (end-to-end, all alignments, <=3 MM, ungapped, PAM-agnostic; first default-preset pass under-reported 2-3MM sites and was superseded)',
    'guides': len(by_id),
    'guides_with_self_hit': len(by_id) - len(no_self),
    'guides_missing_self_hit': no_self,
    'library_offtargets_le3mm_pam_agnostic': lib_ot,
    'library_ot_mm_hist': {str(k): tot_mm[k] for k in sorted(tot_mm)},
    'median_ot_per_guide': sorted(p['offtargets_le3mm_pam_agnostic'] for p in per_guide.values())[len(per_guide)//2],
    'max_ot_per_guide': max(p['offtargets_le3mm_pam_agnostic'] for p in per_guide.values()),
    'leads': lead_rows,
    'note': 'PAM-agnostic counts >= FlashFry PAM-adjacent counts by construction; comparison of record is the bowtie1 sweep (same search space).',
}
json.dump({'summary': summary, 'per_guide': per_guide}, open(OUTJ, 'w'), indent=1)

with open(OUTM, 'w') as m:
    m.write('# bowtie2 2.5.4 whole-genome off-target audit (tool 29/40)\n\n')
    m.write(f"- guides aligned: {summary['guides']} (unique 20-mers), self-hit found for {summary['guides_with_self_hit']}\n")
    m.write(f"- library off-targets <=3 MM (PAM-agnostic): {lib_ot:,}; MM hist {dict(sorted(tot_mm.items()))}\n")
    m.write(f"- median {summary['median_ot_per_guide']}, max {summary['max_ot_per_guide']} off-targets per guide\n\n")
    m.write('## Leads\n\n')
    for name, lr in lead_rows.items():
        m.write(f"- {name}: {lr.get('offtargets_le3mm_pam_agnostic')} OTs <=3MM (hist {lr.get('ot_mm_hist')}), "
                f"self-hits {lr.get('self_hits_mm0')}; FlashFry PAM-adjacent otCount {lr.get('flashfry_otcount')}, "
                f"FF MM0-4 {lr.get('flashfry_mm_0_4')}\n")
print(json.dumps(summary, indent=1)[:2000])
