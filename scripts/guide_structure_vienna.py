"""Guide-spacer secondary-structure audit with ViennaRNA 2.7.2.

Question: do the three v3 lead spacers carry unusually stable self-structure
(MFE) relative to a genomic null? Stable spacer hairpins are a documented
Cas9-efficacy confound, so a lead far above the null would be flagged.
"""
import csv, json, random, statistics
import RNA

random.seed(27)

def mfe(seq):
    struct, energy = RNA.fold(seq)
    return struct, energy

# Leads + Kyrou control
leads = []
with open('results/ranked_designs_v3.csv') as fh:
    for row in csv.DictReader(fh):
        if row['rank_v3'] and int(row['rank_v3']) <= 3:
            leads.append({'name': f"v3-{row['rank_v3']}", 'spacer': row['protospacer'],
                          'efficacy': float(row['efficacy'])})
kyrou = {'name': 'kyrou2018-control', 'spacer': 'GTTTAACACAGGTCAAGCGG', 'efficacy': None}

# Genomic null: 500 random 20-mers from chromosome 2R (same background as sieve)
seq = []
with open('/tmp/chr2R.fasta') as fh:
    for line in fh:
        if not line.startswith('>'):
            seq.append(line.strip())
chr2r = ''.join(seq).upper()
null_mfes = []
for _ in range(500):
    i = random.randrange(0, len(chr2r) - 20)
    kmer = chr2r[i:i+20]
    if 'N' in kmer:
        continue
    _, e = mfe(kmer)
    null_mfes.append(e)

mu = statistics.mean(null_mfes)
sd = statistics.pstdev(null_mfes)

records = []
for rec in leads + [kyrou]:
    struct, e = mfe(rec['spacer'])
    z = (e - mu) / sd if sd > 0 else 0.0
    records.append({**rec, 'mfe_structure': struct, 'mfe_kcal_mol': e,
                    'null_z': round(z, 3),
                    'structured_flag': bool(z < -2.0)})

out = {
    'tool': 'ViennaRNA', 'tool_version': '2.7.2', 'method': 'RNA.fold MFE (DNA alphabet spacer treated as RNA)',
    'null': {'source': '500 random 20-mers from NC_064601.1 (2R)', 'n': len(null_mfes),
             'mean_mfe_kcal_mol': round(mu, 3), 'sd_kcal_mol': round(sd, 3)},
    'records': records,
    'verdict': 'no lead spacer is a z<-2 outlier vs the genomic null'
               if not any(r['structured_flag'] for r in records[:3])
               else 'at least one lead spacer is unusually structured',
}
with open('results/guide_structure_vienna.json', 'w') as fh:
    json.dump(out, fh, indent=2)

lines = ['# Guide-spacer secondary-structure audit (ViennaRNA 2.7.2)', '',
         f"Null: {len(null_mfes)} random 2R 20-mers, mean MFE {mu:.3f} kcal/mol, sd {sd:.3f}.", '',
         '| guide | spacer | MFE (kcal/mol) | null z | structured (z<-2) |',
         '|---|---|---|---|---|']
for r in records:
    lines.append(f"| {r['name']} | {r['spacer']} | {r['mfe_kcal_mol']:.2f} | {r['null_z']} | {r['structured_flag']} |")
lines += ['', f"Verdict: {out['verdict']}.",
          'Method: RNA.fold MFE on the 20-nt spacer; spacers treated as RNA (the species loaded into Cas9).']
with open('results/guide_structure_vienna.md', 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
print('\n'.join(lines))
