"""Independent re-score of archived Ag1000G modern VCF-slice JSONL.

Uses the archived AgamP4 FASTA and AgamP4 reannotation positions, not the
allele_freq implementation's reference constants or position/status helpers.
This is a second computational interpretation of the *same* archived sample
calls, not another population or independent VCF extraction.
"""
import collections
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEQ = ''.join(line.strip().upper() for line in
              (ROOT / 'data/agamp4/AgamP4_2R_dsx_region.fasta').read_text().splitlines()
              if not line.startswith('>'))
ANN = json.loads((ROOT / 'results/agamp4_reannotation.json').read_text())['guides']
PAMS = {'kyrou_2018': (48714637, 48714639),
        'dsx-v3-1': (48711498, 48711500),
        'dsx-v3-2': (48712762, 48712764)}
NAMES = {'kyrou_2018': 'kyrou', 'dsx-v3-1': 'dsx-v3-1', 'dsx-v3-2': 'dsx-v3-2'}
COMPLEMENT = dict(zip('ACGT', 'TGCA'))


def evaluate(row):
    results = {}
    calls = {h['pos']: h for h in row['target_hits']}
    for original_name, info in ANN.items():
        if original_name not in PAMS:
            continue
        name = NAMES[original_name]
        left, right = info['assembly_position_2R']
        pam_left, pam_right = PAMS[original_name]
        positions = list(range(left, right + 1)) + [pam_left, pam_left + 1]
        guide_plus = info['guide'].translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        states = []
        for pos in positions:
            ref = SEQ[pos - 48700000]
            expected = guide_plus[pos - left] if pos >= left else 'C'
            h = calls.get(pos)
            if h is None:
                alleles = (ref, ref)  # all-sites VCF omitted this hom-ref row
            else:
                assert h['ref'] == ref, (row['sample'], pos, h['ref'], ref)
                if h['gq'] is None or h['dp'] is None or h['gq'] < 20 or h['dp'] < 5:
                    states.append(None)
                    continue
                genotypes = h['gt'].replace('|', '/').split('/')
                if len(genotypes) != 2 or '.' in genotypes:
                    states.append(None)
                    continue
                palette = [ref] + h['alt'].split(',')
                alleles = tuple(palette[int(i)] for i in genotypes)
            states.append(sum(base != expected for base in alleles))
        results[name] = 2 if 2 in states else None if None in states else max(states)
    return results


def main():
    tally = {n: collections.Counter() for n in NAMES.values()}
    n = 0
    for line in (ROOT / 'results/ag1000g/on_target.jsonl').open():
        row = json.loads(line)
        if row.get('n_rows') != 3250:
            continue
        n += 1
        for name, status in evaluate(row).items():
            tally[name][status] += 1
    assert n == 4106
    summary = {}
    for name, c in sorted(tally.items()):
        callable_n = sum(c[i] for i in (0, 1, 2))
        summary[name] = {'intact': c[0], 'het': c[1], 'hom_compromised': c[2],
                         'unknown': c[None], 'n_callable': callable_n,
                         'compromised_allele_fraction': (2*c[2]+c[1])/(2*callable_n)}
    expected = json.loads((ROOT / 'results/ag1000g/allele_freq.json').read_text())['guide_summary']
    for name, fields in summary.items():
        for field in ('intact', 'het', 'hom_compromised', 'n_callable', 'compromised_allele_fraction'):
            assert fields[field] == expected[name][field], (name, field, fields[field], expected[name][field])
    out = {'source': 'same archived sample JSONL, second interpretation using archived FASTA',
           'not_independent_extraction_or_population': True, 'n_modern': n,
           'guide_summary': summary}
    p = ROOT / 'results/ag1000g/independent_af_check.json'
    p.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))

if __name__ == '__main__':
    main()
