"""CRISPRoff/CRISPRspec energy-based specificity audit (tool 27/40).

Runs the RTH-tools/crisproff v1.1.2 pipeline (Alkan et al. 2018, Genome Biol
19:177) on the four dsx lead guides, using the CRISPOR AgamP4 off-target sets
(committed lineage: results/offtarget_crispor.json) as the near-cognate input.
Pipeline was validated against the repo's own test.sh (diff-clean vs
test_data.example.out) before use. RNAfold CLI requirement is satisfied by a
thin shim over the ViennaRNA 2.7.2 python bindings already used for
results/guide_structure_vienna.json (same library, same fold algorithm).
"""
import csv, json, subprocess, sys, os

CRISPROFF_DIR = os.environ.get('CRISPROFF_DIR', '/tmp/crisproff-master')
RUN_DIR = os.environ.get('CRISPROFF_RUN_DIR', '/tmp/crisproff_run')

LEADS = {
    'dsx-v3-1': {'gid': '341rev', 'spacer': 'TGGGCAGTATGCGTTAGGGT', 'win': 'w1'},
    'dsx-v3-2': {'gid': '1605rev', 'spacer': 'CATTAAGACCTACGAAGCGC', 'win': 'w1'},
    'dsx-v3-3': {'gid': '797rev', 'spacer': 'GAAGCGAGCCCAATGGCTGT', 'win': 'w2'},
    'kyrou':    {'gid': '175rev', 'spacer': 'GTTTAACACAGGTCAAGCGG', 'win': 'w2'},
}

def build_inputs():
    summary = {}
    for name, L in LEADS.items():
        tsv = f"/tmp/crispor_{L['win']}_ot.tsv"
        gseq, ots = None, []
        with open(tsv) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                if row['guideId'] != L['gid']:
                    continue
                gseq = row['guideSeq'].upper()
                oseq = row['offtargetSeq'].upper()
                if oseq != gseq:
                    ots.append((oseq, int(row['mismatchCount'])))
        assert gseq and gseq[:20] == L['spacer'], f'spacer verify failed {name}'
        seqs = [gseq] + sorted({s for s, _ in ots})
        fn = os.path.join(RUN_DIR, f'{name}_offtargets.txt')
        with open(fn, 'w') as out:
            out.write('\n'.join(seqs) + '\n')
        mmhist = {}
        for _, m in ots:
            mmhist[str(m)] = mmhist.get(str(m), 0) + 1
        summary[name] = {'on_target_23mer': gseq, 'n_offtargets': len(seqs) - 1,
                         'mm_hist': dict(sorted(mmhist.items())), 'file': fn}
    return summary

def run_pipeline(summary):
    out = {}
    for name, info in summary.items():
        g = info['on_target_23mer']
        spec = os.path.join(RUN_DIR, f'{name}_spec.tsv')
        cmd = ['python3', os.path.join(CRISPROFF_DIR, 'CRISPRspec_CRISPRoff_pipeline.py'),
               '--guide', g, '--offtargets', info['file'], '--no_azimuth',
               '--duplex_energy_params', os.path.join(CRISPROFF_DIR, 'energy_dics.pkl'),
               '--specificity_report', spec,
               '--CRISPRoff_scores_folder', RUN_DIR]
        env = dict(os.environ, PATH=CRISPROFF_DIR + ':' + os.environ['PATH'])
        r = subprocess.run(cmd, cwd=RUN_DIR, env=env, capture_output=True, text=True)
        assert r.returncode == 0, f'{name}: {r.stderr[-400:]}'
        row = list(csv.DictReader(open(spec), delimiter='\t'))[0]
        ot_tsv = os.path.join(RUN_DIR, f'{g}.CRISPRoff.tsv')
        per_ot = []
        for line in open(ot_tsv):
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) >= 5:
                per_ot.append({'seq': c[3], 'crisproff_score': float(c[4])})
        out[name] = {
            'on_target_23mer': g,
            'n_offtargets_input': info['n_offtargets'],
            'mm_hist_input': info['mm_hist'],
            'crisprspec_specificity_score': float(row['CRISPRspec_specificity_score']),
            'mm_counts': row['MM_counts'],
            'per_offtarget_crisproff': per_ot,
        }
    return out

if __name__ == '__main__':
    s = build_inputs()
    res = run_pipeline(s)
    json.dump({'tool': 'CRISPRoff/CRISPRspec pipeline', 'tool_version': '1.1.2',
               'source': 'github.com/RTH-tools/crisproff',
               'citation': 'Alkan et al. 2018, Genome Biol 19:177',
               'validation': "repo test.sh reproduces test_data.example.out diff-clean "
                             "with RNAfold shim over ViennaRNA 2.7.2 python bindings",
               'inputs': 'CRISPOR AgamP4 off-target sets per lead (results/offtarget_crispor.json lineage)',
               'leads': res}, open(sys.argv[1], 'w'), indent=1)
    for k, v in res.items():
        print(k, 'CRISPRspec=', round(v['crisprspec_specificity_score'], 3),
              'n_ot=', v['n_offtargets_input'])
