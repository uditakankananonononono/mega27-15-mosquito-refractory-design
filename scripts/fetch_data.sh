#!/usr/bin/env bash
# Live NCBI pulls - run outside CI. Anopheles gambiae doublesex (AGAP004050,
# NCBI Gene 1270904) RefSeq transcripts via elink.
set -euo pipefail
cd "$(dirname "$0")/../data"
IDS=$(curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=nuccore&id=1270904&retmode=json" \
  | python3 -c "import json,sys; d=json.load(sys.stdin); print(','.join(l for ls in d['linksets'] for ldb in ls['linksetdbs'] if ldb['linkname'].startswith('gene_nuccore') for l in ldb['links']))")
curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${IDS}&rettype=fasta&retmode=text" -o dsx_all.fasta
python3 - <<'PY'
seqs, name, buf = {}, None, []
for line in open('dsx_all.fasta'):
    line = line.rstrip()
    if line.startswith('>'):
        if name: seqs[name] = ''.join(buf)
        name, buf = line, []
    else: buf.append(line)
if name: seqs[name] = ''.join(buf)
mrna = {k: v for k, v in seqs.items() if k.startswith(('>XM_', '>NM_'))}
with open('dsx_transcripts.fasta', 'w') as f:
    for k, v in mrna.items():
        f.write(k + '\n')
        for i in range(0, len(v), 60): f.write(v[i:i+60] + '\n')
print(f"kept {len(mrna)} transcripts")
PY
rm -f dsx_all.fasta
echo "Fetched dsx transcripts."
