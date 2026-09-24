#!/usr/bin/env bash
# Fetch AgamP4 evidence for the dsx re-annotation (the ONLY network step for this analysis).
# Sources (fetched 24 Sep 2026):
#   Ensembl Metazoa FTP, current release, AgamP4.63 chromosome 2R GFF3 (VectorBase community models)
#   Ensembl Metazoa FTP, AgamP4 soft-masked chromosome 2R DNA
set -euo pipefail
BASE=https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/current
curl -sL "$BASE/gff3/anopheles_gambiae/Anopheles_gambiae.AgamP4.63.chromosome.2R.gff3.gz" -o /tmp/agamp4_2R.gff3.gz
curl -sL "$BASE/fasta/anopheles_gambiae/dna/Anopheles_gambiae.AgamP4.dna_sm.chromosome.2R.fa.gz" -o /tmp/agamp4_2R.fa.gz
zcat /tmp/agamp4_2R.gff3.gz | grep AGAP004050 > data/agamp4/AGAP004050_AgamP4.gff3
python3 - << 'PY'
import gzip
seq=''.join(l.strip() for l in gzip.open('/tmp/agamp4_2R.fa.gz','rt') if not l.startswith('>')).upper()
sub=seq[48699999:48790000]
with open('data/agamp4/AgamP4_2R_dsx_region.fasta','w') as f:
    f.write('>2R:48700000-48790000 Anopheles gambiae AgamP4 (GCA_000005575.1) dsx region, Ensembl Metazoa dna_sm\n')
    for i in range(0,len(sub),60): f.write(sub[i:i+60]+'\n')
PY
