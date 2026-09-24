#!/usr/bin/env bash
# Live NCBI pull - run outside CI. AgamP5 chromosome 2 (NC_064601.1, ~118 Mb),
# home chromosome of the dsx locus (2R).
set -euo pipefail
OUT="${1:-/tmp/chr2R.fasta}"
curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_064601.1&rettype=fasta&retmode=text" -o "$OUT"
echo "wrote $OUT"
