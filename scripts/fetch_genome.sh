#!/usr/bin/env bash
# Live NCBI pulls - run outside CI. Full AgamP5 reference (GCF_943734735.2,
# idAnoGambNW_F1_1): chromosome 2RL carries the dsx locus; X, 3RL and MT
# complete the nuclear genome. Output: one multi-record FASTA.
set -euo pipefail
OUT="${1:-/tmp/agam_genome.fasta}"
TMP="$(mktemp -d)"
for acc in NC_064601.1 OX030909.1 OX030908.1 OX030910.2; do
  curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" -o "$TMP/$acc.fasta"
done
cat "$TMP"/NC_064601.1.fasta "$TMP"/OX030909.1.fasta "$TMP"/OX030908.1.fasta "$TMP"/OX030910.2.fasta > "$OUT"
rm -rf "$TMP"
echo "wrote $OUT"
