#!/usr/bin/env bash
set -euo pipefail

# Build a splici (spliced + intronic) transcriptome index for Alevin-fry.
# Requires:
#   data/chr5.fa
#   data/annotations.gtf
#
# Produces:
#   index/splici_index/

echo "[build_index] Building splici reference..."

mkdir -p index
mkdir -p splici

FA=data/chr5.fa
GTF=data/annotations.gtf

if [[ ! -f "$FA" ]]; then
    echo "[ERROR] Missing reference FASTA: $FA"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "[ERROR] Missing GTF: $GTF"
    exit 1
fi

# Create splici reference using alevin-fry helper script
# This generates:
#   splici/transcriptome_spliced.fa
#   splici/intronic.fa
#   splici/splici.fa
echo "[build_index] Generating splici transcriptome..."
alevin-fry create-splici \
    -f "$FA" \
    -g "$GTF" \
    -o splici \
    --flank-trim-length 5

echo "[build_index] Building Salmon index..."
salmon index \
    -t splici/splici.fa \
    -i index/splici_index \
    --type quasi \
    -k 31

echo "[build_index] Done."
