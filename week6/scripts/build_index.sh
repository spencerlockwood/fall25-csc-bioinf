#!/usr/bin/env bash
set -euo pipefail

# Build splici reference & Salmon index using alevin-fry generate-splici

echo "[build_index] Building splici reference..."

FA="genome.fa"
GTF="genes.gtf"

if [[ ! -f "$FA" ]]; then echo "[ERROR] Missing $FA"; exit 1; fi
if [[ ! -f "$GTF" ]]; then echo "[ERROR] Missing $GTF"; exit 1; fi

mkdir -p splici
mkdir -p index

echo "[build_index] Generating splici transcriptome with alevin-fry..."
alevin-fry generate-splici \
    -f "$FA" \
    -g "$GTF" \
    -o splici \
    --flank-trim-length 10

echo "[build_index] Building Salmon index..."
salmon index \
    -t splici/splici.fa \
    -i index/splici_index \
    -k 31

echo "[build_index] Done."
