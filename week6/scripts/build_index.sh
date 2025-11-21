#!/usr/bin/env bash
set -euo pipefail

echo "[build_index] Building splici reference..."

mkdir -p index
mkdir -p splici

FA=data/genome.fa
GTF=data/genes.gtf

if [[ ! -f "$FA" ]]; then
    echo "[ERROR] Missing reference FASTA: $FA"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "[ERROR] Missing GTF: $GTF"
    exit 1
fi

echo "[build_index] Generating splici transcriptome with salmon..."
salmon splici \
    -r "$FA" \
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
