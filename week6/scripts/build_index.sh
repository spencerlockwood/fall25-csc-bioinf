#!/usr/bin/env bash
set -euo pipefail

echo "[build_index] Building splici reference..."

GENOME=data/genome.fa
GTF=data/genes.gtf
SPLICI_DIR=splici
INDEX_DIR=index

if [[ ! -f "$GENOME" ]]; then
    echo "Error:  Missing reference FASTA: $GENOME"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "Error:  Missing GTF: $GTF"
    exit 1
fi

mkdir -p $SPLICI_DIR
mkdir -p $INDEX_DIR

echo "[build_index] Generating splici transcriptome with salmon..."

# Salmon 1.10.x syntax (no --transcriptome)
salmon alevin generate-splici \
    -r $GENOME \
    -g $GTF \
    -o $SPLICI_DIR \
    -f 91

echo "[build_index] Building Salmon index with pufferfish..."

salmon index \
    -t $SPLICI_DIR/splici.fa \
    -d $SPLICI_DIR/splici.dfa \
    -i $INDEX_DIR \
    -p 4

echo "[build_index] Done."
