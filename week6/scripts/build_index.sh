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

# Modern splici reference construction
salmon alevin generate-splici \
    --transcriptome $GENOME \
    --gtf $GTF \
    --output $SPLICI_DIR \
    --flank-length 91

echo "[build_index] Building Salmon index with pufferfish..."

salmon index \
    -t $SPLICI_DIR/splici.fa \
    -d $SPLICI_DIR/splici.dfa \
    -i $INDEX_DIR \
    -p 4

echo "[build_index] Done."
