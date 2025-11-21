#!/usr/bin/env bash
set -euo pipefail

R1=data/selected_R1_reads.fastq.gz
R2=data/selected_R2_reads.fastq.gz
WL=data/whitelist.txt.gz
INDEX=index/splici_index

if [[ ! -f $WL ]]; then
    echo "[ERROR] Whitelist not found: $WL"
    exit 1
fi

if [[ ! -d $INDEX ]]; then
    echo "[ERROR] Splici index not found: $INDEX"
    exit 1
fi

mkdir -p af_output

echo "[alevin-fry] Generating RAD file with Salmon..."
salmon alevin \
    -l ISR \
    -i "$INDEX" \
    --chromium \
    -1 "$R1" \
    -2 "$R2" \
    --whitelist "$WL" \
    -p 8 \
    -o af_output/salmon_alevin

echo "[alevin-fry] Generating permit list..."
alevin-fry generate-permit-list \
    -i af_output/salmon_alevin \
    -o af_output/permit \
    --inspect \
    --expect-cells 3000

echo "[alevin-fry] Collating RAD..."
alevin-fry collate \
    -i af_output/salmon_alevin \
    -r af_output/permit \
    -o af_output/collated

echo "[alevin-fry] Quantifying..."
alevin-fry quant \
    -i "$INDEX" \
    -r af_output/collated \
    -o af_output/quant \
    --use-mtx \
    --resolution cr-like

echo "[alevin-fry] Done."
