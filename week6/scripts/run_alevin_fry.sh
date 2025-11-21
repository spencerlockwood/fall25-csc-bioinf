#!/usr/bin/env bash
set -euo pipefail

echo "[alevin-fry] Starting quantification..."

R1="reads_R1.fastq"
R2="reads_R2.fastq"
INDEX="index/splici_index"

if [[ ! -f "$R1" ]]; then echo "[ERROR] Missing $R1"; exit 1; fi
if [[ ! -f "$R2" ]]; then echo "[ERROR] Missing $R2"; exit 1; fi
if [[ ! -d "$INDEX" ]]; then echo "[ERROR] Missing index: $INDEX"; exit 1; fi

mkdir -p af_output

echo "[alevin-fry] Running Salmon alevin (RAD generation)..."
salmon alevin \
    -l ISR \
    -i "$INDEX" \
    --chromium \
    -1 "$R1" \
    -2 "$R2" \
    -p 4 \
    -o af_output/salmon_alevin

echo "[alevin-fry] Permit-list generation..."
alevin-fry generate-permit-list \
    -i af_output/salmon_alevin \
    -o af_output/permit \
    --expect-cells 300

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
