#!/usr/bin/env bash
set -euo pipefail

# This script fetches all required raw inputs for the Week 6 scRNA-seq workflow.
# It expects two environment variables:
#   BOX_FASTQ_URL     - Direct download link to the FASTQ + reference bundle
#   WHITELIST_URL     - Direct link to whitelist_barcodes.txt.gz

mkdir -p data
cd data

echo "[fetch_data] Starting download..."

if [[ -z "${BOX_FASTQ_URL:-}" ]]; then
    echo "[ERROR] Environment variable BOX_FASTQ_URL is not set."
    echo "Provide a direct-download URL for the Box archive."
    exit 1
fi

if [[ -z "${WHITELIST_URL:-}" ]]; then
    echo "[ERROR] Environment variable WHITELIST_URL is not set."
    echo "Provide a direct-download URL for the whitelist barcodes."
    exit 1
fi

echo "[fetch_data] Downloading FASTQ + reference bundle..."
wget -O sample_bundle.tar.gz "$BOX_FASTQ_URL"

echo "[fetch_data] Downloading whitelist barcodes..."
wget -O whitelist.txt.gz "$WHITELIST_URL"

echo "[fetch_data] Extracting FASTQs and reference files..."
tar -xzf sample_bundle.tar.gz

echo "[fetch_data] Done."
