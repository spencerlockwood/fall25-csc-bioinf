#!/usr/bin/env bash
set -euo pipefail

# This version no longer downloads from Box.
# Instead, it uses the toy_read_ref_set.tar.gz already stored in week6/data/.

cd "$(dirname "$0")/../data"

echo "[fetch_data] Using local dataset: toy_read_ref_set.tar.gz"

if [[ ! -f toy_read_ref_set.tar.gz ]]; then
    echo "[ERROR] Dataset not found: data/toy_read_ref_set.tar.gz"
    exit 1
fi

echo "[fetch_data] Extracting local dataset..."
tar -xzf toy_read_ref_set.tar.gz

# After extraction, dataset is under toy_ref_read/
# Move + rename into expected filenames for build_index.sh and run_alevin_fry.sh

echo "[fetch_data] Organizing reference files..."

REF_DIR="toy_ref_read/toy_human_ref"

# Reference
cp "$REF_DIR/fasta/genome.fa" chr5.fa
cp "$REF_DIR/genes/genes.gtf" annotations.gtf

echo "[fetch_data] Organizing FASTQs..."

FASTQ_DIR="toy_ref_read/toy_read_fastq"

# Gzip FASTQs so run_alevin_fry.sh matches R1=*.fastq.gz
gzip -c "$FASTQ_DIR/selected_R1_reads.fastq" > selected_R1_reads.fastq.gz
gzip -c "$FASTQ_DIR/selected_R2_reads.fastq" > selected_R2_reads.fastq.gz

# Download whitelist
echo "[fetch_data] Downloading whitelist..."
wget -O whitelist.txt.gz "$WHITELIST_URL"

echo "[fetch_data] Done."
