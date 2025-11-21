#!/usr/bin/env bash
set -euo pipefail

echo "[fetch_data] Using local dataset: toy_read_ref_set.tar.gz"

mkdir -p data

# Extract
echo "[fetch_data] Extracting dataset..."
tar -xzf data/toy_read_ref_set.tar.gz -C data

# Move reference files into expected locations
REF_DIR="data/toy_ref_read/toy_human_ref"

FASTA_SRC="$REF_DIR/fasta/genome.fa"
GTF_SRC="$REF_DIR/genes/genes.gtf"

FASTA_DST="data/genome.fa"
GTF_DST="data/genes.gtf"

echo "[fetch_data] Organizing reference files..."
mkdir -p data
mkdir -p data

# Move genome
if [[ -f "$FASTA_SRC" ]]; then
    cp "$FASTA_SRC" "$FASTA_DST"
else
    echo "[fetch_data] ERROR: fasta not found at $FASTA_SRC"
    exit 1
fi

# Move genes.gtf
if [[ -f "$GTF_SRC" ]]; then
    cp "$GTF_SRC" "$GTF_DST"
else
    echo "[fetch_data] ERROR: gtf not found at $GTF_SRC"
    exit 1
fi

# Organize FASTQs
echo "[fetch_data] Organizing FASTQs..."
FASTQ_DIR="data/toy_ref_read/toy_read_fastq"
mkdir -p fastq

cp "$FASTQ_DIR/selected_R1_reads.fastq" fastq/reads_1.fastq
cp "$FASTQ_DIR/selected_R2_reads.fastq" fastq/reads_2.fastq

echo "[fetch_data] Done."
