#!/usr/bin/env bash
set -euo pipefail

# This script unpacks the local dataset and places files into stable locations.

echo "[fetch_data] Using local dataset: data/toy_read_ref_set.tar.gz"

cd data

echo "[fetch_data] Extracting dataset..."
tar -xzf toy_read_ref_set.tar.gz

echo "[fetch_data] Organizing reference files..."
cp toy_ref_read/toy_human_ref/fasta/genome.fa ../genome.fa
cp toy_ref_read/toy_human_ref/genes/genes.gtf ../genes.gtf

echo "[fetch_data] Organizing FASTQs..."
cp toy_ref_read/toy_read_fastq/selected_R1_reads.fastq ../reads_R1.fastq
cp toy_ref_read/toy_read_fastq/selected_R2_reads.fastq ../reads_R2.fastq

cd ..
echo "[fetch_data] Done."
