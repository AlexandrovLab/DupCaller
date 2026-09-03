#!/usr/bin/env bash
# Runs the full DupCaller pipeline (index -> trim -> bwa mem -> gatk
# MarkDuplicates -> call -> estimate) against the synthetic dataset in
# data/, entirely inside OUTDIR. Used both by test_mock_pipeline.py (to
# validate a fresh install against the premade expected/ outputs) and to
# regenerate expected/ after a deliberate change to data/ or to DupCaller
# itself.
#
# Usage: run_pipeline.sh OUTDIR
#
# Requires DupCaller.py, bwa, samtools, and gatk on PATH (or overridden via
# the DUPCALLER/BWA/SAMTOOLS/GATK env vars).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTDIR="${1:?usage: run_pipeline.sh OUTDIR}"

DUPCALLER="${DUPCALLER:-DupCaller.py}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
GATK="${GATK:-gatk}"

mkdir -p "$OUTDIR"
cp "$DATA_DIR/reference.fa" "$DATA_DIR/repeats.tsv" "$DATA_DIR/mock_1.fastq" "$DATA_DIR/mock_2.fastq" "$OUTDIR/"
cd "$OUTDIR"

echo "[1/6] index"
"$SAMTOOLS" faidx reference.fa
"$DUPCALLER" index -f reference.fa -rt repeats.tsv

echo "[2/6] trim"
"$DUPCALLER" trim -i mock_1.fastq -i2 mock_2.fastq -p NNNXXXX -o mock_trm

echo "[3/6] align"
"$BWA" index reference.fa
"$BWA" mem -C -R "@RG\tID:mock\tSM:mock\tPL:ILLUMINA" reference.fa mock_trm_1.fastq mock_trm_2.fastq \
    | "$SAMTOOLS" sort -o mock.bam -
"$SAMTOOLS" index mock.bam

echo "[4/6] mark duplicates"
"$GATK" MarkDuplicates \
    -I mock.bam -O mock.mkdped.bam -M mock.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
"$SAMTOOLS" index mock.mkdped.bam

echo "[5/6] call"
"$DUPCALLER" call -b mock.mkdped.bam -f reference.fa -o result/result \
    -r mockchr1 -p 1 -w 5000 --seed 1

echo "[6/6] estimate"
"$DUPCALLER" estimate -i result/result -f reference.fa -r mockchr1

echo "done"
