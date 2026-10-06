#!/usr/bin/env bash
# Runs the full DupCaller pipeline (index -> trim -> bwa-mem2 mem -> gatk
# MarkDuplicates -> call -> estimate) against the synthetic dataset in
# data/, entirely inside OUTDIR. Used both by test_mock_pipeline.py (to
# validate a fresh install against the premade expected/ outputs) and to
# regenerate expected/ after a deliberate change to data/ or to DupCaller
# itself.
#
# Usage: run_pipeline.sh OUTDIR
#
# R2_SHORT=N (optional) cuts N bases off the 3' end of every read 2 with
# make_short_read2.py, so read 2 is N bases shorter than read 1 after trim;
# the outputs must still reproduce expected/ exactly.
#
# BARCODE_LIST=1 (optional) rewrites the reads with make_barcode_list_reads.py
# (variable-length listed codes + Tn5 ME instead of NNNXXXX) and trims with
# -p B + 19 X -bl barcodes.txt; outputs must match expected/ once the codes
# are mapped back to the original barcodes (barcode_map.tsv).
#
# MATCHED_NORMAL=1 (optional) runs call/estimate exactly as the Nextflow
# pipeline's CALL_VARIANTS/ESTIMATE_BURDEN do for the mock data in
# run_nextflow_pipeline.sh's "default" case (the same reads as tumor and
# matched normal, Nextflow's call options, --seed 1), writing OUTDIR/mock.
# Its outputs are compared against expected_matched_normal/, the same files
# the Nextflow test compares the containerized pipeline against.
#
# Requires DupCaller.py, bwa-mem2, samtools, and gatk on PATH (or overridden
# via the DUPCALLER/BWA_MEM2/SAMTOOLS/GATK env vars).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTDIR="${1:?usage: run_pipeline.sh OUTDIR}"

DUPCALLER="${DUPCALLER:-DupCaller.py}"
BWA_MEM2="${BWA_MEM2:-bwa-mem2}"
SAMTOOLS="${SAMTOOLS:-samtools}"
GATK="${GATK:-gatk}"

mkdir -p "$OUTDIR"
cp "$DATA_DIR/reference.fa" "$DATA_DIR/repeats.tsv" "$DATA_DIR/mock_1.fastq" "$DATA_DIR/mock_2.fastq" "$OUTDIR/"
cd "$OUTDIR"
if [ -n "${R2_SHORT:-}" ]; then
    "${PYTHON:-python3}" "${SCRIPT_DIR}/make_short_read2.py" "$DATA_DIR/mock_2.fastq" mock_2.fastq "$R2_SHORT"
fi

echo "[1/6] index"
"$SAMTOOLS" faidx reference.fa
"$DUPCALLER" index -f reference.fa -rt repeats.tsv

echo "[2/6] trim"
if [ -n "${BARCODE_LIST:-}" ]; then
    "${PYTHON:-python3}" "${SCRIPT_DIR}/make_barcode_list_reads.py" mock_1.fastq mock_2.fastq \
        mock_bl_1.fastq mock_bl_2.fastq barcodes.txt barcode_map.tsv
    "$DUPCALLER" trim -i mock_bl_1.fastq -i2 mock_bl_2.fastq -p "BXXXXXXXXXXXXXXXXXXX" \
        -bl barcodes.txt -o mock_trm
else
    "$DUPCALLER" trim -i mock_1.fastq -i2 mock_2.fastq -p NNNXXXX -o mock_trm
fi

echo "[3/6] align"
"$BWA_MEM2" index reference.fa
"$BWA_MEM2" mem -C -T 0 -R "@RG\tID:mock\tSM:mock\tPL:ILLUMINA" reference.fa mock_trm_1.fastq mock_trm_2.fastq \
    | "$SAMTOOLS" sort -o mock.bam -
"$SAMTOOLS" index mock.bam

echo "[4/6] mark duplicates"
"$GATK" MarkDuplicates \
    -I mock.bam -O mock.mkdped.bam -M mock.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
"$SAMTOOLS" index mock.mkdped.bam

if [ -n "${MATCHED_NORMAL:-}" ]; then
    # Same options as nextflow/DupCaller.nf CALL_VARIANTS + ESTIMATE_BURDEN
    # with tests/mock_pipeline/nextflow_test.config.
    echo "[5/6] call (matched normal, Nextflow options)"
    "$DUPCALLER" call -b mock.mkdped.bam -n mock.mkdped.bam -f reference.fa -o mock \
        -p 1 -r mockchr1 --seed 1 -maf 1.0 -gaf 0.001 -d 1 -tt 7 -tr 7 -mq 0

    echo "[6/6] estimate"
    "$DUPCALLER" estimate -i mock -f reference.fa -r mockchr1
else
    echo "[5/6] call"
    "$DUPCALLER" call -b mock.mkdped.bam -f reference.fa -o result/result \
        -r mockchr1 -p 1 -w 5000 --seed 1

    echo "[6/6] estimate"
    "$DUPCALLER" estimate -i result/result -f reference.fa -r mockchr1
fi

echo "done"
