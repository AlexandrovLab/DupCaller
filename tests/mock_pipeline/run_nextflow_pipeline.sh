#!/usr/bin/env bash
# Runs nextflow/DupCaller.nf (index -> trim -> bwa mem -> gatk
# MarkDuplicates -> call -> estimate, containerized) against the synthetic
# dataset in data/, entirely inside OUTDIR. Companion to run_pipeline.sh,
# which runs the same tools directly on PATH instead of through Nextflow
# and containers -- this one instead exercises the Nextflow DAG itself
# (staging, container images, per-process resource directives) using the
# published yuhecheng62/dupcaller image.
#
# The pipeline's CALL_VARIANTS step always wants a matched normal (no
# tumor-only mode), unlike run_pipeline.sh's data/ usage -- so mock_1/2 are
# staged as both tumor and normal fastqs here. That means SBS/indel calls
# will differ from run_pipeline.sh's expected/ output (the "normal" now
# looks identical to the tumor); this script only exercises that the
# Nextflow pipeline itself runs end-to-end on the published image, not
# numeric parity with the direct-CLI run.
#
# Usage: run_nextflow_pipeline.sh OUTDIR [PROFILE]
#   PROFILE defaults to "docker"; pass "singularity" to use that instead.
#
# Requires nextflow, samtools, and bwa on PATH, plus a working docker (or
# singularity) able to pull yuhecheng62/dupcaller, biocontainers/bwa, and
# broadinstitute/gatk.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTDIR="${1:?usage: run_nextflow_pipeline.sh OUTDIR [PROFILE]}"
PROFILE="${2:-docker}"

NEXTFLOW="${NEXTFLOW:-nextflow}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BWA="${BWA:-bwa}"

mkdir -p "$OUTDIR"
cp "$DATA_DIR/reference.fa" "$DATA_DIR/repeats.tsv" "$OUTDIR/"
cp "$DATA_DIR/mock_1.fastq" "$OUTDIR/tumor_1.fastq"
cp "$DATA_DIR/mock_2.fastq" "$OUTDIR/tumor_2.fastq"
cp "$DATA_DIR/mock_1.fastq" "$OUTDIR/normal_1.fastq"
cp "$DATA_DIR/mock_2.fastq" "$OUTDIR/normal_2.fastq"
cd "$OUTDIR"

echo "[1/2] pre-index (bwa index + faidx -- DupCaller.nf expects these already present)"
"$SAMTOOLS" faidx reference.fa
"$BWA" index reference.fa

printf 'sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\nmock\ttumor_1.fastq\ttumor_2.fastq\tnormal_1.fastq\tnormal_2.fastq\n' > sample.map

echo "[2/2] nextflow run (profile: ${PROFILE})"
"$NEXTFLOW" run "${REPO_ROOT}/nextflow/DupCaller.nf" \
    -c "${REPO_ROOT}/nextflow/nextflow.config" \
    -c "${SCRIPT_DIR}/nextflow_test.config" \
    -profile "${PROFILE}",local \
    -w work

echo "done"
