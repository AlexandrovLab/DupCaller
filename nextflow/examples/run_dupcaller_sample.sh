#!/usr/bin/env bash
# run_dupcaller_sample.sh -- Run nextflow/DupCaller.nf end-to-end (fastq ->
# results) for one tumor/normal sample. Scheduler-agnostic: runs on a
# workstation or HPC node with Docker or Singularity. On a cluster, wrap it
# in a submission script and pass the granted CPU count via -p (see
# run_sample.slurm.sh).
#
# Requires: nextflow, and either Docker (running) or Singularity, on PATH.
# The reference must already have a samtools faidx index (.fai), a
# bwa-mem2 index (.0123/.amb/.ann/.bwt.2bit.64/.pac) and a DupCaller.py
# index (.ref.h5/.tn.h5/.hp.h5/.str.h5/.dbs.h5) alongside it.
#
# Usage:
#   run_dupcaller_sample.sh -s SAMPLE_ID \
#       -1 TUMOR_FASTQ_R1 -2 TUMOR_FASTQ_R2 \
#       -3 NORMAL_FASTQ_R1 -4 NORMAL_FASTQ_R2 \
#       -f REFERENCE_FASTA \
#       [-m MASK1.bed.gz[,MASK2.bed.gz,...]] \
#       [-g GERMLINE_VCF] \
#       [-r "chr1 chr2 ..."] \
#       [-b BARCODE_PATTERN] \
#       [-B NORMAL_BAM] \
#       [-p THREADS] \
#       [-M MAX_MEMORY] \
#       [-o OUTDIR] \
#       [-P PROFILE] \
#       [-R REPO_ROOT]
#
# Required:  -s -1 -2 -f, plus either (-3 and -4) or -B
# -B NORMAL_BAM: reuse an already-aligned, indexed (.bai) normal BAM
#   instead of aligning -3/-4 (skips trim/align/markdup for the normal).
#   -3/-4 are ignored and not required when set.
# -M MAX_MEMORY: cap every task's memory request (Nextflow memory string,
#   e.g. 24.GB); default pipeline.config's 96.GB. Lower it on smaller
#   machines -- the local executor refuses tasks that ask for more than
#   the machine has (bwa-mem2 mem and DupCaller call ask for 32 GB).
# Defaults:  -p $(nproc), -o ./results/SAMPLE_ID, -P docker,local,
#            -r (DupCaller.py's own default: chr1-22,chrX),
#            -b pipeline.config's default (NNNXXXX) -- CHECK THIS MATCHES
#               YOUR BARCODE SCHEME (N=barcode base, X=skipped base); a
#               wrong pattern silently breaks duplex family formation.
#               Compare a raw fastq read's length with its aligned length
#               in an existing BAM (the difference is the full pattern
#               length); the DB tag alone shows only the N count,
#            -R this script's own repo checkout
set -euo pipefail

usage() { grep '^#' "${BASH_SOURCE[0]}" | sed 's/^#//' | head -n 30; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROFILE="docker,local"
THREADS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
REGIONS=""
NOISE_MASKS=""
GERMLINE_VCF=""
BARCODE_PATTERN=""
NORMAL_BAM=""
MAX_MEMORY=""
OUTDIR=""
SAMPLE_ID=""

while getopts "s:1:2:3:4:f:m:g:r:b:B:p:M:o:P:R:h" opt; do
    case "$opt" in
        s) SAMPLE_ID="$OPTARG" ;;
        1) TUMOR_FASTQ_1="$OPTARG" ;;
        2) TUMOR_FASTQ_2="$OPTARG" ;;
        3) NORMAL_FASTQ_1="$OPTARG" ;;
        4) NORMAL_FASTQ_2="$OPTARG" ;;
        f) REFERENCE="$OPTARG" ;;
        m) NOISE_MASKS="$OPTARG" ;;
        g) GERMLINE_VCF="$OPTARG" ;;
        r) REGIONS="$OPTARG" ;;
        b) BARCODE_PATTERN="$OPTARG" ;;
        B) NORMAL_BAM="$OPTARG" ;;
        p) THREADS="$OPTARG" ;;
        M) MAX_MEMORY="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        P) PROFILE="$OPTARG" ;;
        R) REPO_ROOT="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

: "${SAMPLE_ID:?-s SAMPLE_ID is required}"
: "${TUMOR_FASTQ_1:?-1 TUMOR_FASTQ_R1 is required}"
: "${TUMOR_FASTQ_2:?-2 TUMOR_FASTQ_R2 is required}"
: "${REFERENCE:?-f REFERENCE_FASTA is required}"
if [ -z "$NORMAL_BAM" ]; then
    : "${NORMAL_FASTQ_1:?-3 NORMAL_FASTQ_R1 is required (or pass -B NORMAL_BAM)}"
    : "${NORMAL_FASTQ_2:?-4 NORMAL_FASTQ_R2 is required (or pass -B NORMAL_BAM)}"
fi

OUTDIR="${OUTDIR:-$(pwd)/results/${SAMPLE_ID}}"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"   # absolute, nextflow work/publishDir dislike relative paths across -resume

# DupCaller.nf always passes -r, so spell out DupCaller.py's default.
if [ -z "$REGIONS" ]; then
    REGIONS="$(printf 'chr%d ' $(seq 1 22))chrX"
fi

if [ -n "$NOISE_MASKS" ]; then
    IFS=',' read -ra _masks <<< "$NOISE_MASKS"
    mask_literal="["
    for m in "${_masks[@]}"; do mask_literal+="\"${m}\", "; done
    mask_literal="${mask_literal%, }]"
else
    mask_literal="null"
fi
if [ -n "$GERMLINE_VCF" ]; then
    germline_literal="\"${GERMLINE_VCF}\""
else
    germline_literal="null"
fi

# -resume history is keyed to the launch directory, so use a stable
# per-sample launch dir for retries to actually resume.
WORKDIR="${OUTDIR}/.nextflow_launch"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

if [ -n "$NORMAL_BAM" ]; then
    printf 'sample_id\ttumor_fastq_1\ttumor_fastq_2\n%s\t%s\t%s\n' \
        "$SAMPLE_ID" "$TUMOR_FASTQ_1" "$TUMOR_FASTQ_2" > sample.map
else
    printf 'sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\n%s\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE_ID" "$TUMOR_FASTQ_1" "$TUMOR_FASTQ_2" "$NORMAL_FASTQ_1" "$NORMAL_FASTQ_2" > sample.map
fi

# Override only per-run settings; calling parameters come from
# pipeline.config. Optional lines are written only when their flag is given,
# so pipeline.config's defaults otherwise apply.
barcode_config_line=""
if [ -n "$BARCODE_PATTERN" ]; then
    barcode_config_line="    barcode_pattern = \"${BARCODE_PATTERN}\""
fi
max_memory_config_line=""
if [ -n "$MAX_MEMORY" ]; then
    max_memory_config_line="    max_memory   = \"${MAX_MEMORY}\""
fi
normal_bam_config_line=""
if [ -n "$NORMAL_BAM" ]; then
    normal_bam_config_line="    normal_bam = \"${NORMAL_BAM}\""
fi
cat > run.config <<EOF
params {
    sample_map   = "sample.map"
    outdir       = "${OUTDIR}"
    reference    = "${REFERENCE}"
    skip_index   = true
    noise_mask   = ${mask_literal}
    germline_vcf = ${germline_literal}
    regions      = "${REGIONS}"
    threads      = ${THREADS}
    max_cpus     = ${THREADS}
${max_memory_config_line}
${barcode_config_line}
${normal_bam_config_line}
}
EOF

nextflow run "${REPO_ROOT}/nextflow/DupCaller.nf" \
    -c "${REPO_ROOT}/nextflow/nextflow.config" \
    -c "${REPO_ROOT}/nextflow/pipeline.config" \
    -c run.config \
    -profile "${PROFILE}" \
    -w "${OUTDIR}/work" \
    -resume

echo "Results: ${OUTDIR}/${SAMPLE_ID}"
