#!/usr/bin/env bash
# run_dupcaller_sample.sh -- Run nextflow/DupCaller.nf end-to-end (fastq ->
# results) for one tumor/normal sample. Platform-agnostic: makes no
# assumptions about a scheduler. Works standing alone on a laptop/
# workstation with Docker installed, or on any HPC node with Docker or
# Singularity available. For cluster resource requests (SLURM/SGE/PBS
# etc.), wrap this script in your scheduler's own submission script and
# pass through whatever CPU count it grants you via -p -- see
# run_sample.slurm.sh for a SLURM example.
#
# Requires: nextflow, and either Docker (running) or Singularity, on PATH.
# The reference must already have a bwa index (.bwt/.pac/.amb/.ann/.sa)
# and, unless you pass -I, a DupCaller.py index (.ref.h5/.tn.h5/.hp.h5/
# .str.h5/.dbs.h5) alongside it.
#
# Usage:
#   run_dupcaller_sample.sh -s SAMPLE_ID \
#       -1 TUMOR_FASTQ_R1 -2 TUMOR_FASTQ_R2 \
#       -3 NORMAL_FASTQ_R1 -4 NORMAL_FASTQ_R2 \
#       -f REFERENCE_FASTA \
#       [-m MASK1.bed.gz[,MASK2.bed.gz,...]] \
#       [-g GERMLINE_VCF] \
#       [-r "chr1 chr2 ..."] \
#       [-p THREADS] \
#       [-o OUTDIR] \
#       [-P PROFILE] \
#       [-R REPO_ROOT]
#
# Required:  -s -1 -2 -3 -4 -f
# Defaults:  -p $(nproc), -o ./results/SAMPLE_ID, -P docker,local,
#            -r (DupCaller.py's own default: chr1-22,chrX),
#            -R this script's own repo checkout
set -euo pipefail

usage() { grep '^#' "${BASH_SOURCE[0]}" | sed 's/^#//' | head -n 30; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROFILE="docker,local"
THREADS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
REGIONS=""
NOISE_MASKS=""
GERMLINE_VCF=""
OUTDIR=""
SAMPLE_ID=""

while getopts "s:1:2:3:4:f:m:g:r:p:o:P:R:h" opt; do
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
        p) THREADS="$OPTARG" ;;
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
: "${NORMAL_FASTQ_1:?-3 NORMAL_FASTQ_R1 is required}"
: "${NORMAL_FASTQ_2:?-4 NORMAL_FASTQ_R2 is required}"
: "${REFERENCE:?-f REFERENCE_FASTA is required}"

OUTDIR="${OUTDIR:-$(pwd)/results/${SAMPLE_ID}}"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"   # absolute, nextflow work/publishDir dislike relative paths across -resume

# DupCaller.py call/estimate's own default region set when -r is omitted
# entirely -- DupCaller.nf always passes -r explicitly, so spell it out
# unless the caller overrides it.
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

WORKDIR="$(mktemp -d)"
cd "$WORKDIR"

printf 'sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\n%s\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE_ID" "$TUMOR_FASTQ_1" "$TUMOR_FASTQ_2" "$NORMAL_FASTQ_1" "$NORMAL_FASTQ_2" > sample.map

# Only overrides pipeline.config's placeholders that must be per-run --
# every calling-parameter default already in pipeline.config (-maf, -gaf,
# -d, -tt, -tr, -mq) matches DupCaller.py's own CLI defaults, so it's
# inherited unchanged from there.
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
