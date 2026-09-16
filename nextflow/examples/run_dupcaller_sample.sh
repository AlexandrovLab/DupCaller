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
#       [-b BARCODE_PATTERN] \
#       [-B NORMAL_BAM] \
#       [-p THREADS] \
#       [-o OUTDIR] \
#       [-P PROFILE] \
#       [-R REPO_ROOT]
#
# Required:  -s -1 -2 -f, plus either (-3 and -4) or -B
# -B NORMAL_BAM: reuse an already-aligned, already-indexed (.bai) normal
#   BAM instead of aligning -3/-4 -- e.g. one matched normal shared across
#   many tumor-only benchmark/mock samples. Skips trim/align/markdup for
#   the normal entirely. -3/-4 are ignored (and not required) when set.
# Defaults:  -p $(nproc), -o ./results/SAMPLE_ID, -P docker,local,
#            -r (DupCaller.py's own default: chr1-22,chrX),
#            -b pipeline.config's default (NNNXXXX) -- CHECK THIS MATCHES
#               YOUR ACTUAL BARCODE SCHEME (N=barcode base, X=skipped
#               constant base) before trusting results; a mismatch here
#               corrupts barcode extraction and every read's sequence
#               offset, silently collapsing duplex family formation
#               without erroring. Verify by comparing a raw fastq read's
#               length to its aligned length in an existing BAM from the
#               same data (the difference is the true total pattern
#               length) -- a barcode tag's own length (e.g. DB:Z:xxx-yyy)
#               only reveals the N-count, not any trailing skipped bases,
#               and checking that alone is NOT sufficient,
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
OUTDIR=""
SAMPLE_ID=""

while getopts "s:1:2:3:4:f:m:g:r:b:B:p:o:P:R:h" opt; do
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

# nextflow's -resume history (.nextflow/cache, .nextflow.log) is keyed to
# the launch directory, not to -w's work dir -- a fresh mktemp dir here
# would give every invocation a clean launch dir with no prior session to
# resume from, silently turning -resume into a full rerun. Reuse a stable,
# per-sample launch dir under OUTDIR instead so a retry after a failure
# actually resumes.
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

# Only overrides pipeline.config's placeholders that must be per-run --
# every calling-parameter default already in pipeline.config (-maf, -gaf,
# -d, -tt, -tr, -mq) matches DupCaller.py's own CLI defaults, so it's
# inherited unchanged from there. barcode_pattern is only written here if
# -b was actually passed, so omitting -b still falls through to
# pipeline.config's own default rather than an empty string.
barcode_config_line=""
if [ -n "$BARCODE_PATTERN" ]; then
    barcode_config_line="    barcode_pattern = \"${BARCODE_PATTERN}\""
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
