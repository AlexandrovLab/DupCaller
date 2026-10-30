#!/usr/bin/env bash
# run_dupcaller_sample.sh -- Run nextflow/DupCaller.nf end-to-end (fastq ->
# results) for one tumor/normal sample. Scheduler-agnostic: runs on a
# workstation or HPC node with Docker or Singularity. On a cluster, wrap it
# in a submission script and pass the granted CPU count via -p (see
# run_sample.slurm.sh).
#
# Requires: nextflow, and either Docker (running) or Singularity, on PATH.
# The reference (-f) must already have a samtools faidx index (.fai), a
# bwa-mem2 index (.0123/.amb/.ann/.bwt.2bit.64/.pac) and a DupCaller.py
# index (.ref.h5/.tn.h5/.hp.h5/.str.h5/.dbs.h5) alongside it. Or pass
# -G hg38|mm39 instead of -f to download a pre-made genome (see -G below).
#
# Usage:
#   run_dupcaller_sample.sh -s SAMPLE_ID \
#       -1 TUMOR_FASTQ_R1 -2 TUMOR_FASTQ_R2 \
#       -3 NORMAL_FASTQ_R1 -4 NORMAL_FASTQ_R2 \
#       -f REFERENCE_FASTA | -G hg38|mm39 [-S STRAIN] [-D PREMADE_DIR] \
#       [-m MASK1.bed.gz[,MASK2.bed.gz,...]] \
#       [-g GERMLINE_VCF] \
#       [-r "chr1 chr2 ..."] \
#       [-c PROTOCOL [-l BARCODE_LIST] | -b BARCODE_PATTERN] \
#       [-B NORMAL_BAM | -T] \
#       [-K] \
#       [-p THREADS] \
#       [-M MAX_MEMORY] \
#       [-o OUTDIR] \
#       [-P PROFILE] \
#       [-R REPO_ROOT]
#
# Required:  -s -1 -2, -f or -G, plus one of (-3 and -4), -B or -T
# -G GENOME: hg38 or mm39 (pipeline param premade). Downloads that genome's
#   FASTA, DupCaller index, germline VCF and noise mask(s), and builds the
#   bwa-mem2 index (~90 GB RAM for hg38, ~80 GB for mm39) into -D on the
#   first run; later runs reuse -D. -m/-g, if given, replace the premade
#   mask(s)/germline VCF. -S STRAIN (mm39 only, e.g. C57BL_6NJ) picks the
#   mouse germline VCF. -D defaults to ./dupcaller_premade.
# -B NORMAL_BAM: reuse an already-aligned, indexed (.bai) normal BAM
#   instead of aligning -3/-4 (skips trim/align/markdup for the normal).
#   -3/-4 are ignored and not required when set.
# -c PROTOCOL: renanoseq, nanoseqv2, udseq or metacs (pipeline param
#   protocol). Sets the barcode pattern and end trimming (trim_template/
#   trim_read) for that library protocol; see the main README's Protocol
#   Settings. metacs also needs -l BARCODE_LIST, the file of allowed Tn5
#   barcodes, one per line. Cannot be combined with -b.
# -T: tumor only, no matched normal (pipeline param tumor_only). Calls
#   with -maf 0.1 unless pipeline.config sets max_af below 1, so variants
#   above 10% allele fraction (mostly germline) are not called. Cannot be
#   combined with -3/-4 or -B.
# -K: keep the trimmed FASTQs in the work directory (pipeline param
#   keep_trimmed_fastq); by default they are cleared once aligned.
# -M MAX_MEMORY: cap every task's memory request (Nextflow memory string,
#   e.g. 24.GB); default pipeline.config's 96.GB. Lower it on smaller
#   machines -- the local executor refuses tasks that ask for more than
#   the machine has (bwa-mem2 mem asks for 32 GB, DupCaller call for
#   6 GB per CPU).
# Defaults:  -p $(nproc), -o ./results/SAMPLE_ID, -P docker,local,
#            -r (DupCaller.py's own default: chr1-22,chrX; chr1-19,chrX
#               with -G mm39),
#            -b pipeline.config's default (NNNXXXX), unless -c is given -- CHECK THIS MATCHES
#               YOUR BARCODE SCHEME (N=barcode base, X=skipped base); a
#               wrong pattern silently breaks duplex family formation.
#               Compare a raw fastq read's length with its aligned length
#               in an existing BAM (the difference is the full pattern
#               length); the DB tag alone shows only the N count,
#            -R this script's own repo checkout
set -euo pipefail

usage() { sed -n '2,/^set /s/^#//p' "${BASH_SOURCE[0]}"; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROFILE="docker,local"
THREADS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
REGIONS=""
NOISE_MASKS=""
GERMLINE_VCF=""
BARCODE_PATTERN=""
PROTOCOL=""
BARCODE_LIST=""
NORMAL_BAM=""
TUMOR_ONLY=""
KEEP_TRIMMED_FASTQ=""
MAX_MEMORY=""
OUTDIR=""
SAMPLE_ID=""
REFERENCE=""
PREMADE=""
PREMADE_STRAIN=""
PREMADE_DIR="$(pwd)/dupcaller_premade"

while getopts "s:1:2:3:4:f:G:S:D:m:g:r:b:c:l:B:TKp:M:o:P:R:h" opt; do
    case "$opt" in
        s) SAMPLE_ID="$OPTARG" ;;
        1) TUMOR_FASTQ_1="$OPTARG" ;;
        2) TUMOR_FASTQ_2="$OPTARG" ;;
        3) NORMAL_FASTQ_1="$OPTARG" ;;
        4) NORMAL_FASTQ_2="$OPTARG" ;;
        f) REFERENCE="$OPTARG" ;;
        G) PREMADE="$OPTARG" ;;
        S) PREMADE_STRAIN="$OPTARG" ;;
        D) PREMADE_DIR="$OPTARG" ;;
        m) NOISE_MASKS="$OPTARG" ;;
        g) GERMLINE_VCF="$OPTARG" ;;
        r) REGIONS="$OPTARG" ;;
        b) BARCODE_PATTERN="$OPTARG" ;;
        c) PROTOCOL="$OPTARG" ;;
        l) BARCODE_LIST="$OPTARG" ;;
        B) NORMAL_BAM="$OPTARG" ;;
        T) TUMOR_ONLY=1 ;;
        K) KEEP_TRIMMED_FASTQ=1 ;;
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
if [ -z "$REFERENCE" ] && [ -z "$PREMADE" ]; then
    echo "-f REFERENCE_FASTA or -G hg38|mm39 is required" >&2; exit 1
fi
if [ -n "$REFERENCE" ] && [ -n "$PREMADE" ]; then
    echo "-f and -G are mutually exclusive" >&2; exit 1
fi
if [ -n "$PROTOCOL" ] && [ -n "$BARCODE_PATTERN" ]; then
    echo "-c PROTOCOL and -b BARCODE_PATTERN are mutually exclusive (the protocol sets the pattern)" >&2; exit 1
fi
if [ -n "$BARCODE_LIST" ]; then
    [ -f "$BARCODE_LIST" ] || { echo "-l ${BARCODE_LIST}: no such file" >&2; exit 1; }
    BARCODE_LIST="$(cd "$(dirname "$BARCODE_LIST")" && pwd)/$(basename "$BARCODE_LIST")"
fi
if [ -n "$TUMOR_ONLY" ]; then
    if [ -n "$NORMAL_BAM" ] || [ -n "${NORMAL_FASTQ_1:-}" ] || [ -n "${NORMAL_FASTQ_2:-}" ]; then
        echo "-T (tumor only) cannot be combined with -3/-4 or -B" >&2; exit 1
    fi
elif [ -z "$NORMAL_BAM" ]; then
    : "${NORMAL_FASTQ_1:?-3 NORMAL_FASTQ_R1 is required (or pass -B NORMAL_BAM, or -T for tumor only)}"
    : "${NORMAL_FASTQ_2:?-4 NORMAL_FASTQ_R2 is required (or pass -B NORMAL_BAM, or -T for tumor only)}"
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

if [ -n "$NORMAL_BAM" ] || [ -n "$TUMOR_ONLY" ]; then
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
elif [ -n "$PROTOCOL" ]; then
    barcode_config_line="    protocol = \"${PROTOCOL}\""
fi
if [ -n "$BARCODE_LIST" ]; then
    barcode_config_line+="
    barcode_list = \"${BARCODE_LIST}\""
fi
max_memory_config_line=""
executor_memory_line=""
if [ -n "$MAX_MEMORY" ]; then
    max_memory_config_line="    max_memory   = \"${MAX_MEMORY}\""
    executor_memory_line="        memory = '${MAX_MEMORY}'"
fi
if [ -n "$PREMADE" ]; then
    mkdir -p "$PREMADE_DIR"
    PREMADE_DIR="$(cd "$PREMADE_DIR" && pwd)"
    reference_config_lines="    premade      = \"${PREMADE}\"
    premade_dir  = \"${PREMADE_DIR}\""
    if [ -n "$PREMADE_STRAIN" ]; then
        reference_config_lines+="
    premade_strain = \"${PREMADE_STRAIN}\""
    fi
else
    reference_config_lines="    reference    = \"${REFERENCE}\""
fi
normal_bam_config_line=""
if [ -n "$NORMAL_BAM" ]; then
    normal_bam_config_line="    normal_bam = \"${NORMAL_BAM}\""
elif [ -n "$TUMOR_ONLY" ]; then
    normal_bam_config_line="    tumor_only = true"
fi
keep_fastq_config_line=""
[ -n "$KEEP_TRIMMED_FASTQ" ] && keep_fastq_config_line="    keep_trimmed_fastq = true"
cat > run.config <<EOF
params {
    sample_map   = "sample.map"
    outdir       = "${OUTDIR}"
${reference_config_lines}
    skip_index   = true
    noise_mask   = ${mask_literal}
    germline_vcf = ${germline_literal}
    regions      = "${REGIONS}"
    threads      = ${THREADS}
    max_cpus     = ${THREADS}
${max_memory_config_line}
${barcode_config_line}
${normal_bam_config_line}
${keep_fastq_config_line}
}
// The local executor otherwise schedules against the whole machine's CPUs
// and memory, not what -p/-M (or a scheduler allocation) grant this run.
executor {
    \$local {
        cpus = ${THREADS}
${executor_memory_line}
    }
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
