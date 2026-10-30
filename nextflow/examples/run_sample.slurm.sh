#!/usr/bin/env bash
#SBATCH --job-name=dc_nf_SAMPLE_ID
#SBATCH --partition=YOUR_PARTITION
#SBATCH --qos=YOUR_QOS
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=7-00:00:00
#SBATCH --output=SAMPLE_ID_nextflow.%j.out
#
# SLURM wrapper around run_dupcaller_sample.sh: adds resource requests and
# site paths only. Per sample, only SAMPLE_ID and the 4 fastq paths change.
#
# EDIT BEFORE THE FIRST sbatch (the placeholders fail as shipped):
#   - the #SBATCH --partition/--qos/--account lines above (YOUR_*), and
#     --cpus-per-task/--mem/--time for your cluster;
#   - REPO_ROOT and OUTDIR below (or pass them in --export);
#   - REFERENCE, SNP_MASK and NOISE_MASK below, or instead pass PREMADE=hg38
#     (or mm39) with PREMADE_DIR=/shared/dupcaller_premade in --export to
#     download a pre-made genome (run_dupcaller_sample.sh -G/-D; for mm39,
#     PREMADE_STRAIN=C57BL_6NJ etc. picks the germline VCF, -S).
#
# Usage:
#   sbatch --export=ALL,SAMPLE_ID=SAMPLE1,\
#     TUMOR_FASTQ_1=...,TUMOR_FASTQ_2=...,\
#     NORMAL_FASTQ_1=...,NORMAL_FASTQ_2=... \
#     run_sample.slurm.sh
#
# Tumor only (no matched normal): pass TUMOR_ONLY=1 instead of
# NORMAL_FASTQ_1/2 (run_dupcaller_sample.sh -T).
# KEEP_TRIMMED_FASTQ=1 keeps the trimmed FASTQs in the work directory
# (run_dupcaller_sample.sh -K).
# MAX_MEMORY caps every task's memory (run_dupcaller_sample.sh -M, a
# Nextflow memory string, e.g. MAX_MEMORY=192.GB). Default: this job's
# Slurm memory allocation (--mem), so DupCaller.py call gets up to 6 GB per
# CPU instead of pipeline.config's 96.GB cap.
#
# PROTOCOL=renanoseq|nanoseqv2|udseq|metacs sets the barcode pattern and end
# trimming for that protocol (metacs also needs BARCODE_LIST=file of allowed
# barcodes); without it, BARCODE_PATTERN is used.
# BARCODE_PATTERN defaults to NNNXXXX (pipeline.config's default). Check it
# against your own data by comparing raw-FASTQ read length with the trimmed
# BAM's: the difference is the full pattern length (e.g. 151bp -> 144bp is 7
# bases, len("NNNXXXX")). A DB tag (e.g. DB:Z:TCC-CTG) only shows the number
# of N (barcode) bases, not how many trailing X (skipped) bases were also
# removed, so it alone cannot determine the pattern.

set -euo pipefail

SAMPLE_ID="${SAMPLE_ID:?set SAMPLE_ID}"
TUMOR_FASTQ_1="${TUMOR_FASTQ_1:?set TUMOR_FASTQ_1}"
TUMOR_FASTQ_2="${TUMOR_FASTQ_2:?set TUMOR_FASTQ_2}"
if [ -n "${TUMOR_ONLY:-}" ]; then
    NORMAL_ARGS=(-T)
else
    NORMAL_ARGS=(-3 "${NORMAL_FASTQ_1:?set NORMAL_FASTQ_1 (or TUMOR_ONLY=1)}"
                 -4 "${NORMAL_FASTQ_2:?set NORMAL_FASTQ_2 (or TUMOR_ONLY=1)}")
fi

REPO_ROOT="${REPO_ROOT:-/path/to/DupCaller}"
OUTDIR="${OUTDIR:-/path/to/results/${SAMPLE_ID}}"

# PREMADE=hg38|mm39 replaces REFERENCE and the masks with a pre-made genome
# kept in PREMADE_DIR (needs --mem of at least ~90 GB for the first run's
# bwa-mem2 index).
if [ -n "${PREMADE:-}" ]; then
    REFERENCE_ARGS=(-G "${PREMADE}" -D "${PREMADE_DIR:?set PREMADE_DIR (shared, ~50 GB for hg38) with PREMADE}")
    [ -n "${PREMADE_STRAIN:-}" ] && REFERENCE_ARGS+=(-S "${PREMADE_STRAIN}")
else
    REFERENCE=/path/to/dupcaller_index/reference.fa
    SNP_MASK=/path/to/snp_mask.bed.gz
    NOISE_MASK=/path/to/noise_mask.bed.gz
    REFERENCE_ARGS=(-f "${REFERENCE}" -m "${SNP_MASK},${NOISE_MASK}")
fi
# SLURM_MEM_PER_NODE / SLURM_MEM_PER_CPU are in MB.
if [ -z "${MAX_MEMORY:-}" ]; then
    mem_mb="${SLURM_MEM_PER_NODE:-0}"
    [ "$mem_mb" -gt 0 ] || mem_mb=$(( ${SLURM_MEM_PER_CPU:-0} * ${SLURM_CPUS_PER_TASK:-1} ))
    if [ $(( mem_mb / 1024 )) -gt 0 ]; then MAX_MEMORY="$(( mem_mb / 1024 )).GB"; fi
fi
MEMORY_ARGS=()
[ -n "${MAX_MEMORY:-}" ] && MEMORY_ARGS=(-M "${MAX_MEMORY}")
if [ -n "${PROTOCOL:-}" ]; then
    BARCODE_ARGS=(-c "${PROTOCOL}")
    [ -n "${BARCODE_LIST:-}" ] && BARCODE_ARGS+=(-l "${BARCODE_LIST}")
else
    BARCODE_ARGS=(-b "${BARCODE_PATTERN:-NNNXXXX}")
fi
KEEP_ARGS=()
[ -n "${KEEP_TRIMMED_FASTQ:-}" ] && KEEP_ARGS=(-K)

module load singularity  # or your cluster's equivalent module name
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity_cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

"${REPO_ROOT}/nextflow/examples/run_dupcaller_sample.sh" \
    -s "${SAMPLE_ID}" \
    -1 "${TUMOR_FASTQ_1}" -2 "${TUMOR_FASTQ_2}" \
    "${NORMAL_ARGS[@]}" \
    "${REFERENCE_ARGS[@]}" \
    "${BARCODE_ARGS[@]}" \
    ${KEEP_ARGS[@]+"${KEEP_ARGS[@]}"} \
    ${MEMORY_ARGS[@]+"${MEMORY_ARGS[@]}"} \
    -p "${SLURM_CPUS_PER_TASK}" \
    -o "${OUTDIR}" \
    -P singularity,local \
    -R "${REPO_ROOT}"
