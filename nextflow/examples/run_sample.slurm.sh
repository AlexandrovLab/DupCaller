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
# SLURM-specific wrapper around the platform-agnostic run_dupcaller_sample.sh
# -- this file only adds scheduler resource requests and your cluster's
# real paths/account; all the actual pipeline logic lives in the generic
# script. Adjust the #SBATCH lines above for your own scheduler/account/
# partition; everything below only needs SAMPLE_ID and the 4 fastq paths
# to change per sample.
#
# Usage:
#   sbatch --export=ALL,SAMPLE_ID=SAMPLE1,\
#     TUMOR_FASTQ_1=...,TUMOR_FASTQ_2=...,\
#     NORMAL_FASTQ_1=...,NORMAL_FASTQ_2=... \
#     run_sample.slurm.sh
#
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
NORMAL_FASTQ_1="${NORMAL_FASTQ_1:?set NORMAL_FASTQ_1}"
NORMAL_FASTQ_2="${NORMAL_FASTQ_2:?set NORMAL_FASTQ_2}"

REPO_ROOT="${REPO_ROOT:-/path/to/DupCaller}"
OUTDIR="${OUTDIR:-/path/to/results/${SAMPLE_ID}}"

REFERENCE=/path/to/dupcaller_index/reference.fa
SNP_MASK=/path/to/snp_mask.bed.gz
NOISE_MASK=/path/to/noise_mask.bed.gz
BARCODE_PATTERN="${BARCODE_PATTERN:-NNNXXXX}"

module load singularity  # or your cluster's equivalent module name
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity_cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

"${REPO_ROOT}/nextflow/examples/run_dupcaller_sample.sh" \
    -s "${SAMPLE_ID}" \
    -1 "${TUMOR_FASTQ_1}" -2 "${TUMOR_FASTQ_2}" \
    -3 "${NORMAL_FASTQ_1}" -4 "${NORMAL_FASTQ_2}" \
    -f "${REFERENCE}" \
    -m "${SNP_MASK},${NOISE_MASK}" \
    -b "${BARCODE_PATTERN}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -o "${OUTDIR}" \
    -P singularity,local \
    -R "${REPO_ROOT}"
