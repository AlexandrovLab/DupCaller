#!/usr/bin/env bash
#SBATCH --job-name=dc_nf_SAMPLE_ID
#SBATCH --partition=platinum
#SBATCH --qos=hcp-ddp302
#SBATCH --account=ddp302
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=7-00:00:00
#SBATCH --output=SAMPLE_ID_nextflow.%j.out
#
# SLURM-specific wrapper around the platform-agnostic run_dupcaller_sample.sh
# -- this file only adds scheduler resource requests and this cluster's
# real paths/account; all the actual pipeline logic lives in the generic
# script. Adjust the #SBATCH lines above for your own scheduler/account/
# partition; everything below only needs SAMPLE_ID and the 4 fastq paths
# to change per sample.
#
# Usage:
#   sbatch --export=ALL,SAMPLE_ID=PD43276,\
#     TUMOR_FASTQ_1=...,TUMOR_FASTQ_2=...,\
#     NORMAL_FASTQ_1=...,NORMAL_FASTQ_2=... \
#     run_sample.slurm.sh
#
# The reference, mask pair, and per-flag defaults below were validated
# (2026-09-08) against a real production sample's actual call command
# (runtime_bechmarking/DC4.call.rerun/2pass_64core/PD*.2pass.sl) to
# confirm they resolve to the exact same DupCaller.py call invocation.

set -euo pipefail

SAMPLE_ID="${SAMPLE_ID:?set SAMPLE_ID}"
TUMOR_FASTQ_1="${TUMOR_FASTQ_1:?set TUMOR_FASTQ_1}"
TUMOR_FASTQ_2="${TUMOR_FASTQ_2:?set TUMOR_FASTQ_2}"
NORMAL_FASTQ_1="${NORMAL_FASTQ_1:?set NORMAL_FASTQ_1}"
NORMAL_FASTQ_2="${NORMAL_FASTQ_2:?set NORMAL_FASTQ_2}"

REPO_ROOT="${REPO_ROOT:-/tscc/nfs/home/yuc211/DupCaller_backup_20251215}"
OUTDIR="${OUTDIR:-/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/NanoSeq_Sanger_Analysis/runtime_bechmarking/nextflow_results/${SAMPLE_ID}}"

REFERENCE=/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/reference/hg38_gdc/GRCh38.d1.vd1.fa
SNP_MASK=/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/NanoSeq_Sanger_Analysis/runtime_bechmarking/snp_mask_dc_benchmark.bed.gz
NOISE_MASK=/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/NanoSeq_Sanger_Analysis/runtime_bechmarking/noise_mask_dc_benchmark.bed.gz

module load singularitypro/3.11
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity_cache}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

"${REPO_ROOT}/nextflow/examples/run_dupcaller_sample.sh" \
    -s "${SAMPLE_ID}" \
    -1 "${TUMOR_FASTQ_1}" -2 "${TUMOR_FASTQ_2}" \
    -3 "${NORMAL_FASTQ_1}" -4 "${NORMAL_FASTQ_2}" \
    -f "${REFERENCE}" \
    -m "${SNP_MASK},${NOISE_MASK}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -o "${OUTDIR}" \
    -P singularity,local \
    -R "${REPO_ROOT}"
