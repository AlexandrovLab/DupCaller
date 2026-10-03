#!/usr/bin/env bash
# Runs nextflow/DupCaller.nf (DupCaller index + bwa-mem2 index -> trim ->
# bwa-mem2 mem -> gatk MarkDuplicates -> call -> estimate, containerized)
# against the synthetic dataset in data/, entirely inside OUTDIR. Companion
# to run_pipeline.sh, which runs the same tools directly on PATH instead of
# through Nextflow and containers -- this one instead exercises the
# Nextflow DAG itself (staging, container images, per-process resource
# directives) using the published yuhecheng62/dupcaller image.
#
# The pipeline's CALL_VARIANTS step always wants a matched normal (no
# tumor-only mode), unlike run_pipeline.sh's data/ usage -- so mock_1/2 are
# staged as both tumor and normal fastqs here. That means SBS/indel calls
# will differ from run_pipeline.sh's expected/ output (the "normal" now
# looks identical to the tumor); this script only exercises that the
# Nextflow pipeline itself runs end-to-end on the published image, not
# numeric parity with the direct-CLI run.
#
# Usage: run_nextflow_pipeline.sh OUTDIR [PROFILE] [CASE ...]
#   PROFILE defaults to "docker"; pass "singularity" to use that instead.
#   CASE (default: "default") selects which optional inputs are set; each
#   case publishes to OUTDIR/results/CASE/mock and reuses earlier cases'
#   cached tasks (-resume, shared work dir):
#     default      every optional resource unset (germline_vcf, noise_mask,
#                  target_bed, indel_bed, gene_bed all null)
#     germline     germline_vcf
#     noise_mask   noise_mask (two masks)
#     target_bed   target_bed
#     indel_bed    indel_bed
#     gene_bed     gene_bed
#     dilute       estimate_dilute = true
#     normal_bam   normal_bam (the default case's own normal BAM; run
#                  after "default")
#     all          every case above, in order
#
# Requires nextflow and python (with pysam, a DupCaller dependency) on
# PATH, plus a working docker (or singularity) able to pull
# yuhecheng62/dupcaller, quay.io/biocontainers/bwa-mem2, and
# broadinstitute/gatk.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTDIR="${1:?usage: run_nextflow_pipeline.sh OUTDIR [PROFILE] [CASE ...]}"
PROFILE="${2:-docker}"
shift $(( $# >= 2 ? 2 : $# ))
CASES=("$@")
[ ${#CASES[@]} -eq 0 ] && CASES=(default)
if [ "${CASES[0]}" = "all" ]; then
    CASES=(default germline noise_mask target_bed indel_bed gene_bed dilute normal_bam)
fi

NEXTFLOW="${NEXTFLOW:-nextflow}"
PYTHON="${PYTHON:-python}"

mkdir -p "$OUTDIR"
cp "$DATA_DIR/reference.fa" "$DATA_DIR/repeats.tsv" "$OUTDIR/"
cp "$DATA_DIR/mock_1.fastq" "$OUTDIR/tumor_1.fastq"
cp "$DATA_DIR/mock_2.fastq" "$OUTDIR/tumor_2.fastq"
cp "$DATA_DIR/mock_1.fastq" "$OUTDIR/normal_1.fastq"
cp "$DATA_DIR/mock_2.fastq" "$OUTDIR/normal_2.fastq"
cd "$OUTDIR"
OUTDIR="$(pwd)"

echo "[setup] faidx + optional test resources (bwa-mem2/DupCaller indexes are built by the pipeline)"
"$PYTHON" -c "import pysam; pysam.faidx('reference.fa')"
"$PYTHON" "${SCRIPT_DIR}/make_test_resources.py" resources reference.fa

printf 'sample_id\ttumor_fastq_1\ttumor_fastq_2\tnormal_fastq_1\tnormal_fastq_2\nmock\ttumor_1.fastq\ttumor_2.fastq\tnormal_1.fastq\tnormal_2.fastq\n' > sample.map

R="${OUTDIR}/resources"
for case in "${CASES[@]}"; do
    case "$case" in
        default)    extra='' ;;
        germline)   extra="germline_vcf = \"${R}/germline.vcf.gz\"" ;;
        noise_mask) extra="noise_mask = [\"${R}/snp_mask.bed.gz\", \"${R}/noise_mask.bed.gz\"]" ;;
        target_bed) extra="target_bed = \"${R}/target.bed.gz\"" ;;
        indel_bed)  extra="indel_bed = \"${R}/indel_pon.bed.gz\"" ;;
        gene_bed)   extra="gene_bed = \"${R}/genes.bed.gz\"" ;;
        dilute)     extra='estimate_dilute = true' ;;
        normal_bam)
            # The default case's normal BAM (aligned from the synthetic
            # normal FASTQs), copied out of work/ so the file is new to
            # Nextflow's cache and CALL_VARIANTS really reruns with it.
            nbam="$(find "${OUTDIR}/work" -name 'mock_normal.mkdped.bam' | head -1)"
            [ -n "$nbam" ] || { echo "normal_bam case needs the default case's normal BAM -- run 'default' first" >&2; exit 1; }
            mkdir -p "${R}/normal_bam"
            cp "$nbam" "${R}/normal_bam/mock_normal.mkdped.bam"
            cp "${nbam}.bai" "${R}/normal_bam/mock_normal.mkdped.bam.bai"
            extra="normal_bam = \"${R}/normal_bam/mock_normal.mkdped.bam\"" ;;
        *) echo "unknown case: $case" >&2; exit 1 ;;
    esac
    printf 'params {\n    outdir = "%s/results/%s"\n    %s\n}\n' "$OUTDIR" "$case" "$extra" > "case_${case}.config"

    echo "[run] case ${case} (profile: ${PROFILE})"
    "$NEXTFLOW" run "${REPO_ROOT}/nextflow/DupCaller.nf" \
        -c "${REPO_ROOT}/nextflow/nextflow.config" \
        -c "${SCRIPT_DIR}/nextflow_test.config" \
        -c "case_${case}.config" \
        -profile "${PROFILE}",local \
        -w work \
        -resume
done

echo "done"
