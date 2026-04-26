#!/usr/bin/env nextflow

/*
 * DupCaller Nextflow Pipeline
 *
 * A pipeline for calling somatic mutations from barcoded error-corrected NGS data
 */

nextflow.enable.dsl=2

/*
 * Pipeline parameters with default values
 */
params.help = false
params.outdir = "./results"

// Input files
params.sample_name = null
params.read1 = null
params.read2 = null
params.normal_read1 = null
params.normal_read2 = null
params.normal_bam = null
params.tumor_bam = null         // pre-aligned tumor BAM (when skipping trim/align)

// Reference files
params.reference = null
params.germline_vcf = null
params.noise_mask = null
params.indel_bed = null

// Barcode trimming parameters
params.barcode_pattern = "NNNXXXX"

// Alignment parameters
params.threads = 4
params.read_group_id = null
params.platform = "ILLUMINA"

// DupCaller calling parameters
params.regions = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
params.max_af = 1.0
params.germline_af_cutoff = 0.001
params.min_n_depth = 10
params.max_zero_qual_frac = 0.5
params.trim_template = 7
params.trim_read = 7
params.mapq = 40
params.window_size = 100000

// Burden estimation parameters
params.estimate_clonal = false
params.estimate_dilute = false
params.gene_bed = null

// Skip options
params.skip_trim = false
params.skip_align = false
params.skip_normal = false

def helpMessage() {
    log.info"""
    ================================================================
    DupCaller Pipeline
    ================================================================

    Usage:
    nextflow run DupCaller.nf --sample_name SAMPLE --read1 R1.fq.gz --read2 R2.fq.gz --reference ref.fa [options]

    Required arguments:
      --sample_name         Sample name for output files
      --reference           Reference genome fasta file (must be indexed with DupCaller.py index)

    Input options (choose one):
      --read1               Read 1 fastq file (start from FASTQ)
      --read2               Read 2 fastq file (start from FASTQ)
      --tumor_bam           Pre-aligned/markdup tumor BAM file (skip trim and align)

    Normal sample (choose one, or use --skip_normal):
      --normal_read1        Normal sample read 1 fastq
      --normal_read2        Normal sample read 2 fastq
      --normal_bam          Pre-aligned normal BAM file (must have .bai index)
      --skip_normal         Run without matched normal (set --max_af appropriately)

    Reference resources:
      --germline_vcf        Indexed germline VCF with AF field
      --noise_mask          BED file for noise masking
      --indel_bed           Enhanced panel of normal for indels

    Barcode trimming:
      --barcode_pattern     Barcode pattern (default: NNNXXXX for NanoSeq)

    Alignment:
      --threads             Number of threads (default: 4)
      --platform            Sequencing platform (default: ILLUMINA)

    Variant calling:
      --regions             Contigs for variant calling (default: chr1-22,X,Y)
      --max_af              Maximum allele frequency (default: 1.0)
      --germline_af_cutoff  Germline AF cutoff (default: 0.001)
      --min_n_depth         Minimum normal depth (default: 10)
      --max_zero_qual_frac  Maximum zero quality fraction (default: 0.5)
      --trim_template       Template end trim distance (default: 7)
      --trim_read           Read end trim distance (default: 7)
      --mapq                Minimum mapping quality (default: 40)
      --window_size         Genomic window size (default: 100000)

    Burden estimation:
      --estimate_clonal     Consider clonal mutations (default: false)
      --estimate_dilute     Dilute mode for same starting material (default: false)
      --gene_bed            Gene BED file for coverage calculation

    Output:
      --outdir              Output directory (default: ./results)

    Example:
      nextflow run DupCaller.nf \\
        --sample_name sample1 \\
        --read1 sample1_R1.fq.gz \\
        --read2 sample1_R2.fq.gz \\
        --normal_bam normal.bam \\
        --reference hg38.fa \\
        --germline_vcf gnomad.vcf.gz \\
        --noise_mask noise.bed.gz \\
        --threads 8

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.sample_name) {
    error "Error: --sample_name is required"
}
if (!params.reference) {
    error "Error: --reference is required"
}
if (!params.read1 && !params.read2 && !params.tumor_bam) {
    error "Error: provide --read1/--read2 (start from FASTQ) or --tumor_bam (pre-aligned BAM)"
}

/*
 * Process: Trim barcodes from fastq files
 */
process TRIM_BARCODES {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/trimmed", mode: 'copy'

    input:
    tuple val(sample_name), path(read1), path(read2)

    output:
    tuple val(sample_name), path("${sample_name}_1.fastq"), path("${sample_name}_2.fastq")

    stub:
    """
    touch ${sample_name}_1.fastq ${sample_name}_2.fastq
    """

    script:
    """
    DupCaller.py trim \\
        -i ${read1} \\
        -i2 ${read2} \\
        -p ${params.barcode_pattern} \\
        -o ${sample_name}
    """
}

/*
 * Process: Align reads with BWA-MEM
 */
process BWA_ALIGN {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/aligned", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sample_name), path(read1), path(read2), path(reference)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")

    stub:
    """
    touch ${sample_name}.bam ${sample_name}.bam.bai
    """

    script:
    def rg_id = params.read_group_id ?: sample_name
    """
    bwa mem -C \\
        -t ${params.threads} \\
        -R "@RG\\tID:${rg_id}\\tSM:${sample_name}\\tPL:${params.platform}" \\
        ${reference} \\
        ${read1} ${read2} | \\
    samtools sort -@ ${params.threads} -o ${sample_name}.bam

    samtools index -@ ${params.threads} ${sample_name}.bam
    """
}

/*
 * Process: Mark duplicates with GATK
 */
process MARK_DUPLICATES {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/markdup", mode: 'copy'

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("${sample_name}.mkdped.bam"), path("${sample_name}.mkdped.bam.bai"), path("${sample_name}.mkdp_metrics.txt")

    stub:
    """
    touch ${sample_name}.mkdped.bam ${sample_name}.mkdped.bam.bai ${sample_name}.mkdp_metrics.txt
    """

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_name}.mkdped.bam \\
        -M ${sample_name}.mkdp_metrics.txt \\
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \\
        --DUPLEX_UMI \\
        --TAGGING_POLICY OpticalOnly \\
        --BARCODE_TAG DB

    samtools index ${sample_name}.mkdped.bam
    """
}

/*
 * Process: Index a pre-existing BAM file
 */
process INDEX_BAM {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/bam", mode: 'copy'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${bam}"), path("${bam}.bai"), path("${sample_name}.metrics_placeholder.txt")

    stub:
    """
    touch ${bam}.bai ${sample_name}.metrics_placeholder.txt
    """

    script:
    """
    samtools index ${bam}
    touch ${sample_name}.metrics_placeholder.txt
    """
}

/*
 * Process: Call variants with DupCaller
 */
process CALL_VARIANTS {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/variants", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sample_name), path(bam), path(bai), path(metrics),
          path(reference), path(ref_h5), path(tn_h5), path(hp_h5)
    tuple path(normal_bam), path(normal_bai)
    path(germline_vcf)
    path(noise_mask)
    path(indel_bed)

    output:
    tuple val(sample_name),
          path("${sample_name}_snv.vcf"),
          path("${sample_name}_indel.vcf"),
          path("${sample_name}_coverage.bed.gz"),
          path("${sample_name}_coverage.bed.gz.tbi"),
          path("${sample_name}_trinuc_by_duplex_group.txt"),
          path("${sample_name}_duplex_group_stats.txt"),
          path("${sample_name}_stats.txt"),
          path("${sample_name}.amp.tn.txt"),
          path("${sample_name}.amp.id.txt"),
          path("${sample_name}.dmg.tn.txt"),
          path("${sample_name}.dmg.id.txt"),
          path("${sample_name}_call_params.log")

    stub:
    """
    touch ${sample_name}_snv.vcf ${sample_name}_indel.vcf
    touch ${sample_name}_coverage.bed.gz ${sample_name}_coverage.bed.gz.tbi
    touch ${sample_name}_trinuc_by_duplex_group.txt ${sample_name}_duplex_group_stats.txt ${sample_name}_stats.txt
    touch ${sample_name}.amp.tn.txt ${sample_name}.amp.id.txt
    touch ${sample_name}.dmg.tn.txt ${sample_name}.dmg.id.txt
    touch ${sample_name}_call_params.log
    """

    script:
    def normal_arg = normal_bam.name != 'NO_NORMAL_BAM' ? "-n ${normal_bam}" : ""
    def germline_arg = germline_vcf.name != 'NO_GERMLINE_VCF' ? "-g ${germline_vcf}" : ""
    def noise_arg = noise_mask.name != 'NO_NOISE_MASK' ? "-m ${noise_mask}" : ""
    def indel_arg = indel_bed.name != 'NO_INDEL_BED' ? "-id ${indel_bed}" : ""

    """
    DupCaller.py call \\
        -b ${bam} \\
        -f ${reference} \\
        -o ${sample_name} \\
        -p ${params.threads} \\
        -r ${params.regions} \\
        ${normal_arg} \\
        ${germline_arg} \\
        ${noise_arg} \\
        ${indel_arg} \\
        -maf ${params.max_af} \\
        -gaf ${params.germline_af_cutoff} \\
        -d ${params.min_n_depth} \\
        -z ${params.max_zero_qual_frac} \\
        -tt ${params.trim_template} \\
        -tr ${params.trim_read} \\
        -mq ${params.mapq} \\
        -w ${params.window_size}
    mv ${sample_name}/* .
    """
}

/*
 * Process: Estimate mutation burden
 */
process ESTIMATE_BURDEN {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/burden", mode: 'copy'

    input:
    tuple val(sample_name),
          path(snv_vcf),
          path(indel_vcf),
          path(coverage),
          path(coverage_tbi),
          path(trinuc),
          path(duplex_stats),
          path(stats),
          path(amp_tn),
          path(amp_id),
          path(dmg_tn),
          path(dmg_id),
          path(call_params_log),
          path(reference),
          path(ref_h5),
          path(tn_h5),
          path(hp_h5)
    path(gene_bed)

    output:
    tuple val(sample_name),
          path("${sample_name}_sbs_burden.txt"),
          path("${sample_name}_indel_burden.txt"),
          path("${sample_name}_sbs_96_corrected.txt"),
          path("${sample_name}_sbs_96_corrected.png"),
          path("${sample_name}_sbs_burden_by_min_read_group_size.txt"),
          path("${sample_name}_sbs_burden_by_min_read_group_size.png"),
          path("${sample_name}_duplex_allele_counts.txt")

    stub:
    """
    touch ${sample_name}_sbs_burden.txt ${sample_name}_indel_burden.txt
    touch ${sample_name}_sbs_96_corrected.txt ${sample_name}_sbs_96_corrected.png
    touch ${sample_name}_sbs_burden_by_min_read_group_size.txt ${sample_name}_sbs_burden_by_min_read_group_size.png
    touch ${sample_name}_duplex_allele_counts.txt
    """

    script:
    def gene_arg = gene_bed.name != 'NO_GENE_BED' ? "-gb ${gene_bed}" : ""
    def clonal_arg = params.estimate_clonal ? "-c true" : ""
    def dilute_arg = params.estimate_dilute ? "-d true" : ""

    """
    mkdir -p ${sample_name}
    for f in \$(ls ${sample_name}_* ${sample_name}.* 2>/dev/null); do
        ln -sf \$(realpath \$f) ${sample_name}/
    done
    DupCaller.py estimate \\
        -i ${sample_name} \\
        -f ${reference} \\
        -r ${params.regions} \\
        ${gene_arg} \\
        ${clonal_arg} \\
        ${dilute_arg}
    mv ${sample_name}/${sample_name}_sbs_burden.txt \\
       ${sample_name}/${sample_name}_indel_burden.txt \\
       ${sample_name}/${sample_name}_sbs_96_corrected.txt \\
       ${sample_name}/${sample_name}_sbs_96_corrected.png \\
       ${sample_name}/${sample_name}_sbs_burden_by_min_read_group_size.txt \\
       ${sample_name}/${sample_name}_sbs_burden_by_min_read_group_size.png \\
       ${sample_name}/${sample_name}_duplex_allele_counts.txt \\
       .
    """
}

/*
 * Main workflow
 */
workflow {

    // Reference file channels — all four files must exist (created by DupCaller.py index)
    ref_ch    = Channel.fromPath(params.reference,              checkIfExists: true)
    ref_h5_ch = Channel.fromPath(params.reference + ".ref.h5", checkIfExists: true)
    tn_h5_ch  = Channel.fromPath(params.reference + ".tn.h5",  checkIfExists: true)
    hp_h5_ch  = Channel.fromPath(params.reference + ".hp.h5",  checkIfExists: true)

    // Optional input channels
    germline_ch = params.germline_vcf
        ? Channel.fromPath(params.germline_vcf, checkIfExists: true)
        : Channel.value(file('NO_GERMLINE_VCF'))
    noise_ch = params.noise_mask
        ? Channel.fromPath(params.noise_mask, checkIfExists: true)
        : Channel.value(file('NO_NOISE_MASK'))
    indel_ch = params.indel_bed
        ? Channel.fromPath(params.indel_bed, checkIfExists: true)
        : Channel.value(file('NO_INDEL_BED'))
    gene_ch = params.gene_bed
        ? Channel.fromPath(params.gene_bed, checkIfExists: true)
        : Channel.value(file('NO_GENE_BED'))

    // -----------------------------------------------------------------------
    // Tumor sample: produce markdup channel
    //   emits: tuple val(sample_name), path(bam), path(bai), path(metrics)
    // -----------------------------------------------------------------------
    if (params.tumor_bam) {
        // Pre-aligned BAM provided — index it and create a placeholder metrics file
        bam_ch = Channel.of([params.sample_name, file(params.tumor_bam)])
        markdup = INDEX_BAM(bam_ch)
    } else {
        // Start from FASTQ
        reads_ch = Channel.of([params.sample_name, file(params.read1), file(params.read2)])
        trimmed  = TRIM_BARCODES(reads_ch)
        aligned  = BWA_ALIGN(trimmed.combine(ref_ch))
        markdup  = MARK_DUPLICATES(aligned)
    }

    // -----------------------------------------------------------------------
    // Normal sample channel
    //   emits: tuple path(normal_bam), path(normal_bai)
    // -----------------------------------------------------------------------
    if (params.normal_bam) {
        normal_ch = Channel.value([
            file(params.normal_bam),
            file(params.normal_bam + ".bai")
        ])
    } else if (params.normal_read1 && params.normal_read2 && !params.skip_normal) {
        normal_reads   = Channel.of(["normal_${params.sample_name}", file(params.normal_read1), file(params.normal_read2)])
        normal_trimmed = TRIM_BARCODES(normal_reads)
        normal_aligned = BWA_ALIGN(normal_trimmed.combine(ref_ch))
        normal_markdup = MARK_DUPLICATES(normal_aligned)
        normal_ch      = normal_markdup.map { _sn, bam, bai, _m -> [bam, bai] }
    } else {
        normal_ch = Channel.value([file('NO_NORMAL_BAM'), file('NO_NORMAL_BAI')])
    }

    // -----------------------------------------------------------------------
    // Variant calling
    // call_input: tuple(sample_name, bam, bai, metrics, ref, ref_h5, tn_h5, hp_h5)
    // -----------------------------------------------------------------------
    call_input = markdup
        .combine(ref_ch)
        .combine(ref_h5_ch)
        .combine(tn_h5_ch)
        .combine(hp_h5_ch)

    variants = CALL_VARIANTS(
        call_input,
        normal_ch,
        germline_ch,
        noise_ch,
        indel_ch
    )

    // -----------------------------------------------------------------------
    // Burden estimation
    // burden_input: variants tuple + ref + ref_h5 + tn_h5 + hp_h5
    // -----------------------------------------------------------------------
    burden_input = variants
        .combine(ref_ch)
        .combine(ref_h5_ch)
        .combine(tn_h5_ch)
        .combine(hp_h5_ch)

    ESTIMATE_BURDEN(burden_input, gene_ch)
}

workflow.onComplete {
    log.info """
    ================================================================
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output directory: ${params.outdir}
    ================================================================
    """.stripIndent()
}
