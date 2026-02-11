#!/usr/bin/env nextflow

/*
 * DupCaller Nextflow Pipeline
 *
 * A pipeline for calling somatic mutations from barcoded error-corrected NGS data
 *
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

// DupCaller parameters
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
      --read1               Read 1 fastq file (paired with --read2)
      --read2               Read 2 fastq file (paired with --read1)
      --skip_trim           Skip trimming step if BAM already exists
      --skip_align          Skip alignment step if BAM already exists

    Normal sample:
      --normal_read1        Normal sample read 1 fastq
      --normal_read2        Normal sample read 2 fastq
      --normal_bam          Pre-aligned normal BAM file
      --skip_normal         Run without matched normal (requires --max_af)

    Reference resources:
      --germline_vcf        Indexed germline VCF with AF field
      --noise_mask          BED file(s) for noise masking (space-separated if multiple)
      --indel_bed           Enhanced panel of normal for indels

    Barcode trimming:
      --barcode_pattern     Barcode pattern (default: NNNXXXX for NanoSeq)

    Alignment:
      --threads             Number of threads (default: 4)
      --platform            Sequencing platform (default: ILLUMINA)

    Variant calling:
      --regions             Contigs for variant calling (default: chr1-22,X,Y)
      --max_af              Maximum allele frequency (default: 1.0, set to 0.1 for panel without normal)
      --germline_af_cutoff  Germline AF cutoff (default: 0.001)
      --min_n_depth         Minimum normal depth (default: 10)
      --max_zero_qual_frac  Maximum zero quality fraction (default: 0.5)
      --trim_template       Template end trim distance (default: 7)
      --trim_read           Read end trim distance (default: 7)
      --mapq                Minimum mapping quality (default: 40)

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
 * Process: Call variants with DupCaller
 */
process CALL_VARIANTS {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/variants", mode: 'copy'
    cpus params.threads

    input:
    tuple val(sample_name), path(bam), path(bai), path(metrics), path(reference), path(ref_h5), path(tn_h5), path(hp_h5)
    path(normal_bam)
    path(germline_vcf)
    path(noise_mask)
    path(indel_bed)

    output:
    tuple val(sample_name),
          path("${sample_name}_snv.vcf"),
          path("${sample_name}_indel.vcf"),
          path("${sample_name}_coverage.bed.gz"),
          path("${sample_name}_trinuc_by_duplex_group.txt"),
          path("${sample_name}_duplex_group_stats.txt"),
          path("${sample_name}_stats.txt"),
          path("${sample_name}.amp.tn.txt"),
          path("${sample_name}.amp.id.txt"),
          path("${sample_name}.dmg.tn.txt"),
          path("${sample_name}.dmg.id.txt")

    script:
    def normal_arg = normal_bam.name != 'NO_FILE' ? "-n ${normal_bam}" : ""
    def germline_arg = germline_vcf.name != 'NO_FILE' ? "-g ${germline_vcf}" : ""
    def noise_arg = noise_mask.name != 'NO_FILE' ? "-m ${noise_mask}" : ""
    def indel_arg = indel_bed.name != 'NO_FILE' ? "-id ${indel_bed}" : ""

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
          path(trinuc),
          path(duplex_stats),
          path(stats),
          path(amp_tn),
          path(amp_id),
          path(dmg_tn),
          path(dmg_id),
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

    script:
    def gene_arg = gene_bed.name != 'NO_FILE' ? "-gb ${gene_bed}" : ""
    def clonal_arg = params.estimate_clonal ? "-c true" : ""
    def dilute_arg = params.estimate_dilute ? "-d true" : ""

    """
    # Copy variant call outputs to current directory for estimate command
    ln -s ${snv_vcf} .
    ln -s ${indel_vcf} .
    ln -s ${coverage} .
    ln -s ${trinuc} .
    ln -s ${duplex_stats} .
    ln -s ${stats} .
    ln -s ${amp_tn} .
    ln -s ${amp_id} .
    ln -s ${dmg_tn} .
    ln -s ${dmg_id} .

    DupCaller.py estimate \\
        -i ${sample_name} \\
        -f ${reference} \\
        -r ${params.regions} \\
        ${gene_arg} \\
        ${clonal_arg} \\
        ${dilute_arg}
    """
}

/*
 * Main workflow
 */
workflow {
    // Check reference files
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
    ref_h5 = params.reference+".ref.h5"
    tn_h5 = params.reference+".tn.h5"
    hp_h5 = params.reference+".hp.h5"


    // Handle optional files
    germline_ch = params.germline_vcf ? Channel.fromPath(params.germline_vcf, checkIfExists: true) : Channel.value(file('NO_FILE'))
    noise_ch = params.noise_mask ? Channel.fromPath(params.noise_mask, checkIfExists: true) : Channel.value(file('NO_FILE'))
    indel_ch = params.indel_bed ? Channel.fromPath(params.indel_bed, checkIfExists: true) : Channel.value(file('NO_FILE'))
    gene_ch = params.gene_bed ? Channel.fromPath(params.gene_bed, checkIfExists: true) : Channel.value(file('NO_FILE'))

    // Tumor sample processing
    if (!params.skip_trim && params.read1 && params.read2) {
        // Start from FASTQ
        reads_ch = Channel.of([params.sample_name, file(params.read1), file(params.read2)])
        trimmed = TRIM_BARCODES(reads_ch)

        // Combine with reference
        align_input = trimmed.combine(ref_ch)
        aligned = BWA_ALIGN(align_input)
    } else if (params.skip_trim) {
        error "BAM file input not yet implemented. Please provide --read1 and --read2"
    }

    // Mark duplicates
    markdup = MARK_DUPLICATES(aligned)

    // Normal sample handling
    if (params.normal_bam) {
        normal_ch = Channel.fromPath(params.normal_bam, checkIfExists: true)
    } else if (params.normal_read1 && params.normal_read2 && !params.skip_normal) {
        normal_reads = Channel.of(["normal_${params.sample_name}", file(params.normal_read1), file(params.normal_read2)])
        normal_trimmed = TRIM_BARCODES(normal_reads)
        normal_align_input = normal_trimmed.combine(indexed_ref)
        normal_aligned = BWA_ALIGN(normal_align_input)
        normal_markdup = MARK_DUPLICATES(normal_aligned)
        normal_ch = normal_markdup.map { it[1] }
    } else {
        normal_ch = Channel.value(file('NO_FILE'))
    }

    // Variant calling
    call_input = markdup.combine(ref_ch)
    variants = CALL_VARIANTS(
        call_input,
        normal_ch,
        germline_ch,
        noise_ch,
        indel_ch
    )

    // Burden estimation
    burden_input = variants.combine(ref_ch)
    burden = ESTIMATE_BURDEN(burden_input, gene_ch)
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
