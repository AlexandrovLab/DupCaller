#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (!params.sample_map) error "params.sample_map is required"
if (!params.reference)  error "params.reference is required"

// ─────────────────────────────────────────────────────────────────────────────
// Step 1: Index reference genome
// ─────────────────────────────────────────────────────────────────────────────
process INDEX_REFERENCE {
    container 'yuhecheng62/dupcaller:1.1.0-amd64'

    input:
    path reference
    path repeat_tsv

    output:
    path "${reference}.ref.h5", emit: ref_h5
    path "${reference}.tn.h5",  emit: tn_h5
    path "${reference}.hp.h5",  emit: hp_h5
    path "${reference}.str.h5", emit: str_h5
    path "${reference}.dbs.h5", emit: dbs_h5

    script:
    """
    DupCaller.py index -f ${reference} -rt ${repeat_tsv}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 2: Trim barcodes
// ─────────────────────────────────────────────────────────────────────────────
process TRIM_BARCODES {
    tag "${sample_id}:${type}"
    container 'yuhecheng62/dupcaller:1.1.0-amd64'

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}_trm_1.fastq"),
          path("${sample_id}_${type}_trm_2.fastq")

    script:
    """
    DupCaller.py trim \
        -i  ${read1} \
        -i2 ${read2} \
        -p  ${params.barcode_pattern} \
        -o  ${sample_id}_${type}_trm
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 3: Align reads
// ─────────────────────────────────────────────────────────────────────────────
process BWA_MEM {
    tag "${sample_id}:${type}"
    container 'biocontainers/bwa:v0.7.17_cv1'
    cpus params.threads

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)
    path bwa_index  // reference FASTA + all BWA index files staged together

    output:
    tuple val(sample_id), val(type), path("${sample_id}_${type}.sam")

    script:
    def ref_name = file(params.reference).name
    """
    bwa mem -C \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}_${type}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${ref_name} ${read1} ${read2} \
        > ${sample_id}_${type}.sam
    """
}

// Step 3b: sort + index (samtools lives in the GATK image)
process SAMTOOLS_SORT {
    tag "${sample_id}:${type}"
    container 'broadinstitute/gatk:4.2.6.0'
    cpus params.threads

    input:
    tuple val(sample_id), val(type), path(sam)

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.bam"),
          path("${sample_id}_${type}.bam.bai")

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}_${type}.bam ${sam}
    samtools index ${sample_id}_${type}.bam
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 4: Mark duplicates
// ─────────────────────────────────────────────────────────────────────────────
process MARK_DUPLICATES {
    tag "${sample_id}:${type}"
    container 'broadinstitute/gatk:4.2.6.0'

    input:
    tuple val(sample_id), val(type), path(bam), path(bai)

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.mkdped.bam"),
          path("${sample_id}_${type}.mkdped.bam.bai")

    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}_${type}.mkdped.bam \
        -M ${sample_id}_${type}.mkdp_metrics.txt \
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \
        --DUPLEX_UMI \
        --TAGGING_POLICY OpticalOnly \
        --BARCODE_TAG DB

    samtools index ${sample_id}_${type}.mkdped.bam
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 5: Call variants
// ─────────────────────────────────────────────────────────────────────────────
process CALL_VARIANTS {
    tag "${sample_id}"
    container 'yuhecheng62/dupcaller:1.1.0-amd64'
    cpus params.threads
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(tumor_bam),  path(tumor_bai),
          path(normal_bam), path(normal_bai)
    path dc_ref             // reference FASTA + .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(germline_vcf), path(germline_tbi)
    tuple path(noise_mask),   path(noise_tbi)
    path target_bed
    path indel_bed

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')

    script:
    def ref_name     = file(params.reference).name
    def normal_arg   = (normal_bam.name   != 'NO_FILE')        ? "-n ${normal_bam}"    : ""
    def germline_arg = (germline_vcf.name != 'NO_GERMLINE_VCF') ? "-g ${germline_vcf}" : ""
    def noise_arg    = (noise_mask.name   != 'NO_NOISE_MASK')   ? "-m ${noise_mask}"   : ""
    def target_arg   = (target_bed.name   != 'NO_TARGET_BED')   ? "-R ${target_bed}"   : ""
    def indel_arg    = (indel_bed.name    != 'NO_INDEL_BED')    ? "-id ${indel_bed}"   : ""
    """
    DupCaller.py call \
        -b  ${tumor_bam} \
        -f  ${ref_name} \
        -o  ${sample_id} \
        -p  ${task.cpus} \
        -r  ${params.regions} \
        ${normal_arg} \
        ${germline_arg} \
        ${noise_arg} \
        ${target_arg} \
        ${indel_arg} \
        -maf ${params.max_af} \
        -gaf ${params.germline_af_cutoff} \
        -d   ${params.min_n_depth} \
        -tt  ${params.trim_template} \
        -tr  ${params.trim_read} \
        -mq  ${params.mapq}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 6: Estimate mutational burden
// ─────────────────────────────────────────────────────────────────────────────
process ESTIMATE_BURDEN {
    tag "${sample_id}"
    container 'yuhecheng62/dupcaller:1.1.0-amd64'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(call_dir)
    path dc_ref        // reference FASTA + .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    path gene_bed

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')

    script:
    def ref_name   = file(params.reference).name
    def gene_arg   = (gene_bed.name != 'NO_GENE_BED') ? "-gb ${gene_bed}" : ""
    def clonal_arg = params.estimate_clonal ? "-c" : ""
    def dilute_arg = params.estimate_dilute ? "-d" : ""
    """
    DupCaller.py estimate \
        -i ${call_dir} \
        -f ${ref_name} \
        -r ${params.regions} \
        ${gene_arg} \
        ${clonal_arg} \
        ${dilute_arg}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Main workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Reference file channels ──────────────────────────────────────────────

    // All BWA index files staged together so 'bwa mem' finds them in the work dir
    bwa_ref_ch = Channel.fromPath([
        params.reference,
        "${params.reference}.fai",
        "${params.reference}.bwt",
        "${params.reference}.pac",
        "${params.reference}.ann",
        "${params.reference}.amb",
        "${params.reference}.sa"
    ], checkIfExists: true).collect()

    // DupCaller h5 index files — may be produced by INDEX_REFERENCE or pre-existing
    if (!params.skip_index) {
        if (!params.repeat_tsv) error "params.repeat_tsv is required when skip_index = false"
        repeat_tsv_ch = Channel.fromPath(params.repeat_tsv, checkIfExists: true)
        ref_fa  = Channel.fromPath(params.reference, checkIfExists: true)
        idx     = INDEX_REFERENCE(ref_fa, repeat_tsv_ch)
        dc_ref_ch = Channel.fromPath([
            params.reference,
            "${params.reference}.fai"
        ], checkIfExists: true)
            .concat(idx.ref_h5)
            .concat(idx.tn_h5)
            .concat(idx.hp_h5)
            .concat(idx.str_h5)
            .concat(idx.dbs_h5)
            .collect()
    } else {
        dc_ref_ch = Channel.fromPath([
            params.reference,
            "${params.reference}.fai",
            "${params.reference}.ref.h5",
            "${params.reference}.tn.h5",
            "${params.reference}.hp.h5",
            "${params.reference}.str.h5",
            "${params.reference}.dbs.h5"
        ], checkIfExists: true).collect()
    }

    // ── Optional resource file channels ─────────────────────────────────────

    germline_ch = params.germline_vcf
        ? Channel.value([
            file(params.germline_vcf,           checkIfExists: true),
            file("${params.germline_vcf}.tbi",  checkIfExists: true)
          ])
        : Channel.value([file('NO_GERMLINE_VCF'), file('NO_GERMLINE_TBI')])

    noise_ch = params.noise_mask
        ? Channel.value([
            file(params.noise_mask,             checkIfExists: true),
            file("${params.noise_mask}.tbi",    checkIfExists: true)
          ])
        : Channel.value([file('NO_NOISE_MASK'), file('NO_NOISE_TBI')])

    target_ch = params.target_bed
        ? Channel.value(file(params.target_bed, checkIfExists: true))
        : Channel.value(file('NO_TARGET_BED'))

    indel_ch = params.indel_bed
        ? Channel.value(file(params.indel_bed,  checkIfExists: true))
        : Channel.value(file('NO_INDEL_BED'))

    gene_ch = params.gene_bed
        ? Channel.value(file(params.gene_bed,   checkIfExists: true))
        : Channel.value(file('NO_GENE_BED'))

    // ── Parse sample map ─────────────────────────────────────────────────────
    // Emits (sample_id, type, fq1, fq2) for both tumor and normal per row

    Channel.fromPath(params.sample_map, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .flatMap { row ->
            [
                [row.sample_id, 'tumor',
                 file(row.tumor_fastq_1,  checkIfExists: true),
                 file(row.tumor_fastq_2,  checkIfExists: true)],
                [row.sample_id, 'normal',
                 file(row.normal_fastq_1, checkIfExists: true),
                 file(row.normal_fastq_2, checkIfExists: true)]
            ]
        }
        .set { reads_ch }

    // ── Steps 2–4: trim → align → mark-dup (tumor and normal in parallel) ───

    trimmed_ch = TRIM_BARCODES(reads_ch)
    sam_ch     = BWA_MEM(trimmed_ch, bwa_ref_ch)
    aligned_ch = SAMTOOLS_SORT(sam_ch)
    markdup_ch = MARK_DUPLICATES(aligned_ch)

    // ── Rejoin tumor and normal by sample_id ─────────────────────────────────

    tumor_md  = markdup_ch
        .filter { it[1] == 'tumor' }
        .map    { sid, _type, bam, bai -> [sid, bam, bai] }

    normal_md = markdup_ch
        .filter { it[1] == 'normal' }
        .map    { sid, _type, bam, bai -> [sid, bam, bai] }

    // join emits: [sample_id, tumor_bam, tumor_bai, normal_bam, normal_bai]
    call_input_ch = tumor_md.join(normal_md)

    // ── Step 5: Call variants ────────────────────────────────────────────────

    variants_ch = CALL_VARIANTS(
        call_input_ch,
        dc_ref_ch,
        germline_ch,
        noise_ch,
        target_ch,
        indel_ch
    )

    // ── Step 6: Estimate burden ──────────────────────────────────────────────

    ESTIMATE_BURDEN(variants_ch, dc_ref_ch, gene_ch)
}

workflow.onComplete {
    log.info """
    ================================================================
    DupCaller pipeline ${workflow.success ? 'completed' : 'FAILED'}
    Duration : ${workflow.duration}
    Output   : ${params.outdir}
    ================================================================
    """.stripIndent()
}
