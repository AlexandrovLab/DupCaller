#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Defaults for params a config may not set (pipeline.config, or a -c file,
// overrides any of these). skip_bwa_index defaults to true so a config
// written before it existed never starts an unrequested bwa-mem2 index;
// likewise skip_index for the DupCaller index (Step 1a).
params.skip_index      = true
params.skip_bwa_index  = true
params.normal_bam      = null
params.germline_vcf    = null
params.noise_mask      = null
params.target_bed      = null
params.indel_bed       = null
params.gene_bed        = null
params.seed            = null
params.p_threshold     = null
params.estimate_dilute = false
params.estimate_clonal = null

if (!params.sample_map) error "params.sample_map is required"
if (!params.reference)  error "params.reference is required"
if (params.estimate_clonal) {
    error "params.estimate_clonal was removed: DupCaller.py estimate has no clonal option"
}

// bwa-mem2 index files expected next to params.reference (bwa-mem2 index).
def BWA_MEM2_INDEX_SUFFIXES = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac']

// ─────────────────────────────────────────────────────────────────────────────
// Optional inputs
//
// Nextflow `path` inputs are staged files, so an unset optional resource is
// bound to a real placeholder file shipped in assets/ (one distinctly named
// file per input, so two placeholders never collide in one task's work dir).
// DupCaller opens every optional VCF/BED with tabix, so each one is staged
// together with its .tbi.
// ─────────────────────────────────────────────────────────────────────────────
def placeholder(name) {
    file("${projectDir}/assets/${name}", checkIfExists: true)
}

// Exact names only, so a real input that happens to start with NO_ is used.
def isPlaceholder(f) {
    f.name in ['NO_GERMLINE_VCF', 'NO_TARGET_BED', 'NO_INDEL_BED', 'NO_GENE_BED', 'NO_NOISE_MASK']
}

def indexedResource(p, param_name) {
    if (!p.toString().endsWith('.gz')) {
        error "params.${param_name} = ${p}: must be bgzip-compressed (.gz) with a tabix index (.gz.tbi) next to it"
    }
    [file(p, checkIfExists: true), file("${p}.tbi", checkIfExists: true)]
}

def optionalIndexedResource(p, param_name, placeholder_name) {
    p ? indexedResource(p, param_name)
      : [placeholder(placeholder_name), placeholder("${placeholder_name}.tbi")]
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 1a: DupCaller reference index (optional)
// ─────────────────────────────────────────────────────────────────────────────
process INDEX_REFERENCE {
    container 'yuhecheng62/dupcaller:1.2.7'

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
// Step 1b: bwa-mem2 reference index (optional)
// ─────────────────────────────────────────────────────────────────────────────
process BWA_MEM2_INDEX {
    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0'

    input:
    path reference

    output:
    path "${reference}.{0123,amb,ann,bwt.2bit.64,pac}"

    script:
    """
    bwa-mem2 index ${reference}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 2: Trim barcodes
// ─────────────────────────────────────────────────────────────────────────────
process TRIM_BARCODES {
    tag "${sample_id}:${type}"
    container 'yuhecheng62/dupcaller:1.2.7'

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
// Step 3: Align reads (bwa-mem2)
// ─────────────────────────────────────────────────────────────────────────────
process BWA_MEM2 {
    tag "${sample_id}:${type}"
    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0'
    cpus params.threads

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)
    path bwa_index  // reference FASTA + all bwa-mem2 index files staged together

    output:
    tuple val(sample_id), val(type), path("${sample_id}_${type}.sam")

    script:
    def ref_name = file(params.reference).name
    """
    bwa-mem2 mem -C -T 0 \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}_${type}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${ref_name} ${read1} ${read2} \
        > ${sample_id}_${type}.sam
    """
}

// Step 3b: sort + index (samtools lives in the GATK image)
process SAMTOOLS_SORT {
    tag "${sample_id}:${type}"
    container 'broadinstitute/gatk:4.3.0.0'
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
    container 'broadinstitute/gatk:4.3.0.0'

    input:
    tuple val(sample_id), val(type), path(bam), path(bai)

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.mkdped.bam"),
          path("${sample_id}_${type}.mkdped.bam.bai")

    script:
    // MarkDuplicates' external sort spills large intermediate chunks to
    // disk once it exceeds its in-memory buffer -- without an explicit
    // TMP_DIR it falls back to Java's default tmp dir, which is usually
    // the compute node's small local /tmp, not sized for a real BAM's
    // spill files ("No space left on device" on anything but a tiny
    // input). Point it at the task's own work dir instead, which is
    // always on whatever (large) filesystem -w points at.
    def heap_gb = Math.max(1, (task.memory.toGiga() * 0.8) as int)
    """
    gatk --java-options "-Xmx${heap_gb}g -Djava.io.tmpdir=\$PWD" MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}_${type}.mkdped.bam \
        -M ${sample_id}_${type}.mkdp_metrics.txt \
        --TMP_DIR \$PWD \
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
    container 'yuhecheng62/dupcaller:1.2.7'
    cpus params.threads
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(tumor_bam),  path(tumor_bai),
          path(normal_bam), path(normal_bai)
    path dc_ref             // reference FASTA + .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(germline_vcf), path(germline_tbi)
    path noise_mask_files     // noise/snp mask .bed.gz + .bed.gz.tbi files, or the NO_NOISE_MASK placeholder
    val  noise_mask_names     // basenames of just the mask files (not their .tbi), in -m order
    tuple path(target_bed), path(target_tbi)
    tuple path(indel_bed),  path(indel_tbi)

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')

    script:
    def ref_name     = file(params.reference).name
    def germline_arg = isPlaceholder(germline_vcf) ? "" : "-g ${germline_vcf}"
    def noise_arg    = noise_mask_names ? "-m " + noise_mask_names.collect { "'${it}'" }.join(' ') : ""
    def target_arg   = isPlaceholder(target_bed)   ? "" : "-R ${target_bed}"
    def indel_arg    = isPlaceholder(indel_bed)    ? "" : "-id ${indel_bed}"
    // Omitted by default (DupCaller.py call itself then generates a fresh
    // random seed every run) -- set params.seed to pin it, e.g. to
    // reproduce/compare against a specific prior run's exact seed.
    def seed_arg     = params.seed != null         ? "--seed ${params.seed}" : ""
    // Omitted by default (DupCaller.py call's own default, 0.05, applies).
    def pt_arg       = params.p_threshold != null  ? "-pt ${params.p_threshold}" : ""
    """
    DupCaller.py call \
        -b  ${tumor_bam} \
        -n  ${normal_bam} \
        -f  ${ref_name} \
        -o  ${sample_id} \
        -p  ${task.cpus} \
        -r  ${params.regions} \
        ${germline_arg} \
        ${noise_arg} \
        ${target_arg} \
        ${indel_arg} \
        ${seed_arg} \
        ${pt_arg} \
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
    container 'yuhecheng62/dupcaller:1.2.7'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(call_dir, stageAs: 'call_in/*')
    path dc_ref        // reference FASTA + .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(gene_bed), path(gene_tbi)

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')
    path "${sample_id}_estimate_params.log"

    script:
    def ref_name   = file(params.reference).name
    def gene_arg   = isPlaceholder(gene_bed) ? "" : "-gb ${gene_bed}"
    def dilute_arg = params.estimate_dilute ? "-d" : ""
    """
    # sigProfilerPlotting's plotSBS() caches a template pickle inside its own
    # site-packages install dir by default, which fails under a container's
    # read-only root filesystem; redirect it into the (always-writable) task
    # work dir instead.
    export SIGPROFILERPLOTTING_VOLUME=\$PWD/spp_templates

    # estimate writes its outputs into the call directory. The staged one is
    # CALL_VARIANTS' own (cached) output, so work in a mirror of symlinks to
    # it instead: new files land in this task, and the call task's output
    # stays unchanged for -resume and for later runs with other options.
    # estimate appends base coverage to _stats.txt, so that file is a real
    # copy (without coverage lines from any earlier estimate) rather than a
    # link through which the append would reach the call output.
    cp -rs "\$(readlink -f ${call_dir})" ${sample_id}
    rm ${sample_id}/${sample_id}_stats.txt
    grep -v -e '^SBS Base Coverage' -e '^Indel Base Coverage' -e '^DBS Base Coverage' \
        "\$(readlink -f ${call_dir})/${sample_id}_stats.txt" > ${sample_id}/${sample_id}_stats.txt

    DupCaller.py estimate \
        -i ${sample_id} \
        -f ${ref_name} \
        -r ${params.regions} \
        ${gene_arg} \
        ${dilute_arg}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Main workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Reference file channels ──────────────────────────────────────────────

    // bwa-mem2 index files staged together with the FASTA so 'bwa-mem2 mem'
    // finds them in the work dir -- pre-built next to params.reference
    // (default), or built here by BWA_MEM2_INDEX when skip_bwa_index = false.
    if (!params.skip_bwa_index) {
        bwa_idx    = BWA_MEM2_INDEX(Channel.fromPath(params.reference, checkIfExists: true))
        bwa_ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
            .concat(bwa_idx.flatten())
            .collect()
    } else {
        bwa_ref_ch = Channel.fromPath(
            [params.reference] + BWA_MEM2_INDEX_SUFFIXES.collect { "${params.reference}.${it}" },
            checkIfExists: true
        ).collect()
    }

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

    germline_ch = Channel.value(optionalIndexedResource(params.germline_vcf, 'germline_vcf', 'NO_GERMLINE_VCF'))
    target_ch   = Channel.value(optionalIndexedResource(params.target_bed,   'target_bed',   'NO_TARGET_BED'))
    indel_ch    = Channel.value(optionalIndexedResource(params.indel_bed,    'indel_bed',    'NO_INDEL_BED'))
    gene_ch     = Channel.value(optionalIndexedResource(params.gene_bed,     'gene_bed',     'NO_GENE_BED'))

    // params.noise_mask may be a single path or a list of paths (DupCaller.py
    // call's -m/--noise takes nargs="+" -- e.g. a SNP mask and a noise mask
    // together, as the real benchmark runs do); normalize to a list.
    noise_mask_list = params.noise_mask
        ? (params.noise_mask instanceof List ? params.noise_mask : [params.noise_mask])
        : []
    noise_files_ch = Channel.value(
        noise_mask_list
            ? noise_mask_list.collectMany { m -> indexedResource(m, 'noise_mask') }
            : [placeholder('NO_NOISE_MASK')]
    )
    noise_names_ch = Channel.value(noise_mask_list.collect { file(it).name })

    // ── Parse sample map ─────────────────────────────────────────────────────
    // Emits (sample_id, type, fq1, fq2). Normal rows are only emitted (and
    // only need normal_fastq_1/2 columns) when params.normal_bam is unset --
    // see the shared-normal branch below.

    Channel.fromPath(params.sample_map, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .flatMap { row ->
            def rows = [
                [row.sample_id, 'tumor',
                 file(row.tumor_fastq_1,  checkIfExists: true),
                 file(row.tumor_fastq_2,  checkIfExists: true)]
            ]
            if (!params.normal_bam) {
                rows << [row.sample_id, 'normal',
                          file(row.normal_fastq_1, checkIfExists: true),
                          file(row.normal_fastq_2, checkIfExists: true)]
            }
            rows
        }
        .set { reads_ch }

    // ── Steps 2–4: trim → align → mark-dup (tumor, and normal unless shared) ─

    trimmed_ch = TRIM_BARCODES(reads_ch)
    sam_ch     = BWA_MEM2(trimmed_ch, bwa_ref_ch)
    aligned_ch = SAMTOOLS_SORT(sam_ch)
    markdup_ch = MARK_DUPLICATES(aligned_ch)

    // ── Rejoin tumor and normal by sample_id ─────────────────────────────────

    tumor_md  = markdup_ch
        .filter { it[1] == 'tumor' }
        .map    { sid, _type, bam, bai -> [sid, bam, bai] }

    if (params.normal_bam) {
        // A single already-aligned, already-indexed normal BAM shared across
        // every sample_id in the map (e.g. one matched normal reused across
        // many tumor-only mock/benchmark samples) -- skip trim/align/markdup
        // for it entirely rather than reprocessing the same normal reads
        // once per tumor sample.
        normal_bam_file = file(params.normal_bam, checkIfExists: true)
        normal_bai_file = file("${params.normal_bam}.bai", checkIfExists: true)
        call_input_ch = tumor_md.map { sid, bam, bai ->
            [sid, bam, bai, normal_bam_file, normal_bai_file]
        }
    } else {
        normal_md = markdup_ch
            .filter { it[1] == 'normal' }
            .map    { sid, _type, bam, bai -> [sid, bam, bai] }
        // join emits: [sample_id, tumor_bam, tumor_bai, normal_bam, normal_bai]
        call_input_ch = tumor_md.join(normal_md)
    }

    // ── Step 5: Call variants ────────────────────────────────────────────────

    variants_ch = CALL_VARIANTS(
        call_input_ch,
        dc_ref_ch,
        germline_ch,
        noise_files_ch,
        noise_names_ch,
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
