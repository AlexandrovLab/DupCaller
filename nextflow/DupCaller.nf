#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Defaults for params a config may not set (pipeline.config or a -c file
// overrides them). Both index steps default to skipped, i.e. pre-built
// indexes next to params.reference.
params.skip_index      = true
params.skip_bwa_index  = true
params.premade         = null
params.premade_dir     = null
params.premade_strain  = null
params.normal_bam      = null
params.tumor_only      = false
params.protocol        = null
params.barcode_list    = null
params.germline_vcf    = null
params.noise_mask      = null
params.target_bed      = null
params.indel_bed       = null
params.gene_bed        = null
params.seed            = null
params.p_threshold     = null
params.estimate_dilute = false
params.keep_trimmed_fastq = false
params.estimate_clonal = null
params.slurm_queue     = null
params.slurm_qos       = null
params.lsf_queue       = null
params.sge_queue       = null
params.sge_pe          = null
params.pbs_queue       = null
params.awsbatch_queue  = null
params.cluster_options = null
params.threads         = 1
params.max_cpus        = null
params.max_memory      = null
params.max_time        = null

if (!params.sample_map) error "params.sample_map is required"
if (params.premade) {
    if (!(params.premade in ['hg38', 'mm39']))
        error "params.premade = ${params.premade}: must be hg38 or mm39"
    if (params.reference)
        error "params.premade and params.reference are both set: premade downloads its own reference, so set reference = null"
    if (params.premade_strain) {
        if (params.premade != 'mm39') error "params.premade_strain applies only to premade = mm39"
        def strains = premadeManifest()*.file.findAll { it.startsWith('mgp.') && it.endsWith('_reh.vcf.gz') }
            .collect { it - 'mgp.' - '_reh.vcf.gz' }
        if (!(params.premade_strain in strains))
            error "params.premade_strain = ${params.premade_strain}: not in assets/premade/mm39.tsv; available strains: ${strains.join(', ')}"
    } else if (params.premade == 'mm39' && !params.germline_vcf) {
        log.warn "premade = mm39 without premade_strain or germline_vcf: calling runs without a germline VCF (-g)."
    }
} else if (!params.reference) {
    error "params.reference is required (or set params.premade to hg38 or mm39)"
}
// Scheduler params that ship as 'PLEASE CHANGE' in pipeline.config.
def SCHEDULER_PARAMS = [
    slurm:    ['slurm_queue', 'slurm_qos'],
    lsf:      ['lsf_queue'],
    sge:      ['sge_queue', 'sge_pe'],
    pbs:      ['pbs_queue'],
    pbspro:   ['pbs_queue'],
    awsbatch: ['awsbatch_queue'],
]
workflow.profile.tokenize(',').each { prof ->
    SCHEDULER_PARAMS.get(prof, []).each { name ->
        if (params[name] == 'PLEASE CHANGE')
            error "params.${name} is still 'PLEASE CHANGE': set it in pipeline.config (or null to omit) for -profile ${prof}"
    }
}
if (params.estimate_clonal) {
    error "params.estimate_clonal was removed: DupCaller.py estimate has no clonal option"
}
if (params.protocol) {
    if (!protocolSpec())
        error "params.protocol = ${params.protocol}: must be one of ${protocols().keySet().join(', ')}"
    log.info "protocol = ${params.protocol}: barcode_pattern = ${barcodePattern()}, trim_template = ${trimTemplate()}, trim_read = ${trimRead()} (replacing pipeline.config's values)."
}
if (barcodePattern().contains('B') && !params.barcode_list)
    error "barcode pattern ${barcodePattern()} has B (one barcode from a list): set barcode_list to the file of allowed barcodes, one per line"
if (params.barcode_list && !barcodePattern().contains('B'))
    error "params.barcode_list is set but barcode pattern ${barcodePattern()} has no B, so the list would not be used"
if (params.tumor_only) {
    if (params.normal_bam)
        error "params.tumor_only and params.normal_bam are both set: tumor_only calls without a matched normal, so set normal_bam = null"
    if ((params.max_af as double) >= 1)
        log.info "tumor_only with max_af = ${params.max_af} (uncapped): DupCaller.py call runs with -maf ${callMaxAf()}, so variants above that allele fraction (mostly germline) are not called."
}
// pipeline.config ships threads = max_cpus = 1; warn so a real sample is
// not aligned and called single-threaded by accident.
def requested_threads = (params.threads ?: 1) as int
def effective_threads = params.max_cpus ? Math.min(requested_threads, params.max_cpus as int) : requested_threads
if (effective_threads <= 1) {
    log.warn "threads = ${params.threads}, max_cpus = ${params.max_cpus}: bwa-mem2, samtools sort and DupCaller.py call will run single-threaded. " +
             "Set both in pipeline.config to the CPUs available for a real dataset."
}
// CALL_VARIANTS asks for 6 GB per CPU (nextflow.config); max_memory caps it.
if (params.max_memory) {
    try {
        def call_need = 6.GB * effective_threads
        def limit = params.max_memory as nextflow.util.MemoryUnit
        if (call_need > limit)
            log.warn "DupCaller.py call asks for 6 GB per CPU (${call_need} at ${effective_threads} CPUs) but max_memory = ${params.max_memory}; " +
                     "it gets ${params.max_memory}. Raise max_memory or lower threads if call runs out of memory."
    } catch (all) { }
}
// bwa-mem2 index needs ~28 bytes per reference base (nextflow.config sizes
// BWA_MEM2_INDEX the same way); a max_memory cap below that fails late.
// premade builds it once, unless premade_dir already holds it.
def bwa_index_fasta_bytes = null
if (params.premade) {
    if (!file("${premadeDir()}/${refName()}.0123").exists())
        bwa_index_fasta_bytes = premadeSpec().fasta_bytes
} else if (!params.skip_bwa_index && file(params.reference).exists()) {
    bwa_index_fasta_bytes = file(params.reference).size()
}
if (bwa_index_fasta_bytes && params.max_memory) {
    def index_need = new nextflow.util.MemoryUnit(28L * bwa_index_fasta_bytes + 2.GB.toBytes())
    try {
        def limit = params.max_memory as nextflow.util.MemoryUnit
        if (index_need > limit)
            log.warn "bwa-mem2 index of ${refName()} needs ~${index_need} but max_memory = ${params.max_memory}; " +
                     (params.premade
                         ? "raise max_memory for the first premade run (the index is then kept in ${premadeDir()})."
                         : "build the index once on a larger machine and set skip_bwa_index = true.")
    } catch (all) { }
}

// bwa-mem2 index files expected next to params.reference.
def BWA_MEM2_INDEX_SUFFIXES = ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac']

// ─────────────────────────────────────────────────────────────────────────────
// premade = 'hg38' / 'mm39'
//
// assets/premade/<genome>.tsv lists every hosted file (file, url, md5; md5
// may be empty). `archive` unpacks to `fasta`; the DupCaller index is
// <fasta>.{ref,tn,hp,str,dbs}.h5. mm39 germline VCFs are per strain
// (mgp.<strain>_reh.vcf.gz), so only params.premade_strain's is used.
// ─────────────────────────────────────────────────────────────────────────────
def premadeSpec() {
    [
        hg38: [archive:     'GRCh38.d1.vd1.fa.tar.gz',
               fasta:       'GRCh38.d1.vd1.fa',
               fasta_bytes: 3152288848L,
               germline:    'af-only-gnomad.hg38.vcf.gz',
               noise_mask:  ['SNP_GRCh38.wgns.bed.gz', 'NOISE_GRCh38.wgns.bed.gz'],
               regions:     null],
        mm39: [archive:     'mm39.fa.gz',
               fasta:       'mm39.fa',
               fasta_bytes: 2782787884L,
               germline:    params.premade_strain ? "mgp.${params.premade_strain}_reh.vcf.gz".toString() : null,
               noise_mask:  ['NOISE.mm39.bed.gz'],
               regions:     (1..19).collect { "chr${it}" }.join(" ") + " chrX"],
    ][params.premade]
}

def premadeManifest() {
    file("${projectDir}/assets/premade/${params.premade}.tsv", checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
}

// Downloads and the bwa-mem2 index are kept here (storeDir) and reused by
// every later run with the same premade_dir.
def premadeDir() {
    file(params.premade_dir ?: "${launchDir}/premade").resolve(params.premade).toString()
}

// Library protocols: barcode pattern (trim -p) and end trimming (call
// -tt/-tr), as in the main README's Protocol Settings table. Names are
// matched case-insensitively, ignoring '-' and '_' (e.g. RE-NanoSeq).
def protocols() {
    [
        renanoseq: [barcode_pattern: 'NNNXXXX',        trim_template: 7,  trim_read: 7],
        nanoseqv2: [barcode_pattern: 'NNNXX',          trim_template: 7,  trim_read: 7],
        udseq:     [barcode_pattern: 'NNNNNNNN',       trim_template: 15, trim_read: 8],
        // Tn5 transposon barcode from barcode_list, then the 19 nt ME.
        metacs:    [barcode_pattern: 'B' + 'X' * 19,   trim_template: 7,  trim_read: 7],
    ]
}

def protocolSpec() {
    params.protocol ? protocols()[params.protocol.toString().toLowerCase().replaceAll(/[-_]/, '')] : null
}

// Process scripts get these (and callRegions/callMaxAf) as val inputs, never
// by calling them: -resume's task hash does not see values returned by a
// function, so a changed protocol/regions/max_af would reuse a stale task.
def barcodePattern() { params.protocol ? protocolSpec().barcode_pattern : params.barcode_pattern }
def trimTemplate()   { params.protocol ? protocolSpec().trim_template   : params.trim_template }
def trimRead()       { params.protocol ? protocolSpec().trim_read       : params.trim_read }
// Read through a function on purpose: it stays out of BWA_MEM2_SORT's task
// hash, so changing keep_trimmed_fastq on -resume does not rerun alignment
// (it applies to alignments that run after the change).
def keepTrimmedFastq() { params.keep_trimmed_fastq.toString() == 'true' }

// Without a matched normal, an uncapped max_af would call germline
// variants, so tumor_only caps it at 0.1 unless set below 1.
def callMaxAf() {
    (params.tumor_only && (params.max_af as double) >= 1) ? 0.1 : params.max_af
}

def refName() {
    params.premade ? premadeSpec().fasta : file(params.reference).name
}

// pipeline.config's regions default is hg38's; premade = 'mm39' replaces
// it (and only it) with the mouse autosomes + X.
def callRegions() {
    def hg38_default = (1..22).collect { "chr${it}" }.join(" ") + " chrX"
    (params.premade && premadeSpec().regions && params.regions == hg38_default)
        ? premadeSpec().regions
        : params.regions
}

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
    f.name in ['NO_GERMLINE_VCF', 'NO_TARGET_BED', 'NO_INDEL_BED', 'NO_GENE_BED', 'NO_NOISE_MASK', 'NO_NORMAL_BAM', 'NO_BARCODE_LIST']
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
// Step 1 (premade only): download, unpack, check and bwa-mem2 index the
// hosted resources. storeDir skips every step whose output premade_dir
// already holds.
// ─────────────────────────────────────────────────────────────────────────────
process PREMADE_DOWNLOAD {
    tag "${name}"
    storeDir "${premadeDir()}"

    input:
    tuple val(name), val(url), val(md5)

    output:
    path "${name}"

    script:
    def md5_check = md5 ? "echo '${md5}  ${name}' | md5sum -c -" : ""
    """
    # gdown (Google Drive and plain https URLs) is bundled in DupCaller
    # images built after 1.2.12. An older image (e.g. one set with
    # dupcaller_container) lacks it, so it is then installed from PyPI into
    # the task dir. HOME may be read-only.
    export HOME=\$PWD
    if ! python -c 'import gdown' 2>/dev/null; then
        echo "gdown is not in this DupCaller image; installing gdown 6.1.0 from PyPI into the task directory" >&2
        PIP_DISABLE_PIP_VERSION_CHECK=1 pip install --quiet --no-cache-dir --target gdown_lib gdown==6.1.0
        export PYTHONPATH=gdown_lib
    fi
    if ! python -m gdown --no-cookies -O ${name} '${url}'; then
        echo "Could not download ${name} from ${url}" >&2
        echo "If Google Drive refused it ('Too many users have viewed or downloaded this file recently'), rerun later with -resume," >&2
        echo "or download it yourself and place it at ${premadeDir()}/${name} (a file already there is used as is)." >&2
        exit 1
    fi
    ${md5_check}
    """
}

process PREMADE_REFERENCE {
    storeDir "${premadeDir()}"

    input:
    path archive

    output:
    tuple path("${refName()}"), path("${refName()}.fai")

    script:
    def fasta  = refName()
    def unpack = archive.name.endsWith('.tar.gz') ? "tar -xzf ${archive} ${fasta}" : "gzip -dc ${archive} > ${fasta}"
    """
    ${unpack}
    python -c "import pysam; pysam.faidx('${fasta}')"
    """
}

// Fails here, before any alignment, if the hosted index was built for a
// different index schema than this DupCaller image reads.
process PREMADE_CHECK_INDEX {

    input:
    path h5s

    output:
    path "*.h5", includeInputs: true

    script:
    """
    python -c 'import sys; from DupCaller_sub.funcs.misc import check_h5_usable as c; bad = [m for ok, m in map(c, sys.argv[1:]) if not ok]; sys.exit("\\n".join(bad) or None)' *.h5
    """
}

process PREMADE_BWA_MEM2_INDEX {
    storeDir "${premadeDir()}"

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

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)
    path barcode_list   // allowed barcodes for a pattern with B, or the NO_BARCODE_LIST placeholder
    val  barcode_pattern

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}_trm_1.fastq"),
          path("${sample_id}_${type}_trm_2.fastq")

    script:
    def bl_arg = isPlaceholder(barcode_list) ? "" : "-bl ${barcode_list}"
    """
    DupCaller.py trim \
        -i  ${read1} \
        -i2 ${read2} \
        -p  ${barcode_pattern} \
        ${bl_arg} \
        -o  ${sample_id}_${type}_trm
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 3: Align (bwa-mem2) piped into samtools sort, then index
// ─────────────────────────────────────────────────────────────────────────────
// Once its output is written, each step below clears its input, which is
// another task's intermediate file (trimmed FASTQs here, the sorted BAM in
// MARK_DUPLICATES): the file is replaced by a sparse file of the same size,
// so it takes no disk space. Its producing task still counts as cached
// (its output exists), and cache 'lenient' hashes inputs by path and size
// only, so -resume still reuses this task. A cleared input starts with a
// zero byte (a FASTQ starts with '@', a BAM with 0x1f) and is refused.
// keep_trimmed_fastq = true leaves the trimmed FASTQs in place.
process BWA_MEM2_SORT {
    tag "${sample_id}:${type}"
    cpus params.threads
    cache 'lenient'

    input:
    tuple val(sample_id), val(type), path(read1), path(read2)
    path bwa_index  // reference FASTA + all bwa-mem2 index files staged together

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.bam"),
          path("${sample_id}_${type}.bam.bai")

    script:
    def ref_name = refName()
    def clear_fastq = keepTrimmedFastq() ? "" : """
    for f in ${read1} ${read2}; do
        real=\$(readlink -f "\$f")
        [ -e "\$(dirname "\$real")/.command.run" ] || continue   # only another task's output
        size=\$(stat -c %s "\$real")
        truncate -s 0 "\$real" && truncate -s "\$size" "\$real" \\
            || echo "WARNING: could not clear \$real" >&2
    done
    """
    """
    set -o pipefail
    for f in ${read1} ${read2}; do
        if [ "\$(head -c 1 "\$f" | od -An -tx1 | tr -d ' ')" = "00" ]; then
            echo "ERROR: \$f was cleared after an earlier run aligned it; delete the TRIM_BARCODES work directory of ${sample_id}:${type} (or run without -resume)" >&2
            exit 1
        fi
    done

    bwa-mem2 mem -C -T 0 \
        -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}_${type}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${ref_name} ${read1} ${read2} \
        | samtools sort -@ ${task.cpus} > ${sample_id}_${type}.bam
    samtools index -@ ${task.cpus} ${sample_id}_${type}.bam
    ${clear_fastq}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 4: Mark duplicates
// ─────────────────────────────────────────────────────────────────────────────
process MARK_DUPLICATES {
    tag "${sample_id}:${type}"
    cache 'lenient'   // see BWA_MEM2_SORT

    input:
    tuple val(sample_id), val(type), path(bam), path(bai)

    output:
    tuple val(sample_id), val(type),
          path("${sample_id}_${type}.mkdped.bam"),
          path("${sample_id}_${type}.mkdped.bam.bai")

    script:
    // Spill MarkDuplicates' sort files to the task work dir rather than
    // Java's default tmp dir (often a small node-local /tmp).
    def heap_gb = Math.max(1, (task.memory.toGiga() * 0.8) as int)
    """
    if [ "\$(head -c 1 ${bam} | od -An -tx1 | tr -d ' ')" = "00" ]; then
        echo "ERROR: ${bam} was cleared after an earlier run marked duplicates on it; delete the BWA_MEM2_SORT work directory of ${sample_id}:${type} (or run without -resume)" >&2
        exit 1
    fi

    gatk --java-options "-Xmx${heap_gb}g -Djava.io.tmpdir=\$PWD" MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}_${type}.mkdped.bam \
        -M ${sample_id}_${type}.mkdp_metrics.txt \
        --TMP_DIR \$PWD \
        --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*\$" \
        --DUPLEX_UMI \
        --TAGGING_POLICY OpticalOnly \
        --BARCODE_TAG DB \
        --CREATE_INDEX true

    # The GATK image has no samtools; MarkDuplicates writes X.mkdped.bai.
    mv ${sample_id}_${type}.mkdped.bai ${sample_id}_${type}.mkdped.bam.bai

    for f in ${bam} ${bai}; do
        real=\$(readlink -f "\$f")
        [ -e "\$(dirname "\$real")/.command.run" ] || continue   # only another task's output
        size=\$(stat -c %s "\$real")
        truncate -s 0 "\$real" && truncate -s "\$size" "\$real" \
            || echo "WARNING: could not clear \$real" >&2
    done
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 5: Call variants
// ─────────────────────────────────────────────────────────────────────────────
process CALL_VARIANTS {
    tag "${sample_id}"
    cpus params.threads
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id),
          path(tumor_bam),  path(tumor_bai),
          path(normal_bam), path(normal_bai)   // or the NO_NORMAL_BAM placeholders (tumor_only)
    path dc_ref             // reference FASTA + .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(germline_vcf), path(germline_tbi)
    path noise_mask_files     // noise/snp mask .bed.gz + .bed.gz.tbi files, or the NO_NOISE_MASK placeholder
    val  noise_mask_names     // basenames of just the mask files (not their .tbi), in -m order
    tuple path(target_bed), path(target_tbi)
    tuple path(indel_bed),  path(indel_tbi)
    val  opts               // [regions, max_af, trim_template, trim_read], resolved from protocol/tumor_only/premade

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')

    script:
    def ref_name     = refName()
    def normal_arg   = isPlaceholder(normal_bam)   ? "" : "-n ${normal_bam}"
    def germline_arg = isPlaceholder(germline_vcf) ? "" : "-g ${germline_vcf}"
    def noise_arg    = noise_mask_names ? "-m " + noise_mask_names.collect { "'${it}'" }.join(' ') : ""
    def target_arg   = isPlaceholder(target_bed)   ? "" : "-R ${target_bed}"
    def indel_arg    = isPlaceholder(indel_bed)    ? "" : "-id ${indel_bed}"
    // Omitted by default (DupCaller then picks and logs a random seed);
    // set params.seed to reproduce a run.
    def seed_arg     = params.seed != null         ? "--seed ${params.seed}" : ""
    // Omitted by default (DupCaller.py call's own default, 0.05, applies).
    def pt_arg       = params.p_threshold != null  ? "-pt ${params.p_threshold}" : ""
    """
    DupCaller.py call \
        -b  ${tumor_bam} \
        ${normal_arg} \
        -f  ${ref_name} \
        -o  ${sample_id} \
        -p  ${task.cpus} \
        -r  ${opts.regions} \
        ${germline_arg} \
        ${noise_arg} \
        ${target_arg} \
        ${indel_arg} \
        ${seed_arg} \
        ${pt_arg} \
        -maf ${opts.max_af} \
        -gaf ${params.germline_af_cutoff} \
        -d   ${params.min_n_depth} \
        -tt  ${opts.trim_template} \
        -tr  ${opts.trim_read} \
        -mq  ${params.mapq}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 6: Estimate mutational burden
// ─────────────────────────────────────────────────────────────────────────────
process ESTIMATE_BURDEN {
    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(call_dir, stageAs: 'call_in/*')
    path dc_ref        // reference FASTA + .fai + .ref.h5 + .tn.h5 + .hp.h5 + .str.h5 + .dbs.h5
    tuple path(gene_bed), path(gene_tbi)
    val  regions

    output:
    tuple val(sample_id), path("${sample_id}", type: 'dir')
    path "${sample_id}_estimate_params.log"

    script:
    def ref_name   = refName()
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
        -r ${regions} \
        ${gene_arg} \
        ${dilute_arg}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Main workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Reference file channels ──────────────────────────────────────────────

    if (params.premade) {
        // Every hosted file except the other strains' germline VCFs.
        def spec = premadeSpec()
        def rows = premadeManifest().findAll { row ->
            !row.file.startsWith('mgp.') || row.file in [spec.germline, "${spec.germline}.tbi".toString()]
        }
        premade_ch = PREMADE_DOWNLOAD(Channel.fromList(rows.collect { [it.file, it.url, it.md5 ?: ''] }))

        ref_fai    = PREMADE_REFERENCE(premade_ch.filter { it.name == spec.archive })
        fasta_ch   = ref_fai.map { fasta, _fai -> fasta }
        bwa_ref_ch = fasta_ch
            .concat(PREMADE_BWA_MEM2_INDEX(fasta_ch).flatten())
            .collect()

        def h5_names = ['ref', 'tn', 'hp', 'str', 'dbs'].collect { "${spec.fasta}.${it}.h5".toString() }
        h5_ch     = PREMADE_CHECK_INDEX(premade_ch.filter { it.name in h5_names }.collect())
        dc_ref_ch = ref_fai.flatten()
            .concat(h5_ch.flatten())
            .collect()
    } else {
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
    }

    // ── Optional resource file channels ─────────────────────────────────────

    // premade supplies germline_vcf and noise_mask unless they are set.
    def premade_germline = (params.premade && !params.germline_vcf) ? premadeSpec().germline : null
    def premade_masks    = (params.premade && !params.noise_mask)   ? premadeSpec().noise_mask : []

    germline_ch = premade_germline
        ? premade_ch.filter { it.name in [premade_germline, "${premade_germline}.tbi".toString()] }
            .toSortedList { a, b -> a.name <=> b.name }
        : Channel.value(optionalIndexedResource(params.germline_vcf, 'germline_vcf', 'NO_GERMLINE_VCF'))
    target_ch   = Channel.value(optionalIndexedResource(params.target_bed,   'target_bed',   'NO_TARGET_BED'))
    indel_ch    = Channel.value(optionalIndexedResource(params.indel_bed,    'indel_bed',    'NO_INDEL_BED'))
    gene_ch     = Channel.value(optionalIndexedResource(params.gene_bed,     'gene_bed',     'NO_GENE_BED'))

    if (premade_masks) {
        noise_files_ch = premade_ch
            .filter { it.name in premade_masks.collectMany { m -> [m, "${m}.tbi".toString()] } }
            .collect()
        noise_names_ch = Channel.value(premade_masks)
    } else {
        // params.noise_mask may be a single path or a list (-m takes several,
        // e.g. a SNP mask and a noise mask); normalize to a list.
        noise_mask_list = params.noise_mask
            ? (params.noise_mask instanceof List ? params.noise_mask : [params.noise_mask])
            : []
        noise_files_ch = Channel.value(
            noise_mask_list
                ? noise_mask_list.collectMany { m -> indexedResource(m, 'noise_mask') }
                : [placeholder('NO_NOISE_MASK')]
        )
        noise_names_ch = Channel.value(noise_mask_list.collect { file(it).name })
    }

    // ── Parse sample map ─────────────────────────────────────────────────────
    // Emits (sample_id, type, fq1, fq2). Normal rows are only emitted (and
    // only need normal_fastq_1/2 columns) when params.normal_bam and
    // params.tumor_only are unset -- see the normal branches below.

    Channel.fromPath(params.sample_map, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .flatMap { row ->
            def rows = [
                [row.sample_id, 'tumor',
                 file(row.tumor_fastq_1,  checkIfExists: true),
                 file(row.tumor_fastq_2,  checkIfExists: true)]
            ]
            if (!params.normal_bam && !params.tumor_only) {
                rows << [row.sample_id, 'normal',
                          file(row.normal_fastq_1, checkIfExists: true),
                          file(row.normal_fastq_2, checkIfExists: true)]
            }
            rows
        }
        .set { reads_ch }

    // ── Steps 2–4: trim → align → mark-dup (tumor, and normal unless shared or tumor_only) ─

    barcode_list_ch = Channel.value(params.barcode_list
        ? file(params.barcode_list, checkIfExists: true)
        : placeholder('NO_BARCODE_LIST'))
    trimmed_ch = TRIM_BARCODES(reads_ch, barcode_list_ch, barcodePattern())
    aligned_ch = BWA_MEM2_SORT(trimmed_ch, bwa_ref_ch)
    markdup_ch = MARK_DUPLICATES(aligned_ch)

    // ── Rejoin tumor and normal by sample_id ─────────────────────────────────

    tumor_md  = markdup_ch
        .filter { it[1] == 'tumor' }
        .map    { sid, _type, bam, bai -> [sid, bam, bai] }

    if (params.tumor_only) {
        // No matched normal: CALL_VARIANTS gets placeholders and omits -n.
        call_input_ch = tumor_md.map { sid, bam, bai ->
            [sid, bam, bai, placeholder('NO_NORMAL_BAM'), placeholder('NO_NORMAL_BAM.bai')]
        }
    } else if (params.normal_bam) {
        // One pre-aligned, indexed normal BAM shared by every sample_id;
        // trim/align/markdup are skipped for it.
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
        indel_ch,
        [regions: callRegions(), max_af: callMaxAf(),
         trim_template: trimTemplate(), trim_read: trimRead()]
    )

    // ── Step 6: Estimate burden ──────────────────────────────────────────────

    ESTIMATE_BURDEN(variants_ch, dc_ref_ch, gene_ch, callRegions())
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
