# DupCaller Nextflow Pipeline

Runs the full DupCaller workflow end-to-end, one or more tumor/normal
(or tumor-only) samples at a time: barcode trimming → bwa-mem2 alignment → GATK duplicate
marking → `DupCaller.py call` → `DupCaller.py estimate`.

Preprocessing uses the same tools and commands as the main README:
bwa-mem2 2.3 (`bwa-mem2 mem -C -T 0`, piped into `samtools sort`),
samtools 1.22.1 and GATK 4.3.0.0 MarkDuplicates (see **Containers**).
(bwa-mem2 v2.3 still reports itself as `2.2.1` in `bwa-mem2 version` and
the BAM `@PG` header.)

## Prerequisites

- Nextflow 23.10.1 to 24.10.4, the versions CI runs the pipeline on
  (`manifest.nextflowVersion`); other versions print a warning at launch
  and are untested
- Docker or Singularity/Apptainer (all processes run in containers — no
  local DupCaller/bwa-mem2/GATK install is required; see **Containers**)
- A reference FASTA with, next to it, its samtools FASTA index (`.fai`),
  its DupCaller index (five files: `.ref.h5`, `.tn.h5`, `.hp.h5`,
  `.str.h5`, `.dbs.h5`, built with `DupCaller.py index -f reference.fa
  -rt repeats.tsv`) and its bwa-mem2 index. Build or download them
  once before the first run (see **Before running a real dataset**); the
  pipeline then reuses them (`skip_index = true`, `skip_bwa_index = true`,
  the defaults). For hg38 or mm39, `--premade` can download and prepare
  all of this instead (see the main README's
  [Quick Start](../README.md#quick-start) and **Pre-made genomes**).

## Quick Start

The Quick Start (hg38/mm39 with `--premade`, FASTQs as the only input)
is in the main README: [Quick Start](../README.md#quick-start).

## Protocol settings

Set `protocol` in `pipeline.config` (`--protocol` on the `nextflow run`
command line, `-c` in `examples/run_dupcaller_sample.sh`) to your library
protocol. It sets the barcode pattern and end trimming, replacing
`barcode_pattern`, `trim_template` and `trim_read`:

| Protocol | `protocol` | barcode pattern | `trim_template` | `trim_read` |
|---|---|---|---|---|
| RE-NanoSeq | `renanoseq` | `NNNXXXX` | 7 | 7 |
| NanoSeq v2 | `nanoseqv2` | `NNNXX` | 7 | 7 |
| UDSeq | `udseq` | `NNNNNNNN` | 15 | 8 |
| META-CS | `metacs` | `BXXXXXXXXXXXXXXXXXXX` (Tn5 barcode from `barcode_list`, then the 19 nt ME) | 7 | 7 |

Names are matched case-insensitively, ignoring `-` and `_` (`RE-NanoSeq`
works). `metacs` also needs `barcode_list` (`-l` in the script): the
allowed Tn5 barcodes, one per line, passed to `DupCaller.py trim -bl`.
Without `protocol`, `barcode_pattern` (default `NNNXXXX`, `-b` in the
script), `trim_template` and `trim_read` (default 7) are used as set; a
custom pattern with `B` also needs `barcode_list`.

`trim_template` ignores mutations within N bp of either end of the
sequenced fragment and `trim_read` within N bp of either end of each read
(`DupCaller.py call -tt`/`-tr`).

**Without a matched normal:** set `tumor_only = true` (`--tumor_only true`,
or `-T` in `examples/run_dupcaller_sample.sh`). The sample map then needs
only `sample_id`, `tumor_fastq_1` and `tumor_fastq_2`, and
`DupCaller.py call` runs without `-n` and with `-maf 0.1` (unless
`max_af` is set below 1), so variants above 10% allele fraction (mostly
germline) are not called. `-maf 0.1` also masks positions with raw tumor
depth below 10x. `tumor_only` cannot be combined with `normal_bam`. The
command-line equivalent is in the main README
([Without Matched Normal](../README.md#without-matched-normal)).

## Before running a real dataset

> **The shipped `pipeline.config` is single-threaded.** `threads = 1` and
> `max_cpus = 1`; `BWA_MEM2_SORT` and `CALL_VARIANTS` take
> their CPU count from these, so a whole-genome sample run with the
> defaults is aligned and called on one core (the pipeline prints a
> warning at launch). Set both to the CPUs you have, e.g.
> `threads = 32` and `max_cpus = 32` (`examples/run_dupcaller_sample.sh`
> sets both from `-p`, default all cores). Also check `max_memory`
> (default `96.GB`): `bwa-mem2 mem | samtools sort` asks for 32 GB + 1 GB
> per CPU and `DupCaller.py call` for 6 GB per CPU (192 GB at 32 CPUs),
> each capped by `max_memory` (`-M` in the example script, which does not
> derive it from `-p`).

**1. Indexes (do this once, outside the pipeline).** This is the
recommended setup:

- **DupCaller index** (`.ref.h5`/`.tn.h5`/`.hp.h5`/`.str.h5`/`.dbs.h5`):
  download the pre-built hg38 or mm39 index from the main README's
  [Resources](../README.md#resources) (place the FASTA and its `.h5` files
  in one directory, and keep the FASTA under the name the `.h5` files are
  prefixed with, e.g. `GRCh38.d1.vd1.fa` for `GRCh38.d1.vd1.fa.ref.h5`: the
  pipeline looks for `<reference>.ref.h5` etc.), or build it with
  `DupCaller.py index` (main README, Step 1).
- **samtools FASTA index**: `samtools faidx reference.fa`.
- **bwa-mem2 index** (`.0123`/`.amb`/`.ann`/`.bwt.2bit.64`/`.pac`): no
  pre-built bwa-mem2 index is hosted, so build it once on a machine with
  enough memory and reuse it for every run:

  ```bash
  bwa-mem2 index reference.fa   # ~28 bytes of RAM per base: ~90 GB and ~15-20 minutes for GRCh38
  ```

  The same container the pipeline uses works for this, e.g.
  `singularity exec docker://quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0 bwa-mem2 index reference.fa`.
  An index from plain `bwa index` (`.bwt`/`.sa`) is not used by bwa-mem2.

Then point `reference` at the FASTA and leave `skip_index = true` and
`skip_bwa_index = true`.

*Advanced / high-memory only:* `skip_bwa_index = false` builds the
bwa-mem2 index inside the run (`BWA_MEM2_INDEX`), and `skip_index = false`
plus a `repeat_tsv` builds the DupCaller index (`INDEX_REFERENCE`). The
`BWA_MEM2_INDEX` memory request is sized from the FASTA (28 bytes per
base + 2 GB, ~90 GB for GRCh38) and capped by `max_memory`, so the job
needs a node, queue and `max_memory` at least that large; the pipeline
warns at launch when `max_memory` is below the estimate. Many schedulers
reject requests that size, which is why the pre-built route is the
default.

**2. Barcode pattern.** Set `protocol` to your library protocol (see
**Protocol settings**), or check `barcode_pattern` (default `NNNXXXX`)
against your library prep; a wrong pattern silently breaks duplex family
formation.

## Pre-made genomes (`--premade`)

For GRCh38 (hg38) and GRCm39 (mm39), the pipeline can fetch everything in
the main README's [Resources](../README.md#resources) table itself, in
place of `reference` and the index steps:

```bash
cd /path/to/DupCaller   # repository root
nextflow run nextflow/DupCaller.nf \
    -c nextflow/nextflow.config -c nextflow/pipeline.config \
    --premade hg38 --premade_dir /shared/dupcaller_premade \
    -profile singularity,local -resume
```

or `premade = "hg38"` in `pipeline.config` (leave `reference = null`).
With `examples/run_dupcaller_sample.sh`, use `-G hg38` instead of `-f`, and
`-D /shared/dupcaller_premade` for the download directory (the script's
default is `./dupcaller_premade` under the directory you run it from, so
without `-D` each new directory downloads ~50 GB again); for mm39, `-S` sets
`premade_strain`.

| | hg38 | mm39 |
|---|---|---|
| Reference FASTA | GDC `GRCh38.d1.vd1.fa` | UCSC `mm39.fa` |
| DupCaller index | lab Google Drive | lab Google Drive |
| Germline VCF | GATK `af-only-gnomad.hg38.vcf.gz` | mouse strain VCF `mgp.<strain>_reh.vcf.gz`, only with `premade_strain` |
| Noise mask(s) | NanoSeq `SNP_GRCh38.wgns.bed.gz` + `NOISE_GRCh38.wgns.bed.gz` | `NOISE.mm39.bed.gz` |
| `regions` default | chr1-chr22, chrX | chr1-chr19, chrX |

On the first run, `PREMADE_DOWNLOAD` downloads each file listed in
`assets/premade/<genome>.tsv` (for mm39, only `premade_strain`'s germline
VCF) with `gdown` and checks its md5
against the TSV (every hg38 file; the mm39 FASTA and DupCaller index), `PREMADE_REFERENCE` unpacks the FASTA and builds its
`.fai` (pysam's `faidx`), `PREMADE_CHECK_INDEX` checks the five DupCaller
index files (`.ref.h5`/`.tn.h5`/`.hp.h5`/`.str.h5`/`.dbs.h5`) and
stops the run, before any alignment, if one cannot be read or was built
for a different index schema than this DupCaller image reads (it runs on
every run, including on files reused from `premade_dir`), and
`PREMADE_BWA_MEM2_INDEX` builds the bwa-mem2 index, which is not hosted.
All of it is kept in `premade_dir/<genome>/` (`storeDir`), so every later
run with the same `premade_dir`, including a different sample or a fresh
work directory, reuses it without downloading or indexing again.

- **Disk and memory:** hg38 needs about 50 GB in `premade_dir` (25 GB of
  DupCaller index, 17 GB bwa-mem2 index, 3 GB gnomAD VCF, 4 GB FASTA);
  mm39 about 40 GB. Building the bwa-mem2 index needs ~90 GB of RAM for
  hg38 and ~80 GB for mm39, so set `max_memory` at least that high for
  the first run (the pipeline warns at launch otherwise).
- **Network:** the download tasks need internet access from wherever they
  run (with `-profile ...,slurm`, the compute nodes). `gdown` is bundled
  in DupCaller images after 1.2.12, so they install nothing from PyPI; the
  1.2.12 image installs `gdown` from PyPI in the task (see **Containers**),
  so it also needs PyPI access. If compute nodes are offline, run the first `--premade` run
  in an interactive job on a node with network access, or pre-download
  the files as described below. The DupCaller index and mm39 files are on Google Drive, which
  blocks a large file for up to a day after many downloads ("Too many
  users have viewed or downloaded this file recently"). A failed download
  is retried twice, a minute or two apart; after that, rerun later with
  `-resume` (files already in `premade_dir` are kept), or download the
  missing file yourself (in a browser, or via **Make a copy** to your own
  Drive) and put it in `premade_dir/<genome>/` under the name in
  `assets/premade/<genome>.tsv`. A file already there is used as is, with
  no md5 check (check it yourself with `md5sum` against the TSV's `md5`
  column); the DupCaller index files are still schema-checked.
- **Overrides:** `germline_vcf` or `noise_mask`, if set, replace the
  premade ones; `regions`, if changed from the hg38 default, is used as
  is. `skip_index`, `skip_bwa_index` and `repeat_tsv` are ignored.
- **mm39 germline VCF:** set `premade_strain` to your mouse strain (e.g.
  `C57BL_6NJ`; the run lists the valid names if it is wrong). Without it,
  calling runs without a germline VCF.

## One sample with your own reference

The easiest way to run a single tumor/normal pair is
`examples/run_dupcaller_sample.sh`, which builds the sample map and a
run-specific config for you:

```bash
nextflow/examples/run_dupcaller_sample.sh \
    -s SAMPLE_ID \
    -1 tumor_R1.fastq.gz  -2 tumor_R2.fastq.gz \
    -3 normal_R1.fastq.gz -4 normal_R2.fastq.gz \
    -f /path/to/reference.fa \
    -m snp_mask.bed.gz,noise_mask.bed.gz \
    -g germline.vcf.gz \
    -p 32 \
    -M 192.GB \
    -o ./results/SAMPLE_ID \
    -P singularity,local
```

Required: `-s -1 -2`, `-f` (or `-G hg38|mm39`, see **Pre-made genomes**),
plus one of `-3 -4` (normal FASTQs), `-B` (an
already-aligned, indexed normal BAM) or `-T` (tumor only, no matched
normal; see **Protocol settings**). Run with `-h` for the full flag list
(masks/germline/regions/profile/output dir/thread count all have sensible
defaults). `-M` sets `max_memory`; without it, `pipeline.config`'s `96.GB`
caps `DupCaller.py call` (6 GB per CPU) below what 32 CPUs ask for. This
script also sets up a stable `-resume` launch directory
under `-o`, so re-running the same command after a failure resumes instead
of restarting from scratch.

On SLURM, wrap it with `examples/run_sample.slurm.sh` (or copy/adapt it).
It adds `#SBATCH` resource lines and calls the script above with
`-P singularity,local` (`PROFILE` overrides it), `-p "$SLURM_CPUS_PER_TASK"`
and `-M` set to the job's Slurm memory allocation (`--mem`; set
`MAX_MEMORY`, e.g. `MAX_MEMORY=192.GB`, to override). It hard-codes no
partition, QOS, account or module: pass the scheduler flags your site
needs on the `sbatch` command line, and set `CONTAINER_MODULE` (e.g.
`singularity` or `apptainer`) if the container runtime must be loaded with
`module load`. **Edit it before the first `sbatch`**: as shipped, its
placeholder paths fail. Set `REPO_ROOT`, `OUTDIR`, `REFERENCE`, `SNP_MASK`
and `NOISE_MASK` in the script (`REPO_ROOT` and `OUTDIR` can instead be
passed in `--export`). Then, from the repository root:

```bash
sbatch --partition=YOUR_PARTITION --account=YOUR_ACCOUNT \
    --cpus-per-task=32 --mem=200G \
    --export=ALL,REPO_ROOT=$PWD,OUTDIR=/path/to/results/SAMPLE1,SAMPLE_ID=SAMPLE1,CONTAINER_MODULE=singularity,\
TUMOR_FASTQ_1=...,TUMOR_FASTQ_2=...,NORMAL_FASTQ_1=...,NORMAL_FASTQ_2=... \
    nextflow/examples/run_sample.slurm.sh
```

For tumor-only samples, pass `TUMOR_ONLY=1` instead of `NORMAL_FASTQ_1`/`NORMAL_FASTQ_2`.
For a pre-made genome instead of `REFERENCE` and the masks, add
`PREMADE=hg38,PREMADE_DIR=/shared/dupcaller_premade` (the script's `-G`/`-D`;
for mm39 also `PREMADE_STRAIN=...`, `-S`) to `--export`; the first run needs
`--mem` of at least ~90 GB to build the bwa-mem2 index.

## Multiple samples: sample map

For more than one sample, run `DupCaller.nf` directly against a sample
map — a tab-separated file with a header row and one row per sample:

```
sample_id	tumor_fastq_1	tumor_fastq_2	normal_fastq_1	normal_fastq_2
SAMPLE1	/data/SAMPLE1_tumor_R1.fastq.gz	/data/SAMPLE1_tumor_R2.fastq.gz	/data/SAMPLE1_normal_R1.fastq.gz	/data/SAMPLE1_normal_R2.fastq.gz
SAMPLE2	...	...	...	...
```

With `normal_bam` or `tumor_only` set, the `normal_fastq_1`/`normal_fastq_2`
columns are not read and may be left out.

Fill in your paths in `nextflow/pipeline.config`, or in a copy passed with
`-c` in its place (see **Parameters** below), then, from the repository
root (relative paths such as `sample_map` resolve from where you launch):

```bash
nextflow run nextflow/DupCaller.nf \
    -c nextflow/nextflow.config -c nextflow/pipeline.config \
    -profile singularity,local \
    -w /fast/scratch/work \
    -resume
```

Every sample's tumor and normal are processed independently through
trim/align/markdup, then rejoined by `sample_id` for calling and burden
estimation — samples run in parallel up to what your profile/executor
allows.

## Parameters

All parameters except `seed` live in `pipeline.config` (documented inline
there); any of them can also be given on the command line as
`--<name> <value>`. The tables below are a quick reference.

### Required

| Parameter | Description |
|---|---|
| `sample_map` | Path to the tab-separated sample map (see above) |
| `reference` | Reference genome FASTA (unless `premade` is set) |

Results are published to `outdir` (default `./results`, relative to the
launch directory; see **Output**).

### Matched normal

By default each sample's normal is aligned from its `normal_fastq_1`/`normal_fastq_2`.

| Parameter | Description | Default |
|---|---|---|
| `normal_bam` | One already-aligned, indexed (`.bai`) normal BAM shared by every sample; skips trim/align/markdup for the normal | `null` |
| `tumor_only` | No matched normal: `call` runs without `-n`, and with `-maf 0.1` unless `max_af` is below 1 (see **Protocol settings**). Cannot be combined with `normal_bam` | `false` |

### Pre-made genome

| Parameter | Description | Default |
|---|---|---|
| `premade` | `hg38` or `mm39`: download and prepare that genome's reference, indexes, germline VCF and noise mask(s) instead of using `reference` (see **Pre-made genomes**) | `null` |
| `premade_dir` | Where premade files are downloaded and kept, reused across runs | `<launch dir>/premade` |
| `premade_strain` | mm39 only: mouse strain whose germline VCF to use, e.g. `C57BL_6NJ` | `null` |

### Reference indexing

| Parameter | Description | Default |
|---|---|---|
| `skip_index` | Skip `DupCaller.py index` (reuse an existing `.ref.h5`/`.tn.h5`/`.hp.h5`/`.str.h5`/`.dbs.h5` set next to `reference`) | `true` |
| `repeat_tsv` | PERF-format repeat TSV, required only when `skip_index = false` | `null` |
| `skip_bwa_index` | Skip `bwa-mem2 index` (reuse an existing `.0123`/`.amb`/`.ann`/`.bwt.2bit.64`/`.pac` set next to `reference`). `false` builds it in the run: high-memory, see **Before running a real dataset** | `true` |

### Optional resource files

Every file here must be bgzip-compressed (`.gz`) with its tabix index
(`.gz.tbi`) next to it: DupCaller opens them all with tabix, and the
pipeline stages each file together with its `.tbi` (a path that does not
end in `.gz`, or a missing `.tbi`, stops the run at launch). Leaving a
parameter `null` omits the matching DupCaller option; the task is then
given an empty placeholder from `assets/` instead.

| Parameter | Description | DupCaller option | Default |
|---|---|---|---|
| `germline_vcf` | Germline VCF with an AF field (`.vcf.gz` + `.tbi`) | `call -g` | `null` |
| `noise_mask` | Noise/SNP mask BED (`.bed.gz` + `.tbi`), or a list of several (e.g. `["snp_mask.bed.gz", "noise_mask.bed.gz"]`) | `call -m` | `null` |
| `target_bed` | Restrict calling to target regions (`.bed.gz` + `.tbi`) | `call -R` | `null` |
| `indel_bed` | Indel-enhanced panel of normals (`.bed.gz` + `.tbi`) | `call -id` | `null` |
| `gene_bed` | Gene BED for per-gene duplex coverage, 4th column `{gene}_{exon}` (`.bed.gz` + `.tbi`; format in the main README's [Per-Gene Coverage](../README.md#per-gene-coverage)). For hg38 exomes, `resources/xgen-exome-hyb-panel-v2-targets-hg38.bed.gz` | `estimate -gb` | `null` |

### Calling / burden options

| Parameter | Description | Default |
|---|---|---|
| `protocol` | Library protocol: `renanoseq`, `nanoseqv2`, `udseq` or `metacs`; replaces `barcode_pattern`, `trim_template` and `trim_read` (see **Protocol settings**) | `null` |
| `barcode_pattern` | Barcode pattern (`N`=barcode base, `X`=skipped, `B`=one barcode from `barcode_list`); ignored when `protocol` is set | `NNNXXXX` |
| `barcode_list` | Allowed barcodes, one per line, for a pattern with `B` (`trim -bl`); required for `protocol = "metacs"` | `null` |
| `regions` | Contigs to call, space-separated | `chr1`...`chr22 chrX` (`chr1`...`chr19 chrX` with `premade = "mm39"`) |
| `threads` | Threads for bwa-mem2, `samtools sort` and `DupCaller.py call`; capped by `max_cpus`. **Raise both for real data** | `1` |
| `max_af` | Max tumor allele fraction to call a mutation (`-maf`/`--maxAF`; not `--naf`, the matched-normal VAF cap). Values below 1 also mask positions with raw tumor depth < 1/`max_af`. With `tumor_only`, a value of 1 is replaced by 0.1 | `1` |
| `germline_af_cutoff` | Skip positions above this population AF (`-gaf`) | `0.001` |
| `min_n_depth` | Minimum normal depth for a called variant (`-d`) | `10` |
| `trim_template` | Ignore mutations within N bp of template ends (`-tt`) | `7` |
| `trim_read` | Ignore mutations within N bp of read ends (`-tr`) | `7` |
| `mapq` | Minimum alignment MAPQ (`-mq`) | `40` |
| `seed` | Pin the Monte Carlo seed (omit to let `DupCaller.py call` pick a fresh one every run). Not in `pipeline.config`: pass `--seed N` or add it to your config | unset |
| `p_threshold` | Strand independence hypothesis SBS/indel filter (`-pt`/`--p_threshold`): per-strand p threshold, not corrected for the number of calls (strand p <= this); failed calls get FILTER `strand_independence`; `0` disables. Unset uses `DupCaller.py call`'s default (`0.05`) | unset |
| `estimate_dilute` | Set when sample and matched normal share starting DNA material (`estimate -d`) | `false` |
| `keep_trimmed_fastq` | Keep the trimmed FASTQs in the work directory instead of clearing them after alignment (see **Disk (work directory)**) | `false` |

`estimate_clonal` was removed: it passed `-c`, which `DupCaller.py
estimate` has never accepted. Setting it now stops the run with an error.

These calling/burden defaults are deliberately kept in sync with
`DupCaller.py call`/`estimate`'s own CLI defaults — if you change one in
the underlying tool, update `pipeline.config` to match.

### Compute limits

`max_cpus`/`max_memory`/`max_time` (set in `pipeline.config`) cap the
per-process resource requests computed in `nextflow.config`; per-process
defaults (memory/time, scaled up automatically on retry) live in that file's `process { withName: ... }`
blocks — edit them there if a step needs more than its default. The
multi-threaded steps (`BWA_MEM2_SORT`, `CALL_VARIANTS`) ask
for `threads` CPUs, capped by `max_cpus`; both default to 1. The defaults
are sized for a whole-genome duplex library (~65 GB fq.gz per read):
`bwa-mem2 mem | samtools sort` asks for 32 GB + 1 GB per CPU (24 h),
MarkDuplicates 64 GB (24 h), trimming 12 h and `DupCaller.py call` 6 GB per CPU (96 h),
so keep `max_memory` (default `96.GB`) and `max_time` (default `120.h`)
at least that high. `BWA_MEM2_INDEX`
(only with `skip_bwa_index = false`) asks for 28 bytes per reference base
+ 2 GB, ~90 GB for GRCh38. On a smaller machine, lower `max_memory` (the
example script's `-M`, e.g. `-M 24.GB`); the local executor refuses any
task that asks for more memory than the machine has.

### Disk (work directory)

Intermediate files are cleared once the next step has written its output:
`BWA_MEM2_SORT` clears the trimmed FASTQs (uncompressed) and
`MARK_DUPLICATES` clears the sorted BAM. A cleared file is a sparse file of
the original size, so it takes no disk space but `-resume` still reuses the
finished tasks. Disk use per library peaks during alignment (trimmed FASTQs,
~5x the input fq.gz, plus the sorted BAM being written) and afterwards is
mostly the duplicate-marked BAM, which stays in the work directory (only `call` and
`estimate` results are published to `outdir`). If a step must rerun after
its input was cleared (e.g. its own work directory was deleted), it stops
with an error naming the upstream task whose work directory to delete, or
run without `-resume`.

`keep_trimmed_fastq = true` (`--keep_trimmed_fastq true`, `-K` in
`examples/run_dupcaller_sample.sh`, `KEEP_TRIMMED_FASTQ=1` in
`examples/run_sample.slurm.sh`) keeps the trimmed FASTQs
(`work/.../{sample}_{tumor|normal}_trm_{1,2}.fastq`). Changing it on
`-resume` does not rerun alignment: it applies only to alignments that
run after the change.

## Execution profiles

Combine one container profile with one executor profile, comma-separated:

| Profile | Purpose |
|---|---|
| `docker` | Run containers via Docker |
| `singularity` | Run containers via Singularity/Apptainer (typical on HPC) |
| `local` | Execute directly on the current node |
| `slurm` | Each process is its own SLURM job; set `slurm_queue` (partition) and `slurm_qos` |
| `lsf` | Each process is its own LSF job; set `lsf_queue` |
| `sge` | Each process is its own SGE/UGE job; set `sge_queue` and `sge_pe` (parallel environment) |
| `pbs` / `pbspro` | Each process is its own Torque/OpenPBS (`pbs`) or PBS Pro (`pbspro`) job; set `pbs_queue` |
| `awsbatch` | Each process is an AWS Batch job; set `awsbatch_queue`. Use on its own, without `docker`/`singularity` |

Scheduler settings live in `pipeline.config` and ship as `PLEASE CHANGE`;
the run stops until the ones for your profile are set (`null` omits a
flag, e.g. no QOS or the default queue). `cluster_options` adds extra
submit flags to any HPC profile, e.g. `"--account=myproj"` for SLURM or
`"-P myproj"` for LSF. Anything else site-specific (LSF memory units,
PBS resource strings) goes in your own config file passed with another
`-c`.

Only `local` and SLURM (via `examples/run_sample.slurm.sh`) have been run
on real data; `lsf`, `sge`, `pbs`, `pbspro` and `awsbatch` are untested
starting points.

`awsbatch` also needs an S3 work directory (`-w s3://bucket/work`), a
region (`aws.region` in a config file, `AWS_REGION`, or `~/.aws/config`),
an AMI with the AWS CLI (or `aws.batch.cliPath` pointing at it), and a
compute environment with instances large enough for the biggest request:
`BWA_MEM2_INDEX` asks for 90 GB, so on AWS keep `skip_bwa_index = true`
and pass a pre-built index unless the instances are that large. Inputs and
`outdir` can be `s3://` paths.

e.g. `-profile singularity,local` (single node, everything runs locally
inside containers) or `-profile singularity,slurm` (each step becomes its
own SLURM job). `examples/run_sample.slurm.sh` uses the first pattern:
one `sbatch` allocation runs the whole per-sample pipeline locally within
it, which is usually simpler to reason about than nesting SLURM inside
SLURM.

## Containers

| Process(es) | Image | Registry |
|---|---|---|
| `INDEX_REFERENCE`, `TRIM_BARCODES`, `CALL_VARIANTS`, `ESTIMATE_BURDEN`, `PREMADE_DOWNLOAD`, `PREMADE_REFERENCE`, `PREMADE_CHECK_INDEX` | `ghcr.io/alexandrovlab/dupcaller:1.2.12` | GitHub Container Registry (AlexandrovLab) |
| `BWA_MEM2_INDEX`, `PREMADE_BWA_MEM2_INDEX` | `quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0` | quay.io |
| `BWA_MEM2_SORT` | `community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113` (bwa-mem2 2.3, samtools 1.22.1; nf-core's `bwamem2/mem` image) | Seqera Community |
| `MARK_DUPLICATES` | `quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0` | quay.io |

The DupCaller image is built from this repository's `Dockerfile` at each
release's tag and published by `.github/workflows/release.yml` after the
Nextflow mock test passes on it; it includes `gdown` for
`PREMADE_DOWNLOAD`. (The 1.2.12 image predates this and lacks `gdown`; with
it, `PREMADE_DOWNLOAD` falls back to installing `gdown` 6.1.0 from PyPI
into its task directory.) `--dupcaller_container IMAGE` (or
`dupcaller_container` in a config file) replaces the image in every
process above, e.g. with a local `.sif` file or a mirror.

The quay.io and Seqera Community images are pinned in `nextflow.config` by `sha256`
digest rather than by tag, so every run pulls exactly the same image even
if a tag is later re-pushed (the tag is noted next to each digest).
Singularity does not accept a reference with both a tag and a digest, so
only the digest is used. The GATK image contains the same GATK 4.3.0.0
(HTSJDK 3.0.1, Picard 2.27.5) as `broadinstitute/gatk:4.3.0.0`, which
earlier versions of this pipeline used.

**Image caching.** With Singularity, set
`NXF_SINGULARITY_CACHEDIR` to a shared directory (as
`examples/run_sample.slurm.sh` does): each image is pulled and converted
once and reused by every later run, including on offline compute nodes. To
use a local copy of the DupCaller image:

```bash
singularity pull dupcaller-1.2.12.sif docker://ghcr.io/alexandrovlab/dupcaller:1.2.12
nextflow run nextflow/DupCaller.nf ... --dupcaller_container "$PWD/dupcaller-1.2.12.sif"
```

## Output

Each sample's `CALL_VARIANTS` and `ESTIMATE_BURDEN` outputs are copied
(via `publishDir`) into `${outdir}/${sample_id}/` — the same directory
`DupCaller.py call -o` / `estimate -i` would produce standalone (VCFs,
coverage beds, error-profile tables, burden/signature outputs) — plus
`${outdir}/${sample_id}_estimate_params.log`. `Mutation number per genome`
in the `_burden.txt` files is per haploid genome (see the main README). Nextflow's
own `work/` directory (set via `-w`) holds all intermediate trim/align/
markdup files and can be deleted once you're happy with the results in
`outdir`.

## Troubleshooting

**"params.reference is required (or set params.premade to hg38 or mm39)"**
— `reference` ships as `null` in `pipeline.config`; set it, or set
`premade`. ("params.sample_map is required" appears only if `sample_map`,
which `pipeline.config` sets to `sample.map`, is cleared.)

**`checkIfExists` error on a reference/mask/VCF file** — every file input
is staged with `checkIfExists: true`, so a typo'd path fails immediately
at pipeline-launch time rather than partway through a run. Check that
every optional resource (`germline_vcf`, `noise_mask`, `target_bed`,
`indel_bed`, `gene_bed`) is bgzipped and has its `.tbi` right next to it
(a sorted BED: `bgzip file.bed && tabix -p bed file.bed.gz`; a VCF:
`bgzip file.vcf && tabix -p vcf file.vcf.gz`).

**Resuming after a failure** — always pass `-resume` and reuse the same
`-w` work directory; `examples/run_dupcaller_sample.sh` does this for you
automatically via a stable per-sample launch directory.

## Tests

`tests/mock_pipeline/run_nextflow_pipeline.sh OUTDIR [docker|singularity] all`
runs the whole pipeline (including `DupCaller.py index` and `bwa-mem2
index`) on a 3.2 kb synthetic contig: once with every optional resource
unset, then once per optional input (`germline_vcf`, `noise_mask`,
`target_bed`, `indel_bed`, `gene_bed`, `estimate_dilute`, `normal_bam`,
`tumor_only`, and the `udseq` and `metacs` protocols; the `metacs` case
must reproduce the matched-normal reference once its barcodes are mapped
back).
`tests/mock_pipeline/test_mock_pipeline_nextflow.py` wraps it in pytest
and checks that each case passed the right option to DupCaller (and that
the default case passed none). `DUPCALLER_IMAGE=IMAGE` runs it with that
DupCaller image instead of the pinned one and checks that every DupCaller
task ran in it. CI (`.github/workflows/test.yml`) builds the image from the
checked-out commit and runs this test on it with Nextflow 23.10.1 and
24.10.4; `.github/workflows/release.yml` does the same for a release tag
before publishing the image, then again with the published image pulled
anonymously. See the main README's
[Verify Installation](../README.md#verify-installation) for all three
test suites.
