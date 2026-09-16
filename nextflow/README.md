# DupCaller Nextflow Pipeline

Runs the full DupCaller workflow end-to-end, one or more tumor/normal
samples at a time: barcode trimming → BWA-MEM alignment → GATK duplicate
marking → `DupCaller.py call` → `DupCaller.py estimate`.

## Prerequisites

- Nextflow >= 21.10.0
- Docker or Singularity (all processes run in containers — no local
  DupCaller/BWA/GATK install is required)
- A reference FASTA with:
  - a BWA index (`.bwt`/`.pac`/`.ann`/`.amb`/`.sa`), and
  - a DupCaller index (`.ref.h5`/`.tn.h5`/`.hp.h5`/`.str.h5`/`.dbs.h5`) —
    either pre-built (point `reference` at it and leave `skip_index = true`,
    the default), or built as part of the run (`skip_index = false` plus
    a `repeat_tsv`; see `DupCaller.py index --help`)

## Quick start: one sample

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
    -o ./results/SAMPLE_ID \
    -P singularity,local
```

Required: `-s -1 -2 -3 -4 -f`. Run with `-h` for the full flag list
(masks/germline/regions/profile/output dir/thread count all have sensible
defaults). This script also sets up a stable `-resume` launch directory
under `-o`, so re-running the same command after a failure resumes instead
of restarting from scratch.

On SLURM, wrap it with `examples/run_sample.slurm.sh` (or copy/adapt it) —
it just adds `#SBATCH` resource lines and calls the script above with
`-P singularity,local` and `-p "$SLURM_CPUS_PER_TASK"`:

```bash
sbatch --cpus-per-task=32 \
    --export=ALL,SAMPLE_ID=SAMPLE1,\
TUMOR_FASTQ_1=...,TUMOR_FASTQ_2=...,NORMAL_FASTQ_1=...,NORMAL_FASTQ_2=... \
    nextflow/examples/run_sample.slurm.sh
```

## Multiple samples: sample map

For more than one sample, run `DupCaller.nf` directly against a sample
map — a tab-separated file with a header row and one row per sample:

```
sample_id	tumor_fastq_1	tumor_fastq_2	normal_fastq_1	normal_fastq_2
SAMPLE1	/data/SAMPLE1_tumor_R1.fastq.gz	/data/SAMPLE1_tumor_R2.fastq.gz	/data/SAMPLE1_normal_R1.fastq.gz	/data/SAMPLE1_normal_R2.fastq.gz
SAMPLE2	...	...	...	...
```

Copy `pipeline.config` and fill in your paths (see **Parameters** below),
then:

```bash
nextflow run DupCaller.nf \
    -c nextflow.config -c pipeline.config \
    -profile singularity,local \
    -w /fast/scratch/work \
    -resume
```

Every sample's tumor and normal are processed independently through
trim/align/markdup, then rejoined by `sample_id` for calling and burden
estimation — samples run in parallel up to what your profile/executor
allows.

## Parameters

All parameters live in `pipeline.config` (documented inline there); the
table below is a quick reference.

### Required

| Parameter | Description |
|---|---|
| `sample_map` | Path to the tab-separated sample map (see above) |
| `reference` | Reference genome FASTA |

### Reference indexing

| Parameter | Description | Default |
|---|---|---|
| `skip_index` | Skip `DupCaller.py index` (reuse an existing `.ref.h5`/`.tn.h5`/`.hp.h5`/`.str.h5`/`.dbs.h5` set next to `reference`) | `true` |
| `repeat_tsv` | PERF-format repeat TSV, required only when `skip_index = false` | `null` |

### Optional resource files

| Parameter | Description | Default |
|---|---|---|
| `germline_vcf` | Tabix-indexed germline VCF with an AF field | `null` |
| `noise_mask` | Tabix-indexed noise/SNP mask BED, or a list of several (e.g. `["snp_mask.bed.gz", "noise_mask.bed.gz"]`) — passed to `DupCaller.py call`'s `-m`/`--noise` (`nargs="+"`) | `null` |
| `target_bed` | Restrict calling to target regions | `null` |
| `indel_bed` | Indel-enhanced panel of normals | `null` |
| `gene_bed` | Gene BED for per-gene duplex coverage in burden estimation | `null` |

### Calling / burden options

| Parameter | Description | Default |
|---|---|---|
| `barcode_pattern` | Barcode pattern (`N`=barcode base, `X`=skipped) | `NNNXXXX` |
| `regions` | Contigs to call, space-separated | `chr1`...`chr22 chrX` |
| `threads` | Threads for BWA and `DupCaller.py call` | `1` |
| `max_af` | Max allele fraction in matched normal (`-maf`/`--naf`) — matches `DupCaller.py call`'s own default | `0.01` |
| `germline_af_cutoff` | Skip positions above this population AF (`-gaf`) | `0.001` |
| `min_n_depth` | Minimum normal depth for a called variant (`-d`) | `10` |
| `trim_template` | Ignore mutations within N bp of template ends (`-tt`) | `7` |
| `trim_read` | Ignore mutations within N bp of read ends (`-tr`) | `7` |
| `mapq` | Minimum alignment MAPQ (`-mq`) | `40` |
| `seed` | Pin the Monte Carlo seed (omit to let `DupCaller.py call` pick a fresh one every run) | unset |
| `estimate_clonal` | Treat multi-molecule mutations as one in burden estimation | `false` |
| `estimate_dilute` | Set when sample and matched normal share starting DNA material | `false` |

These calling/burden defaults are deliberately kept in sync with
`DupCaller.py call`/`estimate`'s own CLI defaults — if you change one in
the underlying tool, update `pipeline.config` to match.

### Compute limits (`nextflow.config`)

`max_cpus`/`max_memory`/`max_time` cap the per-process resource requests
computed in `nextflow.config`; per-process defaults (memory/time, scaled
up automatically on retry) live in that file's `process { withName: ... }`
blocks — edit them there if a step needs more than its default.

## Execution profiles

Combine one container profile with one executor profile, comma-separated:

| Profile | Purpose |
|---|---|
| `docker` | Run containers via Docker |
| `singularity` | Run containers via Singularity/Apptainer (typical on HPC) |
| `local` | Execute directly on the current node |
| `slurm` | Submit each process as its own SLURM job — edit the `slurm` block in `nextflow.config` (`queue`/`clusterOptions`) for your partition/QOS/account first |

e.g. `-profile singularity,local` (single node, everything runs locally
inside containers) or `-profile singularity,slurm` (each step becomes its
own SLURM job). `examples/run_sample.slurm.sh` uses the first pattern:
one `sbatch` allocation runs the whole per-sample pipeline locally within
it, which is usually simpler to reason about than nesting SLURM inside
SLURM.

## Output

Each sample's `CALL_VARIANTS` and `ESTIMATE_BURDEN` outputs are copied
(via `publishDir`) into `${outdir}/${sample_id}/` — the same directory
`DupCaller.py call -o` / `estimate -i` would produce standalone (VCFs,
coverage beds, error-profile tables, burden/signature outputs). Nextflow's
own `work/` directory (set via `-w`) holds all intermediate trim/align/
markdup files and can be deleted once you're happy with the results in
`outdir`.

## Troubleshooting

**"params.reference is required" / "params.sample_map is required"** —
these two are the only params without a default; both must be set.

**`checkIfExists` error on a reference/mask/VCF file** — every file input
is staged with `checkIfExists: true`, so a typo'd path fails immediately
at pipeline-launch time rather than partway through a run. Check that
tabix-indexed inputs (`germline_vcf`, `noise_mask`, `indel_bed`) have
their `.tbi` sitting right next to them.

**Reusing a completed run's error profiles** — `DupCaller.py call` only
re-runs round-0 error-rate learning if the six error files don't already
exist at the resolved output prefix. If you're rerunning the same
`outdir`/`sample_id` after a DupCaller code change that affects error
learning, delete or move that sample's prior output first so learning
actually recomputes instead of silently reusing stale rates.

**Resuming after a failure** — always pass `-resume` and reuse the same
`-w` work directory; `examples/run_dupcaller_sample.sh` does this for you
automatically via a stable per-sample launch directory.

## Known gaps

`test_pipeline.sh` in this directory predates the current `sample_map`-based
parameter interface (it still checks for `--sample_name`/`--read1`-style
flags that no longer exist) and does not currently reflect the pipeline
as implemented above — treat it as unmaintained until it's rewritten
against the current interface.
