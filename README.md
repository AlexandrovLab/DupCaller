# DupCaller

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/AlexandrovLab/DupCaller/blob/main/README.md) [![License](<https://img.shields.io/badge/License-BSD%202--Clause-orange.svg>)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://github.com/AlexandrovLab/DupCaller/actions/workflows/test.yml/badge.svg)](https://github.com/AlexandrovLab/DupCaller/actions)[![Uptime Robot status](https://img.shields.io/uptimerobot/status/m795312784-02766a79f207f67626cef289)](https://stats.uptimerobot.com/jjqW4Ulymx)

DupCaller is a tool for calling somatic mutations and calculating somatic mutational burden from barcoded error-corrected next generation sequencing (ecNGS/duplex sequencing) data with matched normal (e.g. NanoSeq, UDSeq).

---

## Table of Contents

- [Quick Start](#quick-start)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Verify Installation](#verify-installation)
- [Pipeline](#pipeline)
  - [Protocol Settings](#protocol-settings)
  - [Step 1: Index Reference Genome](#step-1-index-reference-genome)
  - [Step 2: Trim Barcodes](#step-2-trim-barcodes)
  - [Step 3: Align Reads](#step-3-align-reads)
  - [Step 4: Mark Duplicates](#step-4-mark-duplicates)
  - [Step 5: Call Variants](#step-5-call-variants)
  - [Step 6: Estimate Mutational Burden](#step-6-estimate-mutational-burden)
  - [Step 7: Summarize Across Samples](#step-7-summarize-across-samples)
  - [Optional: Aggregate Error Profiles Across Samples](#optional-aggregate-error-profiles-across-samples)
- [Results](#results)
- [End-to-End Examples](#end-to-end-examples)
- [Resources](#resources)
- [Citation](#citation)
- [Copyright](#copyright)
- [Contact](#contact)

---

## Quick Start

Use this if you have no preferred reference genome file and want to follow the developers' defaults: GRCh38 (human) or GRCm39 (mouse), with the reference FASTA, DupCaller index, germline VCF and noise/SNP mask(s) that `--premade` downloads (listed in [Resources](#resources)).

**What you need:**

- **This repository**, cloned: `git clone https://github.com/AlexandrovLab/DupCaller.git`
- **Anaconda or equivalent**: if you use another package manager such as mamba or pixi, replace the installation commands with their equivalents.
- **[Nextflow](https://www.nextflow.io/) 23.10.1 to 24.10.4** on your `PATH` (the versions CI tests; the commands below install 24.10.4). Other versions print a warning at launch and are untested.
- **Docker** (for local/cloud) **or Singularity** (for HPC) on your `PATH`.  If you are on HPC and do not have Singularity, contact your admin to install it, or follow [Installation](#installation) and [Pipeline](#pipeline) instead.
- **Paired-end FASTQs** (R1 and R2) for the duplex sample and, if you have one, its matched normal. Without a matched normal, see [Without a matched normal](#quick-start-without-a-matched-normal) below. We do not recommend calling samples <=100x duplex coverage without a matched normal.
- **Internet access** on the first run, from wherever the tasks run, to pull the containers and download the genome (see [Pre-made genome downloads](#pre-made-genome-downloads) if your compute nodes are offline).
- **Disk:** ~50 GB for the downloaded genome (hg38), plus space for the Nextflow work directory and results. Expect ~300GB of space per 30x data.
- **CPUs and memory:** see the box below.

> **Compute requirements — set these for real data.**
>
> - `nextflow/pipeline.config` ships with conservative defaults, `threads = 1`, `max_cpus = 1` and `max_memory = 96.GB`, which align and call a whole genome on **one core**. Always set `--threads`, `--max_cpus` and `--max_memory` to what your machine or cluster job allows, as in the commands below (e.g. `--threads 32 --max_cpus 32 --max_memory 192.GB`). `DupCaller.py call` asks for 6 GB per CPU and `bwa-mem2 mem | samtools sort` for 32 GB + 1 GB per CPU, each capped by `--max_memory`. The pipeline warns at launch when it would run single-threaded, or when `--max_memory` is below what calling or the bwa-mem2 index needs.
> - **The first run builds the bwa-mem2 index locally: ~90 GB of RAM on one node for hg38 (~80 GB for mm39), about 20 minutes, once per `--premade_dir`.** No pre-built bwa-mem2 index is hosted, so set `--max_memory` to at least 90.GB for that first run. Every later run reuses the index and needs only the per-CPU memory above.
> - `--premade` (below) is the recommended route: it downloads the reference and the pre-built DupCaller index, so you never run `DupCaller.py index` yourself. With your own reference, build both indexes once on a large-memory node (bwa-mem2 index needs ~28 bytes of RAM per reference base) and reuse them; see [Before running a real dataset](nextflow/README.md#before-running-a-real-dataset).

**Installation**

With Singularity (HPC):

```bash
# 1. Nextflow and Singularity
module load singularity # Module name varies by cluster. Contact HPC admin if you have problem.
conda create -y -n nextflow -c conda-forge -c bioconda nextflow=24.10.4
conda activate nextflow

# 2. This repository
git clone https://github.com/AlexandrovLab/DupCaller.git
cd DupCaller
```

With Docker (Linux workstation with sudo):

```bash
# 1. Docker, and permission to run it without sudo
curl -fsSL https://get.docker.com | sudo sh
sudo usermod -aG docker "$USER"
newgrp docker                    # or log out and back in
docker run --rm hello-world      # check: prints "Hello from Docker!"
### If using macOS, install and launch Docker Desktop and ignore step 1

# 2. Nextflow
conda create -y -n nextflow -c conda-forge -c bioconda nextflow=24.10.4
conda activate nextflow

# 3. This repository
git clone https://github.com/AlexandrovLab/DupCaller.git
cd DupCaller
```

**Configure your scheduler profile:**

- `-profile singularity,local` runs every step on the current machine (e.g. inside one batch job). To submit each step as its own cluster job instead, use `-profile singularity,slurm` (or `lsf`, `sge`, `pbs`, `pbspro`), and first replace the `"PLEASE CHANGE"` values for that scheduler in the "Scheduler settings" block of `nextflow/pipeline.config` with your cluster's information (nothing site-specific is preset): the run stops with an error until they are set. `null` omits a flag (e.g. no QOS), and `cluster_options` adds extra submit flags such as your account. For SLURM, for example:

  ```groovy
  slurm_queue     = "your_partition"       // sbatch --partition
  slurm_qos       = "your_qos"             // sbatch --qos, or null
  cluster_options = "--account=your_account"
  ```

**Fill in the sample map:** copy the example into the DupCaller directory and edit it:

```bash
cd /path/to/DupCaller   # change to the actual path
cp nextflow/examples/sample.map sample.map                   # with a matched normal
cp nextflow/examples/sample_tumor_only.map sample.map        # or: without a matched normal
```

Replace `SAMPLE1` with your sample name and each `/path/to/...` with the absolute path of that FASTQ. Alternatively, create your own tab-delimited `sample.map` with the same header as the example.

**With a matched normal (Singularity):** replace each `/path/to/...` and adjust `--protocol`, the resource flags and `-profile` as described in the notes below the commands.

```bash
conda activate nextflow
export NXF_SINGULARITY_CACHEDIR=/path/to/sif_cache   # Singularity only
cd /path/to/DupCaller                                # change to the actual path

nextflow run nextflow/DupCaller.nf \
    -c nextflow/nextflow.config -c nextflow/pipeline.config \
    --premade hg38 --premade_dir /path/to/dupcaller_premade \
    --protocol renanoseq \
    --sample_map sample.map --outdir results \
    --threads 32 --max_cpus 32 --max_memory 192.GB \
    -profile singularity,slurm -resume
```

<a id="quick-start-without-a-matched-normal"></a>
**Without a matched normal (Singularity):** use the tumor-only sample map and add `--tumor_only true`:

```bash
conda activate nextflow
export NXF_SINGULARITY_CACHEDIR=/path/to/sif_cache   # Singularity only
cd /path/to/DupCaller                                # change to the actual path

nextflow run nextflow/DupCaller.nf \
    -c nextflow/nextflow.config -c nextflow/pipeline.config \
    --premade hg38 --premade_dir /path/to/dupcaller_premade \
    --protocol renanoseq \
    --tumor_only true \
    --sample_map sample.map --outdir results \
    --threads 32 --max_cpus 32 --max_memory 192.GB \
    -profile singularity,slurm -resume
```

- `-profile`: `singularity,slurm` submits each step to SLURM (configure it as above); replace `slurm` with your scheduler (`lsf`, `sge`, `pbs`, `pbspro`), or use `singularity,local` to run everything on the current machine. With Docker, use `-profile docker,local` and skip `NXF_SINGULARITY_CACHEDIR`.
- `NXF_SINGULARITY_CACHEDIR` (on shared storage) keeps the container images, so they are downloaded once instead of into every run's work directory.
- `--premade hg38` or `--premade mm39`. For mm39, add `--premade_strain STRAIN` (e.g. `C57BL_6NJ`) to use that mouse strain's germline VCF.
- `--premade_dir` is where the genome is downloaded and kept (~50 GB for hg38). The first run downloads it and builds the bwa-mem2 index (~90 GB of RAM, about 20 minutes for hg38); every later run with the same `--premade_dir`, for any sample, reuses it. Put it on shared storage.
- Each downloaded file is checked against the md5 listed in `nextflow/assets/premade/<genome>.tsv` (all hg38 files; for mm39, the FASTA and the DupCaller index). On every run, before any alignment, `PREMADE_CHECK_INDEX` opens each of the five DupCaller index files and compares its index schema version with the one the pipeline's DupCaller image reads. An outdated or incompatible index stops the run at that step with an error naming the file; it is never used silently. (With your own reference, `DupCaller.py call` and `estimate` run the same check when they start.)
- `--protocol` is your library protocol: `renanoseq`, `nanoseqv2`, `udseq` or `metacs`. It sets the barcode pattern and how far from fragment and read ends mutations are ignored (see [Protocol Settings](#protocol-settings)); a wrong barcode pattern silently breaks duplex family formation. `metacs` also needs `--barcode_list barcode_list.txt`, the allowed Tn5 barcodes, one per line. For another protocol, give the barcode pattern with `--barcode_pattern` instead (N = barcode base, X = skipped base).
- `--threads`/`--max_cpus`/`--max_memory`: **required for real data** (the shipped defaults are 1 CPU; see the box above). Set them to your machine or the largest job your cluster allows. `--max_memory` should be at least 6 GB × `--threads` (what `DupCaller.py call` asks for; e.g. `--threads 16 --max_memory 96.GB`), and at least ~90 GB for the first run.
- Results go to `results/SAMPLE1/`. Rerunning the same command in the same directory resumes after a failure (`-resume`).

<a id="pre-made-genome-downloads"></a>
**Pre-made genome downloads.** The first `--premade` run downloads the files listed in `nextflow/assets/premade/<genome>.tsv` (GDC, GATK and Google Drive) from wherever its `PREMADE_DOWNLOAD` tasks run: the current machine with `-profile ...,local`, or the compute nodes with `-profile ...,slurm`. The download tool (`gdown`) is bundled in DupCaller images after 1.2.12, so nothing is installed at run time (the 1.2.12 image still installs `gdown` from PyPI in the download task), but the download itself needs internet access. If your compute nodes have none, either run that first `--premade` run inside an interactive job on a node that does, or download each file in the TSV yourself on such a node and put it in `<premade_dir>/<genome>/` under the name in the TSV's `file` column; files already there are used as is (see [Pre-made genomes](nextflow/README.md#pre-made-genomes---premade)). Later runs need no network access for the genome.

If you have your own reference genome, see [Before running a real dataset](nextflow/README.md#before-running-a-real-dataset) and [One sample with your own reference](nextflow/README.md#one-sample-with-your-own-reference) in the Nextflow README, or run the command-line steps in [Pipeline](#pipeline) below.

---

## Prerequisites

DupCaller requires python>=3.11 to run. Earlier versions are not supported.
The complete DupCaller pipeline also requires the following tools for data preprocessing. The Nextflow pipeline (`nextflow/`) runs exactly these versions and commands. Other versions may work but have not been validated.

- bwa-mem2 2.3 (https://github.com/bwa-mem2/bwa-mem2), run as `bwa-mem2 mem -C -T 0` (Step 3). (The v2.3 release still prints `2.2.1` from `bwa-mem2 version` and in the BAM `@PG` header.)
- GATK 4.3.0.0 (https://github.com/broadinstitute/gatk/releases), for MarkDuplicates (Step 4)
- samtools, for sorting and indexing BAMs
- htslib, which provides `tabix` and `bgzip` for compressed, indexed genomic files (recommended installation: `conda install bioconda::htslib`)
- [PERF](https://github.com/rkmlab/perf) (Pattern-based Exhaustive Repeat Finder) v0.4.6 (MIT license) is bundled with DupCaller, runs on the same Python, and is installed as the `PERF` command, so it needs no separate installation.

## Installation

### pip Installation

The tool uses pip for installing scripts and prerequisites. We recommend creating a new environment that also holds the [Prerequisites](#prerequisites) tools (bwa-mem2, GATK, samtools, htslib). If you do not have conda, install Miniforge first:

```bash
curl -fsSL -o Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/bin/activate"
```

Then create the environment:

```bash
conda create -y -n DupCaller -c conda-forge -c bioconda python=3.12 bwa-mem2=2.3 gatk4=4.3.0.0 samtools htslib
```

To install DupCaller, simply clone this repository and install via pip:

```bash
conda activate DupCaller
git clone https://github.com/AlexandrovLab/DupCaller.git
cd DupCaller
pip install .
```

### Docker / Singularity

A pre-built image is published by the Alexandrov Lab on the GitHub Container Registry at `ghcr.io/alexandrovlab/dupcaller:1.2.12` (public, no login needed). From the release after 1.2.12 on, each release's image is built from that release's tagged source by GitHub Actions (`.github/workflows/release.yml`), which runs the Nextflow mock test on it before publishing it and again on a fresh runner that pulls it anonymously.

**Pull and run with Docker:**

```bash
docker pull ghcr.io/alexandrovlab/dupcaller:1.2.12
docker run --rm ghcr.io/alexandrovlab/dupcaller:1.2.12 DupCaller.py --help
```

**Pull and run with Singularity:**

Pull the image (only needed once):

```bash
singularity pull dupcaller-1.2.12.sif docker://ghcr.io/alexandrovlab/dupcaller:1.2.12
```

For quick verification:

```bash
singularity exec dupcaller-1.2.12.sif DupCaller.py --help
```

For installation-free execution of DupCaller commands, run all DupCaller.py commands with `singularity exec` and binding of current directories:

```bash
singularity exec --bind $(pwd):$(pwd) dupcaller-1.2.12.sif DupCaller.py {your commands}
```

---

## Verify Installation

The repository ships three automated checks, the same ones GitHub Actions runs on every commit (`.github/workflows/test.yml`). Run them from the repository root. The two end-to-end tests run every step on a 3.2 kb synthetic genome with planted mutations and compare every deterministic output (VCFs, burdens, spectra, error profiles, stats) with the expected files in `tests/mock_pipeline/`, so a passing run reproduces the reference results exactly.

`-rs` lists skipped tests with the reason: an end-to-end test **skips** instead of failing when a tool it needs is not on `PATH`, so a skip means "not checked", not "passed".

**1. Unit tests** (DupCaller installed with pip, plus `tabix` from htslib):

```bash
conda activate DupCaller        # the environment from pip Installation
pip install ".[test]"           # adds pytest
pytest tests/unit/ -v -rs
```

**2. Mock end-to-end test, command line** (`index` → `trim` → `bwa-mem2 mem` → `gatk MarkDuplicates` → `call` → `estimate`, run directly; needs bwa-mem2, GATK, samtools and htslib, all in the environment from [pip Installation](#pip-installation)):

```bash
conda activate DupCaller
pytest tests/mock_pipeline/test_mock_pipeline.py -v -rs
```

**3. Mock end-to-end test, Nextflow** (the same pipeline through `nextflow/DupCaller.nf` in containers, in 11 configurations: defaults, each optional input, tumor-only, and the `udseq` and `metacs` protocols). Needs Nextflow and Docker or Singularity on `PATH`, plus pytest and pysam (both in the DupCaller environment):

```bash
conda activate DupCaller
conda install -y -c conda-forge -c bioconda nextflow=24.10.4   # if nextflow is not already on PATH
module load singularity                                         # HPC only; skip with Docker
pytest tests/mock_pipeline/test_mock_pipeline_nextflow.py -v -rs
```

It uses the DupCaller image pinned in `nextflow/nextflow.config`; set `DUPCALLER_IMAGE` to test another one (CI tests an image built from the checked-out commit this way). Without pytest, `tests/mock_pipeline/run_nextflow_pipeline.sh OUTDIR [docker|singularity] all` runs the same cases.

---

## Pipeline

### Protocol Settings

The barcode pattern (`trim -p`, Step 2) and how far from fragment and read ends mutations are ignored (`call -tt`/`-tr`, Step 5) depend on the library protocol:

| Protocol   | Nextflow `protocol` | `trim -p`                                                         | `call -tt` | `call -tr` |
| ---------- | --------------------- | ------------------------------------------------------------------- | ------------ | ------------ |
| RE-NanoSeq | `renanoseq`         | `NNNXXXX`                                                         | 7 (default)  | 7 (default)  |
| NanoSeq v2 | `nanoseqv2`         | `NNNXX`                                                           | 7 (default)  | 7 (default)  |
| UDSeq      | `udseq`             | `NNNNNNNN`                                                        | 15           | 8            |
| META-CS    | `metacs`            | `BXXXXXXXXXXXXXXXXXXX` with `-bl barcode_list.txt` (see Step 2) | 7 (default)  | 7 (default)  |

The Nextflow pipeline sets all three from the protocol name (`--protocol`, or `protocol` in `pipeline.config`; `metacs` also needs `--barcode_list`); see [Quick Start](#quick-start).

`-tt`/`--trimF` ignores mutations within N bp of either end of the sequenced fragment, and `-tr`/`--trimR` within N bp of either end of each read. For example, for UDSeq:

```bash
DupCaller.py trim -i ${SAMPLE}_R1.fq.gz -i2 ${SAMPLE}_R2.fq.gz -p NNNNNNNN -o ${SAMPLE}
...
DupCaller.py call -b ${SAMPLE}.mkdped.bam -n ${NORMAL}.mkdped.bam -f ${REF} -o ${SAMPLE} \
    -g ${GERMLINE} -m ${NOISE} -p ${THREADS} -tt 15 -tr 8
```

**Without a matched normal**, whatever the protocol, omit `-n` and add `-maf 0.1` to `call` so that variants above 10% allele fraction (mostly germline) are not called; see [Without Matched Normal](#without-matched-normal) (the Nextflow pipeline does this with `tumor_only`; see [Quick Start](#quick-start-without-a-matched-normal)). A `-maf` below 1 also masks positions with raw tumor depth below `1/maf` (10x for 0.1).

### Step 1: Index Reference Genome

DupCaller uses a numpyrized reference genome to perform memory-efficient reference fetching, trinucleotide context lookup, and repeat (homopolymer/short tandem repeat, STR) annotation used to improve indel calling near repetitive regions.

STR annotation comes from a single tsv produced by [PERF](https://github.com/rkmlab/perf). First, run PERF against the reference FASTA:

```bash
PERF -m 2 -M 10 -u 2 -i reference.fa -o repeats.tsv
```

- `-m`/`-M` — min/max repeat unit (motif) size to search for. `-m 2 -M 10` covers 2–10bp-unit STRs; going much higher increases PERF's runtime and memory non-trivially. Homopolymers (unit length 1) are not taken from this tsv: `index` derives them directly from the reference sequence and ignores any unit-length-1 rows, so `-m 1` only adds runtime.
- `-u 2` (`--min-units`) — minimum number of repeat copies to report.

Then index the reference, passing PERF's tsv directly:

```bash
DupCaller.py index -f reference.fa -rt repeats.tsv
```

The command will generate five h5 files in the same folder as the reference: `{reference}.ref.h5`, `{reference}.tn.h5`, `{reference}.hp.h5`, `{reference}.str.h5`, and `{reference}.dbs.h5` — numpyrized reference sequences, trinucleotide contexts, homopolymer annotations, STR (unit length >=2) annotations, and per-position dinucleotide (DBS) classes, respectively.

For human reference genome hg38 and mouse reference genome mm39, we provided pre-built indexes and resource files in the [Resources](#resources).

The reference also needs a samtools FASTA index and a bwa-mem2 index (for Step 3), built once:

```bash
samtools faidx reference.fa
bwa-mem2 index reference.fa   # writes reference.fa.{0123,amb,ann,bwt.2bit.64,pac}; ~90 GB RAM for GRCh38 (28 bytes per base)
```

#### Parameters

| Short | Long        | Description                                                                                                                                                                                                                                               |
| ----- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| -f    | --reference | Reference genome fasta file (required)                                                                                                                                                                                                                    |
| -rt   | --repeatTsv | PERF-format repeat tsv (chrom, start, end, motif, length, strand, num_units, motif_repeat) for the reference (required). Repeat unit length and repeat count are read directly from the motif and num_units columns; rows with unit length 1 are ignored. |

### Step 2: Trim Barcodes

`DupCaller.py trim` extracts 5-prime barcodes from paired-end FASTQs:

```bash
DupCaller.py trim -i read1.fq -i2 read2.fq -p barcode_pattern -o sample_name
```

- `read1.fq` / `read2.fq` — FASTQ files from read 1 and read 2. Both unzipped and gzip-compressed files are supported.
- `barcode_pattern` — pattern of barcodes starting from the 5-prime end, with `N` representing a barcode base and `X` representing a skipped base (similar notation to [UMI-tools](https://github.com/CGATOxford/UMI-tools)). For example, NanoSeq uses a 3-base barcode followed by 4 constant bases, so the pattern should be `NNNXXXX`.
- Barcodes from a fixed list, of varying length (e.g. Tn5-based protocols such as META-CS, where each read starts with one of a set of 11–13 nt transposon barcodes followed by the 19 nt ME): put `B` in the pattern for "one barcode from the list" and pass the list with `-bl/--barcode-list` (one barcode per line). Positions after `B` count from the end of the matched barcode, e.g. `-p BXXXXXXXXXXXXXXXXXXX` takes the barcode and skips the 19 nt ME. Each read's barcode may have up to `-mm/--max-mismatch` mismatches (default 1); the barcode with the fewest mismatches wins, then the shortest. Up to 1 mismatch is matched with lookup tables (fast, any list size); larger limits walk all listed barcodes at once, dropping each as soon as it exceeds the limit. The listed (corrected) sequence is written to `DB`. Read pairs whose read 1 or read 2 matches no barcode, or matches two equally well, are dropped and counted in the summary `trim` prints.
- `sample_name` — prefix of output paired FASTQs. After the run completes, `{sample_name}_1.fastq` and `{sample_name}_2.fastq` will be generated. Barcodes are appended to each read name (`{original_read_name}_{read1_barcode}+{read2_barcode}`, for optical-duplicate read-name matching only) and also written as a `DB:Z:{read1_barcode}-{read2_barcode}` FASTQ comment; the comment is promoted to a BAM tag by `bwa-mem2 mem -C` (Step 3) and is what `DupCaller.py call` actually reads to assign duplex family membership — it never parses the read name for barcodes.

If the matched normal is prepared in the same way as the sample, apply trimming with the same scheme to the matched normal FASTQs. For traditional bulk normal, trimming is not needed.

---

### Step 3: Align Reads

Align the trimmed FASTQs of both sample and matched normal with bwa-mem2 (`bwa-mem2 index` from Step 1 must already exist). GATK requires read group fields `ID`, `SM`, and `PL`, so adding those tags during alignment is recommended. **FASTQ tags must be kept — for `bwa-mem2 mem` (as for `bwa mem`) this requires the `-C` option.** The Nextflow pipeline and the mock-pipeline regression test also use `-T 0`, which outputs every alignment regardless of its score, so low-scoring alignments are filtered by DupCaller (MAPQ/AS-XS) rather than dropped by the aligner.

```bash
bwa-mem2 mem -C -T 0 -t {threads} -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" \
    reference.fa {sample_name}_1.fastq {sample_name}_2.fastq \
    | samtools sort -@ {threads} > {sample_name}.bam
samtools index -@ {threads} {sample_name}.bam
```

- `threads` — number of cores used for aligning
- `reference.fa` — reference genome FASTA file, with its bwa-mem2 index next to it
- `{sample_name}_1.fastq` / `{sample_name}_2.fastq` — trimmed FASTQ files from Step 2

---

### Step 4: Mark Duplicates

Run GATK MarkDuplicates on sample and matched-normal BAMs. Optical and PCR duplicates must be treated differently in ecNGS variant calling:

- Set `--TAGGING_POLICY OpticalOnly` to differentiate optical from PCR duplicates.
- Set `--DUPLEX_UMI true` for duplex UMI handling.
- Set `--READ_NAME_REGEX` as shown below, because the read names were modified in Step 2.

> **Note:** The published analyses used GATK 4.3.0.0. Older GATK releases do not have the `--DUPLEX_UMI` flag.

```bash
gatk MarkDuplicates \
    -I sample.bam -O sample.mkdped.bam -M sample.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
samtools index sample.mkdped.bam
```

`DupCaller.py call` needs the BAM index (`.bai`) of every input BAM.

---

### Step 5: Call Variants

After preprocessing, run `DupCaller.py call` to call somatic mutations. Usage depends on your experimental design:

**Whole genome / whole exome / reduced genome (e.g. NanoSeq) with a matched normal:**

```bash
DupCaller.py call -b ${sample}.mkdped.bam -f reference.fa -o {output_prefix} \
    -p {threads} -n {normal.bam} -g germline.vcf.gz -m noise_mask.bed.gz
```

**Mutagenesis panel without a matched normal:**

```bash
DupCaller.py call -b ${sample}.mkdped.bam -f reference.fa -o {output_prefix} \
    -p {threads} -g germline.vcf.gz -m noise_mask.bed.gz -maf 0.1
```

> Note: DupCaller partitions jobs by genomic region; multithreading is less effective for small targeted panels. In this case, use at most one thread per distinct targeted region.

**Non-default barcode tags:** BAMs produced by DupCaller's own `trim` (Step 2) carry a shared `DB:Z:{bc1}-{bc2}` tag per read pair, and this is what `call` reads by default. If instead you are calling directly from a BAM that stores the duplex barcode differently, pass `-bc`/`--barcode TAG,NORMALIZE,SEP`: `TAG` is the bam tag holding the barcode, and `NORMALIZE` says whether it still needs splitting/reordering. For example, the Sanger NanoSeq pipeline tags each read with a single `RB` tag that already identifies the duplex family (e.g. `chrom,fragment_start,fragment_end,bc1,bc2`) independent of read orientation — pass `-bc RB,0,-` to read it directly (`SEP` is ignored in this case). DupCaller's own `DB:Z:{bc1}-{bc2}` tag is the `NORMALIZE=1` case (the default, `-bc DB,1,-`): it's an unnormalized `bc1<SEP>bc2` pair that still needs read-orientation-based reordering to group both strand-copies of a molecule together.

**Input validation:** DupCaller checks for the existence of all required files (BAM, reference, h5 index files, and optional files) before starting analysis, providing clear error messages for missing files.

**Multi-threading:** Coverage files from different threads are automatically merged post-processing with region-aware boundary detection and tabix indexing.

See the [Results](#results) section for descriptions of all output files.

#### Parameters

##### Required

| Short | Long        | Description                                                                                                               |
| ----- | ----------- | ------------------------------------------------------------------------------------------------------------------------- |
| -b    | --bam       | BAM file of ecNGS data                                                                                                    |
| -f    | --reference | Reference genome FASTA file                                                                                               |
| -o    | --output    | Output directory; created if needed. Files inside are named`{sample}_*`, where `{sample}` is the directory's basename |

##### Recommended

These options should be understood and customized accordingly.

| Short | Long             | Description                                                                                                                          | Default         |
| ----- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------ | --------------- |
| -r    | --regions        | Contigs to consider for variant calling. For non-human species, set accordingly (e.g. for mouse:`-r chr{1..19} chrX`)              | chr{1..22} chrX |
| -g    | --germline       | Indexed germline VCF with AF field                                                                                                   | None            |
| -p    | --threads        | Number of threads                                                                                                                    | 1               |
|       | --minChunkLength | Minimum chunk length in bases; short contigs stay unsplit. May use fewer workers than requested. Also accepts`--min-chunk-length`. | 10000           |
| -n    | --normalBams     | BAM file(s) of matched normals. When unavailable, set`-maf` to an appropriate value (e.g. 0.1)                                     | None            |
| -m    | --noise          | BED interval file(s) masking noisy positions                                                                                         | None            |
| -R    | --regionfile     | Inclusive BED file specifying target regions                                                                                         | None            |
| -maf  | --maxAF          | Maximum allele fraction to call a somatic mutation. Must be set when matched normal (`-n`) is unavailable                          | 1               |
| -tt   | --trimF          | Ignore mutations less than n bp from template ends                                                                                   | 7               |
| -tr   | --trimR          | Ignore mutations less than n bp from read ends                                                                                       | 7               |
| -bc   | --barcode        | Molecular/duplex barcode read tag, as`TAG,NORMALIZE,SEP` (see note above)                                                          | DB,1,-          |

##### Optional

The effect of changing these parameters should be evaluated before implementation.

| Short    | Long                    | Description                                                                                                                                                                                                                                                                                                                                                            | Default |
| -------- | ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- |
|          | --naf                   | Maximum VAF in matched normal for a mutation to be called                                                                                                                                                                                                                                                                                                              | 0.01    |
| -res     | --rescue                | Output discarded variants whose LR clears their channel's final threshold, with the reason in the FILTER field (called in the coverage pass; not with `--skipCoveragePass`)                                                                                                                                                                                          | False   |
| -nm      | --nmflt                 | Drop a read family when at least half of its reads have an edit distance (NM) above this value                                                                                                                                                                                                                                                                         | 5       |
| -ax      | --minMeanASXS           | Minimum mean AS-XS alignment score difference for a read group to be considered                                                                                                                                                                                                                                                                                        | 50      |
| -gaf     | --germlineAfCutoff      | Skip positions with germline AF above this threshold                                                                                                                                                                                                                                                                                                                   | 0.001   |
| -d       | --minNdepth             | Minimum coverage in normal for called variants                                                                                                                                                                                                                                                                                                                         | 10      |
| -sc      | --skipCoveragePass      | Skip round 2 (the coverage-only pass), roughly halving total runtime: every output except the per-locus `coverage.bed.gz` is still written, with coverage from round 1's power tallies contracted at the final thresholds. `estimate` then skips the duplex allele counts and per-gene coverage, and `-rb/--reestimatebed` is unavailable. Not with `--rescue` | False   |
| -nms     | --doNotOutputMaskedSite | Skip snp_mask/noise_mask (`-m`) candidates entirely: no depth extraction for them and no `masked` records in the fail VCFs. Calls, mutation rates and corrected burdens are unchanged; `Unmasked Coverage` (`_stats.txt`) and `Unmasked burden` (`_sbs_burden.txt`) are written as N/A                                                                     | False   |

##### Advanced

These are variant calling model parameters; adjustment is unnecessary for general use.

| Short        | Long             | Description                                                                                                                                                                                                                                                                                                                                                                                                              | Default              |
| ------------ | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------- |
| -E           | --errprefix      | Prefix for all six error files ({prefix}.amp.tn.srd.txt, {prefix}.amp.hp.txt, {prefix}.amp.str.txt, {prefix}.dmg.tn.txt, {prefix}.dmg.hp.txt, {prefix}.dmg.str.txt); overrides the default (output prefix)                                                                                                                                                                                                               | None                 |
| -lo          | --learnOnly      | Stop after estimating/writing the error-rate files (ERROR/ dir); skip variant calling entirely                                                                                                                                                                                                                                                                                                                           | False                |
| -lfdr        | --lfdrThreshold  | Maximum local false discovery rate (lFDR) of a PASS call. A call's lFDR is `(1-mu)/(LR*mu+1-mu)`, where `mu` is the estimated mutation rate of its mutation channel. A call passes when its lFDR is at most this value and its log10 LR is at least `--minLR`, i.e. each channel's log10 LR threshold is `max(minLR, log10((1-lfdr)(1-mu)/(lfdr*mu)))`                                                                               | 0.05                 |
| -a           | --pseudocount    | Pseudocount (> 0) of the SBS amplification-error (SRD) fit and of the error-matrix smoothing. The per-channel mutation-rate solve uses none                                                                                                                                                                                                                                                                              | 0.5                  |
| -mlr         | --minLR          | Minimum log10 LR for a PASS call (applied together with `-lfdr`): floor on every channel's log10 LR threshold, used for calling and for the detection-power simulation behind the burden's sensitivity correction                                                                                                                                                                                                                                  | 5                    |
| -mq          | --mapq           | Minimum MAPQ for an alignment to be considered                                                                                                                                                                                                                                                                                                                                                                           | 40                   |
| -w           | --windowSize     | Genomic window size for coverage calculation and BAM partitioning                                                                                                                                                                                                                                                                                                                                                        | 100000               |
| -bq          | --minBq          | Bases with quality below this value are zeroed out and excluded from variant calling/learning (usable iff BQ >= minBq)                                                                                                                                                                                                                                                                                                   | 18                   |
| -aq          | --minAltQual     | Minimum summed per-strand consensus base quality at a position for SBS damage-rate learning to use it                                                                                                                                                                                                                                                                                                                    | 90                   |
|              | --minRef         | Minimum per-position read depth (ref+alt combined) for damage-rate learning to use it; the stricter of`--minRef`/`--minAlt` is applied, since the underlying check doesn't separate ref vs alt counts                                                                                                                                                                                                                | 3                    |
|              | --minAlt         | Minimum per-position read depth (ref+alt combined) for damage-rate learning to use it; the stricter of`--minRef`/`--minAlt` is applied, since the underlying check doesn't separate ref vs alt counts                                                                                                                                                                                                                | 3                    |
|              | --srdMinRead     | Minimum reads on a single strand (F1R2 or F2R1) for that strand to be included in SBS single-read discrepancy (SRD) rate learning; evaluated independently per strand, and also gates the minimum BQ-qualifying bases required for a site to be counted                                                                                                                                                                       | 3                    |
|              | --ssmMinRead     | Minimum reads required on EACH strand (F1R2 and F2R1 both) for a duplex family to be considered for single-strand mismatch (SSM) damage-rate learning                                                                                                                                                                                                                                                                 | 3                    |
| -pt          | --p_threshold    | Strand independence filter threshold for SBS and indel calls: a PASS call fails with FILTER `strand_independence` if either strand's p-value (INFO `MSP`) is <= this. Per strand, not corrected for the number of calls. Must be in [0, 1); 0 disables the filter (MSP is still reported)                                                                                                                              | 0.05                 |
| -id          | --indelbed       | Indel enhanced Panel of Normals (ePoN) for indel calling                                                                                                                                                                                                                                                                                                                                                                 | None                 |
| -rt          | --regionst       | Contigs to consider for error-profile training, if different from`-r`/`--regions`                                                                                                                                                                                                                                                                                                                                    | same as`--regions` |
| -pd          | --maxPileupDepth | Maximum depth for samtools mpileup                                                                                                                                                                                                                                                                                                                                                                                       | 1000000              |
| -mr          | --muterateprefix | Reuse per-channel mutation-rate tables from an earlier run instead of estimating them from this run's own candidates (useful when calling a small region). Pass`{out}/tmp/{sample}` from that run, which holds `{sample}_sbs96_rate_n1.txt` and `{sample}_indel_rate_by_hp_str.txt` (tables written before the per-run indel Eeff change are refused). LR thresholds are still solved against this run's `-lfdr` | None                 |
|              | --seed           | RNG seed for the Monte Carlo detection-power simulation, for reproducible results across runs and`-p` values. If unset, a random seed is generated and recorded in `{sample}_call_params.log`                                                                                                                                                                                                                        | random               |

#### Germline and Noise Masks

The `-m` option accepts BED files and will ignore any mutations overlapping an excluded locus. Any custom BED file can be used as input.

---

### Step 6: Estimate Mutational Burden

After mutation calling, run burden estimation on the `call` output directory (`-i` is the same path that was passed to `call -o`):

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX
```

Adjust the `-r` regions according to the reference genome used.

#### Per-Gene Coverage

For dNdScv coverage correction, DupCaller can output mean duplex depth per gene using the `-gb` option:

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX -gb {target}.bed.gz
```

The gene BED is a tab-separated BED with one row per exon (or target interval). DupCaller reads only the first four columns; any further columns (score, strand, ...) are ignored:

| Column | Content                               | Example     |
| ------ | ------------------------------------- | ----------- |
| 1      | Chromosome, named as in the reference | `chr17`   |
| 2      | Start (0-based)                       | `7669608` |
| 3      | End (exclusive)                       | `7669690` |
| 4      | `{gene_name}_{exon_number}`         | `TP53_12` |

```
chr17	7669608	7669690	TP53_12	0	-
chr17	7670608	7670715	TP53_11	0	-
chr17	7673218	7673266	TP53_10	0	-
```

The gene name is the part of column 4 before the first `_`, and all rows with the same gene name are pooled, so a gene name must not itself contain `_`. Each gene's depth is averaged over the summed length of its rows. The file must be sorted by position within each chromosome, bgzip-compressed and tabix-indexed (`bgzip {target}.bed && tabix -p bed {target}.bed.gz`).

For hg38 exomes, a ready-to-use gene BED of the IDT xGen Exome Hyb Panel v2 targets (197,769 targets, 4th column `{gene}_{exon}`) is included in this repository: `resources/xgen-exome-hyb-panel-v2-targets-hg38.bed.gz` (with its `.tbi`):

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX \
    -gb resources/xgen-exome-hyb-panel-v2-targets-hg38.bed.gz
```

In the Nextflow pipeline, pass the same file as `--gene_bed`.

#### Re-estimation for Specific Regions

To re-estimate trinucleotide-corrected mutational burden in specific regions without re-running variant calling, use `-rb`. This is a separate mode: it writes only the `_re_estimate` files (see [Optional / Conditional Files](#optional--conditional-files)) and the parameter log, not the main burden outputs, so run it after a normal `estimate`. It needs the per-locus coverage, so it is unavailable when `call` was run with `-sc/--skipCoveragePass`:

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX -rb {re_estimate}.bed.gz
```

As with `-gb`, the BED file must be bgzip-compressed and tabix-indexed.

#### Parameters

##### Required

| Short | Long        | Description                                                             |
| ----- | ----------- | ----------------------------------------------------------------------- |
| -i    | --prefix    | Output directory of the`call` command (the `-o` value)              |
| -f    | --reference | FASTA file of reference genome (its h5 index files must sit next to it) |

##### Optional

| Short | Long            | Description                                                                                                                                                                                                                                                                                               | Default         |
| ----- | --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------- |
| -r    | --regions       | Contigs to consider for trinucleotide calculation                                                                                                                                                                                                                                                         | chr{1..22} chrX |
| -ot   | --outTrinuc     | Write the reference trinucleotide composition to this file, for reuse with`-ft`                                                                                                                                                                                                                         | None            |
| -ft   | --refTrinuc     | Load a trinucleotide composition written by`-ot` instead of rescanning the reference (`-f` is still required)                                                                                                                                                                                         | None            |
| -d    | --dilute        | Set when sample and matched normal come from the same starting DNA material: SNVs with a tumor alt allele count (`AC`) above 1 are dropped when their tumor and normal allele counts differ significantly (Barnard's exact test, p <= 0.05), and kept calls are written to `SBS/{sample}_sbs_flt.vcf` | False           |
| -gb   | --genebed       | bgzipped, tabix-indexed gene BED file for per-gene coverage calculation, 4th column `{gene}_{exon}` (format in [Per-Gene Coverage](#per-gene-coverage)). For hg38 exomes, use the included `resources/xgen-exome-hyb-panel-v2-targets-hg38.bed.gz`                                                     | None            |
| -rb   | --reestimatebed | bgzipped, tabix-indexed BED file for burden re-estimation in specific regions. Writes only the `_re_estimate` outputs (see [Re-estimation for Specific Regions](#re-estimation-for-specific-regions))                                                                                                     | None            |

---

### Step 7: Summarize Across Samples

After running `estimate` on all samples, use `DupCaller.py summarize` to collate burden metrics and SBS96 profiles into multi-sample tables:

```bash
DupCaller.py summarize -i sample1 sample2 sample3 -o cohort_summary.txt
```

This reads `{sample}/{sample}_stats.txt`, `{sample}/SBS/{sample}_sbs_burden.txt`, `{sample}/INDEL/{sample}_indel_burden.txt`, `{sample}/DBS/{sample}_dbs_burden.txt`, and `{sample}/SBS/{sample}_sbs_96_corrected.txt` from each sample folder (all five are required) and writes four output files:

| File                                     | Description                                                                       |
| ---------------------------------------- | --------------------------------------------------------------------------------- |
| `cohort_summary.txt`                   | One row per sample with SBS, indel, and DBS burden metrics and library statistics |
| `cohort_summary_SBS96_uncorrected.txt` | 96-context raw mutation counts across all samples (SigProfiler-compatible format) |
| `cohort_summary_SBS96_corrected.txt`   | 96-context trinucleotide-corrected counts across all samples                      |
| `cohort_summary_SBS96_genome.txt`      | 96-context estimated mutations per genome across all samples                      |

All SBS96 file can be directly input into [SigProfilerPlotting](https://github.com/SigProfilerSuite/SigProfilerPlotting)

#### Parameters

| Short | Long     | Description                                                        |
| ----- | -------- | ------------------------------------------------------------------ |
| -i    | --input  | One or more sample folders (space-separated)                       |
| -o    | --output | Output filename for the summary table (e.g.`cohort_summary.txt`) |

### Optional: Aggregate Error Profiles Across Samples

`DupCaller.py aggregate` pools the learned error profiles of several `call` runs into one set, which can then be passed to `call -E` (for example, for a low-coverage sample from the same library prep):

```bash
DupCaller.py aggregate -i sample1 sample2 sample3 -o pooled
DupCaller.py call ... -E pooled
```

It sums each sample's `ERROR/{sample}.*.txt` count tables and re-fits the SRD amplification matrix from the summed base-quality histograms in each sample's `tmp/{sample}.amp.tn.bqhist.npz`, so the `tmp/` folder of each input run must still exist. It writes `pooled.amp.tn.txt`, `pooled.amp.tn.srd.txt`, `pooled.dmg.tn.txt`, and the four `pooled.{amp,dmg}.{hp,str}.txt` files.

| Short | Long          | Description                                                                                                       | Default |
| ----- | ------------- | ----------------------------------------------------------------------------------------------------------------- | ------- |
| -i    | --input       | One or more`call` output directories                                                                            | None    |
| -f    | --input-file  | File with one error-file prefix per line (e.g.`sample1/ERROR/sample1`), as an alternative or addition to `-i` | None    |
| -o    | --output      | Output prefix for the aggregated error files (required)                                                           | None    |
| -a    | --pseudocount | Pseudocount for the SRD matrix re-fit (same as`call -a`)                                                        | 0.5     |

---

## Results

For detailed column-by-column and field-by-field descriptions of every output file, see the [`docs/`](docs/) folder:

- [`docs/call_outputs.md`](docs/call_outputs.md) — all files from `DupCaller.py call`
- [`docs/estimate_outputs.md`](docs/estimate_outputs.md) — all files from `DupCaller.py estimate`
- [`docs/summarize_outputs.md`](docs/summarize_outputs.md) — the cohort summary table from `DupCaller.py summarize`

SBS/indel/DBS-specific outputs (VCFs, corrected-context tables, burden files, and signature plots) are written into per-type `SBS/`, `INDEL/`, and `DBS/` subfolders of the sample's output directory.

### Core Output Files

| File                                                      | Description                                                                                                                                         |
| --------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SBS/{sample}_sbs.vcf`                                  | VCF of detected SNVs (adjacent SNVs from the same read family are reported as DBS instead)                                                          |
| `INDEL/{sample}_indel.vcf`                              | VCF of detected short indel mutations                                                                                                               |
| `DBS/{sample}_dbs.vcf`                                  | VCF of detected dinucleotide substitution (DBS) mutations                                                                                           |
| `{sample}_coverage.bed.gz`                              | Duplex coverage depths across genomic positions, merged across threads and tabix-indexed. Not written with `-sc/--skipCoveragePass`                  |
| `{sample}_coverage.bed.gz.tbi`                          | Tabix index for the coverage BED file                                                                                                               |
| `SBS/{sample}_trinuc_by_duplex_group.txt`               | Trinucleotide context counts grouped by duplex read number, used for SBS burden estimation                                                          |
| `INDEL/{sample}_indel_by_duplex_group.txt`              | Indel context counts grouped by duplex read number, used for indel burden estimation                                                                |
| `DBS/{sample}_dbs_by_duplex_group.txt`                  | 144-class dinucleotide context counts grouped by duplex read number, used for DBS burden estimation                                                 |
| `{sample}_duplex_family_strand_composition.txt`         | Strand composition statistics for duplex read families                                                                                              |
| `{sample}_duplex_family_strand_composition_heatmap.pdf` | Heatmap visualization of duplex family strand composition                                                                                           |
| `{sample}_call_params.log`                              | Full record of all resolved parameters used for the`call` run                                                                                     |
| `{sample}_stats.txt`                                    | Overall sequencing and analysis metrics                                                                                                             |

### Error Profile Files

Written to the `ERROR/` subfolder of the `call` output directory. If a complete set of the six files used for calling already exists there, learning is skipped and they are reused.

| File                              | Description                                                                       |
| --------------------------------- | --------------------------------------------------------------------------------- |
| `ERROR/{sample}.amp.tn.txt`     | Amplification SBS mismatch counts by trinucleotide context (diagnostic only)      |
| `ERROR/{sample}.amp.tn.srd.txt` | SRD-EM-fitted amplification SBS error-rate matrix used for calling                |
| `ERROR/{sample}.dmg.tn.txt`     | Damage SBS error counts by trinucleotide context                                  |
| `ERROR/{sample}.amp.hp.txt`     | Amplification indel error counts for homopolymer contexts                         |
| `ERROR/{sample}.amp.str.txt`    | Amplification indel error counts for short-tandem-repeat contexts                 |
| `ERROR/{sample}.dmg.hp.txt`     | Damage indel error counts for homopolymer contexts                                |
| `ERROR/{sample}.dmg.str.txt`    | Damage indel error counts for short-tandem-repeat contexts                        |

### Burden Estimation Files

| File                                              | Description                                                                                                                                               |
| ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SBS/{sample}_sbs_burden.txt`                   | SBS burden with uncorrected and corrected estimates and 95% confidence intervals                                                                          |
| `INDEL/{sample}_indel_burden.txt`               | Indel burden with uncorrected and corrected estimates and 95% confidence intervals                                                                        |
| `DBS/{sample}_dbs_burden.txt`                   | DBS burden with uncorrected and corrected estimates and 95% confidence intervals                                                                          |
| `SBS/{sample}_sbs_96_corrected.txt`             | Corrected SBS counts across 96 trinucleotide contexts for signature analysis, including a`mutations_per_opportunity` column                             |
| `INDEL/{sample}_indel_83_corrected.txt`         | Corrected indel counts across the 83 ID contexts, including a`mutations_per_opportunity` column                                                         |
| `DBS/{sample}_dbs_78_corrected.txt`             | Corrected DBS counts across the 78 DBS contexts, including a`mutations_per_opportunity` column                                                          |
| `SBS/{sample}_sbs_burden_by_group_size.txt`     | SBS burden stratified by duplex group size, both cumulative ("min group size >= N") and exact ("group size == N"), corrected and uncorrected, with 95% CI |
| `INDEL/{sample}_indel_burden_by_group_size.txt` | Same stratification as above, for indels                                                                                                                  |
| `DBS/{sample}_dbs_burden_by_group_size.txt`     | Same stratification as above, for DBS                                                                                                                     |
| `{sample}_duplex_allele_counts.txt`             | Duplex depths and allele counts for each PASS SNV and indel. Not written when `call` used `-sc/--skipCoveragePass`                                        |
| `{sample}_estimate_params.log`                  | Full record of all resolved parameters used for the`estimate` run. Written next to the sample folder (`{-i}_estimate_params.log`), not inside it            |

`Mutation number per genome` (and its 95% bounds) in each `_burden.txt` is the corrected burden multiplied by `Reference base number`, the number of reference bases considered (one copy of each region in `-r`). It is therefore the expected number of mutations **per haploid genome**; double it for a diploid genome. The same applies to the per-channel `mutation_number_genome` column of `_sbs_96_corrected.txt` and `_indel_83_corrected.txt` and to `summarize`'s `*_mutations_per_genome` columns.

### Visualization Files

| File                                              | Description                                                                  |
| ------------------------------------------------- | ---------------------------------------------------------------------------- |
| `SBS/SBS_96_plots_{sample}.pdf`                 | 96-context SBS mutational signature plots (uncorrected and corrected counts) |
| `INDEL/ID_83_plots_{sample}.pdf`                | 83-context indel (ID83) mutational signature plots                           |
| `DBS/DBS_78_plots_{sample}.pdf`                 | 78-context DBS mutational signature plots                                    |
| `SBS/{sample}_sbs_burden_by_group_size.pdf`     | SBS burden across cumulative and exact duplex group sizes                    |
| `INDEL/{sample}_indel_burden_by_group_size.pdf` | Same, for indels                                                             |
| `DBS/{sample}_dbs_burden_by_group_size.pdf`     | Same, for DBS                                                                |

### Optional / Conditional Files

| File                                              | Condition      | Description                                                   |
| ------------------------------------------------- | -------------- | ------------------------------------------------------------- |
| `{sample}_gene_coverage.txt`                    | `-gb` option | Mean duplex coverage per gene for dNdScv correction           |
| `SBS/{sample}_sbs_burden_re_estimate.txt`       | `-rb` option | Re-estimated SBS burden for specific regions                  |
| `INDEL/{sample}_indel_burden_re_estimate.txt`   | `-rb` option | Re-estimated indel burden for specific regions                |
| `SBS/{sample}_sbs_96_corrected_re_estimate.txt` | `-rb` option | Re-estimated 96-context SBS counts for specific regions       |
| `SBS/SBS_96_plots_{sample}_re_estimate.pdf`     | `-rb` option | Signature plots for re-estimated regions                      |
| `SBS/{sample}_sbs_flt.vcf`                      | `-d` option  | SNV calls kept after the`--dilute` normal-comparison filter |

---

## End-to-End Examples

The following examples assume the reference genome has already been indexed (`DupCaller.py index`, `samtools faidx`, `bwa-mem2 index`; Step 1) and that the germline VCF and noise mask are available. Replace all placeholder paths with your actual files.

### With Matched Normal

```bash
SAMPLE=sample1
NORMAL=normal1
REF=/path/to/hg38.fa
GERMLINE=/path/to/gnomad.hg38.vcf.gz
NOISE=/path/to/noise_mask.bed.gz
THREADS=16

# 1. Trim barcodes
DupCaller.py trim -i ${SAMPLE}_R1.fq.gz -i2 ${SAMPLE}_R2.fq.gz -p NNNXXXX -o ${SAMPLE}
DupCaller.py trim -i ${NORMAL}_R1.fq.gz -i2 ${NORMAL}_R2.fq.gz -p NNNXXXX -o ${NORMAL}

# 2. Align
bwa-mem2 mem -C -T 0 -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    ${REF} ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq | samtools sort -@ ${THREADS} > ${SAMPLE}.bam
samtools index -@ ${THREADS} ${SAMPLE}.bam

bwa-mem2 mem -C -T 0 -t ${THREADS} -R "@RG\tID:${NORMAL}\tSM:${NORMAL}\tPL:ILLUMINA" \
    ${REF} ${NORMAL}_1.fastq ${NORMAL}_2.fastq | samtools sort -@ ${THREADS} > ${NORMAL}.bam
samtools index -@ ${THREADS} ${NORMAL}.bam

# 3. Mark duplicates
gatk MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.mkdped.bam -M ${SAMPLE}.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
samtools index ${SAMPLE}.mkdped.bam

gatk MarkDuplicates -I ${NORMAL}.bam -O ${NORMAL}.mkdped.bam -M ${NORMAL}.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
samtools index ${NORMAL}.mkdped.bam

# 4. Call variants
DupCaller.py call -b ${SAMPLE}.mkdped.bam -f ${REF} -o ${SAMPLE} \
    -n ${NORMAL}.mkdped.bam -g ${GERMLINE} -m ${NOISE} -p ${THREADS}

# 5. Estimate mutational burden
DupCaller.py estimate -i ${SAMPLE} -f ${REF} -r chr{1..22} chrX
```

### Without Matched Normal

When no matched normal is available, omit `-n` and set `-maf` to cap the maximum allele frequency of called variants. A value of 0.1 is appropriate for most somatic applications to exclude common germline variants. The Nextflow pipeline runs the same steps with `--tumor_only true`; see [Quick Start](#quick-start-without-a-matched-normal).

```bash
SAMPLE=sample1
REF=/path/to/hg38.fa
GERMLINE=/path/to/gnomad.hg38.vcf.gz
NOISE=/path/to/noise_mask.bed.gz
THREADS=16

# 1. Trim barcodes
DupCaller.py trim -i ${SAMPLE}_R1.fq.gz -i2 ${SAMPLE}_R2.fq.gz -p NNNXXXX -o ${SAMPLE}

# 2. Align
bwa-mem2 mem -C -T 0 -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    ${REF} ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq | samtools sort -@ ${THREADS} > ${SAMPLE}.bam
samtools index -@ ${THREADS} ${SAMPLE}.bam

# 3. Mark duplicates
gatk MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.mkdped.bam -M ${SAMPLE}.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
samtools index ${SAMPLE}.mkdped.bam

# 4. Call variants (no matched normal; filter by maximum allele frequency)
DupCaller.py call -b ${SAMPLE}.mkdped.bam -f ${REF} -o ${SAMPLE} \
    -g ${GERMLINE} -m ${NOISE} -maf 0.1 -p ${THREADS}

# 5. Estimate mutational burden
DupCaller.py estimate -i ${SAMPLE} -f ${REF} -r chr{1..22} chrX
```

---

## Resources

Pre-built reference indexes, germline VCFs, and noise masks are available for download:

### Human (GRCh38/hg38)

| Resource         | Link/Source                                                                                                                                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Reference genome | [TCGA hg38 reference file](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834)                                                                                      |
| Reference index  | [Pre-built hg38 DupCaller reference](https://drive.google.com/drive/folders/1v8ut5aky01zujXlZ0b-dhCqImWMyMbiK?usp=drive_link)                                                         |
| Germline VCF     | af-only-gnomad.hg38.vcf.gz from the GATK somatic resource bundle ([download](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz), [index](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi)) |
| Noise mask       | NanoSeq `SNP_GRCh38.wgns.bed.gz` and `NOISE_GRCh38.wgns.bed.gz` masks from [NanoSeq](https://github.com/cancerit/NanoSeq); pass both to `call -m` (`--premade hg38` does)               |
| Gene BED (exome) | `resources/xgen-exome-hyb-panel-v2-targets-hg38.bed.gz` in this repository: IDT xGen Exome Hyb Panel v2 targets, for `estimate -gb` (see [Per-Gene Coverage](#per-gene-coverage)) |

### Mouse (GRCm39/mm39)

| Resource         | Link                                                                                                                                                                                       |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Reference genome | [UCSC mm39 reference file](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz)                                                                                              |
| Reference index  | [Pre-built mm39 DupCaller reference](https://drive.google.com/drive/folders/1apvqjyO3OhApQ5cnmOe1LWR2wi7_cMqK?usp=drive_link)                                                               |
| Germline VCF     | [mgp_strains](https://drive.google.com/drive/folders/1apvqjyO3OhApQ5cnmOe1LWR2wi7_cMqK?usp=drive_link) — includes VCF for SNPs in popular mouse strains. Use the vcf file for your strain. |
| Noise mask       | [NOISE.mm39.bed.gz](https://drive.google.com/drive/folders/1apvqjyO3OhApQ5cnmOe1LWR2wi7_cMqK?usp=drive_link), in-house noise mask from mouse duplex sequencing data                         |

## Citation

DupCaller: Cheng, Y. et al. Improved Mutation Detection in Duplex Sequencing Data with Sample-Specific Error Profiles. bioRxiv (2025). https://doi.org/10.1101/2025.07.13.664565

PERF: Avvaru, A. K., Sowpati, D. T. & Mishra, R. K. PERF: an exhaustive algorithm for ultra-fast and efficient identification of microsatellites from large DNA sequences. Bioinformatics 34, 943-948 (2018). https://doi.org/10.1093/bioinformatics/btx721

---

## Copyright

BSD 2-Clause License

Copyright (c) 2024, Alexandrov Lab

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## Contact

Yuhe Cheng (yuc211@ucsd.edu)
