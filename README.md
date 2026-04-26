# DupCaller

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/AlexandrovLab/DupCaller/blob/main/README.md) [![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://github.com/AlexandrovLab/DupCaller/actions/workflows/test.yml/badge.svg)](https://github.com/AlexandrovLab/DupCaller/actions)[![Uptime Robot status](https://img.shields.io/uptimerobot/status/m795312784-02766a79f207f67626cef289)](https://stats.uptimerobot.com/jjqW4Ulymx)

DupCaller is a universal tool for calling somatic mutations and calculating somatic mutational burden from barcoded error-corrected next generation sequencing (ecNGS) data with matched normal (e.x. NanoSeq).

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Pipeline](#pipeline)
  - [Step 1: Index Reference Genome](#step-1-index-reference-genome)
  - [Step 2: Trim Barcodes](#step-2-trim-barcodes)
  - [Step 3: Align Reads](#step-3-align-reads)
  - [Step 4: Mark Duplicates](#step-4-mark-duplicates)
  - [Step 5: Call Variants](#step-5-call-variants)
  - [Step 6: Estimate Mutational Burden](#step-6-estimate-mutational-burden)
- [Results](#results)
- [End-to-End Examples](#end-to-end-examples)
- [Resources](#resources)
- [Citation](#citation)
- [Copyright](#copyright)
- [Contact](#contact)

---

## Prerequisites

DupCaller requires python>=3.10 to run. Earlier versions may be sufficient to run DupCaller but have not been tested.
The complete DupCaller pipeline also requires the following tools for data preprocessing. The versions are used by the developer and other versions may or may not work.

- BWA version 0.7.17 (https://bio-bwa.sourceforge.net)
- GATK version 4.2.6 (https://github.com/broadinstitute/gatk/releases)
- Tabix for indexing compressed genomic files (recommended installation: `conda install bioconda::tabix`)

---

## Installation

### pip Installation

The tool uses pip for installing scripts and prerequisites. We recommend creating a new environment to install DupCaller:

```bash
conda create -n DupCaller python=3.12 bioconda::tabix
```

To install DupCaller, simply clone this repository and install via pip:

```bash
conda activate DupCaller
git clone https://github.com/AlexandrovLab/DupCaller.git
cd DupCaller
pip install .
```

### Docker / Singularity

A pre-built Docker image is available on Docker Hub at `yuhecheng62/dupcaller:1.1.0-amd64`.

**Pull and run with Singularity:**

Pull the image from Docker Hub (only needed once):
```bash
singularity pull dupcaller-1.1.0.sif docker://yuhecheng62/dupcaller:1.1.0-amd64
```

For quick verification:
```bash
singularity exec dupcaller-1.1.0.sif DupCaller.py --help
```

For installation-free execution of DupCaller commands, run all DupCaller.py commands with `singularity exec` and binding of current directories:
```bash
singularity exec dupcaller-1.1.0.sif --bind $(pwd):$(pwd) DupCaller.py {your commands}
```

---

## Pipeline

### Step 1: Index Reference Genome

DupCaller uses a numpyrized reference genome to perform memory-efficient reference fetching and trinucleotide context indexing. Indexing requires a tabix-indexed BED file of short tandem repeat (STR) regions, which is used to annotate STR loci for improved indel calling near repetitive regions. To index the reference genome, run:

```bash
DupCaller.py index -f reference.fa -s str_regions.bed.gz
```

The command will generate three h5 files in the same folder: `ref.h5`, `tn.h5` and `hp.h5`, which are numpyrized reference sequences, trinucleotide contexts, and homopolymer/STR annotations, respectively. Make sure that when running other DupCaller utilities, the three files are within the same folder as the reference genome.

The STR BED file must be tabix-indexed (`.bed.gz` with a `.tbi` index).

For human reference genome hg38 and mouse reference genome mm39, we provided pre-built indexes and resource files in the [Resources](#resources).

#### Parameters

| Short | Long | Description |
| --- | --- | --- |
| -f | --reference | Reference genome fasta file (required) |
| -s | --strbed | Tabix-indexed BED file of STR regions (required) |

#### Pre-built Indexes

We have built indexes for human reference genome GRCh38/hg38 and mouse reference genome GRCm39/mm39, and can be downloaded at:

#### Short Tandem Repeat File

For other reference genomes, a BED file of simple repeats is needed to build the index. The file can be obtained from the UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables) with the following steps:

1. For **Genome**, select the reference genome.
2. For **Group**, select "Variation and Repeats".
3. For **Table**, select "repeat masker" or similar table.
4. Add filter: "repClass Does match Simple_repeat".
5. For **Output format**, select "BED - browser extensible data".
6. Input desired filename (e.g. `mm39_str.bed`) and download the output BED file.

---

### Step 2: Trim Barcodes

`DupCaller.py trim` extracts 5-prime barcodes from paired-end FASTQs:

```bash
DupCaller.py trim -i read1.fq -i2 read2.fq -p barcode_pattern -o sample_name
```

- `read1.fq` / `read2.fq` — FASTQ files from read 1 and read 2. Both unzipped and gzip-compressed files are supported.
- `barcode_pattern` — pattern of barcodes starting from the 5-prime end, with `N` representing a barcode base and `X` representing a skipped base (similar notation to [UMI-tools](https://github.com/CGATOxford/UMI-tools)). For example, NanoSeq uses a 3-base barcode followed by 4 constant bases, so the pattern should be `NNNXXXX`.
- `sample_name` — prefix of output paired FASTQs. After the run completes, `{sample_name}_1.fastq` and `{sample_name}_2.fastq` will be generated. Barcodes will be recorded in each read name as `{original_read_name}:{read1_barcode}+{read2_barcode}`.

If the matched normal is prepared in the same way as the sample, apply trimming with the same scheme to the matched normal FASTQs. For traditional bulk normal, trimming is not needed.

---

### Step 3: Align Reads

Use a DNA NGS aligner such as BWA-MEM to align the trimmed FASTQs of both sample and matched normal. GATK requires read group fields `ID`, `SM`, and `PL`, so adding those tags during BWA alignment is recommended. **FASTQ tags must be kept — for `bwa mem` this requires the `-C` option.**

```bash
bwa mem -C -t {threads} -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" \
    reference.fa {sample_name}_1.fastq {sample_name}_2.fastq \
    | samtools sort -@ {threads} > {sample_name}.bam
samtools index -@ {threads} {sample_name}.bam
```

- `threads` — number of cores used for aligning
- `reference.fa` — reference genome FASTA file
- `{sample_name}_1.fastq` / `{sample_name}_2.fastq` — trimmed FASTQ files from Step 2

---

### Step 4: Mark Duplicates

Run GATK MarkDuplicates on sample and matched-normal BAMs. Optical and PCR duplicates must be treated differently in ecNGS variant calling:

- Set `--TAGGING_POLICY OpticalOnly` to differentiate optical from PCR duplicates.
- Set `--DUPLEX_UMI true` for duplex UMI handling.
- Set `--READ_NAME_REGEX` as shown below, because the read names were modified in Step 2.

> **Note:** Outdated versions of GATK do not have the `--DUPLEX_UMI` flag. Please update GATK to the latest version if this occurs.

```bash
gatk MarkDuplicates \
    -I sample.bam -O sample.mkdped.bam -M sample.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
```

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

**Input validation:** DupCaller checks for the existence of all required files (BAM, reference, h5 index files, and optional files) before starting analysis, providing clear error messages for missing files.

**Multi-threading:** Coverage files from different threads are automatically merged post-processing with region-aware boundary detection and tabix indexing.

See the [Results](#results) section for descriptions of all output files.

#### Parameters

##### Required

| Short | Long | Description |
| --- | --- | --- |
| -b | --bam | BAM file of ecNGS data |
| -f | --reference | Reference genome FASTA file |
| -o | --output | Prefix of the output files |

##### Recommended

These options should be understood and customized accordingly.

| Short | Long | Description | Default |
| --- | --- | --- | --- |
| -r | --regions | Contigs to consider for variant calling. For non-human species, set accordingly (e.g. for mouse: `-r chr{1..19} chrX chrY`) | chr{1..22} chrX chrY |
| -g | --germline | Indexed germline VCF with AF field | None |
| -p | --threads | Number of threads | 1 |
| -n | --normalBam | BAM file(s) of matched normals. When unavailable, set `-maf` to an appropriate value (e.g. 0.1) | None |
| -m | --noise | BED interval file(s) masking noisy positions | None |
| -R | --regionfile | Inclusive BED file specifying target regions | None |
| -maf | --maxAF | Maximum allele fraction to call a somatic mutation. Must be set when matched normal (`-n`) is unavailable | 1 |
| -tt | --trimF | Ignore mutations less than n bp from template ends | 7 |
| -tr | --trimR | Ignore mutations less than n bp from read ends | 7 |

##### Optional

The effect of changing these parameters should be evaluated before implementation.

| Short | Long | Description | Default |
| --- | --- | --- | --- |
| --naf | | Maximum VAF in matched normal for a mutation to be called | 0.01 |
| --rescue | | Output discarded variants with reason in the FILTER field | False |
| -nm | --nmflt | Filter reads with edit distance larger than this value | 5 |
| -ax | --minMeanASXS | Minimum mean AS-XS alignment score difference for a read group to be considered | 50 |
| -gaf | --germlineAfCutoff | Skip positions with germline AF above this threshold | 0.001 |
| -d | --minNdepth | Minimum coverage in normal for called variants | 10 |

##### Advanced

These are variant calling model parameters; adjustment is unnecessary for general use.

| Short | Long | Description | Default |
| --- | --- | --- | --- |
| -AS | --amperrfile | Pre-learned error profile for amplification SBS error | None |
| -AI | --amperrfileindel | Pre-learned error profile for amplification indel error | None |
| -DS | --dmgerrfile | Pre-learned error profile for SBS damage | None |
| -DI | --dmgerrfileindel | Pre-learned error profile for indel damage | None |
| -mr | --mutRate | Prior somatic mutation rate per base | 2.5e-7 |
| -ts | --thresholdSnv | Log likelihood ratio threshold for SNV calls | 0.5 |
| -ti | --thresholdIndel | Log likelihood ratio threshold for indel calls | 0.5 |
| -mq | --mapq | Minimum MAPQ for an alignment to be considered | 40 |
| -w | --windowSize | Genomic window size for coverage calculation and BAM partitioning | 100000 |
| -bq | --minBq | Bases with quality below this value will be set to 6 | 18 |
| -aq | --minAltQual | Minimum consensus quality of alt allele in a read group | 60 |
| --minRef | | Minimum consensus quality of ref allele in a read group | 2 |
| --minAlt | | Minimum consensus quality of alt allele in a read group | 2 |
| -z | --maxZeroQualFrac | Maximum fraction of zero-quality bases in a read family | 0.5 |
| -id | --indelBed | Indel enhanced Panel of Normals (ePoN) for indel calling | None |

#### Germline and Noise Masks

The `-m` option accepts BED files and will ignore any mutations overlapping an excluded locus. Any custom BED file can be used as input.

---

### Step 6: Estimate Mutational Burden

After mutation calling, run burden estimation from the output folder:

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX
```

Adjust the `-r` regions according to the reference genome used.

#### Per-Gene Coverage

For dNdScv coverage correction, DupCaller can output mean duplex depth per gene using the `-gb` option:

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX -gb {target}.bed
```

The gene BED should have the fourth column formatted as `{gene_name}_{exon_number}` (e.g. `tp53_1`).

#### Re-estimation for Specific Regions

To re-estimate trinucleotide-corrected mutational burden in specific regions without re-running variant calling, use `-rb`:

```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX -rb {re_estimate}.bed
```

#### Parameters

##### Required

| Short | Long | Description |
| --- | --- | --- |
| -i | --prefix | Input prefix of results from `call` command |

##### Reference (choose one)

| Short | Long | Description | Default |
| --- | --- | --- | --- |
| -f | --reference | FASTA file of reference genome | None |
| -ft | --refTrinuc | Precomputed trinucleotide composition of reference genome | None |

##### Optional

| Short | Long | Description | Default |
| --- | --- | --- | --- |
| -r | --regions | Contigs to consider for trinucleotide calculation | chr{1..22} chrX |
| -ot | --outTrinuc | Output the computed trinucleotide composition file for future use | None |
| -c | --clonal | Treat mutations detected in more than one molecule as one mutation | False |
| -d | --dilute | Set to true when sample and matched normal are from the same starting DNA | False |
| -gb | --genebed | Gene BED file for per-gene coverage calculation | None |
| -rb | --reestimatebed | BED file for burden re-estimation in specific regions | None |

---

## Results

### Core Output Files

| File | Description |
| --- | --- |
| `{sample}_snv.vcf` | VCF of detected SNVs and MNVs |
| `{sample}_indel.vcf` | VCF of detected short indel mutations |
| `{sample}_coverage.bed.gz` | Duplex coverage depths across genomic positions. For multi-threaded runs, files from different threads are automatically merged with tabix indexing |
| `{sample}_trinuc_by_duplex_group.txt` | Trinucleotide context counts grouped by duplex read number, used for burden estimation |
| `{sample}_duplex_group_stats.txt` | Statistics for duplex groups including read counts and quality metrics |
| `{sample}_stats.txt` | Overall sequencing and analysis metrics |

### Error Profile Files

| File | Description |
| --- | --- |
| `{sample}.amp.tn.txt` | Amplification SBS error profile by trinucleotide context |
| `{sample}.amp.id.txt` | Amplification indel error profile |
| `{sample}.dmg.tn.txt` | Damage SBS error profile by trinucleotide context |
| `{sample}.dmg.id.txt` | Damage indel error profile |

### Burden Estimation Files

| File | Description |
| --- | --- |
| `{sample}_sbs_burden.txt` | SBS burden with uncorrected and corrected estimates and 95% confidence intervals |
| `{sample}_indel_burden.txt` | Indel burden with 95% confidence intervals, including masked and unmasked calculations |
| `{sample}_sbs_96_corrected.txt` | Corrected SBS counts across 96 trinucleotide contexts for signature analysis |
| `{sample}_sbs_burden_by_min_read_group_size.txt` | Burden estimates stratified by minimum read group size |
| `{sample}_duplex_allele_counts.txt` | Duplex depths and allele counts for each unique mutation |

### Visualization Files

| File | Description |
| --- | --- |
| `{sample}_sbs_96_corrected.png` | 96-context mutational signature plot |
| `{sample}_sbs_burden_by_min_read_group_size.png` | Burden estimates across minimum read group sizes |

### Optional / Conditional Files

| File | Condition | Description |
| --- | --- | --- |
| `{sample}_snv_flt.vcf` | Dilute mode | Filtered SNV VCF from dilute analysis |
| `{sample}_gene_coverage.txt` | `-gb` option | Mean duplex coverage per gene for dNdScv correction |
| `{sample}_sbs_burden_re_estimate.txt` | `-rb` option | Re-estimated burden for specific regions |
| `{sample}_sbs_96_corrected_re_estimate.png` | `-rb` option | Signature plot for re-estimated regions |

---

## End-to-End Examples

The following examples assume the reference genome has already been indexed (`DupCaller.py index`) and that the germline VCF and noise mask are available. Replace all placeholder paths with your actual files.

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
bwa mem -C -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    ${REF} ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq | samtools sort -@ ${THREADS} > ${SAMPLE}.bam
samtools index -@ ${THREADS} ${SAMPLE}.bam

bwa mem -C -t ${THREADS} -R "@RG\tID:${NORMAL}\tSM:${NORMAL}\tPL:ILLUMINA" \
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

When no matched normal is available, omit `-n` and set `-maf` to cap the maximum allele frequency of called variants. A value of 0.1 is appropriate for most somatic applications to exclude common germline variants.

```bash
SAMPLE=sample1
REF=/path/to/hg38.fa
GERMLINE=/path/to/gnomad.hg38.vcf.gz
NOISE=/path/to/noise_mask.bed.gz
THREADS=16

# 1. Trim barcodes
DupCaller.py trim -i ${SAMPLE}_R1.fq.gz -i2 ${SAMPLE}_R2.fq.gz -p NNNXXXX -o ${SAMPLE}

# 2. Align
bwa mem -C -t ${THREADS} -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
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
| Resource | Link/Source|
| --- | --- |
| Reference genome | [TCGA hg38 reference file](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834) |
| Reference index | [Pre-built hg38 DupCaller reference](https://drive.google.com/drive/folders/1v8ut5aky01zujXlZ0b-dhCqImWMyMbiK?usp=drive_link) |
| Germline VCF | af-only-gnomad.hg38.vcf.gz file from the legacy GATK resource bundle. A copy of the file can be found [here](https://www.bcgsc.ca/downloads/morinlab/reference/). |
| Noise mask | NanoSeq noise mask from [Abascal et. al.](https://doi.org/10.1038/s41586-021-03477-4). Can be obtained with approved access to the dataset |

### Mouse (GRCm39/mm39)

| Resource | Link |
| --- | --- |
| Reference genome | [UCSC mm39 reference file](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz) |
| Reference index | [Pre-built mm39 DupCaller reference](https://drive.google.com/drive/folders/1apvqjyO3OhApQ5cnmOe1LWR2wi7_cMqK?usp=drive_link) |
| Germline VCF | [mgp_strains](https://drive.google.com/drive/folders/1apvqjyO3OhApQ5cnmOe1LWR2wi7_cMqK?usp=drive_link) — includes VCF for SNPs in popular mouse strains. Use the vcf file for your strain. |
| Noise mask | [NOISE.mm39.bed.gz](https://drive.google.com/drive/folders/1apvqjyO3OhApQ5cnmOe1LWR2wi7_cMqK?usp=drive_link), in-house noise mask from mouse duplex sequencing data |


## Citation

Cheng, Y. et al. Improved Mutation Detection in Duplex Sequencing Data with Sample-Specific Error Profiles. bioRxiv (2025). https://doi.org/10.1101/2025.07.13.664565

---

## Copyright

Copyright (c) 2024, Yuhe Cheng [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## Contact

Yuhe Cheng (yuc211@ucsd.edu)
