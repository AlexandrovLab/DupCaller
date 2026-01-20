# DupCaller Nextflow Pipeline

This directory contains a Nextflow implementation of the DupCaller pipeline for calling somatic mutations from error-corrected next generation sequencing (ecNGS) data.

## Prerequisites

- Nextflow >= 21.04.0
- DupCaller installed and in PATH
- BWA-MEM >= 0.7.17
- GATK >= 4.2.6
- Samtools
- Reference genome indexed with `DupCaller.py index`

## Quick Start

### Basic usage with matched normal

```bash
nextflow run DupCaller.nf \
  --sample_name sample1 \
  --read1 sample1_R1.fq.gz \
  --read2 sample1_R2.fq.gz \
  --normal_bam normal.mkdped.bam \
  --reference /path/to/hg38.fa \
  --germline_vcf /path/to/gnomad.vcf.gz \
  --threads 8
```

### Panel without matched normal

For mutagenesis panels without a matched normal, set `--max_af` and use `--skip_normal`:

```bash
nextflow run DupCaller.nf \
  --sample_name panel_sample \
  --read1 sample_R1.fq.gz \
  --read2 sample_R2.fq.gz \
  --skip_normal \
  --max_af 0.1 \
  --reference /path/to/hg38.fa \
  --germline_vcf /path/to/gnomad.vcf.gz \
  --threads 4
```

### With normal FASTQ files

Process both tumor and normal from FASTQ:

```bash
nextflow run DupCaller.nf \
  --sample_name tumor \
  --read1 tumor_R1.fq.gz \
  --read2 tumor_R2.fq.gz \
  --normal_read1 normal_R1.fq.gz \
  --normal_read2 normal_R2.fq.gz \
  --reference /path/to/hg38.fa \
  --germline_vcf /path/to/gnomad.vcf.gz \
  --noise_mask /path/to/noise_mask.bed.gz \
  --threads 8
```

## Pipeline Steps

The pipeline executes the following steps:

1. **Index Reference** - Creates h5 index files if not present
2. **Trim Barcodes** - Extracts 5-prime barcodes from paired-end FASTQ
3. **BWA Alignment** - Aligns reads with BWA-MEM (with -C flag to keep FASTQ tags)
4. **Mark Duplicates** - GATK MarkDuplicates with duplex UMI settings
5. **Call Variants** - DupCaller variant calling with error profiling
6. **Estimate Burden** - Mutation burden estimation with trinucleotide correction

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--sample_name` | Sample name for output files |
| `--reference` | Reference genome FASTA file |

### Input Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--read1` | Read 1 FASTQ file | - |
| `--read2` | Read 2 FASTQ file | - |
| `--normal_bam` | Pre-aligned normal BAM file | - |
| `--normal_read1` | Normal read 1 FASTQ | - |
| `--normal_read2` | Normal read 2 FASTQ | - |
| `--skip_normal` | Run without matched normal | false |

### Reference Resources

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--germline_vcf` | Indexed germline VCF with AF field | - |
| `--noise_mask` | BED file for noise masking | - |
| `--indel_bed` | Enhanced panel of normal for indels | - |

### Pipeline Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | ./results |
| `--threads` | Number of threads | 4 |
| `--barcode_pattern` | Barcode pattern (N=barcode, X=skip) | NNNXXXX |

### Variant Calling Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--regions` | Contigs for calling | chr1-22,X,Y |
| `--max_af` | Maximum allele frequency | 1.0 |
| `--germline_af_cutoff` | Germline AF cutoff | 0.001 |
| `--min_n_depth` | Minimum normal depth | 10 |
| `--max_zero_qual_frac` | Maximum zero quality fraction | 0.5 |
| `--trim_template` | Template end trim distance | 7 |
| `--trim_read` | Read end trim distance | 7 |
| `--mapq` | Minimum mapping quality | 40 |
| `--window_size` | Genomic window size | 100000 |

### Burden Estimation

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--estimate_clonal` | Consider clonal mutations | false |
| `--estimate_dilute` | Dilute mode (same starting material) | false |
| `--gene_bed` | Gene BED for coverage calculation | - |

## Execution Profiles

### Local execution (default)

```bash
nextflow run DupCaller.nf [options]
```

### SLURM cluster

```bash
nextflow run DupCaller.nf -profile slurm [options]
```

### SGE cluster

```bash
nextflow run DupCaller.nf -profile sge [options]
```

### PBS cluster

```bash
nextflow run DupCaller.nf -profile pbs [options]
```

### Docker

```bash
nextflow run DupCaller.nf -profile docker [options]
```

### Singularity

```bash
nextflow run DupCaller.nf -profile singularity [options]
```

## Output Files

Results are organized in `${params.outdir}/${sample_name}/`:

### Trimmed FASTQ
- `trimmed/` - Barcode-trimmed FASTQ files

### Aligned BAMs
- `aligned/` - BWA-aligned BAM files
- `markdup/` - Duplicate-marked BAM files and metrics

### Variant Calls
- `variants/`
  - `*_snv.vcf` - Single nucleotide variants
  - `*_indel.vcf` - Insertion/deletion variants
  - `*_coverage.bed.gz` - Duplex coverage depths
  - `*_trinuc_by_duplex_group.txt` - Trinucleotide contexts
  - `*.amp.tn.txt`, `*.dmg.tn.txt` - Error profiles (SBS)
  - `*.amp.id.txt`, `*.dmg.id.txt` - Error profiles (indels)

### Burden Estimates
- `burden/`
  - `*_sbs_burden.txt` - SBS burden estimates
  - `*_indel_burden.txt` - Indel burden estimates
  - `*_sbs_96_corrected.txt` - 96 trinucleotide context counts
  - `*_sbs_96_corrected.png` - Mutational signature plot
  - `*_duplex_allele_counts.txt` - Allele counts per mutation

### Pipeline Info
- `pipeline_info/`
  - `execution_timeline.html` - Timeline of process execution
  - `execution_report.html` - Resource usage report
  - `execution_trace.txt` - Detailed trace of all processes
  - `pipeline_dag.svg` - Directed acyclic graph of pipeline

## Example: Complete Workflow

```bash
# 1. Index reference genome (if not already done)
DupCaller.py index -f /path/to/hg38.fa

# 2. Run the Nextflow pipeline
nextflow run DupCaller.nf \
  --sample_name SAMPLE001 \
  --read1 SAMPLE001_R1.fastq.gz \
  --read2 SAMPLE001_R2.fastq.gz \
  --normal_bam NORMAL001.mkdped.bam \
  --reference /path/to/hg38.fa \
  --germline_vcf /path/to/gnomad.hg38.vcf.gz \
  --noise_mask /path/to/noise_mask.bed.gz \
  --indel_bed /path/to/indel_ePoN.bed.gz \
  --threads 16 \
  --outdir ./results \
  -profile slurm

# 3. Results will be in ./results/SAMPLE001/
```

## Customizing Resource Requirements

Edit `nextflow.config` to adjust memory, CPU, and time limits:

```groovy
process {
    withName: CALL_VARIANTS {
        cpus = 32
        memory = 64.GB
        time = 48.h
    }
}
```

## Troubleshooting

### Reference index files not found

Make sure to run `DupCaller.py index -f reference.fa` before running the pipeline, and ensure the `.h5` files are in the same directory as the reference FASTA.

### GATK MarkDuplicates fails with DUPLEX_UMI error

Update GATK to version >= 4.2.6. Older versions do not support the `--DUPLEX_UMI` flag.

### Out of memory errors

Increase memory allocations in `nextflow.config` for the failing process, or reduce `--threads` for parallel processes.

### Pipeline resumes from checkpoint

Nextflow automatically caches completed tasks. To restart from scratch:

```bash
nextflow run DupCaller.nf [options] -resume
```

To force a clean restart:

```bash
rm -rf work/
nextflow run DupCaller.nf [options]
```

## Citation

If you use this pipeline, please cite:

Cheng, Y. et al. Improved Mutation Detection in Duplex Sequencing Data with Sample-Specific Error Profiles. bioRxiv (2025). https://doi.org/10.1101/2025.07.13.664565

## Support

For issues and questions:
- GitHub: https://github.com/YuheCheng62/DupCaller
- Email: yuc211@ucsd.edu
