# DupCaller

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/s93d5/wiki/home/) [![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://travis-ci.com/AlexandrovLab/SigProfilerMatrixGenerator.svg?branch=master)](https://app.travis-ci.com/AlexandrovLab/SigProfilerMatrixGenerator)
[![Uptime Robot status](https://img.shields.io/uptimerobot/status/m795312784-02766a79f207f67626cef289)](https://stats.uptimerobot.com/jjqW4Ulymx)

DupCaller is a universal tool for calling somatic mutations and calculating somatic mutational burden from barcoded error-corrected next generation sequencing (ecNGS) data with matched normal (e.x. NanoSeq).

## Prerequisites
DupCaller requires python>=3.10 to run. Earlier versions may be sufficient to run DupCaller but have not been tested.
The complete DupCaller pipeline also requires the following tools for data preprocessing. The versions are used by the developer and other versions may or may not work.

- BWA version 0.7.17 (https://bio-bwa.sourceforge.net)
- GATK version 4.2.6 (https://github.com/broadinstitute/gatk/releases)

## Installation
The tool uses pip for installing scripts and prerequisites. To install DupCaller, simply clone this repository and install via pip:

```bash
git clone https://github.com/AlexandrovLab/DupCaller.git
cd DupCaller
pip install .
```

## Pipeline

### Trim barcodes from reads:

DupCallerTrim.py is a scripts that can extract 5-prime barcodes from paired-end fastqs. The usage is as follows:

```bash
DupCallerTrim.py -i read1.fq -i2 read2.fq -p barcode_pattern -o sample_name
```

where

'read1.fq' and 'read2.fq' are fastq files from read1 and read2 of the paired-end sequencing data, respectively. Both unzipped and gzip compressed files can be correctly processed.

'barcode_pattern' is the pattern of barcodes starting from the 5-prime end, with N representing a barcode base and X representing a skipped base. The notation is similar to what has been used in UMI-tools(https://github.com/CGATOxford/UMI-tools). For example, NanoSeq uses 3-base barcode followed by 4 constant bases, so the pattern should be NNNXXXX.

'sample_name' is the prefix of output paired fastqs. After run complete, '{sample_name}\_1.fastq' and '{sample_name}\_2.fastq' will be generated. The barcodes will be recorded in the each read name as {original_read_name}:{read1_barcode}+{read2_barcode}

If the matched normal is prepared in the same way as the sample, also apply the trimming with the same scheme to the matched normal fastqs. For traditional bulk normal, trimming is not needed.

### Align trimmed fastqs

Use a DNA NGS aligner, such as BWA-MEM, to align the trimmed fastqs of both sample and matched normal from the last step. Notice that GATK requires read group ID,SM and PL to be set, so adding those tags during bwa alignment is recommended. For example:

```bash
bwa mem -t {threads} -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" reference.fa {sample_name}\_1.fastq {sample_name}\_2.fastq | samtools sort -@ {threads} > {sample_name}.bam
samtools index -@ {threads} {sample_name}.bam
```

where

'threads' is the number of cores used for aligning

'reference.fa' is the reference genome fasta file

'{sample_name}\_1.fastq' and '{sample_name}\_2.fastq' are trimmed fastq file from last step.

### MarkDuplicates with optical duplicates tags and new read name configuration

Run GATK MarkDuplicates on sample and matched-normal bams. Notice that optical duplicates and PCR duplicates should be treated differently in ecNGS variant calling, so the "TAGGING_POLICY" of GATK MarkDuplicates should be set to OpticalOnly to differentiate optical duplicate from PCR duplicate. Also, since the read name of trimmed fastq is non-traditional, the READ_NAME_REGEX option should also be set to "(?:.*:)?([0-9]+)[^:]_:([0-9]+)[^:]_:([0-9]+)[^:]\_$". The MarkDuplicates commands should be looking like this:

```bash
gatk MarkDuplicates -I sample.bam -O sample.mkdped.bam -M sample.mkdp_metrics.txt --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" --TAGGING_POLICY OpticalOnly
```

### Variant Calling

After appropriate data preprocessing, DupCallerCall.py should be used to call somatic mutations. The usage depends on your experimental design, and here are some examples.

1. Whole genome/ whole exome / reduced genome (e.x. NanoSeq) with a matched normal

```bash
DupCallerCall.py -b ${sample}.bam -f reference.fa -o {output_predix} -p {threads} -n {normal.bam} -g germline.vcf.gz -m noise_mask.bed.gz
```

2. Mutagenesis panel without a matched normal

```bash
DupCallerCall.py -b ${sample}.bam -f reference.fa -o {output_predix} -p {threads} -g germline.vcf.gz -m noise_mask.bed.gz -ma 0.3
```

3. Deep sequencing of very small gene panels (e.x. less than 5 genes)

```bash
DupCallerCall.py -b ${sample}.bam -f reference.fa -o {output_predix} -p 2 -n {normal.bam} -g germline.vcf.gz -m noise_mask.bed.gz 
```
Note that since DupCaller partition jobs based on genomic regions, multithreading capability will be significantly compromised for small targeted panel. In this case we suggest to run at most one thread per distince targeted region

Please see "Parameters" section for explanation of all parameters. See "Results" section for descriptions of all result files in the output folder

#### Parameters

**Required**

These options are required for each run.

| short option | long option | description |
| --- | --- | --- |
| -b | --bam | bam file of ecNGS data |
| -f | --reference | reference genome fasta file |
| -o | --output | prefix of the output files |

**Recommended**

These options should be understood by user and customized accordingly. Some of them involve resources that should be used when available. All resources for GRCh38/hg38 and GRCm39/mm39 are provided at ... and should be used for the matching reference genome. 

| short option | long option | description | default |
| --- | --- | --- | --- |
| -r | --regions | contigs to consider for variant calling. The default is set for human. For any other species, please set the contigs accordingly. For example, for mouse, please set to "-r chr{1..19} chrX chrY" | default: chr{1..22} chrX chrY |
| -g | --germline | indexed germline vcf with AF field. | None |
| -p | --threads | number of threads | 1 |
| -n | --normalBam | bam file of matched normal. When matched normal is not available, set the maximum allele frequency (-ma) to an appropriate value (e.x. 0.3) | None |
| -m | --noise | a bed interval file that masked noisy positions | None |
| -maf | --maxAF | maximum allele fraction to call a somatic mutation. Must be set to appropriate value when a matched normal (-n) is not available | 1 |
| -mac | --maxAltCount | maximum allele count of alt allele | 100000 |
| -tt | --trimF | ignore mutation if it is less than n bps from ends of template | 30 |
| -tr | --trimR | ignore mutation if it is less than n bps from ends of read | 15 |
| -id | --indelBed | an indel enhanced panel of normal (ePoN) used for indel calling | None | 

**Advanced**

These are variant calling parameters and adjustment is unnecessary for general use.

| short option | long option | description | default |
| --- | --- | --- | --- |
| -aes | --amperrs | prior polymerase substitutionerror rate | 1e-5 |
| -aei | --amperri | prior polymerase indel error rate | 3e-7 |
| -mr | --mutRate | prior somatic mutation rate per base | 2.5e-7 |
| -t | --threshold | log10 likelihood ratio threshold of making a mutation call | 2 |
| -mq | --mapq |minumum mapq for an alignment to be considered | 30 |
| -d | --minNdepth | minumum coverage in normal for called variants | 10 |
| -nad | --maxAltCount | maximum allele count of alt allele in matched-normal | 0 |
| -mnv | --maxMNVlen | maximum length of MNV for mutation calls | 2 |
| -gaf | --germlineAfCutoff | locations at which there is a germline mutation with population af larger than this threshold will be skipped | 0.001 |

#### Results

snv.vcf:
vcf of detected single nucleotide mutations in the sample. The vcf also includes multiple nucleotide mutations (MNVs).

snv_burden.txt:
naive burden and least-square snv burden estimation of the samples with 95% confidence interval.

indel.vcf:
vcf of detected short insertion/deletion (indel) mutations in the sample.

indel_burden.txt:
Indel burden estimation in sample.

## Citation

## Copyright

Copyright (c) 2024, Yuhe Cheng [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Contact Information
