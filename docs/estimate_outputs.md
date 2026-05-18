# `DupCaller.py estimate` — Output File Reference

All files are written to the sample directory specified by `-i / --prefix`.

---

## `{sample}_sbs_burden.txt`

Tab-separated key–value file with SNV mutational burden estimates for the sample (using all reads, i.e. minimum duplex group size = 1).

| Field | Description |
| --- | --- |
| `Uncorrected burden` | Raw SNV burden: total mutation count / total effective coverage |
| `Uncorrected burden 95% lower` | Poisson 95% confidence interval lower bound on uncorrected burden |
| `Uncorrected burden 95% upper` | Poisson 95% confidence interval upper bound on uncorrected burden |
| `Uncorrected mutation number` | Total count of called SNVs used in the estimate |
| `Corrected burden` | Trinucleotide-corrected SNV burden (see correction method below) |
| `Corrected burden 95% lower` | Poisson 95% CI lower bound on corrected burden |
| `Corrected burden 95% upper` | Poisson 95% CI upper bound on corrected burden |
| `Corrected mutation number` | Sum of correction-ratio-weighted per-context mutation counts |
| `Mutation number per genome` | Estimated number of mutations per haploid genome: sum over all 96 SBS contexts of (per-context mutation rate × reference genome trinucleotide count) |
| `Genome coverage` | Total trinucleotide count of the reference genome (used to compute per-genome mutation number) |
| `Unmasked burden` | SNV burden including masked positions (before noise mask is applied) |
| `Unmasked burden 95% lower` | Poisson 95% CI lower bound on unmasked burden |
| `Unmasked burden 95% upper` | Poisson 95% CI upper bound on unmasked burden |

### Trinucleotide correction method

The uncorrected burden treats all genomic positions equally. Because the trinucleotide composition of the sequenced regions may differ from the genome-wide average, a correction is applied:

```
correction_ratio[context] = (ref_genome_freq[context]) / (observed_coverage_freq[context])
corrected_mutation_count[context] = uncorrected_count[context] × correction_ratio[context]
corrected_burden = sum(corrected_mutation_count) / total_effective_coverage
```

This reweights mutation counts so that the burden estimate reflects what would be observed if the trinucleotide composition matched the genome average, enabling comparison across samples with different target regions.

---

## `{sample}_indel_burden.txt`

Tab-separated key–value file with indel mutational burden estimates.

| Field | Description |
| --- | --- |
| `Indel burden` | Raw indel burden: indel count / effective indel coverage |
| `Indel burden 95% lower` | Poisson 95% CI lower bound |
| `Indel burden 95% upper` | Poisson 95% CI upper bound |
| `Indel number` | Total count of called indels |
| `Unmasked Indel burden` | Indel burden before noise masking |
| `Unmasked Indel burden 95% lower` | Poisson 95% CI lower bound on unmasked indel burden |
| `Unmasked Indel burden 95% upper` | Poisson 95% CI upper bound on unmasked indel burden |

---

## `{sample}_sbs_96_corrected.txt`

Tab-separated matrix with one row per SBS96 trinucleotide context (e.g., `A[C>A]A`, `A[C>A]T`, …). Row labels follow the standard `flanking[REF>ALT]flanking` notation.

### Columns

| Column | Description |
| --- | --- |
| `mutation_number_uncorrected` | Raw observed mutation count for this context |
| `mutation_number_corrected` | Correction-ratio-weighted mutation count: `uncorrected × correction_ratio` |
| `correction_ratio` | `(ref_genome_context_freq) / (observed_coverage_context_freq)` for this context |
| `mutation_number_genome` | Estimated mutations per haploid genome for this context: `mutation_rate × ref_genome_trinuc_count` |
| `trinuc_mumber_genome` | Reference genome trinucleotide count for this context (the denominator used to compute mutation_number_genome) |

This file can be used directly for mutational signature analysis (e.g., SigProfiler, MutationalPatterns).

---

## `{sample}_sbs_96.pdf`

Three-page PDF with SBS96 bar-chart profiles for the most inclusive duplex group (minimum group size = 1). Each page shows one of the three count types from `_sbs_96_corrected.txt`:

1. `mutation_number_uncorrected` — raw observed counts
2. `mutation_number_corrected` — trinucleotide-corrected counts
3. `mutation_number_genome` — estimated counts per haploid genome

Bars are coloured by substitution type (C>A, C>G, C>T, T>A, T>C, T>G) following COSMIC convention.

---

## `{sample}_sbs_burden_by_min_read_group_size.txt`

Tab-separated table with burden estimates stratified by minimum duplex group size (1–5). Each row corresponds to a threshold: only read families with both F1R2 ≥ threshold and F2R1 ≥ threshold are included.

### Columns

| Column | Description |
| --- | --- |
| `read number` | Minimum duplex group size threshold (1 = most inclusive, 5 = most restrictive) |
| `Uncorrected_burden` | Uncorrected SNV burden at this threshold |
| `Uncorrected_burden_lower` | Poisson 95% CI lower bound |
| `Uncorrected_burden_upper` | Poisson 95% CI upper bound |
| `Coverage base` | Total effective coverage at this threshold |
| `Corrrected_burden` | Trinucleotide-corrected burden at this threshold |
| `Corrected_burden_lower` | Poisson 95% CI lower bound |
| `Corrected_burden_upper` | Poisson 95% CI upper bound |

Higher minimum group sizes produce more accurate (lower error rate) but less sensitive estimates due to reduced coverage. Comparing rows helps assess whether the burden estimate is stable across group sizes.

---

## `{sample}_sbs_burden_by_min_read_group_size.png`

Line plot of the uncorrected burden (with 95% CI bands) across minimum duplex group sizes 1–5 on a log y-axis. A flat line across group sizes indicates robustness of the burden estimate.

---

## `{sample}_duplex_allele_counts.txt`

See [`call_outputs.md`](call_outputs.md#sampleduplex_allele_countstxt) for the column definitions. This file is written by `estimate` (not `call`) because it requires the coverage BED produced by `call`.

---

## `{sample}_gene_coverage.txt`

**Condition:** only written when `--genebed / -gb` is provided.

Tab-separated file with two columns (no header): gene name and mean duplex coverage across the gene's CDS exons. Coverage is the average per-base effective SNV coverage (`_coverage.bed.gz` column 4) across all exon bases.

Intended as input to dNdScv for correction of mutation rate by gene-level sequencing depth.

---

## Re-estimation Files (`-rb` option)

When `--reestimatebed / -rb` is provided, burden is re-estimated using only genomic positions within the supplied BED file (without re-running variant calling). The following files are written.

### `{sample}_sbs_burden_re_estimate.txt`

Same structure as `_sbs_burden.txt` but computed over the re-estimation region only. Fields:

| Field | Description |
| --- | --- |
| `Uncorrected burden` | SNV burden in the re-estimation region |
| `Uncorrected burden 95% lower/upper` | Poisson 95% CI |
| `Uncorrected mutation number` | SNV count in the re-estimation region |
| `Corrected burden` | Trinucleotide-corrected burden |
| `Corrected burden 95% lower/upper` | Poisson 95% CI on corrected burden |
| `Corrected mutation number` | Array of correction-ratio-weighted per-context counts |
| `mutation number per genome` | Estimated mutations per haploid genome |
| `genome coverage` | Reference genome trinucleotide total |

### `{sample}_indel_burden_re_estimate.txt`

Same structure as `_indel_burden.txt` for the re-estimation region.

| Field | Description |
| --- | --- |
| `Indel burden` | Indel burden in the re-estimation region |
| `Indel burden 95% lower/upper` | Poisson 95% CI |
| `Indel number` | Indel count in the re-estimation region |
| `Indel coverage` | Effective indel coverage in the re-estimation region |

### `{sample}_sbs_96_corrected_re_estimate.txt`

Same structure as `_sbs_96_corrected.txt` but computed over the re-estimation region.

### `{sample}_sbs_96_corrected_re_estimate.pdf`

Same structure as `_sbs_96.pdf` (three pages: uncorrected, corrected, per-genome) but for the re-estimation region.
