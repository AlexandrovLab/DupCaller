# `DupCaller.py estimate` — Output File Reference

All files except the run log are written under the sample directory specified by `-i / --prefix` (call it `{prefix}`, with sample name `{sample}` = the basename of `-i`), following the same `SBS/`, `INDEL/`, `DBS/` subfolder layout `call` uses. Files that are shared across mutation types (`_duplex_allele_counts.txt`, `_gene_coverage.txt`, and the appended `_stats.txt` lines) stay at the top level of `{prefix}`. The run log is written to the path `{prefix}_estimate_params.log`, i.e. **next to** the sample directory, not inside it (for `-i results/sample1` it is `results/sample1_estimate_params.log`; with a trailing `/` on `-i` it lands inside the directory as `_estimate_params.log`).

`estimate` runs in one of two modes. Without `-rb/--reestimatebed` it writes the main outputs:

```
{prefix}_estimate_params.log            (beside {prefix}/, see above)
{prefix}/
├── {sample}_duplex_allele_counts.txt   (not if call used --skipCoveragePass)
├── {sample}_gene_coverage.txt          (conditional, -gb; not if call used --skipCoveragePass)
├── {sample}_stats.txt                  (three lines appended by this step — see call_outputs.md)
├── SBS/   _sbs_burden.txt, _sbs_96_corrected.txt, SBS_96_plots_{sample}.pdf,
│          _sbs_burden_by_group_size.txt/.pdf,
│          _sbs_flt.vcf   (conditional, -d)
├── INDEL/ _indel_burden.txt, _indel_83_corrected.txt, ID_83_plots_{sample}.pdf,
│          _indel_burden_by_group_size.txt/.pdf
└── DBS/   _dbs_burden.txt, _dbs_78_corrected.txt, DBS_78_plots_{sample}.pdf,
           _dbs_burden_by_group_size.txt/.pdf
```

With `-rb/--reestimatebed` it writes **only** the re-estimation files (plus the run log) — none of the files above, and nothing is appended to `_stats.txt`:

```
{prefix}_estimate_params.log
{prefix}/
├── SBS/   _sbs_burden_re_estimate.txt, _sbs_96_corrected_re_estimate.txt,
│          SBS_96_plots_{sample}_re_estimate.pdf
└── INDEL/ _indel_burden_re_estimate.txt
```

`-ot/--outTrinuc FILE` (either mode) additionally writes the reference trinucleotide composition used by the run to `FILE`: a `# regions: …` comment line (the `-r` contigs), then a tab-separated table with a `trinuc` index (the 64 raw contexts) and a `count` column. Pass it back with `-ft/--refTrinuc` to skip rescanning the reference.

The `SBS_96_plots_*.pdf`/`ID_83_plots_*.pdf`/`DBS_78_plots_*.pdf` filenames are chosen by the `sigProfilerPlotting` library itself (not by DupCaller), from the sample name passed in.

---

## Burden files (`_sbs_burden.txt`, `_indel_burden.txt`, `_dbs_burden.txt`)

Tab-separated key–value files (one metric per line, no header), each computed at the most inclusive stratum: minimum duplex group size = 1 (i.e. `min(F1R2, F2R1) >= 1`, all groups combined). All three files share the same field order and meaning; differences are called out below.

Each file's `Corrected mutation number` line is immediately followed by three genome-wide-extrapolated figures — `Mutation number per genome` and its 95% lower/upper bounds — computed as `Corrected burden × Reference base number` (and the corresponding CI bounds scaled the same way): a genuine extrapolation to the full considered reference footprint, not just the effective duplex-covered subset, so it reads larger than `Corrected mutation number`. `Reference base number` counts each considered reference base once, so this is the number of mutations **per haploid genome** (one copy of the `-r` regions); multiply by 2 for a diploid genome. This is distinct from the per-context `mutation_number_genome` column in `_sbs_96_corrected.txt`/`_indel_83_corrected.txt` (one value per SBS96/ID83 channel rather than a single aggregate; `_dbs_78_corrected.txt` has no such column).

**Confidence interval method:** `Uncorrected burden`'s CI (and `Unmasked burden`'s) is the standard exact Poisson interval (chi-square-based). `Corrected burden`'s CI (and `Mutation number per genome`'s, scaled from it) instead uses the Fay-Feuer method (Fay & Feuer 1997) — the same one behind age-standardized-rate CIs (e.g. SEER*Stat) — because a correction-ratio-weighted sum of per-channel counts isn't itself Poisson-distributed, so the closed-form Poisson interval doesn't apply to it. Fay-Feuer models the weighted sum as a gamma distribution matched to its first two moments (mean and variance under a per-channel Poisson-MLE assumption), with a built-in conservative correction for the upper bound (one hypothetical extra event in the single largest-weight channel) that keeps it well-defined and non-degenerate even when zero mutations are observed. Deterministic — no RNG/seed involved.

`Duplex coverage` is the opportunity coverage at `min_group_size=1`, normalized to true per-base units:
- **SBS:** the raw per-alt-base coverage sum, divided by 3 (the 3 non-ref alt bases per locus).
- **INDEL:** the raw per-opportunity-column coverage sum, divided by `indel_locus_multiplier = genome_indel_cov / reference_base_number` — a genome-wide, sequence-composition-only ratio of total ID83 opportunity-channel count to total reference base count, since one locus can count toward several overlapping ID83 channels at once (see `call_outputs.md`'s `_coverage.bed.gz` section).
- **DBS:** the raw DBS-opportunity sum, divided by 9 (the 9 non-ref alt pairs per reference dinucleotide).

These same three normalized values are also appended to `_stats.txt` as `SBS Base Coverage`/`Indel Base Coverage`/`DBS Base Coverage` (see `call_outputs.md`).

### `SBS/{sample}_sbs_burden.txt`

| Field | Description |
| --- | --- |
| `Uncorrected burden` | Raw SNV burden: total observed mutation count / total effective coverage, at `min(F1R2,F2R1)>=1`. |
| `Uncorrected burden 95% lower/upper` | Poisson 95% confidence interval on the uncorrected burden. |
| `Uncorrected mutation number` | Total count of `PASS` SNVs used in the estimate (only those whose SBS96 class has nonzero observed coverage). |
| `Corrected burden` | Trinucleotide-corrected SNV burden (see method below). |
| `Corrected burden 95% lower/upper` | Fay-Feuer 95% CI on corrected burden (see CI method above). |
| `Corrected mutation number` | Sum of correction-ratio-weighted per-SBS96-class mutation counts. |
| `Mutation number per genome` | Genome-wide extrapolation per **haploid** genome (see above); ×2 for diploid. |
| `Mutation number per genome 95% lower/upper` | Fay-Feuer CI on `Corrected burden`, scaled to genome-wide units. |
| `Duplex coverage` | See normalization above. |
| `Unmasked burden` | SNV burden counting every `_sbs.vcf` record plus the `_sbs_fail.vcf` records with FILTER `masked`, `SNPM` or `NOISEM` set and tumor `DP` > 0 (snp_mask/noise_mask candidates that cleared LR and got real depth extracted), divided by `_stats.txt`'s `Unmasked Coverage` — i.e. the burden before the snp/noise masks are applied as filters. Not affected by `-d`. `N/A` (with its CI) when the call used `--doNotOutputMaskedSite`. |
| `Unmasked burden 95% lower/upper` | Poisson 95% CI on unmasked burden. |
| `Reference base number` | Total trinucleotide count of the reference genome across the considered regions (`ref_trinuc_64.sum()`). |

**Trinucleotide correction method:** because the trinucleotide composition of the sequenced/covered regions may differ from the genome-wide average, each of the 96 SBS classes gets its own correction ratio:

```
correction_ratio[class] = (ref_genome_trinuc_count[class] / sum(ref_genome_trinuc_count))
                        / (observed_coverage[class] / sum(observed_coverage))
corrected_mutation_count[class] = uncorrected_count[class] × correction_ratio[class]
corrected_burden = sum(corrected_mutation_count) / total_effective_coverage
```

(A class with zero observed coverage gets ratio 0.) This reweights mutation counts so the burden estimate reflects what would be observed if the trinucleotide composition matched the genome-wide average, enabling comparison across samples with different target regions.

**`-d/--dilute`:** a PASS SNV with tumor `AC` > 1 whose tumor and matched-normal allele counts differ significantly (Barnard's exact test, p ≤ 0.05) is left out of all SBS counts — the `_sbs_burden.txt` figures other than `Unmasked burden`, `_sbs_96_corrected.txt`, `SBS_96_plots` and `_sbs_burden_by_group_size` — as well as out of `_sbs_flt.vcf` (see below). It still counts in `_duplex_allele_counts.txt`.

### `INDEL/{sample}_indel_burden.txt`

Same field list as `_sbs_burden.txt`, minus the three `Unmasked burden*` lines (same omission and same reason as `_dbs_burden.txt` below: normalizing an indel-opportunity sum to true per-base units needs `indel_locus_multiplier`, which doesn't exist at `call` time, so there's no `call`-time "Unmasked Indel Coverage" denominator in `_stats.txt` to build it from): `Uncorrected burden`, `..95% lower/upper`, `Uncorrected mutation number`, `Corrected burden`, `..95% lower/upper`, `Corrected mutation number`, `Mutation number per genome`, `..95% lower/upper`, `Duplex coverage`, `Reference base number`. Correction here is per-ID83-channel instead of per-SBS96-class. The burdens and their CIs are multiplied by `indel_locus_multiplier`, `Mutation number per genome` (and bounds) is `Corrected burden × Reference base number` on those rescaled burdens, and `Duplex coverage` is divided by it, so they read in the same per-reference-base units as SBS/DBS. `Uncorrected mutation number` and `Corrected mutation number` are counts and are not rescaled.

### `DBS/{sample}_dbs_burden.txt`

Same field list, minus the three `Unmasked burden*` lines (DBS calling has no masked tier — a DBS is only formed from two unmasked PASS-candidate SBS, so there's no unmasked-vs-masked distinction to report): `Uncorrected burden`, `..95% lower/upper`, `Uncorrected mutation number`, `Corrected burden`, `..95% lower/upper`, `Corrected mutation number`, `Mutation number per genome`, `..95% lower/upper`, `Duplex coverage`, `Reference base number`.

---

## Per-context corrected count tables

### `SBS/{sample}_sbs_96_corrected.txt`

One row per SBS96 trinucleotide context (standard `flanking[REF>ALT]flanking` notation, e.g. `A[C>A]A`).

| Column | Description |
| --- | --- |
| `mutation_number_uncorrected` | Raw observed mutation count for this context. |
| `mutation_number_corrected` | `uncorrected × correction_ratio` for this context. |
| `correction_ratio` | `ref_genome_context_freq / observed_coverage_context_freq` for this context. |
| `mutation_number_genome` | Estimated mutation count across the whole (haploid) reference genome for this context: `mutation_rate[context] × ref_genome_trinuc_count[context]`, where `mutation_rate` is the uncorrected count over this context's observed coverage. |
| `trinuc_number_genome` | Reference genome trinucleotide count for this context (the multiplier used to compute `mutation_number_genome`). |
| `mutations_per_opportunity` | `mutation_number_corrected / opportunity[context]`, at the `min(F1R2,F2R1)>=1` stratum (all duplex groups, the same stratum as every other column); `opportunity` is the context's alt-base-specific coverage (not divided by 3). |

Directly usable for mutational-signature analysis (SigProfiler, MutationalPatterns, etc.).

### `SBS_96_plots_{sample}.pdf`

Two-page PDF (one page per column) rendering `mutation_number_uncorrected` and `mutation_number_corrected` from the file above as COSMIC-style SBS96 bar charts, coloured by substitution type (C>A, C>G, C>T, T>A, T>C, T>G). `mutation_number_genome` is **not** plotted (it stays in the `.txt` table only).

### `INDEL/{sample}_indel_83_corrected.txt`

Same idea, at ID83 resolution (83 rows — homopolymer, STR, and microhomology indel channels; our own label syntax, e.g. `Cdelhp1`, `2delstr3`, `5+delMH2`).

| Column | Description |
| --- | --- |
| `mutation_number_uncorrected` | Raw observed indel count for this ID83 channel. |
| `mutation_number_corrected` | Correction-ratio-weighted count. |
| `correction_ratio` | Per-channel correction ratio (grouped across microhomology-length sub-channels). |
| `mutation_number_genome` | Genome-extrapolated indel count for this channel, per haploid genome. |
| `indel83_number_genome` | Reference genome ID83-opportunity count for this channel. |
| `mutations_per_opportunity` | Corrected count divided by observed opportunity, at the `min(F1R2,F2R1)>=1` stratum. |

### `ID_83_plots_{sample}.pdf`

Two-page PDF, same structure as `SBS_96_plots`, using COSMIC ID83 label translations for the plot axis (our internal labels are relabeled via `INDEL83_TO_SIGPROFILER_LABELS` for this plot only — the `.txt` table keeps our own label syntax).

### `DBS/{sample}_dbs_78_corrected.txt`

Same idea, at DBS78 resolution (78 rows, e.g. `AC>CA` — already in standard COSMIC DBS78 label form).

| Column | Description |
| --- | --- |
| `mutation_number_uncorrected` | Raw observed DBS count for this channel. |
| `mutation_number_corrected` | Correction-ratio-weighted count. |
| `correction_ratio` | Per-channel correction ratio. |
| `dbs78_number_genome` | Reference genome DBS78-opportunity count for this channel. |
| `mutations_per_opportunity` | Corrected count divided by observed opportunity (summed over all duplex groups, not divided by 9). |

There is no `mutation_number_genome` column for DBS.

### `DBS_78_plots_{sample}.pdf`

Two-page PDF, same structure, using our own DBS78 labels directly (they already match SigProfilerPlotting's convention, so no relabeling is needed).

---

## Burden-by-group-size files

**`{outdir}/{sample}_{sbs|indel|dbs}_burden_by_group_size.txt`** (written to the matching `SBS/`/`INDEL/`/`DBS/` subfolder)

One wide-format table stratifying burden by duplex group size, with two parallel column blocks sharing a `read_number` axis (1 through 5, where 5 means "5 or more"):

| Column prefix | Meaning |
| --- | --- |
| `min_group_size_*` | Cumulative: "at least this group size" (`min(F1R2,F2R1) >= N`) — the original, most-inclusive-first behavior. |
| `group_size_*` | Exact: "at exactly this group size" (`min(F1R2,F2R1) == N`, with N=5 meaning `>=5` pooled) — non-cumulative, independent bins. |

Each block has the same 7 sub-columns: `Uncorrected_burden`, `Uncorrected_burden_lower`, `Uncorrected_burden_upper`, `Coverage_base`, `Corrected_burden`, `Corrected_burden_lower`, `Corrected_burden_upper`. For indel, all burden/coverage values are already rescaled by `indel_locus_multiplier` into per-reference-base units, matching SBS/DBS.

Higher minimum group sizes give more accurate (lower error rate) but less sensitive estimates due to reduced coverage; comparing rows helps assess whether the burden estimate is stable across group sizes.

**`{outdir}/{sample}_{sbs|indel|dbs}_burden_by_group_size.pdf`**

Two-panel figure (left: cumulative `min_group_size` curve, right: exact `group_size` curve) plotting uncorrected and corrected burden (solid lines with markers) and their 95% CI bounds (dashed lines) across read numbers 1–5 on a log y-axis. A flat line across group sizes indicates the burden estimate is robust to duplex depth.

---

## `{sample}_duplex_allele_counts.txt`

Tab-separated table, one row per unique `PASS` mutation in `_sbs.vcf` or `_indel.vcf` (DBS events, and therefore their `dbs_member` SBS, are not included), used to assess clonal expansion and duplex depth at mutation sites. Written by `estimate` (not `call`), since it requires the merged `_coverage.bed.gz`; not written for a call run with `--skipCoveragePass` (estimate prints "Skipping duplex allele count due to no coverage file."), nor in `-rb` mode. SNVs removed by `-d/--dilute` are still included.

| Column | Description |
| --- | --- |
| `chromosome` | Chromosome of the mutation. |
| `position_start` | 1-based genomic position (VCF `POS`). |
| `ref` | Reference allele. |
| `alt` | Alt allele. |
| `count` | Number of distinct read families carrying this exact mutation (the number of PASS VCF records with this chromosome/position/ref/alt). |
| `duplex_depth` | Power-weighted duplex coverage at this position, for this specific mutation's own alt-base/indel-category column (looked up in `_coverage.bed.gz`) — not a sum across all alt bases or indel categories. For an SNV the lookup is the `_coverage.bed.gz` row whose 0-based `start` is `POS - 1` (the substituted base). For an indel it is the row whose `start` is `POS`, i.e. the base **after** the VCF anchor base, because the indel coverage columns are keyed to the first changed/potentially-run-extending base rather than the left-aligned anchor. Rows where this comes out to 0 are dropped (logged as a warning) rather than emitting a divide-by-zero `duplex_vaf`. |
| `bam_alt_count` | Raw alt allele count from the first-seen record's `TUMOR` `AC` FORMAT field. |
| `bam_depth` | Raw total depth from the first-seen record's `TUMOR` `DP` FORMAT field. |
| `duplex_vaf` | `count / duplex_depth`. |
| `bam_vaf` | `bam_alt_count / bam_depth`. |
| `gene` | Gene name looked up in the `-gb/--genebed` file at this position, or `.` if `-gb` wasn't given (or no gene-bed interval covers the position). |

---

## `{prefix}_estimate_params.log`

**Location:** `-i` value + `_estimate_params.log`, i.e. beside the sample directory rather than inside it (for `-i results/sample1`: `results/sample1_estimate_params.log`). Written in both modes, before anything else.

Plain-text run log: DupCaller version, run timestamp, the exact command line invoked, and the fully resolved value of every `estimate` argument. `estimate` has no RNG-dependent step (all confidence intervals are computed by closed-form/deterministic methods — see the CI method note above) and takes no `--seed`.

---

## `{sample}_gene_coverage.txt`

**Condition:** only written when `-gb/--genebed` is provided (a bgzip-compressed, tabix-indexed BED file — it is read with `pysam.TabixFile`, so a plain `.bed` without a `.tbi` fails: chrom, start, end, and a 4th column of the form `GENE_exonN` — everything before the first `_` is taken as the gene name), the call run wrote `_coverage.bed.gz` (not with `--skipCoveragePass`), and `-rb` is not given.

Tab-separated, **with a header row**, one row per gene:

| Column | Description |
| --- | --- |
| `gene` | Gene name (from the gene-bed's 4th column, split on `_`). |
| `sbs_duplex_depth` | Mean per-base SNV duplex depth across the gene's annotated exon bases: summed `_coverage.bed.gz` columns 4–7 (A/T/C/G), divided by 3 (the same opportunity-to-base normalization as `_sbs_burden.txt`'s `Duplex coverage`), divided by total exon length (`end - start` summed across all of the gene's exon rows; overlapping rows are counted twice, in both the coverage sum and the length). |
| `indel_duplex_depth` | Mean per-base indel duplex depth: summed `_coverage.bed.gz` columns 8–23 (all 16 indel-opportunity categories), divided by the genome-wide `indel_locus_multiplier`, divided by the same exon length. |

Intended as gene-level sequencing-depth input for tools like dNdScv.

---

## Re-estimation Files (`-rb/--reestimatebed`)

**Condition:** only written when `-rb/--reestimatebed` is provided, and then **instead of** every main output above (burden, corrected-table, plot, by-group-size, allele-count and gene-coverage files; no `_stats.txt` lines are appended). Burden is re-estimated using only genomic positions within the supplied BED file, without re-running variant calling — no re-estimate output exists for DBS. The BED must be bgzip-compressed and tabix-indexed (read with `pysam.TabixFile`), its intervals should not overlap (each interval is counted separately), and `-rb` needs the per-locus `_coverage.bed.gz`, so it exits with an error for a call run with `--skipCoveragePass`.

Only the observed side (coverage and PASS calls) is restricted to the region. The reference side — `Reference base number`, `trinuc_number_genome`, the reference composition in the correction ratios and the indel `indel_locus_multiplier` — still comes from the whole `-r` contig set, as in the main mode.

### `SBS/{sample}_sbs_burden_re_estimate.txt`

Same field list as `_sbs_burden.txt`'s uncorrected/corrected block, computed over the re-estimation region only, with **no** `Unmasked burden*` lines:

| Field | Description |
| --- | --- |
| `Uncorrected burden`, `..95% lower/upper` | SNV burden and Poisson CI within the re-estimation region. |
| `Uncorrected mutation number` | SNV count within the region. |
| `Corrected burden`, `..95% lower/upper` | Trinucleotide-corrected burden and Fay-Feuer CI. |
| `Corrected mutation number` | Correction-ratio-weighted mutation count. |
| `Mutation number per genome`, `..95% lower/upper` | Genome-wide extrapolation per haploid genome (`Corrected burden × Reference base number`, CI scaled the same way), matching the main `_sbs_burden.txt` convention. |
| `Duplex coverage` | Real per-locus-equivalent duplex coverage within the region (`trinuc_cov_96.sum()/3`) — not the reference genome's raw trinucleotide total. |
| `Reference base number` | Reference trinucleotide total over the whole `-r` contig set (not just the region) — the same value as in the main `_sbs_burden.txt`. |

### `SBS/{sample}_sbs_96_corrected_re_estimate.txt`

Same structure as `_sbs_96_corrected.txt`, computed over the re-estimation region, but with **no** `mutations_per_opportunity` column (5 columns: `mutation_number_uncorrected`, `mutation_number_corrected`, `correction_ratio`, `mutation_number_genome`, `trinuc_number_genome`).

### `SBS_96_plots_{sample}_re_estimate.pdf`

Same two-page structure as `SBS_96_plots_{sample}.pdf` (uncorrected, corrected), for the re-estimation region.

### `INDEL/{sample}_indel_burden_re_estimate.txt`

A separate, shorter field set (note the field *names* differ slightly from the main `_indel_burden.txt` — "Indel" rather than a bare noun):

| Field | Description |
| --- | --- |
| `Uncorrected indel burden`, `..95% lower/upper` | Indel burden and Poisson CI within the region, per reference base (multiplied by `indel_locus_multiplier`, as in the main `_indel_burden.txt`). |
| `Uncorrected indel number` | Indel count within the region. |
| `Corrected indel burden`, `..95% lower/upper` | Corrected burden and Fay-Feuer CI, per reference base (multiplied by `indel_locus_multiplier`). |
| `Corrected indel number` | Correction-ratio-weighted indel count (not rescaled). |
| `Mutation number per genome`, `..95% lower/upper` | Genome-wide extrapolation per haploid genome (`Corrected indel burden × ` the `-r` reference base count, CI scaled the same way), matching the main `_sbs_burden.txt` convention. This file has no `Reference base number` line. |
| `Indel coverage` | Total ID83-resolution opportunity coverage within the region, in **raw opportunity units** (not divided by `indel_locus_multiplier`, unlike the main `_indel_burden.txt`'s `Duplex coverage`). Because the burdens above are rescaled and this is not, `Uncorrected indel burden` = `Uncorrected indel number / Indel coverage × indel_locus_multiplier`, not the plain ratio. |

---

## `SBS/{sample}_sbs_flt.vcf`

**Condition:** only written when `-d/--dilute` is set, and not by `-rb/--reestimatebed` (which applies the same filter to its counts only).

A filtered subset of `_sbs.vcf`: SNVs with a tumor alt allele count (`AC`) above 1 are dropped if their tumor and matched-normal allele counts differ significantly (Barnard's exact test, p ≤ 0.05). The same SNVs are also left out of the SBS burden, `_sbs_96_corrected.txt` and by-group-size counts (see the `-d/--dilute` note under `_sbs_burden.txt`). Intended for the case where the sample and matched normal come from the same starting DNA material.
