# `DupCaller.py estimate` — Output File Reference

All files are written under the sample directory specified by `-i / --prefix` (call it `{prefix}`, with sample name `{sample}` = the basename of `-i`), following the same `SBS/`, `INDEL/`, `DBS/` subfolder layout `call` uses. Files that are shared across mutation types (`_duplex_allele_counts.txt`, `_gene_coverage.txt`, `_estimate_params.log`, and the appended `_stats.txt` lines) stay at the top level of `{prefix}`.

```
{prefix}/
├── {sample}_stats.txt                  (three lines appended by this step — see call_outputs.md)
├── {sample}_duplex_allele_counts.txt
├── {sample}_gene_coverage.txt          (conditional, -gb)
├── {sample}_estimate_params.log
├── SBS/   _sbs_burden.txt, _sbs_96_corrected.txt, SBS_96_plots_{sample}.pdf,
│          _sbs_burden_by_group_size.txt/.pdf,
│          _sbs_burden_re_estimate.txt, _sbs_96_corrected_re_estimate.txt,
│          SBS_96_plots_{sample}_re_estimate.pdf   (conditional, -rb)
├── INDEL/ _indel_burden.txt, _indel_83_corrected.txt, ID_83_plots_{sample}.pdf,
│          _indel_burden_by_group_size.txt/.pdf,
│          _indel_burden_re_estimate.txt   (conditional, -rb)
└── DBS/   _dbs_burden.txt, _dbs_78_corrected.txt, DBS_78_plots_{sample}.pdf,
           _dbs_burden_by_group_size.txt/.pdf
```

The `SBS_96_plots_*.pdf`/`ID_83_plots_*.pdf`/`DBS_78_plots_*.pdf` filenames are chosen by the `sigProfilerPlotting` library itself (not by DupCaller), from the sample name passed in.

---

## Burden files (`_sbs_burden.txt`, `_indel_burden.txt`, `_dbs_burden.txt`)

Tab-separated key–value files (one metric per line, no header), each computed at the most inclusive stratum: minimum duplex group size = 1 (i.e. `min(F1R2, F2R1) >= 1`, all groups combined). All three files share the same field order and meaning; differences are called out below.

**`Corrected mutation number` intentionally appears twice in a row** in all three files. This is a recent redesign: that second line used to report a different, genome-wide-extrapolated `Mutation number per genome` figure, but as of the change that added `Duplex coverage`/`Base Coverage` reporting, it was changed to just repeat the corrected-mutation-number value above it, since the true per-context genome-extrapolated figure remains available separately (in `_sbs_96_corrected.txt`'s/`_indel_83_corrected.txt`'s/`_dbs_78_corrected.txt`'s own `mutation_number_genome` column, per-context rather than a single aggregate). Do not rely on the second occurrence being a distinct value — parse it by line position (index 8, 0-based) if you need to, but expect it to equal line 7.

`Duplex coverage` (recently renamed from `Genome coverage`) is the opportunity coverage at `min_group_size=1`, normalized to true per-base units:
- **SBS:** the raw per-alt-base coverage sum, divided by 3 (the 3 non-ref alt bases per locus).
- **INDEL:** the raw per-opportunity-column coverage sum, divided by `indel_locus_multiplier = genome_indel_cov / reference_base_number` — a genome-wide, sequence-composition-only ratio of total ID83 opportunity-channel count to total reference base count, since one locus can count toward several overlapping ID83 channels at once (see `call_outputs.md`'s `_coverage.bed.gz` section).
- **DBS:** already per-locus; no further normalization needed.

These same three normalized values are also appended to `_stats.txt` as `SBS Base Coverage`/`Indel Base Coverage`/`DBS Base Coverage` (see `call_outputs.md`).

### `SBS/{sample}_sbs_burden.txt`

| Field | Description |
| --- | --- |
| `Uncorrected burden` | Raw SNV burden: total observed mutation count / total effective coverage, at `min(F1R2,F2R1)>=1`. |
| `Uncorrected burden 95% lower/upper` | Poisson 95% confidence interval on the uncorrected burden. |
| `Uncorrected mutation number` | Total count of `PASS` SNVs used in the estimate (only those whose SBS96 class has nonzero observed coverage). |
| `Corrected burden` | Trinucleotide-corrected SNV burden (see method below). |
| `Corrected burden 95% lower/upper` | Poisson-bootstrap 95% CI on corrected burden. |
| `Corrected mutation number` (×2, see above) | Sum of correction-ratio-weighted per-SBS96-class mutation counts. |
| `Duplex coverage` | See normalization above. |
| `Unmasked burden` | SNV burden including PASS calls plus noise-masked candidates that cleared LR and got real depth extracted (i.e. before the noise mask is applied as a filter). |
| `Unmasked burden 95% lower/upper` | Poisson 95% CI on unmasked burden. |
| `Reference base number` | Total trinucleotide count of the reference genome across the considered regions (`ref_trinuc_64.sum()`). |

**Trinucleotide correction method:** because the trinucleotide composition of the sequenced/covered regions may differ from the genome-wide average, each of the 96 SBS classes gets its own correction ratio:

```
correction_ratio[class] = ref_genome_trinuc_count[class] / observed_coverage[class]
corrected_mutation_count[class] = uncorrected_count[class] × correction_ratio[class]
corrected_burden = sum(corrected_mutation_count) / total_effective_coverage
```

This reweights mutation counts so the burden estimate reflects what would be observed if the trinucleotide composition matched the genome-wide average, enabling comparison across samples with different target regions.

### `INDEL/{sample}_indel_burden.txt`

Same field list and order as `_sbs_burden.txt` (`Uncorrected burden`, `..95% lower/upper`, `Uncorrected mutation number`, `Corrected burden`, `..95% lower/upper`, `Corrected mutation number` ×2, `Duplex coverage`, `Unmasked burden`, `..95% lower/upper`, `Reference base number`), correction here being per-ID83-channel instead of per-SBS96-class. All burden/mutation-number values are rescaled by `indel_locus_multiplier` so they read in the same per-reference-base units as SBS/DBS.

### `DBS/{sample}_dbs_burden.txt`

Same field list, minus the three `Unmasked burden*` lines (DBS calling has no masked tier — `_detect_dbs_pairs` only ever emits `filter=="PASS"` events, so there's no unmasked-vs-masked distinction to report): `Uncorrected burden`, `..95% lower/upper`, `Uncorrected mutation number`, `Corrected burden`, `..95% lower/upper`, `Corrected mutation number` ×2, `Duplex coverage`, `Reference base number`.

---

## Per-context corrected count tables

### `SBS/{sample}_sbs_96_corrected.txt`

One row per SBS96 trinucleotide context (standard `flanking[REF>ALT]flanking` notation, e.g. `A[C>A]A`).

| Column | Description |
| --- | --- |
| `mutation_number_uncorrected` | Raw observed mutation count for this context. |
| `mutation_number_corrected` | `uncorrected × correction_ratio` for this context. |
| `correction_ratio` | `ref_genome_context_freq / observed_coverage_context_freq` for this context. |
| `mutation_number_genome` | Genuine estimated mutation count across the whole reference genome for this context: `mutation_rate[context] × ref_genome_trinuc_count[context]` (this is the real per-context version of the figure the aggregate `_sbs_burden.txt` used to report before that line was changed — see the burden-file section above). |
| `trinuc_number_genome` | Reference genome trinucleotide count for this context (the multiplier used to compute `mutation_number_genome`). |
| `mutations_per_opportunity` | `mutation_number_corrected / opportunity[context]`, at the `min(F1R2,F2R1)>=5` (highest-confidence) stratum. |

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
| `mutation_number_genome` | Genome-extrapolated indel count for this channel. |
| `indel83_number_genome` | Reference genome ID83-opportunity count for this channel. |
| `mutations_per_opportunity` | Corrected count divided by observed opportunity. |

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
| `mutations_per_opportunity` | Corrected count divided by observed opportunity. |

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

Two-panel figure (left: cumulative `min_group_size` curve, right: exact `group_size` curve) plotting uncorrected and corrected burden with 95% CI bands across read numbers 1–5 on a log y-axis. A flat line across group sizes indicates the burden estimate is robust to duplex depth.

(These replace an older, SBS-only `_sbs_burden_by_min_read_group_size.txt`/`.png` pair with the same purpose — the current names, shared code path, and PDF-not-PNG format apply uniformly to SBS, INDEL, and DBS.)

---

## `{sample}_duplex_allele_counts.txt`

Tab-separated table, one row per unique `PASS` mutation (SNV or indel) detected across all read families, used to assess clonal expansion and duplex depth at mutation sites. Written by `estimate` (not `call`), since it requires the merged `_coverage.bed.gz`.

| Column | Description |
| --- | --- |
| `chromosome` | Chromosome of the mutation. |
| `position_start` | 1-based genomic position (VCF `POS`). |
| `ref` | Reference allele. |
| `alt` | Alt allele. |
| `count` | Number of distinct read families carrying this exact mutation. |
| `duplex_depth` | Power-weighted duplex coverage at this position, for this specific mutation's own alt-base/indel-category column (looked up in `_coverage.bed.gz`) — not a sum across all alt bases or indel categories. For an indel, the lookup position is the VCF anchor position itself (0-based), one base after where you might expect from the 1-based `POS`, because the coverage columns are keyed to the first changed/potentially-run-extending base rather than the left-aligned anchor. Rows where this comes out to 0 are dropped (logged as a warning) rather than emitting a divide-by-zero `duplex_vaf`. |
| `bam_alt_count` | Raw alt allele count from the first-seen record's `TUMOR` `AC` FORMAT field. |
| `bam_depth` | Raw total depth from the first-seen record's `TUMOR` `DP` FORMAT field. |
| `duplex_vaf` | `count / duplex_depth`. |
| `bam_vaf` | `bam_alt_count / bam_depth`. |
| `gene` | Gene name looked up in the `-gb/--genebed` file at this position, or `.` if `-gb` wasn't given. |

---

## `{sample}_gene_coverage.txt`

**Condition:** only written when `-gb/--genebed` is provided (a BED file: chrom, start, end, and a 4th column of the form `GENE_exonN` — everything before the first `_` is taken as the gene name).

Tab-separated, **with a header row**, one row per gene:

| Column | Description |
| --- | --- |
| `gene` | Gene name (from the gene-bed's 4th column, split on `_`). |
| `sbs_duplex_depth` | Mean per-base SNV duplex depth across the gene's annotated exon bases: summed `_coverage.bed.gz` columns 4–7 (A/T/C/G), divided by 3 (the same opportunity-to-base normalization as `_sbs_burden.txt`'s `Duplex coverage`), divided by total exon length (CDS length, summed across all of the gene's exon rows). |
| `indel_duplex_depth` | Mean per-base indel duplex depth: summed `_coverage.bed.gz` columns 8–23 (all 16 indel-opportunity categories), divided by the genome-wide `indel_locus_multiplier`, divided by the same exon length. |

Intended as gene-level sequencing-depth input for tools like dNdScv. (Older versions of this file had no header and only a single, un-normalized SBS-only coverage column — if you have scripts parsing the old 2-column format, they will need updating for the new 3-column, header-bearing format.)

---

## `{sample}_estimate_params.log`

Plain-text run log: run timestamp, the exact command line invoked, and the fully resolved value of every `estimate` argument (including the RNG `--seed` used for the parametric-bootstrap confidence intervals on corrected burden, generated automatically and recorded here if not passed explicitly).

---

## Re-estimation Files (`-rb/--reestimatebed`)

**Condition:** only written when `-rb/--reestimatebed` is provided. Burden is re-estimated using only genomic positions within the supplied BED file, without re-running variant calling — no re-estimate output exists for DBS.

### `SBS/{sample}_sbs_burden_re_estimate.txt`

Same field list as `_sbs_burden.txt`'s uncorrected/corrected block, computed over the re-estimation region only, with **no** `Unmasked burden*` lines:

| Field | Description |
| --- | --- |
| `Uncorrected burden`, `..95% lower/upper` | SNV burden and Poisson CI within the re-estimation region. |
| `Uncorrected mutation number` | SNV count within the region. |
| `Corrected burden`, `..95% lower/upper` | Trinucleotide-corrected burden and bootstrap CI. |
| `Corrected mutation number` (×2) | Correction-ratio-weighted mutation count (same value repeated, matching the main burden file's convention above). |
| `Duplex coverage` | Real per-locus-equivalent duplex coverage within the region (`trinuc_cov_96.sum()/3`) — not the reference genome's raw trinucleotide total. |
| `Reference base number` | Reference genome trinucleotide total within the region. |

### `INDEL/{sample}_indel_burden_re_estimate.txt`

A separate, shorter field set (note the field *names* differ slightly from the main `_indel_burden.txt` — "Indel" rather than a bare noun):

| Field | Description |
| --- | --- |
| `Uncorrected indel burden`, `..95% lower/upper` | Indel burden and Poisson CI within the region. |
| `Uncorrected indel number` | Indel count within the region. |
| `Corrected indel burden`, `..95% lower/upper` | Corrected burden and bootstrap CI. |
| `Corrected indel number` | Correction-ratio-weighted indel count. |
| `Indel coverage` | Total ID83-resolution opportunity coverage within the region (not rescaled by `indel_locus_multiplier` — unlike the main `_indel_burden.txt`, this file's burden values are already computed directly against this same raw denominator, so no rescaling is needed here). |

### `SBS/{sample}_sbs_96_corrected_re_estimate.txt`

Same structure as `_sbs_96_corrected.txt`, computed over the re-estimation region, but with **no** `mutations_per_opportunity` column (5 columns: `mutation_number_uncorrected`, `mutation_number_corrected`, `correction_ratio`, `mutation_number_genome`, `trinuc_number_genome`).

### `SBS_96_plots_{sample}_re_estimate.pdf`

Same two-page structure as `SBS_96_plots_{sample}.pdf` (uncorrected, corrected), for the re-estimation region.
