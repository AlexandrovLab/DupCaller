# `DupCaller.py call` — Output File Reference

All files are written under the directory specified by `-o / --output` (call it `{out}`; the sample name `{sample}` is the basename of the `-o` value). Per-mutation-type outputs live in `{out}/SBS/`, `{out}/INDEL/`, and `{out}/DBS/` subfolders; learned error-profile matrices live in `{out}/ERROR/`; everything else (coverage, stats, family-composition, logs) stays at the top level of `{out}`.

```
{out}/
├── {sample}_call_params.log
├── {sample}_coverage.bed.gz(.tbi)
├── {sample}_duplex_family_strand_composition.txt
├── {sample}_duplex_family_strand_composition_heatmap.pdf
├── {sample}_stats.txt
├── SBS/{sample}_sbs.vcf, _sbs_fail.vcf, _sbs_flt.vcf, _trinuc_by_duplex_group.txt
├── INDEL/{sample}_indel.vcf, _indel_fail.vcf, _indel_by_duplex_group.txt
├── DBS/{sample}_dbs.vcf, _dbs_fail.vcf, _dbs_by_duplex_group.txt
├── ERROR/{sample}.amp.tn.txt, .amp.tn.srd.txt, .amp.hp.txt, .amp.str.txt,
│         .dmg.tn.txt, .dmg.hp.txt, .dmg.str.txt
└── tmp/  (per-worker intermediate files, safe to delete after a successful run)
```

`{sample}_duplex_allele_counts.txt` is documented in [`estimate_outputs.md`](estimate_outputs.md#sampleduplex_allele_countstxt) — it requires `_coverage.bed.gz` and is actually written by `estimate`, not `call`.

---

## `SBS/{sample}_sbs.vcf` and `SBS/{sample}_sbs_fail.vcf`

Called single-base substitutions and MNVs. `_sbs.vcf` contains only `FILTER=PASS` records; `_sbs_fail.vcf` contains everything else (unless `--rescue` is set — see below). The two files together are a complete partition of all evaluated SBS candidates.

## `INDEL/{sample}_indel.vcf` and `INDEL/{sample}_indel_fail.vcf`

Same PASS/fail split, for short indels.

## `DBS/{sample}_dbs.vcf` and `DBS/{sample}_dbs_fail.vcf`

Same PASS/fail split, for dinucleotide substitutions (DBS) — two adjacent SBS calls that form a single dinucleotide-substitution event. A PASS DBS event's two constituent positions are always individually PASS SBS calls too (so a DBS event double-counts as 2 SBS calls unless you exclude its positions when combining SBS and DBS mutation counts).

### FILTER values (shared across all three types)

| Value | Meaning |
| --- | --- |
| `PASS` | All filters passed, including snp_mask/noise_mask (`SNPM`/`NOISEM` INFO fields both 0) — a fully unmasked call. |
| `masked` | Blocked only by `snp_mask`/`noise_mask` (see `SNPM`/`NOISEM` INFO fields for which). Still fully evaluated (LR + depth) and feeds the unmasked-burden numerator, but not `PASS`. A masked candidate that never got real depth extracted (LR below the channel threshold, or `--skipCoveragePass`) is dropped entirely and never appears in `_fail.vcf` unless `--rescue` is set. |
| `underpowered` | Passed the default calling threshold but failed the per-channel FDR-refined threshold. |

The remaining values are only ever emitted when `--rescue` is set (otherwise those candidates are dropped instead of appearing in `_fail.vcf`):

| Value | Meaning |
| --- | --- |
| `high_nm` | Read family failed the NM/blacklist filter. |
| `low_mapq` | Read family failed the mapq filter. |
| `low_ASXS` | Read family failed the AS-XS filter. |
| `cov_mask` | Blocked by the coverage-depth mask (`n_cov_mask`); no coverage/depth extracted. |
| `nm_mask` | Blocked by the per-family NM mask at this position; no coverage/depth extracted. |
| `trim_mask` | Falls in the read-end trim zone; no coverage/depth extracted. |
| `indel_mask` | Indel blocked by `indel_mask`; no coverage/depth extracted. |
| `other_mask` | Blocked by a mask type not covered by a more specific reason above; no coverage/depth extracted. |
| `no_good_alt_read` | Real depth-extraction found zero alt-supporting reads in the tumor BAM pileup, even though LR cleared its channel's threshold on family-consensus evidence alone — real AC/RC/DP attached. |
| `duplex_vaf` | The extracted tumor allele fraction (AC/DP) exceeds `--maxAF` — real AC/RC/DP attached. |
| `normal_vaf` | The extracted matched-normal allele fraction exceeds `--naf` (likely germline or a systematic artifact) — real AC/RC/DP attached. |
| `n_cov_mask` | Matched-normal depth at this position fell below `--minNdepth` once real depth was extracted — real AC/RC/DP attached. |

### INFO fields — SBS and INDEL records

| Tag | Type | Description |
| --- | --- | --- |
| `F1R2` | Integer | Number of F1R2 reads in the read family (top-strand read count). |
| `F2R1` | Integer | Number of F2R1 reads in the read family (bottom-strand read count). |
| `LR` | Float | Log-likelihood ratio of the major base over the minor base. |
| `LM` | Float | Maximum log-likelihood ratio of major base over minor base across alt bases. |
| `TC` | Integer×4 | Top-strand base counts, order A, T, C, G. |
| `BC` | Float×4 | Bottom-strand base counts, order A, T, C, G. |
| `DF` | Integer | Distance of the variant from the fragment (template) end. |
| `DR` | Integer | Distance of the variant from the read end. |
| `TAG1` | String | 5′ barcode of the top-strand read. |
| `TAG2` | String | 5′ barcode of the bottom-strand read. |
| `SP` | Integer | Read family's reference start position. |
| `TN` | String | Trinucleotide context of the variant (reference base ± 1 flanking base). |
| `HP` | Integer | Homopolymer length at the site. Always 0 for SBS. |
| `STR` | Integer | Reference-allele-length bin of short tandem repeats at the site: 0 = no STR or <10bp, 1 = 10–24bp, 2 = 25–39bp, 3 = 40bp+. Always 0 for SBS. |
| `SNPM` | Integer | 1 if this position falls in `snp_mask`. SBS only, always 0 for indel. Never set on a `PASS` record. |
| `NOISEM` | Integer | 1 if this position (or, for indels, any base of its span) falls in `noise_mask`. Never set on a `PASS` record. |
| `FDR` | Float | Local false discovery rate for this call: `(1-mu)/(rawLR*mu + 1-mu)`, where `mu` is this call's channel's final FDR-refined mutation-rate estimate. Reported regardless of filter, so a fail-VCF record's own `FDR` explains why it didn't clear its channel's threshold. `1.0` if the channel has no determinable `mu` (e.g. an indel with no HP/STR channel at all). |

### INFO fields — DBS records (a separate, smaller INFO set)

| Tag | Type | Description |
| --- | --- | --- |
| `F1R2`, `F2R1`, `TAG1`, `TAG2`, `SP` | — | Same meaning as the SBS/INDEL fields above. |
| `TL` | Integer | Template length of the read pair. |
| `FDR` | Float | Same meaning as above, combined across the DBS event's two constituent SBS positions: `1 - (1-fdr1)*(1-fdr2)`. |

DBS records do **not** carry `LR`/`LM`/`TC`/`BC`/`DF`/`DR`/`TN`/`HP`/`STR`/`SNPM`/`NOISEM` — those are only meaningful for a single-base call.

### FORMAT fields (per sample, all three VCF types)

| Tag | Type | Description |
| --- | --- | --- |
| `AC` | Integer | Count of the alt allele in the read family. |
| `RC` | Integer | Count of the ref allele in the read family. |
| `DP` | Integer | Total read depth at the position. |

Two samples are always written, `TUMOR` then `NORMAL` (matching the VCF header's `...FORMAT\tTUMOR\tNORMAL` column order). If no matched normal BAM was given (`-n`/`--normalBams` omitted), the `NORMAL` sample's `AC`/`RC`/`DP` are all `0`.

---

## `SBS/{sample}_sbs_flt.vcf`

**Condition:** only written when `--dilute` is set.

A filtered subset of `_sbs.vcf` that excludes variants with a statistically significant allele-fraction difference between the tumor and matched normal (Barnard's exact test, p ≤ 0.05). Intended for the case where the sample and matched normal come from the same starting DNA material, so real somatic variants should look the same in both.

---

## `SBS/{sample}_trinuc_by_duplex_group.txt`

Tab-separated matrix, one row per **raw, un-folded** (trinucleotide context, alt base) combination — 64 contexts × 3 non-ref alt bases = **192 rows**, labeled `{trinuc}>{alt}` (e.g. `ACC>T`). Reverse-complement pairs are only combined into the 96 canonical (pyrimidine-only) SBS96 classes later, at estimation time.

- **Columns:** duplex-group labels in `F1R2+F2R1` format (e.g. `1+1`, `1+2`, `2+2`), restricted to groups that actually occurred (`duplex_group_number > 0` in the strand-composition table in `_duplex_family_strand_composition.txt` below).
- **Values:** number of callable genomic positions in that (trinuc, alt) context, observed in read families belonging to that duplex group.

---

## `INDEL/{sample}_indel_by_duplex_group.txt`

Same shape and purpose as `SBS/{sample}_trinuc_by_duplex_group.txt` above, but for indel opportunity: 100 rows (the raw indel100 resolution, `build_indel100_labels` in `funcs/misc.py`) × the same duplex-group columns.

---

## `DBS/{sample}_dbs_by_duplex_group.txt`

Same shape and purpose, for DBS opportunity: 144 rows (raw dinucleotide-substitution resolution, `build_dbs_raw144_labels`) × the same duplex-group columns.

---

## `{sample}_call_params.log`

Plain-text run log: DupCaller version, run timestamp, the exact command line invoked, the fully resolved value of every `call` argument (including the RNG `--seed`, generated automatically and recorded here if not passed explicitly), and the resolved paths of the four error-profile files actually used for calling (`amperr_file`, `amperri_file`, `dmgerr_file`, `dmgerri_file`).

---

## `{sample}_coverage.bed.gz`

**Condition:** only written unless `--skipCoveragePass` is set.

A bgzip-compressed, tabix-indexed BED file with per-position, per-alt-base, power-weighted duplex coverage — this is the "effective coverage" that accounts for how likely a real mutation at this depth/context would have been detected, not a raw read count. Only positions where at least one of the 20 coverage values below is nonzero are written; there is no header line.

### Columns (tab-separated, 0-based BED coordinates)

| Column | Name | Description |
| --- | --- | --- |
| 1 | `chrom` | Chromosome name |
| 2 | `start` | 0-based start coordinate |
| 3 | `end` | `start + 1` |
| 4–7 | `cov_A`, `cov_T`, `cov_C`, `cov_G` | L-weighted duplex SNV coverage for calling that specific alt base at this position (float). Exactly 3 of these 4 are ever nonzero at a given locus — whichever alt bases differ from the locus's own reference base; the self/ref-base column is always 0. |
| 8–23 | 16 indel-opportunity categories (see below) | Power-weighted indel-calling opportunity, one column per category (float). |

The 16 indel columns, in order (from `INDEL_COVERAGE_CATEGORY_LABELS` in `funcs/misc.py`):

1. Deletion of Repeat Unit
2. Insertion of Repeat Unit
3. Deletion Length 2
4. Deletion Length 3
5. Deletion Length 4
6. Deletion Length 5+
7. Insertion A
8. Insertion T
9. Insertion C
10. Insertion G
11. Insertion Length 2
12. Insertion Length 3
13. Insertion Length 4
14. Insertion Length 5+
15. Deletion of Homopolymer
16. Insertion of Homopolymer

These 16 columns are not mutually exclusive with each other in the sense of "which category applies" — several can be nonzero at the same locus simultaneously, because a single position can count toward more than one ID83-style opportunity channel at once (e.g. it may sit inside both a homopolymer and an annotated STR unit). Columns 1/2 ("...of Repeat Unit") are STR-only (0 wherever no STR ≥2bp unit is annotated); columns 15/16 ("...of Homopolymer") carry the unconditional 1bp-slip opportunity that exists at every position regardless of STR status. Because of this overlap, summing all 16 columns is **not** a literal per-locus depth the way summing the 3 nonzero SNV columns is — `Estimate.py` corrects for this with a genome-wide `indel_locus_multiplier` (see [`estimate_outputs.md`](estimate_outputs.md)).

---

## `{sample}_duplex_family_strand_composition.txt`

Tab-separated table describing how read families are distributed across strand-count combinations.

| Column | Description |
| --- | --- |
| `duplex_group_strand_composition` | Strand composition label as `F1R2+F2R1` (e.g. `0+1`, `1+1`, `2+3`). |
| `read_family_number` | Total number of read families with this strand composition. |
| `duplex_group_number` | Number of **effective** duplex families (F1R2 ≥ 1 and F2R1 ≥ 1) in this group that cover at least one callable position. |
| `effective_coverage` | Total power-weighted coverage contributed by families in this group. |

Families with `0` on either strand (e.g. `0+2`) appear in the table but have `duplex_group_number = 0`, since they aren't true duplex (one strand is missing).

---

## `{sample}_duplex_family_strand_composition_heatmap.pdf`

Heatmap visualization of the strand-composition table above. The x-axis is F1R2 count, the y-axis is F2R1 count, and colour intensity is the proportion of all read families in each cell. The axis range is capped automatically so that at most ~1% of families fall outside the plotted grid (rather than always showing the single largest observed count, which could be a rare outlier).

---

## Error Profile Files (`{out}/ERROR/`)

These seven files capture the sample-specific error model learned automatically before variant calling, unless `-E/--errprefix` points at an already-learned set (in which case learning is skipped and those files are read from the given prefix instead — six of the seven are required: `.amp.tn.srd.txt`, `.amp.hp.txt`, `.amp.str.txt`, `.dmg.tn.txt`, `.dmg.hp.txt`, `.dmg.str.txt`). Learning is also skipped, even without `-E`, if a prior run already left a complete set of six at the default `ERROR/{sample}` prefix.

| File | Description |
| --- | --- |
| `{sample}.amp.tn.txt` | Raw amplification (PCR) SBS mismatch profile: rows are the 32 pyrimidine-folded trinucleotide contexts, columns are the four alt bases A/T/C/G. A diagnostic/intermediate byproduct of learning — not itself used for calling. |
| `{sample}.amp.tn.srd.txt` | The single-read-damage (SRD) EM-fitted SBS amplification-error rate matrix actually used for calling (`amperr_file`). Same row/column shape as `.amp.tn.txt`. |
| `{sample}.dmg.tn.txt` | Damage SBS error rate matrix (`dmgerr_file`) — errors from DNA damage (e.g. oxidative damage producing C→A artefacts), symmetric on both strands. Same shape as `.amp.tn.txt`. |
| `{sample}.amp.hp.txt` | Amplification indel error rates for homopolymer contexts. |
| `{sample}.dmg.hp.txt` | Damage indel error rates for homopolymer contexts. |
| `{sample}.amp.str.txt` | Amplification indel error rates for short-tandem-repeat (STR) contexts. |
| `{sample}.dmg.str.txt` | Damage indel error rates for short-tandem-repeat (STR) contexts. |

---

## `{sample}_stats.txt`

Tab-separated key–value file (one metric per line) with library and calling statistics. **Written in two stages**: most lines are written immediately by `call`; the last three (`SBS Base Coverage`, `Indel Base Coverage`, `DBS Base Coverage`) are appended later by `estimate`, once it has computed the opportunity-normalized coverage figures — a `_stats.txt` from a `call`-only run (before `estimate` has ever been run on it) will not yet have those three lines.

| Field | Formula / Source | Description |
| --- | --- | --- |
| Number of Read Families | `unique_read_num` | Total distinct read families (barcode + fragment-start groups) identified. Each represents one original DNA molecule, regardless of whether both strands were captured. |
| Number of Pass-filter Reads | `pass_read_num` | Total individual reads passing alignment quality filters (properly paired, not supplementary/secondary/QC-failed, no synthetic duplex-tag marking, not 5′ soft-clipped). |
| Number of Effective Read Families | `duplex_count` | Read families where **both** strand orientations are captured (F1R2 ≥ 1 **and** F2R1 ≥ 1) **and** the family covers at least one callable genomic position. Only these families contribute to variant calling. |
| Unmasked Coverage | `sum(unmasked_coverage) / 3` | Total power-weighted SNV coverage from effective read families, computed **before** the noise mask is applied, already normalized to true per-base units (divided by 3, the number of possible alt bases per locus — same normalization `SBS Base Coverage` below applies). Used to report the pre-masking (unmasked) SNV burden in `_sbs_burden.txt`. |
| Per Read Family Coverage | `(sum(coverage) / 3) / Number of Effective Read Families` | Mean number of callable positions covered per effective read family, in true per-base units. Reflects how much of the target each duplex molecule spans. |
| Pass-filter Duplication Rate | `1 − (Number of Read Families / Number of Pass-filter Reads)` | Fraction of pass-filter reads that are PCR duplicates within their family. Near 0 = little amplification; higher = more copies per molecule. Only pass-filter reads are in the denominator, so this can read higher than duplication estimates from other tools. |
| Efficiency | `Number of Effective Read Families / total alignment records read from the BAM` | Considers read *number* only, not base number; modulated by duplication rate and by other conditions (shallow matched normal, large read1/read2 overlap, etc.) that reduce the number of duplex bases. |
| Total SBS FDR | mean local `FDR` (INFO field) over `PASS` SBS calls whose channel has a determinable `mu` | Mean per-call local false discovery rate actually realized among the SBS calls that passed. `NaN` if there are no such calls. |
| Total Indel FDR | same, over `PASS` indel calls | Mean realized local FDR among passing indel calls. |
| Total DBS FDR | same, over `PASS` DBS calls | Mean realized local FDR among passing DBS calls (each combined from its two constituent SBS positions' FDRs). |
| SBS Base Coverage | *(appended by `estimate`)* | SNV opportunity coverage normalized to true per-base units (divided by 3) at `min_group_size=1` — same value as `_sbs_burden.txt`'s `Duplex coverage` line. |
| Indel Base Coverage | *(appended by `estimate`)* | Indel opportunity coverage normalized to true per-base units (divided by the genome-wide `indel_locus_multiplier`) at `min_group_size=1` — same value as `_indel_burden.txt`'s `Duplex coverage` line. |
| DBS Base Coverage | *(appended by `estimate`)* | DBS opportunity coverage (already per-locus, no further normalization needed) — same value as `_dbs_burden.txt`'s `Duplex coverage` line. |

`Unmasked Coverage` and `Per Read Family Coverage` are both SNV-only and only written when the coverage pass ran (i.e. `--skipCoveragePass` was **not** set). There is no indel equivalent of `Unmasked Coverage`: unlike the SNV `/3` correction above, normalizing an indel-opportunity sum to true per-base units needs the genome-wide `indel_locus_multiplier`, which only exists once `estimate` has run — so a call-time "Unmasked Indel Coverage" field could only ever be raw, not-yet-comparable opportunity units, and isn't written. `_indel_burden.txt` correspondingly has no `Unmasked burden` row (see [`estimate_outputs.md`](estimate_outputs.md)).
