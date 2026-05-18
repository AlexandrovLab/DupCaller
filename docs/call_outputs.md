# `DupCaller.py call` — Output File Reference

All files are written to the directory specified by `-o / --output`.

---

## `{sample}_snv.vcf` and `{sample}_indel.vcf`

Standard VCF files containing called somatic SNVs/MNVs and short indels respectively. Both files include all variants regardless of filter status; use the FILTER column to select PASS-only calls.

### FILTER values

| Value | Meaning |
| --- | --- |
| `PASS` | Variant passed all filters |
| `masked` | Overlaps a position in the noise mask BED (`-m`) |
| `underpowered` | Log-likelihood ratio below the calling threshold |
| `high_nm` | Read family filtered by the high edit-distance filter (`-nm`) |

### INFO fields

| Tag | Type | Description |
| --- | --- | --- |
| `F1R2` | Integer | Number of F1R2 reads in the read bundle (top-strand read count) |
| `F2R1` | Integer | Number of F2R1 reads in the read bundle (bottom-strand read count) |
| `CS` | Float | Confidence score: log ratio of major to minor consensus base probability |
| `LR` | Float | Log-likelihood ratio of major base over minor base |
| `LM` | Float | Maximum log-likelihood ratio across all alt bases |
| `TC` | Integer×4 | Top-strand base counts in order A, T, C, G |
| `BC` | Float×4 | Bottom-strand base counts in order A, T, C, G |
| `DF` | Integer | Distance of the variant from the fragment end (used for end-trimming filter) |
| `DR` | Integer | Distance of the variant from the read end |
| `TAG1` | String | 5′ barcode of the top-strand read |
| `TAG2` | String | 5′ barcode of the bottom-strand read |
| `SP` | Integer | Reference start position of the read family |
| `TN` | String | Trinucleotide context of the variant (reference base ±1 flanking base) |
| `HP` | Integer | Homopolymer length at the site (always 0 for SNVs) |
| `STR` | Integer | STR length bin at the site: 0 = no STR or <10 bp, 1 = 10–24 bp, 2 = 25–39 bp, 3 = ≥40 bp (always 0 for SNVs) |

### FORMAT fields (per sample)

| Tag | Type | Description |
| --- | --- | --- |
| `AC` | Integer | Count of the alt allele in the read family |
| `RC` | Integer | Count of the ref allele in the read family |
| `DP` | Integer | Total read depth at the position |

---

## `{sample}_snv_flt.vcf`

**Condition:** only written when `--dilute` is set.

A filtered subset of `_snv.vcf` that excludes variants with statistically significant allele-fraction differences between the tumor and normal samples (Barnard's exact test p ≤ 0.05). Used when the sample and matched normal originate from the same starting DNA material.

---

## `{sample}_coverage.bed.gz`

A bgzip-compressed, tabix-indexed BED file with per-position duplex coverage. Each line represents a single genomic base.

### Columns (tab-separated, no header)

| Column | Description |
| --- | --- |
| 1 `chrom` | Chromosome name |
| 2 `start` | 0-based start coordinate |
| 3 `end` | 0-based end coordinate (start + 1 for each base) |
| 4 `snv_coverage` | Number of effective duplex read families covering this position for SNV calling |
| 5 `indel_coverage` | Number of effective duplex read families covering this position for indel calling |

SNV and indel coverage differ because indel calling applies additional site filters (e.g., homopolymer and STR exclusion zones).

Only positions with non-zero coverage in either column are written.

---

## `{sample}_stats.txt`

Tab-separated key–value file with overall library and calling statistics. One metric per line.

| Field | Formula / Source | Description |
| --- | --- | --- |
| Number of Read Families | `unique_read_num` | Total distinct read families (barcode + fragment-start groups) identified. Each family represents one original DNA molecule, regardless of whether both strands were captured. |
| Number of Pass-filter Reads | `pass_read_num` | Total individual reads that pass alignment quality filters: properly paired, not supplementary or secondary, not QC-failed, no duplex-tag marking as synthetic, and not 5′ soft-clipped. |
| Number of Effective Read Families | `duplex_count` | Read families where **both** strand orientations are captured (F1R2 ≥ 1 **and** F2R1 ≥ 1) **and** the family covers at least one callable genomic position. Only these families contribute to variant calling. |
| Effective Coverage | `sum(coverage)` | Total base-level coverage from effective read families at SNV-callable positions. Each position covered by one effective family adds 1. Used as the denominator for SNV burden. |
| Unmasked Coverage | `sum(unmasked_coverage)` | As above, but computed before the noise mask is applied. Used to report the pre-masking (unmasked) SNV burden. |
| Effective Indel Coverage | `sum(coverage_indel)` | Total base-level coverage from effective read families at indel-callable positions. Used as the denominator for indel burden. |
| Unmasked Indel Coverage | `sum(unmasked_coverage_indel)` | Indel coverage before noise masking. |
| Per Read Family Coverage | `Effective Coverage / Number of Effective Read Families` | Mean number of callable positions covered per effective read family. Reflects how much of the target each duplex molecule spans. |
| Pass-filter Duplication Rate | `1 − (Number of Read Families / Number of Pass-filter Reads)` | Fraction of pass-filter reads that are PCR duplicates within their family. Near 0 = little amplification; higher values = more copies per molecule. The denominator only considers pass-filter reads, so the duplication may be higher than the estimation from other tools.|
| Efficiency | `Number of Effective Read Families / total reads processed` | The efficiency only considers read number, not base number, and are generally modulated by duplication rate. Other suboptimal conditions, such as shallow matched normal or large read1-read2 overlap, may further reduce number of duplex bases |

---

## `{sample}_trinuc_by_duplex_group.txt`

Tab-separated matrix used by `DupCaller.py estimate` to compute trinucleotide-corrected mutational burden.

- **Rows:** 32 trinucleotide contexts (pyrimidine reference base only: C or T, with all combinations of ±1 flanking bases, e.g., `ACC`, `ACT`, ..., `TTT`). Row labels are in the first column.
- **Columns:** Duplex group labels in `F1R2+F2R1` format (e.g., `1+1`, `1+2`, `2+2`). Each column represents reads from families with exactly that strand-count combination.
- **Values:** Number of callable genomic positions in that trinucleotide context observed in read families belonging to that duplex group.

---

## `{sample}_duplex_family_strand_composition.txt`

Tab-separated table describing how read families are distributed across strand-count combinations.

### Columns

| Column | Description |
| --- | --- |
| `duplex_group_strand_composition` | Strand composition label as `F1R2+F2R1` (e.g., `0+1`, `1+1`, `2+3`) |
| `read_family_number` | Total number of read families with this strand composition |
| `duplex_group_number` | Number of **effective** duplex families (F1R2 ≥ 1 and F2R1 ≥ 1) in this group that cover at least one callable position |
| `effective_coverage` | Total callable base positions covered by families in this group |

Families with `0` on either strand (e.g., `0+2`) are present in the table but have `duplex_group_number = 0` because they are not true duplex (one strand missing).

---

## `{sample}_duplex_family_strand_composition_heatmap.pdf`

Heatmap visualization of `_duplex_family_strand_composition.txt`. The x-axis is F1R2 count, the y-axis is F2R1 count, and the colour intensity represents the proportion of read families in each cell. Only cells up to the 99th-percentile family-size are shown.

---

## `{sample}_duplex_allele_counts.txt`

Tab-separated table with one row per unique mutation (SNV or indel) detected across all read families. Used to assess clonal expansion and duplex depth at mutation sites.

### Columns

| Column | Description |
| --- | --- |
| `chromosome` | Chromosome of the mutation |
| `position_start` | 1-based genomic position |
| `ref` | Reference allele |
| `alt` | Alt allele |
| `count` | Number of read families carrying the mutation |
| `duplex_depth` | Total duplex depth at this position (from `_coverage.bed.gz`) |
| `bam_alt_count` | Raw alt allele count from the BAM (AC FORMAT field) |
| `bam_depth` | Raw total depth from the BAM (DP FORMAT field) |
| `duplex_vaf` | `count / duplex_depth` |
| `bam_vaf` | `bam_alt_count / bam_depth` |
| `gene` | Gene name if `-gb` was provided, otherwise `.` |

---

## Error Profile Files

These four files capture the sample-specific error model learned during the error-estimation phase (run automatically before variant calling unless pre-learned profiles are provided via `-AS`, `-AI`, `-DS`, `-DI`).

### `{sample}.amp.tn.txt` and `{sample}.dmg.tn.txt`

Tab-separated matrices of SBS error rates. Rows are 32 trinucleotide contexts (same as `_trinuc_by_duplex_group.txt`). Columns are the four possible alt bases: `A`, `T`, `C`, `G`. Values are floating-point error rates (errors per callable site per read family).

- `.amp.tn.txt` — Amplification (PCR) error profile: asymmetric errors arising during PCR amplification of the library.
- `.dmg.tn.txt` — Damage error profile: errors from DNA damage (e.g., oxidative damage producing C→A artefacts), which appear symmetrically on both strands.

### `{sample}.amp.id.txt` and `{sample}.dmg.id.txt`

Tab-separated integer count matrices for indel errors. Used internally for indel variant calling. Format mirrors the SBS error files but for indel-length bins.

- `.amp.id.txt` — Amplification indel error counts.
- `.dmg.id.txt` — Damage indel error counts.
