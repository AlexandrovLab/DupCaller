# `DupCaller.py summarize` — Output File Reference

`summarize` collects each input sample's `estimate` outputs (`-i` takes one or more sample result directories, each of which must already have been through both `call` and `estimate`) into a single cross-sample table, plus six per-sample-column matrices (SBS96 and ID83, three count bases each). `-o / --output` gives the main table's path; the matrix files are derived from it (see below).

For each input directory `{input}` (trailing `/` stripped; `{sample}` = its basename), `summarize` reads six files: `{input}/{sample}_stats.txt`, `{input}/SBS/{sample}_sbs_burden.txt`, `{input}/INDEL/{sample}_indel_burden.txt`, `{input}/DBS/{sample}_dbs_burden.txt`, `{input}/SBS/{sample}_sbs_96_corrected.txt` and `{input}/INDEL/{sample}_indel_83_corrected.txt`. It raises `FileNotFoundError` if any of them is missing. The check runs per sample, inside the loop that writes the table: the header and the rows of all earlier samples have already been written to `{output}` by then, and none of the matrix files is written. It reads `_stats.txt` and the burden files by key (via `parse_stats_file`), not by fixed line position, so a field being added/removed/reordered upstream in `Estimate.py` doesn't silently break it — only a key being renamed or removed would.

---

## `{output}` — main table

Tab-separated, one header row plus one row per sample, in the order given to `-i`. Columns, in order:

| Column | Source | Description |
| --- | --- | --- |
| `sample` | basename of the sample's input path | Sample name. |
| `pass_filter_reads` | `_stats.txt`: `Number of Pass-filter Reads` | |
| `unique_reads` | `_stats.txt`: `Number of Read Families` | |
| `read_families` | `_stats.txt`: `Number of Effective Read Families` | |
| `duplication_rate` | `_stats.txt`: `Pass-filter Duplication Rate` | |
| `read_family_efficiency` | `_stats.txt`: `Efficiency` | |
| `sbs_base_coverage` | `_stats.txt`: `SBS Base Coverage` | Requires `estimate` to have already run for this sample (see `call_outputs.md` — this line is appended by `estimate`, not written by `call`). |
| `sbs_uncorrected_mutations` | `_sbs_burden.txt`: `Uncorrected mutation number` | |
| `sbs_uncorrected_burden` | `_sbs_burden.txt`: `Uncorrected burden` | |
| `sbs_uncorrected_burden_upper_ci` | `_sbs_burden.txt`: `Uncorrected burden 95% upper` | |
| `sbs_uncorrected_burden_lower_ci` | `_sbs_burden.txt`: `Uncorrected burden 95% lower` | |
| `sbs_corrected_mutations` | `_sbs_burden.txt`: `Corrected mutation number` | |
| `sbs_mutations_per_genome` | `_sbs_burden.txt`: `Mutation number per genome` | Distinct from `sbs_corrected_mutations`: this is the corrected burden extrapolated to the full considered reference footprint (not just the effective duplex-covered subset), so it is numerically larger. It is per **haploid** genome (one copy of each considered reference base); ×2 for diploid. See `write_genome_extrapolation` in `Estimate.py` and `estimate_outputs.md`. |
| `sbs_mutations_per_genome_upper_ci` | `_sbs_burden.txt`: `Mutation number per genome 95% upper` | |
| `sbs_mutations_per_genome_lower_ci` | `_sbs_burden.txt`: `Mutation number per genome 95% lower` | |
| `genome_length` | `_sbs_burden.txt`: `Reference base number` | Written as an integer. |
| `sbs_corrected_burden` | `_sbs_burden.txt`: `Corrected burden` | |
| `sbs_corrected_burden_upper_ci` | `_sbs_burden.txt`: `Corrected burden 95% upper` | |
| `sbs_corrected_burden_lower_ci` | `_sbs_burden.txt`: `Corrected burden 95% lower` | |
| `indel_base_coverage` | `_stats.txt`: `Indel Base Coverage` | Same "requires `estimate` to have run" caveat as `sbs_base_coverage`. |
| `indel_uncorrected_mutations` | `_indel_burden.txt`: `Uncorrected mutation number` | |
| `indel_uncorrected_burden` | `_indel_burden.txt`: `Uncorrected burden` | |
| `indel_uncorrected_burden_upper_ci` | `_indel_burden.txt`: `Uncorrected burden 95% upper` | |
| `indel_uncorrected_burden_lower_ci` | `_indel_burden.txt`: `Uncorrected burden 95% lower` | |
| `indel_corrected_mutations` | `_indel_burden.txt`: `Corrected mutation number` | |
| `indel_corrected_burden` | `_indel_burden.txt`: `Corrected burden` | |
| `indel_corrected_burden_upper_ci` | `_indel_burden.txt`: `Corrected burden 95% upper` | |
| `indel_corrected_burden_lower_ci` | `_indel_burden.txt`: `Corrected burden 95% lower` | |
| `indel_mutations_per_genome` | `_indel_burden.txt`: `Mutation number per genome` | Per haploid genome, as for SBS. |
| `indel_mutations_per_genome_upper_ci` | `_indel_burden.txt`: `Mutation number per genome 95% upper` | |
| `indel_mutations_per_genome_lower_ci` | `_indel_burden.txt`: `Mutation number per genome 95% lower` | |
| `dbs_base_coverage` | `_stats.txt`: `DBS Base Coverage` | Same caveat again. |
| `dbs_mutations_per_genome` | `_dbs_burden.txt`: `Mutation number per genome` | Per haploid genome, as for SBS. |
| `dbs_mutations_per_genome_upper_ci` | `_dbs_burden.txt`: `Mutation number per genome 95% upper` | |
| `dbs_mutations_per_genome_lower_ci` | `_dbs_burden.txt`: `Mutation number per genome 95% lower` | |

Every CI pair is written upper bound first, then lower. Only the per-genome figures are taken from `_dbs_burden.txt`; read that file directly for the DBS burden and mutation counts.

---

## Matrix files

Derived from `{output}` by stripping a trailing `.txt` (`{output_prefix}`) and appending a suffix; written after all samples have been read. Each is a tab-separated matrix with a leading `MutationType` column (the row label) and one column per input sample (column header = sample name), rows sorted by label.

| File | Rows | Column values pulled from each sample's table |
| --- | --- | --- |
| `{output_prefix}_SBS96_uncorrected.txt` | 96 SBS96 contexts (e.g. `A[C>A]A`) | `_sbs_96_corrected.txt`: `mutation_number_uncorrected` |
| `{output_prefix}_SBS96_corrected.txt` | same | `_sbs_96_corrected.txt`: `mutation_number_corrected` |
| `{output_prefix}_SBS96_genome.txt` | same | `_sbs_96_corrected.txt`: `mutation_number_genome` (the per-context genome-extrapolated count — see `estimate_outputs.md`) |
| `{output_prefix}_ID83_uncorrected.txt` | 83 ID83 channels, relabeled to SigProfiler/COSMIC syntax (e.g. `1:Del:C:0`) via `INDEL83_TO_SIGPROFILER_LABELS` | `_indel_83_corrected.txt`: `mutation_number_uncorrected` |
| `{output_prefix}_ID83_corrected.txt` | same | `_indel_83_corrected.txt`: `mutation_number_corrected` |
| `{output_prefix}_ID83_genome.txt` | same | `_indel_83_corrected.txt`: `mutation_number_genome` |

These six files are directly usable as multi-sample input to mutational-signature tools (SigProfiler, MutationalPatterns, etc.), one file per mutation type and count basis (raw observed, context-corrected, or genome-extrapolated). No DBS78 matrix is written.
