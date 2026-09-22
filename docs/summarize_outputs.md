# `DupCaller.py summarize` — Output File Reference

`summarize` collects each input sample's `estimate` outputs (`-i` takes one or more sample result directories, each of which must already have been through both `call` and `estimate`) into a single cross-sample table, plus three per-sample-column SBS96 matrices. `-o / --output` gives the main table's path; the matrix files are derived from it (see below).

For each input sample, `summarize` reads `{sample}/{sample}_stats.txt`, `{sample}/SBS/{sample}_sbs_burden.txt`, `{sample}/INDEL/{sample}_indel_burden.txt`, and `{sample}/SBS/{sample}_sbs_96_corrected.txt`, and errors out immediately (before writing anything) if any of the four is missing for any sample. It reads `_sbs_burden.txt`/`_indel_burden.txt` by key (via `parse_stats_file`), not by fixed line position, so a field being added/removed/reordered upstream in `Estimate.py` doesn't silently break it — only a key being renamed or removed would.

---

## `{output}` — main table

Tab-separated, one header row plus one row per sample, in the order given to `-i`.

| Column | Source | Description |
| --- | --- | --- |
| `sample` | basename of the sample's input path | Sample name. |
| `pass_filter_reads` | `_stats.txt`: `Number of Pass-filter Reads` | |
| `unique_reads` | `_stats.txt`: `Number of Read Families` | |
| `read_families` | `_stats.txt`: `Number of Effective Read Families` | |
| `duplication_rate` | `_stats.txt`: `Pass-filter Duplication Rate` | |
| `read_family_efficiency` | `_stats.txt`: `Efficiency` | |
| `sbs_base_coverage` | `_stats.txt`: `SBS Base Coverage` | Requires `estimate` to have already run for this sample (see `call_outputs.md` — this line is appended by `estimate`, not written by `call`). |
| `uncorrected_mutations` | `_sbs_burden.txt` key `Uncorrected mutation number` | |
| `uncorrected_burden` | `_sbs_burden.txt` key `Uncorrected burden` | |
| `uncorrected_burden_upper_ci` | `_sbs_burden.txt` key `Uncorrected burden 95% upper` | |
| `uncorrected_burden_lower_ci` | `_sbs_burden.txt` key `Uncorrected burden 95% lower` | |
| `corrected_mutations` | `_sbs_burden.txt` key `Corrected mutation number` | |
| `mutations_per_genome` | `_sbs_burden.txt` key `Mutation number per genome` | Distinct from `corrected_mutations`: this is `corrected_mutations` extrapolated to the full considered reference footprint (not just the effective duplex-covered subset), so it is numerically larger. See `write_genome_extrapolation` in `Estimate.py` and `estimate_outputs.md`. |
| `genome_length` | `_sbs_burden.txt` key `Reference base number` | |
| `corrected_burden` | `_sbs_burden.txt` key `Corrected burden` | |
| `corrected_burden_upper_ci` | `_sbs_burden.txt` key `Corrected burden 95% upper` | |
| `corrected_burden_lower_ci` | `_sbs_burden.txt` key `Corrected burden 95% lower` | |
| `indel_base_coverage` | `_stats.txt`: `Indel Base Coverage` | Same "requires `estimate` to have run" caveat as `sbs_base_coverage`. |
| `indel_number` | `_indel_burden.txt` key `Uncorrected mutation number` | |
| `indel_burden` | `_indel_burden.txt` key `Uncorrected burden` | |
| `indel_burden_upper_ci` | `_indel_burden.txt` key `Uncorrected burden 95% upper` | |
| `indel_burden_lower_ci` | `_indel_burden.txt` key `Uncorrected burden 95% lower` | |
| `dbs_base_coverage` | `_stats.txt`: `DBS Base Coverage` | Same caveat again. |

Two things worth noting about what this table does **not** include:
- The indel columns are all **uncorrected** (`indel_number`/`indel_burden*`) — there is no `corrected_indel_burden`/`corrected_indel_number` column here, even though `_indel_burden.txt` has trinucleotide/ID83-corrected equivalents. Read `_indel_burden.txt` directly for those.
- There is no DBS mutation-count or DBS-burden column at all, only `dbs_base_coverage`. Read `_dbs_burden.txt` directly for DBS burden/mutation-count figures.

---

## SBS96 matrix files

Derived from `{output}` by stripping a trailing `.txt` and appending a suffix; written unconditionally alongside the main table. Each is a tab-separated matrix with one row per SBS96 trinucleotide context (sorted by context label) and one column per input sample (column header = sample name), plus a leading `MutationType` column repeating the row label.

| File | Column values pulled from each sample's `_sbs_96_corrected.txt` |
| --- | --- |
| `{output_prefix}_SBS96_uncorrected.txt` | `mutation_number_uncorrected` |
| `{output_prefix}_SBS96_corrected.txt` | `mutation_number_corrected` |
| `{output_prefix}_SBS96_genome.txt` | `mutation_number_genome` (the real, per-context genome-extrapolated count — see `estimate_outputs.md`) |

These three files are directly usable as multi-sample input to mutational-signature tools (SigProfiler, MutationalPatterns, etc.), one file per count basis (raw observed, trinucleotide-corrected, or genome-extrapolated).
