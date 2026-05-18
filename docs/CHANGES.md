# Changelog

## 1.1.2 (2026-05-18)

### Bug fixes

- **`call` output directory creation:** `DupCaller.py call -o path/to/sample` now consistently creates the full output directory for any path. Previously, directory creation was silently skipped when the `-o` value contained a `/`, causing output files to land in the wrong location and breaking downstream `estimate` and `summarize` calls.
- **`estimate` internal path resolution:** several file paths inside `estimate` (`_stats.txt`, `_indel.vcf`, `_coverage.bed.gz`) were constructed using only the sample basename, causing `FileNotFoundError` when the `-i` prefix contained a directory component (e.g. `-i results/sample1`). All internal paths now use the full prefix.
- **`estimate` SBS output overhaul:**
  - `{sample}_sbs_96_corrected.txt` now contains five columns — `mutation_number_uncorrected`, `mutation_number_corrected`, `correction_ratio`, `mutation_number_genome`, `trinuc_number_genome` — instead of the previous single `number` column. This makes the file directly usable for signature analysis without post-processing.
  - The SBS96 visualization is now saved as a multi-page PDF (`{sample}_sbs_96.pdf`, 3 pages: uncorrected counts, corrected counts, estimated mutations per genome) instead of a single-panel PNG. The same change applies to the re-estimation output (`{sample}_sbs_96_corrected_re_estimate.pdf`).
  - `{sample}_sbs_burden.txt` field names updated: `mutation number per genome` is now split into `Corrected mutation number` and `Mutation number per genome`; `genome coverage` is now `Genome coverage`.
  - The re-estimation burden calculation (`-rb` option) is now computed with the same trinucleotide-corrected method as the main estimate path, making the two results directly comparable.
- **`summarize` internal path resolution:** same class of bug — file lookup and sample label in the output table now correctly use the basename of the `-i` path rather than the full path string.
- **`genotypeDSSnv` log-domain guard:** changed `LR_diff < 0` to `LR_diff <= 0` before `log10`, preventing a divide-by-zero when the likelihood ratio difference is exactly zero.

### New features

- **Reference index validation:** `call` and `estimate` now verify that the h5 index files (`*.ref.h5`, `*.tn.h5`, `*.hp.h5`) are readable and structurally valid before starting analysis. If a file is corrupt or was built by an incompatible version of DupCaller, the tool exits immediately with a clear error message and the exact `DupCaller.py index` command needed to fix the problem.
