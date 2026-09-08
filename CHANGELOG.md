# Changelog

All notable changes to DupCaller are recorded here, most recent first.

## [1.2.0-dev] - Unreleased

### Fixed
- `call.py`: Skip germline BED/VCF records with `ALT="."` (non-variant rows, e.g. reference-confirmation rows some multi-strain-merged VCFs carry) during germline masking, instead of crashing on `len(rec.alts)`.
- `learn.py`: Fix a crash in `profileTriNucMismatches`'s per-read SBS amp-error BQ-histogram accumulation when a read's only "mismatching" position (relative to its antimask-surviving span) is an N/deletion/off-read call that happens to carry a non-zero base quality (survives the `minBq` zeroing). Its base index falls outside the 4-base ATCG range the histogram's `bincount` is reshaped into, so `.reshape([4, 96, NUM_BQ])` raised `ValueError: cannot reshape array of size ... `. Now excluded from that accumulation the same way a zero-quality N already was, instead of vetoing the whole position for every read in the family (the F1R2/F2R1 antimask's `seq_mat == 4` exclusion, disabled since the 2024 source-layout refactor, stays disabled -- this fix is scoped to the actual crash instead).

### Changed
- `nextflow/DupCaller.nf`, `nextflow/nextflow.config`: Bump the pinned `yuhecheng62/dupcaller` container image from `1.1.0-amd64` to `1.2.0-dev` for all DupCaller-running processes (INDEX_REFERENCE, TRIM_BARCODES, CALL_VARIANTS, ESTIMATE_BURDEN).

### Notes
- Added a regression test (`tests/unit/test_profile_trinuc_mismatches.py`) reproducing the crash above.
