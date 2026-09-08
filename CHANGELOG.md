# Changelog

All notable changes to DupCaller are recorded here, most recent first.

## [1.2.0-dev] - Unreleased

### Fixed
- `call.py`: Skip germline BED/VCF records with `ALT="."` (non-variant rows, e.g. reference-confirmation rows some multi-strain-merged VCFs carry) during germline masking, instead of crashing on `len(rec.alts)`.
- `learn.py`: Restore exclusion of N-base/deletion/off-read positions from the per-strand (F1R2/F2R1) antimask used to learn amp/damage error rates. This check was accidentally disabled (commented out) during the 2024 source-layout refactor and has been dormant since.
