# Changelog

All notable changes to DupCaller are recorded here, most recent first.

## [1.2.0-dev] - Unreleased

### Fixed
- `call.py`: Skip germline BED/VCF records with `ALT="."` (non-variant rows, e.g. reference-confirmation rows some multi-strain-merged VCFs carry) during germline masking, instead of crashing on `len(rec.alts)`.
- `learn.py`: Fix a crash in `profileTriNucMismatches`'s per-read SBS amp-error BQ-histogram accumulation when a read's only "mismatching" position (relative to its antimask-surviving span) is an N/deletion/off-read call that happens to carry a non-zero base quality (survives the `minBq` zeroing). Its base index falls outside the 4-base ATCG range the histogram's `bincount` is reshaped into, so `.reshape([4, 96, NUM_BQ])` raised `ValueError: cannot reshape array of size ... `. Now excluded from that accumulation the same way a zero-quality N already was, instead of vetoing the whole position for every read in the family (the F1R2/F2R1 antimask's `seq_mat == 4` exclusion, disabled since the 2024 source-layout refactor, stays disabled -- this fix is scoped to the actual crash instead).
- `nextflow/DupCaller.nf`: `ESTIMATE_BURDEN` now sets `SIGPROFILERPLOTTING_VOLUME` to a task-local directory before calling `DupCaller.py estimate`. `sigProfilerPlotting`'s `plotSBS()` otherwise tries to cache a template pickle inside its own site-packages install directory, which fails under a container's read-only root filesystem (confirmed end-to-end against the published `yuhecheng62/dupcaller:1.2.0-dev` image via singularity: `PermissionError: ... sigProfilerPlotting/templates/`).
- `nextflow/DupCaller.nf`: `CALL_VARIANTS`'s noise-mask input (`params.noise_mask`) now accepts a single path or a list of paths, matching `DupCaller.py call`'s own `-m`/`--noise` (`nargs="+"`). Previously only a single file could be wired through, so the pipeline couldn't express a real invocation that passes both a SNP mask and a noise mask together (as the production benchmark runs do) -- confirmed end-to-end against a real sample (see Notes).

### Changed
- `nextflow/DupCaller.nf`, `nextflow/nextflow.config`: Bump the pinned `yuhecheng62/dupcaller` container image from `1.1.0-amd64` to `1.2.0-dev` for all DupCaller-running processes (INDEX_REFERENCE, TRIM_BARCODES, CALL_VARIANTS, ESTIMATE_BURDEN).

### Added
- `nextflow/pipeline.config` is now tracked in git (it previously wasn't, despite `nextflow/README.md`'s Quick Start depending on it -- a fresh clone had no way to get it).
- `nextflow/examples/run_dupcaller_sample.sh`: a platform-agnostic script that runs `DupCaller.nf` end-to-end for one tumor/normal sample given just a sample ID, 4 fastq paths, and a reference -- no scheduler assumptions, works with `-P docker,local` on a workstation or `-P singularity,local` on any HPC node.
- `nextflow/examples/run_sample.slurm.sh`: a thin SLURM wrapper around the above carrying only this cluster's scheduler/account/reference/mask specifics (`platinum`/`hcp-ddp302`/`ddp302`, matching the real `PD*.2pass.sl` resource convention).
- Both validated end-to-end (2026-09-08) against a 100k-read-pair subsample of a real production sample (PD43276), real hg38 reference, and the real production mask pair -- the resolved `DupCaller.py call` command matched the actual `PD43276.2pass.sl` production script argument-for-argument (aside from the deliberately-scoped `-r`/`-p` for a fast validation run). Also confirmed the trim step is byte-identical between the containerized and direct-CLI invocation on the same real fastq data.

### Notes
- Added a regression test (`tests/unit/test_profile_trinuc_mismatches.py`) reproducing the crash above.
- Added `tests/mock_pipeline/test_mock_pipeline_nextflow.py` (with `run_nextflow_pipeline.sh` and `nextflow_test.config`): runs `nextflow/DupCaller.nf` end-to-end against the synthetic dataset in `tests/mock_pipeline/data/` through Nextflow + a container runtime (docker if reachable, else singularity), verifying the containerized pipeline itself -- staging, image pulls, per-process resource directives -- works from fastq through SBS/indel/DBS calls and burden estimation. Skips automatically if nextflow/bwa/samtools or a usable docker/singularity aren't available. Reused mock_1/2.fastq as both tumor and normal (the pipeline has no tumor-only mode), so its calls aren't expected to numerically match `test_mock_pipeline.py`'s tumor-only run.
