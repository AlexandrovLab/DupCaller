# Changelog

All notable changes to DupCaller are recorded here, most recent first.

## [1.2.1-dev] - Unreleased

### Fixed
- `learn.py`: `estimate_sbs_srd_rates`'s EM M-step denominator changed from `N + 3*pseudocount` to `N + 4*pseudocount`, matching a proper symmetric Dirichlet(a,a,a,a) prior over all 4 categories (reference + 3 alt bases) instead of smoothing only the 3 alt categories and leaving the residual reference-rate unsmoothed (it could drift toward exactly 0 for a high-base-quality/near-zero-error context). With the corrected denominator, a zero-observation trinuc context now yields the uniform 1/4 directly from the general formula.
- `misc.py`/`Caller.py`: Removed `regularizeErrorMat`, a flat post-hoc additive floor (`+1e-6`/`+1e-8`) applied to every cell of the SBS/indel error matrices after normalization. It was silently dominating/overriding genuinely small EM-fitted rates. Replaced everywhere by Dirichlet pseudocount regularization applied at the count level (`(count + a) / (total + n*a)`), consistent with how `estimate_sbs_srd_rates` already regularizes internally. `-a`/`--pseudocount` now also controls this smoothing (previously hardcoded `1e-6`/`1e-8`).
- `learn.py`: `profileTriNucMismatches`'s SRD (single-read-damage) site-inclusion antimask now uses BQ-qualifying counts (`F1R2_hq_count_mat`/`F1R2_hq_ref_count`) consistently across every check, instead of mixing them with the raw, BQ-blind `F1R2_count_mat` in the "≥2 distinct alleles, else exclude if ref_count==0" fallback. A low-base-quality base is now invisible to every step of the antimask, not just some.
- `funcs/call.py`: `callBam`'s NM-blacklist family filtering restored to the original fractional/aggregate check (skip a family only if ≥50% of its reads, or one entire strand, are individually NM-blacklisted) after briefly trying a per-read filter; the per-strand `F1R2`/`F2R1` counts used to build `duplex_no` are now recomputed from `readSet` immediately after the `all_dup` check instead of reusing a stale count from earlier in the loop.
- `Caller.py`: `merge_and_combine_coverage_files`'s boundary-merge step (peeling overlap rows between adjacent workers' coverage files) is now wrapped in the same try/except shield already used for the BGZF-combine and tabix-index steps -- a corrupt, truncated, or unexpectedly multi-chromosome boundary file now degrades gracefully (warns and leaves per-worker files in place for diagnosis/retry) instead of crashing the whole run and losing completed calling work.
- `Caller.py`: Replaced a hand-copied BGZF EOF magic-byte hex literal in `_concatenate_bgzf_files` with Biopython's own `bgzf._bgzf_eof` constant (verified byte-identical), removing a duplicate, driftable copy of the same value.
- `misc.py`: Extracted a single `_cigar_ref_consumed` helper (M/D/N/=/X reference-consuming ops) shared by `_mate_reference_end` and `determineTrimLength`; the latter previously hand-rolled a narrower M/D-only regex that silently missed N/=/X CIGAR operations.
- `Caller.py`/`Learn.py`: Fixed a `-p 1` vs `-p >1` dtype divergence -- the single-thread branches were returning `callBam`'s raw float64 count matrices unconverted, while the multi-thread branches cast the same matrices to `int` after summing. Both single-thread branches now match their file's own multi-thread branch exactly (`Caller.py`'s embedded auto-learn pass casts 7 matrices, `Learn.py`'s standalone learn command casts all 7 consistently in both branches).
- `misc.py`: `check_mate_cigar_tags` now samples only the first 100,000 fetched BAM records instead of scanning the whole file with `until_eof=True`, avoiding a slow/unnecessary full-file pass on runs that only touch a small region.
- `misc.py`: `splitBamRegions` gained a `min_chunk_length` parameter (default 10000) rejecting a candidate chunk-boundary cut too close to the previous cut or either contig end, bounding how small a worker's chunk (and therefore its risk of MC-tag overflow reaching past the next worker) can get.

### Added
- `-lo`/`--learnOnly`: stop after estimating/writing the error-rate files, skipping variant calling entirely.
- `--minChunkLength`/`--min-chunk-length` (default 10000): minimum chunk length in bases for BAM-splitting across worker processes.
- `--srdMinRead` (default 3, non-negative): minimum reads on a single strand (F1R2 or F2R1) for that strand to be included in SBS single-read-damage (SRD) rate learning, evaluated independently per strand.
- `--ssmMinRead` (default 3, non-negative): minimum reads required on both strands for a duplex family to be considered for single-strand-mutation (SSM)/damage-rate learning.
- `tests/unit/test_caller_mc_tags.py`, `tests/unit/test_learn_n_handling.py`, `tests/unit/test_learn_only.py`.

### Notes
- `--srdMinRead`/`--ssmMinRead` now documented in README.md's parameter table (previously present in `--help` only).
- Removed `test_learn_n_handling.py::test_two_ns_in_one_read_drops_the_whole_read`, which asserted the now-removed per-read `n_mismatch > 1` exclusion gate.

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
