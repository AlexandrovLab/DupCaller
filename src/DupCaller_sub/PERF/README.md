# Vendored PERF

This directory is a copy of the `PERF/` package from
[PERF](https://github.com/rkmlab/perf) (Pattern-based Exhaustive Repeat
Finder) by Divya Tej Sowpati & Akshay Kumar Avvaru, Lab of Dr. Rakesh Mishra,
distributed under the MIT license (see `LICENSE` in this directory).

- Upstream commit: `e343a1fa437033afce5b5a794079230530619983` (PERF v0.4.6)
- The Python sources and `lib/` assets are unmodified.
- The only change is in packaging: upstream's `setup.py` pins
  `biopython==1.69`, which conflicts with DupCaller's own biopython pin. This
  copy is installed as part of DupCaller instead, using DupCaller's biopython
  version, and exposed as the same `PERF` command.

Citation: Avvaru, A. K., Sowpati, D. T. & Mishra, R. K. PERF: an exhaustive
algorithm for ultra-fast and efficient identification of microsatellites from
large DNA sequences. Bioinformatics 34, 943-948 (2018).
https://doi.org/10.1093/bioinformatics/btx721
