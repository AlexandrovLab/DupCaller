# Changes Summary for Dev Branch

## Overview
Total changes: 11 files modified/added (2013 insertions, 215 deletions)

## New Files Added

### 1. nextflow/DupCaller.nf (436 lines)
- Complete Nextflow DSL2 pipeline implementation for DupCaller
- Processes for reference indexing, barcode trimming, alignment, and variant calling
- Parameters for all DupCaller features including burden estimation
- Comprehensive help documentation

### 2. nextflow/README.md (280 lines)
- Complete documentation for the Nextflow pipeline
- Usage instructions and examples
- Parameter descriptions

### 3. nextflow/nextflow.config (200 lines)
- Configuration file for Nextflow pipeline
- Process resource allocations
- Default parameter settings

### 4. src/subcommands/funcs/indels_old.py (65 lines)
- Backup of previous indels.py implementation

### 5. src/subcommands/funcs/prob_old.py (465 lines)
- Backup of previous prob.py implementation

## Modified Files

### 1. src/DupCaller.py (1 deletion)
- Removed empty line after args parsing (line 467)
- Minor formatting cleanup

### 2. src/subcommands/funcs/call.py (410 modifications)
- Major refactoring of variant calling logic
- Updated filtering and quality control mechanisms
- Improved code organization

### 3. src/subcommands/funcs/depth.py (177 additions)
- Added new depth calculation functionality
- Enhanced coverage metrics

### 4. src/subcommands/funcs/indels.py (26 modifications)
- Complete refactoring of indel detection logic
- Changed from separate alt/ref quality and count arrays to unified seqArr/qualArr
- New return signature: `return seqArr, qualArr` instead of `return altQualArr, refQualArr, altCountArr, refCountArr`
- Simplified logic for insertion and deletion detection
- Added handling for ambiguous cases with seqArr = -1

Key changes:
- Line 28-32: Commented out old array definitions, added new seqArr and qualArr
- Line 45-50: Changed insertion handling to set seqArr and qualArr
- Line 60-70: Changed deletion handling to set seqArr and qualArr
- Line 72: New return statement

### 5. src/subcommands/funcs/learn.py (13 modifications)
- Updated to use new indels.py interface
- Added F1R2_trinuc_seq_err_count_mat initialization (line 203)
- Changed getIndelArr() usage to handle new return signature:
  - Old: `aq, rq, ac, rc = getIndelArr(seq, indels_masked)`
  - New: `seqArr, qualArr = getIndelArr(seq, indels_masked)`
- Added calculations to derive counts from seqArr:
  - `ac = np.count_nonzero(seqArr==1)` (alt count)
  - `rc = np.count_nonzero(seqArr==0)` (ref count)
  - `aq = np.sum(qualArr[seqArr==1])` (alt quality sum)
  - `rq = np.sum(qualArr[seqArr==0])` (ref quality sum)
- Applied to both F1R2 and F2R1 read processing (lines 295-310)

### 6. src/subcommands/funcs/prob.py (155 modifications)
- Updated probability calculation logic
- Improved statistical modeling
- Enhanced filtering mechanisms

## Summary of Changes by Category

### Pipeline Infrastructure
- Complete Nextflow pipeline implementation (3 new files)
- Production-ready workflow automation

### Indel Processing Refactor
- Simplified indel detection interface
- Changed from 4 return values to 2 (seqArr, qualArr)
- More efficient memory usage
- Better handling of ambiguous indel cases

### Code Organization
- Created backup files for major refactors (indels_old.py, prob_old.py)
- Maintained code history while implementing new versions

### Quality Improvements
- Enhanced variant calling logic (call.py)
- Improved depth calculations (depth.py)
- Updated probability models (prob.py)
- Minor code cleanup (DupCaller.py)

## Git Commit History Context
Recent commits show focus on:
- a868bc5: bypass zero duplex depth variants
- 11a46f0: Update README with maxZeroQualFrac parameter documentation
- c9e7084: remove unnecessary base quality filters
- cd0c447: remove unnecessary base quality filters
- 73d1d4c: Revert "Improve code formatting and VAF imbalance filtering"
