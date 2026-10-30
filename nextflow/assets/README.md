Empty placeholder files bound to DupCaller.nf's optional `path` inputs when
the matching parameter is unset (Nextflow can only stage files that exist).
DupCaller.nf detects them by their exact names and leaves the matching
DupCaller.py option off (`NO_NOISE_MASK` instead by the empty list of mask
names, so `-m` is omitted). Do not delete or rename them. `NO_NORMAL_BAM`
(and `.bai`) stands in for the matched normal with `tumor_only = true`, and
`NO_BARCODE_LIST` for an unset `barcode_list`.

`premade/<genome>.tsv` lists the files `premade = "<genome>"` downloads
(columns `file`, `url`, `md5`; `md5` empty where the source publishes
none). File names must stay as DupCaller.nf's `premadeSpec()` expects.
