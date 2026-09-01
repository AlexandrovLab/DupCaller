# DupCaller Output File Reference

This folder documents every file produced by DupCaller in detail.

| Document | Contents |
| --- | --- |
| [call_outputs.md](call_outputs.md) | All files written by `DupCaller.py call` |
| [estimate_outputs.md](estimate_outputs.md) | All files written by `DupCaller.py estimate` |
| [summarize_outputs.md](summarize_outputs.md) | The cross-sample table written by `DupCaller.py summarize` |

---

## Quick File Index

Paths are relative to the sample's output directory (`-o` for `call`, `-i` for `estimate`); `{sample}` is that directory's basename.

### From `DupCaller.py call`

| File | Description |
| --- | --- |
| `SBS/{sample}_sbs.vcf`, `_sbs_fail.vcf` | Called SNVs/MNVs, PASS-only and everything-else |
| `INDEL/{sample}_indel.vcf`, `_indel_fail.vcf` | Called short indels, PASS-only and everything-else |
| `DBS/{sample}_dbs.vcf`, `_dbs_fail.vcf` | Called dinucleotide substitutions, PASS-only and everything-else |
| `SBS/{sample}_sbs_flt.vcf` | Dilute-mode filtered SNVs (conditional, `--dilute`) |
| `{sample}_coverage.bed.gz` | Per-position, per-alt-base/indel-category duplex coverage (23 columns) |
| `{sample}_stats.txt` | Library and calling quality metrics (three lines appended later by `estimate`) |
| `{sample}_call_params.log` | Resolved `call` parameters and error-file paths for this run |
| `SBS/{sample}_trinuc_by_duplex_group.txt` | Raw 192-row trinucleotide-context coverage by duplex group |
| `INDEL/{sample}_indel_by_duplex_group.txt` | Raw 100-row indel-opportunity coverage by duplex group |
| `DBS/{sample}_dbs_by_duplex_group.txt` | Raw 144-row DBS-opportunity coverage by duplex group |
| `{sample}_duplex_family_strand_composition.txt` | Read-family strand-composition table |
| `{sample}_duplex_family_strand_composition_heatmap.pdf` | Strand-composition heatmap |
| `ERROR/{sample}.amp.tn.txt` | Raw learned amplification SBS mismatch profile (diagnostic) |
| `ERROR/{sample}.amp.tn.srd.txt` | Learned amplification SBS error rates actually used for calling |
| `ERROR/{sample}.dmg.tn.txt` | Learned damage SBS error rates |
| `ERROR/{sample}.amp.hp.txt` / `.amp.str.txt` | Learned amplification indel error rates (homopolymer / STR) |
| `ERROR/{sample}.dmg.hp.txt` / `.dmg.str.txt` | Learned damage indel error rates (homopolymer / STR) |

### From `DupCaller.py estimate`

| File | Description |
| --- | --- |
| `SBS/{sample}_sbs_burden.txt` | SNV mutational burden with confidence intervals |
| `INDEL/{sample}_indel_burden.txt` | Indel mutational burden with confidence intervals |
| `DBS/{sample}_dbs_burden.txt` | DBS mutational burden with confidence intervals |
| `SBS/{sample}_sbs_96_corrected.txt`, `SBS_96_plots_{sample}.pdf` | Per-context SBS96 mutation counts and profile plot |
| `INDEL/{sample}_indel_83_corrected.txt`, `ID_83_plots_{sample}.pdf` | Per-channel ID83 mutation counts and profile plot |
| `DBS/{sample}_dbs_78_corrected.txt`, `DBS_78_plots_{sample}.pdf` | Per-channel DBS78 mutation counts and profile plot |
| `{type}/{sample}_{sbs\|indel\|dbs}_burden_by_group_size.txt`/`.pdf` | Burden stratified by minimum/exact duplex group size, per mutation type |
| `{sample}_duplex_allele_counts.txt` | Duplex depth and allele counts per unique mutation |
| `{sample}_gene_coverage.txt` | Mean per-base SBS and indel duplex depth per gene (conditional, `-gb`) |
| `{sample}_estimate_params.log` | Resolved `estimate` parameters for this run |
| `SBS/{sample}_sbs_burden_re_estimate.txt` | Re-estimated SNV burden for a sub-region (conditional, `-rb`) |
| `INDEL/{sample}_indel_burden_re_estimate.txt` | Re-estimated indel burden for a sub-region (conditional, `-rb`) |
| `SBS/{sample}_sbs_96_corrected_re_estimate.txt`, `SBS_96_plots_{sample}_re_estimate.pdf` | Re-estimated SBS96 counts and profile plot (conditional, `-rb`) |

### From `DupCaller.py summarize`

| File | Description |
| --- | --- |
| `{output}` | Cross-sample table: read/coverage/burden summary, one row per input sample |
| `{output_prefix}_SBS96_uncorrected.txt` | Multi-sample matrix of raw SBS96 counts |
| `{output_prefix}_SBS96_corrected.txt` | Multi-sample matrix of trinucleotide-corrected SBS96 counts |
| `{output_prefix}_SBS96_genome.txt` | Multi-sample matrix of genome-extrapolated SBS96 counts |
