# DupCaller Output File Reference

This folder documents every file produced by DupCaller in detail.

| Document | Contents |
| --- | --- |
| [call_outputs.md](call_outputs.md) | All files written by `DupCaller.py call` |
| [estimate_outputs.md](estimate_outputs.md) | All files written by `DupCaller.py estimate` |

---

## Quick File Index

### From `DupCaller.py call`

| File | Description |
| --- | --- |
| `{sample}_snv.vcf` | Called SNVs and MNVs |
| `{sample}_indel.vcf` | Called short indels |
| `{sample}_snv_flt.vcf` | Dilute-mode filtered SNVs (conditional) |
| `{sample}_coverage.bed.gz` | Per-position duplex coverage (SNV and indel) |
| `{sample}_stats.txt` | Library and sequencing quality metrics |
| `{sample}_trinuc_by_duplex_group.txt` | Trinucleotide context coverage by duplex group |
| `{sample}_duplex_family_strand_composition.txt` | Read family strand composition table |
| `{sample}_duplex_family_strand_composition_heatmap.pdf` | Strand composition heatmap |
| `{sample}_duplex_allele_counts.txt` | Per-mutation duplex and BAM allele counts |
| `{sample}.amp.tn.txt` | Learned amplification SBS error rates |
| `{sample}.amp.id.txt` | Learned amplification indel error counts |
| `{sample}.dmg.tn.txt` | Learned damage SBS error rates |
| `{sample}.dmg.id.txt` | Learned damage indel error counts |

### From `DupCaller.py estimate`

| File | Description |
| --- | --- |
| `{sample}_sbs_burden.txt` | SNV mutational burden with confidence intervals |
| `{sample}_indel_burden.txt` | Indel mutational burden with confidence intervals |
| `{sample}_sbs_96_corrected.txt` | Per-context SBS96 mutation counts |
| `{sample}_sbs_96.pdf` | SBS96 profile plots (3 pages) |
| `{sample}_sbs_burden_by_min_read_group_size.txt` | Burden stratified by minimum duplex group size |
| `{sample}_sbs_burden_by_min_read_group_size.png` | Burden-vs-min-group-size plot |
| `{sample}_duplex_allele_counts.txt` | Duplex depth and allele counts per mutation |
| `{sample}_gene_coverage.txt` | Mean duplex depth per gene (conditional) |
| `{sample}_sbs_burden_re_estimate.txt` | Re-estimated SNV burden for a sub-region (conditional) |
| `{sample}_indel_burden_re_estimate.txt` | Re-estimated indel burden for a sub-region (conditional) |
| `{sample}_sbs_96_corrected_re_estimate.txt` | Re-estimated SBS96 counts (conditional) |
| `{sample}_sbs_96_corrected_re_estimate.pdf` | Re-estimated SBS96 profile plots (conditional) |
