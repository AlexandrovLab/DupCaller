"""Write small bgzip-compressed, tabix-indexed optional resources for the
mock contig (data/reference.fa), so run_nextflow_pipeline.sh can exercise
every optional DupCaller.nf input: germline VCF, noise/SNP masks, target
BED, indel BED, gene BED.

Usage: python make_test_resources.py OUTDIR REFERENCE_FASTA

REFERENCE_FASTA is the run's faidx-indexed copy of data/reference.fa, so
nothing is written into the source tree.

Intervals stay clear of the truth mutations (data/truth.txt: SBS at 1000,
deletion at 2053) except the target BED, which spans the whole contig, and
the gene BED, whose exons cover the truth sites so per-gene coverage is
non-zero.
"""

import os
import sys

import pysam

CHROM = "mockchr1"
LENGTH = 3200

BEDS = {
    "target.bed": [(0, LENGTH, "target")],
    "indel_pon.bed": [(100, 200, "pon")],
    "genes.bed": [(900, 1100, "GENEA_exon1"), (2000, 2100, "GENEB_exon1")],
    "snp_mask.bed": [(300, 310, "snp")],
    "noise_mask.bed": [(400, 410, "noise")],
}


def main():
    outdir, reference = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    for name, rows in BEDS.items():
        path = os.path.join(outdir, name)
        with open(path, "w") as f:
            for start, end, label in rows:
                f.write(f"{CHROM}\t{start}\t{end}\t{label}\n")
        pysam.tabix_index(path, preset="bed", force=True)

    ref = pysam.FastaFile(reference)
    pos = 500  # 1-based
    ref_base = ref.fetch(CHROM, pos - 1, pos)
    alt = "T" if ref_base != "T" else "G"
    vcf = os.path.join(outdir, "germline.vcf")
    with open(vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CHROM},length={LENGTH}>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write(f"{CHROM}\t{pos}\t.\t{ref_base}\t{alt}\t.\tPASS\tAF=0.4\n")
    pysam.tabix_index(vcf, preset="vcf", force=True)


if __name__ == "__main__":
    main()
