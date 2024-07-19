import numpy as np
import pysam


def extractDepthSnv(bam, chrom, pos, ref, alt, params, minbq=1):
    altAlleleCount = 0
    refAlleleCount = 0
    otherAlleleCount = 0
    indelAlleleCount = 0
    processed_read_names = dict()
    for pileupcolumn in bam.pileup(
        chrom, pos - 1, pos, min_base_quality=minbq, truncated=True
    ):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    or processed_read_names.get(pileupread.alignment.query_name)
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                ):
                    continue
                processed_read_names[pileupread.alignment.query_name] = 1
                if pileupread.alignment.is_duplicate:
                    continue
                if pileupread.is_del or pileupread.indel != 0:
                    indelAlleleCount += 1
                    otherAlleleCount += 1
                elif (
                    pileupread.alignment.query_sequence[pileupread.query_position]
                    == alt
                ):
                    altAlleleCount += 1
                elif (
                    pileupread.alignment.query_sequence[pileupread.query_position]
                    == ref
                ):
                    refAlleleCount += 1
                else:
                    otherAlleleCount += 1
    depth = refAlleleCount + altAlleleCount + otherAlleleCount
    return altAlleleCount, refAlleleCount, indelAlleleCount, depth


def extractDepthIndel(bam, chrom, pos, ref, alt, params, minbq=1):
    indel_size = len(alt) - len(ref)
    altAlleleCount = 0
    refAlleleCount = 0
    otherAlleleCount = 0
    otherIndelCount = 0
    processed_read_names = dict()
    for pileupcolumn in bam.pileup(
        chrom, pos - 1, pos, min_base_quality=minbq, truncate=True
    ):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_del
                    or pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    # or pileupread.alignment.is_duplicate
                    or processed_read_names.get(pileupread.alignment.query_name)
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                ):
                    continue
                processed_read_names[pileupread.alignment.query_name] = 1
                if pileupread.alignment.is_duplicate:
                    continue
                if pileupread.indel == indel_size:
                    altAlleleCount += 1
                elif pileupread.indel == 0:
                    refAlleleCount += 1
                else:
                    otherAlleleCount += 1
                    if pileupread.indel != 0:
                        otherIndelCount += 1
            break
    depth = refAlleleCount + altAlleleCount + otherAlleleCount
    return altAlleleCount, refAlleleCount, otherIndelCount, depth


def extractDepthRegion(bam, chrom, start, end, params):
    depth = np.zeros(end - start)
    indelmask = np.zeros(end - start, dtype=bool)
    # processed_read_names = {}
    mapq = params["mapq"]
    depth = np.zeros(end - start)
    max_depth = params["minNdepth"]
    # for line in pysam.depth("-q","30","-Q",f"{mapq}","-J","-r",f"{chrom}:{start}-{end}",bam).split("\n"):
    for line in pysam.mpileup(
        "-Q", "30", "-q", f"{mapq}", "-r", f"{chrom}:{start}-{end}", bam
    ).split("\n"):
        if line:
            try:
                pos = int(line.split("\t")[1]) - 1  # 0-based position
                depth[pos - start] = int(line.split("\t")[3])
            except:
                1
    return depth, indelmask


def detectOverlapDiscord(bam, chrom, pos, ref, alt, params, bc1, bc2, start):
    discord_num = 0
    for pileupcolumn in bam.pileup(
        chrom,
        pos - 1,
        pos,
        min_base_quality=0,
        truncated=True,
        stepper="samtools",
        flag_filter=2828,
    ):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    # or pileupread.alignment.is_duplicate
                    or pileupread.alignment.has_tag("DT")
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                    or pileupread.is_del
                ):
                    continue
                read_name = pileupread.alignment.query_name
                read_bc1 = (read_name.split("_")[1].split("+"))[0]
                read_bc2 = (read_name.split("_")[1].split("+"))[1]
                qual = pileupread.alignment.query_qualities[pileupread.query_position]
                if (
                    (
                        (read_bc1 == bc1 and read_bc2 == bc2)
                        or (read_bc1 == bc2 and read_bc2 == bc1)
                    )
                    and qual == 0
                    and pileupread.alignment.reference_start == start
                ):
                    discord_num += 1
                    if discord_num >= 2:
                        return True
    return False
