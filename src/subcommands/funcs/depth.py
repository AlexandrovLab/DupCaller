import numpy as np
import pysam
from .misc import getAlignmentObject as BAM


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
                # if pileupread.alignment.query_name
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
    # print(f"calculate depth time:{(time.time()-start_time)/60}")
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

def prepareAlignMask(bam, chrom, start, end, params):
    """
    sum_nm = np.zeros(end - start)
    count_nm = np.zeros(end - start,dtype=int)
    bam=BAM(bam,"rb")
    for rec in bam.fetch(chrom,start,end):
        if (
            rec.is_supplementary
            or rec.is_secondary
            or rec.is_duplicate
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        # If 5 prime is soft clipped
        if (rec.is_forward and rec.cigartuples[0][0] == 4) or (
            rec.is_reverse and rec.cigartuples[-1][0] == 4
        ):
            continue
        #if rec.mapq <=40 : continue
        #if rec.mapping_quality <= params["mapq"]: continue
        id_length = 0
        id_num = 0
        for cigar in rec.cigartuples:
            if cigar[0] == 1 or cigar[0] == 2:
                id_length += cigar[1]
                id_num += 1
        NM_no_id = rec.get_tag("NM") - id_length + id_num
        sum_nm[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += NM_no_id
        count_nm[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += 1

    
    sum_nm[count_nm == 0] = 200
    count_nm[count_nm == 0] = -1

    avg_nm = sum_nm/count_nm
    avg_nm[avg_nm < 0] = np.inf
    for nn in range(avg_nm.size):
        if avg_nm[nn] != np.inf:
            print(avg_nm[nn],nn+start)
    """
    sum_nm_f = np.zeros(end - start)
    count_f = np.zeros(end - start,dtype=int)
    sum_nm_r = np.zeros(end - start)
    count_r = np.zeros(end - start,dtype=int)
    
    sum_asxs_f = np.zeros(end - start)
    sum_asxs_r = np.zeros(end - start)
    
    max_nm_f = np.zeros(end - start)
    max_nm_r = np.zeros(end - start)
    min_asxs_f = np.zeros(end - start)
    min_asxs_r = np.zeros(end - start)

    bam=BAM(bam,"rb")
    for rec in bam.fetch(chrom,start,end):
        if (
            rec.is_supplementary
            or rec.is_secondary
            or rec.is_duplicate
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        if (rec.is_forward and rec.cigartuples[0][0] == 4) or (
            rec.is_reverse and rec.cigartuples[-1][0] == 4
        ):
            continue
        ##Check covered base
        #qualities_pass_with_indel = np.zeros(rec.reference_end-rec.reference_start,dtype=bool)
        #qualities_pass = (np.array(rec.query_alignment_qualities) >= params["minBq"])



        
        ##NM average
        id_length = 0
        id_num = 0
        """
        current_seq_ind = 0
        reference_ind = 0
        for ct in rec.cigartuples:
            if ct[0] == 0:
                qualities_pass_with_indel[
                    reference_ind : reference_ind + ct[1]
                ] = qualities_pass[current_seq_ind : current_seq_ind + ct[1]]
                current_seq_ind += ct[1]
                reference_ind += ct[1]
            elif ct[0] == 1:
                current_seq_ind += ct[1]
            elif ct[0] == 2:
                qualities_pass_with_indel[
                    reference_ind : reference_ind + ct[1]
                ] = False
                reference_ind += ct[1]
            else:
                current_seq_ind += ct[1]

        qualities_pass_trimmed = qualities_pass_with_indel[max(start - rec.reference_start,0):min(end - rec.reference_start,end-start)]
        """    
        
        NM_no_id = rec.get_tag("NM") - id_length + id_num
        asxs = rec.get_tag("AS") - rec.get_tag("XS")
        if rec.is_forward:
            """
            sum_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += NM_no_id
            sum_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += asxs
            count_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += 1
            
            max_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.maximum(max_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],NM_no_id)
            min_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.minimum(min_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],asxs)
            """
            sum_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += NM_no_id
            sum_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += asxs
            count_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += 1
            
            max_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] = \
                np.maximum(max_nm_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)],NM_no_id)
            min_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] = \
                np.minimum(min_asxs_f[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)],asxs)
        
        if rec.is_reverse:
            """
            sum_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += NM_no_id
            sum_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += asxs
            count_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] += 1

            max_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.maximum(max_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],NM_no_id)
            min_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed] = \
                np.minimum(min_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)][qualities_pass_trimmed],asxs)
            """
            sum_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += NM_no_id
            sum_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += asxs
            count_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] += 1
            
            max_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] = \
                np.maximum(max_nm_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)],NM_no_id)
            min_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)] = \
                np.minimum(min_asxs_r[max(rec.reference_start - start,0):min(rec.reference_end - start,end-start)],asxs)
    sum_nm_f[count_f > 1] = sum_nm_f[count_f > 1] - max_nm_f[count_f > 1]
    sum_nm_r[count_r > 1] = sum_nm_r[count_r > 1] - max_nm_r[count_r > 1]
    sum_asxs_f[count_f > 1] = sum_asxs_f[count_f > 1] - min_asxs_f[count_f > 1]
    sum_asxs_r[count_r > 1] = sum_asxs_r[count_r > 1] - min_asxs_r[count_r > 1]
    count_f[count_f > 1] -= 1
    count_r[count_r > 1] -= 1
    
    sum_nm_f[count_f == 0] = 200
    sum_asxs_f[count_f == 0] = -np.inf
    count_f[count_f == 0] = -1
    sum_nm_r[count_r == 0] = 200
    sum_asxs_r[count_r == 0] = -np.inf
    count_r[count_r == 0] = -1

    avg_nm_f = sum_nm_f/count_f
    avg_nm_r = sum_nm_r/count_r
    avg_nm = np.max(np.vstack([avg_nm_f,avg_nm_r]),axis=0)
    avg_nm[avg_nm < 0] = np.inf 

    avg_asxs_f = sum_asxs_f/count_f
    avg_asxs_r = sum_asxs_r/count_r
    avg_asxs = np.min(np.vstack([avg_asxs_f,avg_asxs_r]),axis=0)
    avg_asxs[avg_asxs == np.inf] = 0
    """
    for nn in range(avg_asxs.size):
        if avg_asxs[nn] != 0:
            print(avg_asxs[nn],avg_asxs_f[nn],avg_asxs_r[nn],nn+start)
        if avg_nm[nn] != np.inf:
            print(avg_nm[nn],avg_nm_f[nn],avg_nm_r[nn],nn+start)
    """

    align_mask = np.logical_or(avg_nm >= params["maxNM"]/2,avg_asxs <= params["minMeanASXS"])
    return align_mask


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
