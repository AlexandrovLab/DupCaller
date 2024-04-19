#!/usr/bin/env python3

import gzip

from pysam import FastaFile as FASTA
import time

import numpy as np
from pysam import AlignmentFile as BAM
from pysam import TabixFile as BED
from pysam import VariantFile as VCF
from Bio import SeqIO

from DCutils.funcs import *


def bamIterateMultipleRegion(bam, regions, refpath):
    bamObject = getAlignmentObject(bam, refpath=refpath)
    for region in regions:
        for rec in bamObject.fetch(*region):
            if len(region) >= 2:
                if rec.reference_start < region[1]:
                    continue
            yield rec, region


def callBam(params, processNo, chunkSize):
    # Get parameters
    bam = params["tumorBam"]
    nbam = params["normalBam"]
    regions = params["regions"]
    refpath = params["reference"]
    if params["germline"]:
        germline = VCF(params["germline"], "rb")
    else:
        germline = None
    all_chroms = [_[0] for _ in regions]
    # fasta = dict()
    # fasta = SeqIO.index(params["reference"], "fasta")
    fasta = FASTA(params["reference"])
    start_time = time.time()
    ##for ch in all_chroms:
    # fasta[ch] = whole_fasta[ch]
    # Add this record to our list
    minMapq = params["mapq"]
    mutRate = params["mutRate"]
    pcut = params["pcutoff"]
    isLearn = params.get("isLearn", False)
    nn = processNo
    output = "tmp/" + params["output"] + "_" + str(nn)
    if params["noise"]:
        noise = BED(params["noise"])
    else:
        noise = None
    if params["indel_bed"]:
        indel_bed = BED(params["indel_bed"])
    else:
        indel_bed = None
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    muts = []
    muts_dict = dict()
    muts_indels = []
    duplex_read_num_dict = dict()
    unique_read_num = 0
    pass_read_num = 0
    FPs = []
    RPs = []
    indel_dict = dict()
    mismatch_dict = dict()
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                mismatch_dict[minus_base + ref_base + plus_base] = [0, 0, 0, 0]

    # Initialize
    total_coverage = 0
    total_coverage_indel = 0
    starttime = time.time()
    tumorBam = getAlignmentObject(bam, refpath)
    if nbam:
        normalBam = getAlignmentObject(nbam, refpath)
    else:
        normalBam = None
    currentReadSet = []
    currentBc = ""
    currentStart = -1
    currentReadDict = {}
    # chromPast = 'startChrom'
    # fasta = reference
    recCount = 0
    currentCheckPoint = 1000000
    lastTime = 0
    duplex_count = 0
    reference_mat_chrom = "anyChrom"
    reference_mat_start = 0
    locus_bed = gzip.open(output + "_coverage.bed.gz", "wb")
    processed_read_names = set()
    if len(regions[0]) == 1:
        region_start = regions[0][0] + ":0"
    else:
        region_start = regions[0][0] + ":" + str(regions[0][1])
    if len(regions[-1]) != 3:
        regions_end = (
            regions[-1][0] + ":" + str(tumorBam.get_reference_length(regions[-1][0]))
        )
    else:
        regions_end = regions[-1][0] + ":" + str(regions[-1][2])

    print(f"Process {str(processNo)}: Initiated. Regions:{region_start}-{regions_end}")
    # for region in regions:
    for rec, region in bamIterateMultipleRegion(bam, regions, refpath):
        recCount += 1
        if recCount == currentCheckPoint:
            currentTime = (time.time() - starttime) / 60
            usedTime = currentTime - lastTime
            lastTime = currentTime
            print(
                f"Process {str(processNo)}: processed {str(recCount)} reads in {currentTime : .2f} minutes. Time for process last 1000000 reads:{usedTime : .2f} minutes. Current position:{rec.reference_name}:{rec.reference_start}. End Position:{regions_end}"
            )

            currentCheckPoint += 1000000
        if (
            rec.is_supplementary
            or rec.is_secondary
            or rec.has_tag("DT")
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        # if rec.get_tag("AS") - rec.get_tag("XS") <= 50:
        # continue
        pass_read_num += 1
        start = rec.reference_start
        bc = rec.query_name.split("_")[1]
        bcsplit = bc.split("+")
        bc1 = bcsplit[0]
        bc2 = bcsplit[1]
        chrom = tumorBam.get_reference_name(rec.reference_id)
        if currentStart == -1:
            currentStart = start
        if start == currentStart:
            if currentReadDict.get(bc1 + "+" + bc2) != None:
                if not currentReadDict[bc1 + "+" + bc2]["names"].get(rec.query_name):
                    currentReadDict[bc1 + "+" + bc2]["seqs"].append(rec)
                    currentReadDict[bc1 + "+" + bc2]["names"][rec.query_name] = (
                        len(currentReadDict[bc1 + "+" + bc2]["seqs"]) - 1
                    )
                    if (rec.is_forward and rec.is_read1) or (
                        rec.is_reverse and rec.is_read2
                    ):
                        currentReadDict[bc1 + "+" + bc2]["F1R2"] += 1
                    else:
                        currentReadDict[bc1 + "+" + bc2]["F2R1"] += 1
                else:
                    seq_ind = currentReadDict[bc1 + "+" + bc2]["names"][rec.query_name]
                    seq1_Bq = currentReadDict[bc1 + "+" + bc2]["seqs"][
                        seq_ind
                    ].query_alignment_qualities
                    seq1_meanBq = sum(seq1_Bq) / len(seq1_Bq)
                    seq2_Bq = rec.query_alignment_qualities
                    seq2_meanBq = sum(seq2_Bq) / len(seq2_Bq)
                    if seq2_meanBq >= seq1_meanBq:
                        currentReadDict[bc1 + "+" + bc2]["seqs"][seq_ind] = rec
            elif currentReadDict.get(bc2 + "+" + bc1) != None:
                if not currentReadDict[bc2 + "+" + bc1]["names"].get(rec.query_name):
                    currentReadDict[bc2 + "+" + bc1]["seqs"].append(rec)
                    currentReadDict[bc2 + "+" + bc1]["names"][rec.query_name] = (
                        len(currentReadDict[bc2 + "+" + bc1]["seqs"]) - 1
                    )
                    if (rec.is_forward and rec.is_read1) or (
                        rec.is_reverse and rec.is_read2
                    ):
                        currentReadDict[bc2 + "+" + bc1]["F1R2"] += 1
                    else:
                        currentReadDict[bc2 + "+" + bc1]["F2R1"] += 1
                else:
                    seq_ind = currentReadDict[bc2 + "+" + bc1]["names"][rec.query_name]
                    seq1_Bq = currentReadDict[bc2 + "+" + bc1]["seqs"][
                        seq_ind
                    ].query_alignment_qualities
                    seq1_meanBq = sum(seq1_Bq) / len(seq1_Bq)
                    seq2_Bq = rec.query_alignment_qualities
                    seq2_meanBq = sum(seq2_Bq) / len(seq2_Bq)
                    if seq2_meanBq >= seq1_meanBq:
                        currentReadDict[bc2 + "+" + bc1]["seqs"][seq_ind] = rec
            else:
                currentReadDict.update(
                    {
                        bc: {
                            "seqs": [rec],
                            "F1R2": 0,
                            "F2R1": 0,
                            "names": {rec.query_name: 0},
                        }
                    }
                )
                if (rec.is_forward and rec.is_read1) or (
                    rec.is_reverse and rec.is_read2
                ):
                    currentReadDict[bc1 + "+" + bc2]["F1R2"] += 1
                else:
                    currentReadDict[bc1 + "+" + bc2]["F2R1"] += 1

        else:
            """
            Calling block starts
            """
            for key in currentReadDict.keys():
                readSet = currentReadDict[key]["seqs"]
                mean_mapq = sum([seq.mapping_quality for seq in readSet]) / len(readSet)
                if mean_mapq < params["mapq"]:
                    continue
                meanASXS = sum(
                    [seq.get_tag("AS") - seq.get_tag("XS") for seq in readSet]
                ) / len(readSet)
                if meanASXS < 50:
                    continue
                setBc = key.split("+")
                setBc1 = setBc[0]
                setBc2 = setBc[1]
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]
                duplex_no = f"{min([F1R2,F2R1])}+{max([F1R2,F2R1])}"
                if duplex_read_num_dict.get(duplex_no) is None:
                    duplex_read_num_dict[duplex_no] = [0, 0]
                unique_read_num += 1
                if setBc1 != setBc2 and F2R1 >= 1 and F1R2 >= 1:
                    rs_reference_end = max([r.reference_end for r in readSet])
                    chromNow = readSet[0].reference_name
                    if (
                        chromNow != reference_mat_chrom
                        or rs_reference_end >= reference_mat_end
                    ):
                        ### Output coverage
                        if "coverage" in locals():
                            non_zero_positions = np.nonzero(coverage + coverage_indel)
                            for pos in non_zero_positions[0].tolist():
                                locus_bed.write(
                                    (
                                        "\t".join(
                                            [
                                                reference_mat_chrom,
                                                str(pos + reference_mat_start),
                                                str(pos + 1 + reference_mat_start),
                                                str(coverage[pos]),
                                                str(coverage_indel[pos]),
                                            ]
                                        )
                                        + "\n"
                                    ).encode("utf-8")
                                )
                            total_coverage += np.sum(coverage)
                            total_coverage_indel += np.sum(coverage_indel)
                        # if chromNow != reference_mat_chrom:
                        reference_mat_chrom = chromNow
                        # current_reference = str(fasta[reference_mat_chrom].seq)
                        reference_mat_start = readSet[0].reference_start
                        try:
                            region_end = region[2]
                        except:
                            region_end = 10e10
                        contig_end = tumorBam.get_reference_length(chromNow)
                        reference_mat_end = min(
                            readSet[0].reference_start + 100000,
                            max(
                                region_end, max([seq.reference_end for seq in readSet])
                            ),
                            contig_end,
                        )
                        current_fasta = fasta.fetch(
                            reference_mat_chrom, reference_mat_start, reference_mat_end
                        )
                        (
                            prior_mat,
                            snp_mask,
                            indel_mask,
                            noise_mask,
                            n_cov_mask,
                            ref_np,
                        ) = prepare_reference_mats(
                            reference_mat_chrom,
                            reference_mat_start,
                            reference_mat_end,
                            current_fasta,
                            germline,
                            noise,
                            indel_bed,
                            normalBam,
                            tumorBam,
                            params,
                        )
                        coverage = np.zeros(100000, dtype=int)
                        coverage_indel = np.zeros(100000, dtype=int)

                    ### Record read names to check if mate has been processed
                    processed_flag = 0
                    for seq in readSet:
                        if seq.query_name in processed_read_names:
                            processed_read_names.remove(seq.query_name)
                            processed_flag = 1
                            break
                    if processed_flag == 0:
                        processed_read_names.add(readSet[0].query_name)

                    start_ind = readSet[0].reference_start - reference_mat_start
                    reference_length_min = min(
                        [read.reference_length for read in readSet]
                    )
                    end_ind = (
                        readSet[0].reference_start
                        + reference_length_min
                        - reference_mat_start
                    )

                    reference_length_max = max(
                        [read.reference_length for read in readSet]
                    )
                    end_ind_max = (
                        readSet[0].reference_start
                        + reference_length_max
                        - reference_mat_start
                    )
                    masks = np.zeros([4, end_ind - start_ind], dtype=bool)
                    masks[0, :] = snp_mask[start_ind:end_ind]
                    masks[1, :] = noise_mask[start_ind:end_ind]
                    masks[2, :] = n_cov_mask[start_ind:end_ind]
                    left, right = determineTrimLength(
                        readSet[0], params=params, processed_flag=processed_flag
                    )
                    masks[3, :left] = True
                    masks[3, -right:] = True
                    antimask = np.all(~masks, axis=0)
                    ### If the whole reads are masked:
                    if not np.any(antimask):
                        continue
                    indel_bool = [
                        ("I" in seq.cigarstring or "D" in seq.cigarstring)
                        for seq in readSet
                    ]
                    if any(indel_bool):
                        masks_indel = np.zeros([4, end_ind_max - start_ind], dtype=bool)
                        masks_indel[0, :] = indel_mask[start_ind:end_ind_max]
                        masks_indel[1, :] = noise_mask[start_ind:end_ind_max]
                        masks_indel[2, :] = n_cov_mask[start_ind:end_ind_max]
                        left, right = determineTrimLength(
                            readSet[0], params=params, processed_flag=processed_flag
                        )
                        masks_indel[3, :left] = True
                        masks_indel[3, -right:] = True
                        antimask_indel = np.all(~masks_indel, axis=0)
                        (
                            F1R2_ARLR,
                            F2R1_ARLR,
                            indels,
                            F1R2_ref_count,
                            F1R2_alt_count,
                            F2R1_ref_count,
                            F2R1_alt_count,
                        ) = genotypeDSIndel(readSet, tumorBam, antimask_indel, params)
                        DCS = np.logical_and(
                            F1R2_ARLR >= params["pcutoff"],
                            F2R1_ARLR >= params["pcutoff"],
                        )
                        pass_inds = np.nonzero(DCS)[0].tolist()
                        indels_pass = [indels[_] for _ in pass_inds]
                        coverage_indel[start_ind:end_ind_max][antimask_indel] += 1
                        if len(indels_pass) == 1:
                            indel = indels_pass[0]
                            indel_chrom = chromNow
                            indel_pos = int(indel.split(":")[0])
                            indel_size = int(indel.split(":")[1])
                            NMs = [seq.get_tag("NM") for seq in readSet]
                            averageNM = sum(NMs) / len(NMs)
                            if averageNM - abs(indel_size) > 1:
                                continue
                            if indel_size < 0:
                                indel_ref = str(
                                    current_fasta[
                                        indel_pos
                                        - reference_mat_start : indel_pos
                                        - reference_mat_start
                                        - indel_size
                                        + 1
                                    ]
                                ).upper()
                                indel_alt = current_fasta[
                                    indel_pos - reference_mat_start
                                ].upper()
                            else:
                                indel_ref = current_fasta[
                                    indel_pos - reference_mat_start
                                ].upper()
                                indel_alt = indel_ref + indel.split(":")[2]
                            indel_str = (
                                str(indel_chrom)
                                + ":"
                                + str(indel_pos)
                                + str(indel_ref)
                                + ":"
                                + str(indel_alt)
                            )
                            # if indel_dict.get(indel_str):
                            # continue
                            """
                            ta, tr, tdp = extractDepthIndel(
                                tumorBam,
                                indel_chrom,
                                indel_pos + 1,
                                indel_ref,
                                indel_alt,
                                params,
                            )
                            if ta + tr != tdp:
                                continue
                            if ta == 0:
                                continue
                            if ta > params["maxAltCount"]:
                                continue
                            if IndelFilterByWindows(tumorBam,indel_chrom,indel_pos+1,5,params):
                                continue

                            if normalBam:
                                na, nr, ndp = extractDepthIndel(
                                    normalBam,
                                    indel_chrom,
                                    indel_pos + 1,
                                    indel_ref,
                                    indel_alt,
                                    params,
                                )
                                if ndp - nr > 0:
                                    continue
                                if ndp < params["minNdepth"]:
                                    continue
                                if IndelFilterByWindows(normalBam,indel_chrom,indel_pos+1,5,params):
                                    continue
                            else:
                                na, nr, ndp = (0, 0, 0)
                                if ta / tdp > params["maxAF"]:
                                    continue
                            """
                            indel_rec = {
                                "chrom": chromNow,
                                "pos": indel_pos + 1,
                                "ref": indel_ref,
                                "alt": indel_alt,
                                "infos": {
                                    "F1R2": int(
                                        F1R2_alt_count[pass_inds[0]]
                                        + F1R2_ref_count[pass_inds[0]]
                                    ),
                                    "F2R1": int(
                                        F2R1_alt_count[pass_inds[0]]
                                        + F2R1_ref_count[pass_inds[0]]
                                    ),
                                    "TG": F1R2_ARLR[pass_inds[0]],
                                    "BG": F2R1_ARLR[pass_inds[0]],
                                    "TC": ",".join(
                                        [
                                            str(F1R2_alt_count[pass_inds[0]]),
                                            str(F1R2_ref_count[pass_inds[0]]),
                                        ]
                                    ),
                                    "BC": ",".join(
                                        [
                                            str(F2R1_alt_count[pass_inds[0]]),
                                            str(F2R1_ref_count[pass_inds[0]]),
                                        ]
                                    ),
                                },
                                "formats": ["AC", "RC", "DP"],
                                # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_indels.append(indel_rec)
                            indel_dict[indel_str] = 1
                    else:
                        if isLearn:
                            if sum([int(seq.get_tag("NM")) for seq in readSet]) == 0:
                                continue
                        ### Calculate genotype probability
                        (
                            F1R2_ARLR,
                            F2R1_ARLR,
                            ref_int,
                            alt_int,
                            antimask,
                            F1R2_count,
                            F2R1_count,
                        ) = genotypeDSSnv(
                            readSet,
                            ref_np[start_ind:end_ind],
                            prior_mat[start_ind:end_ind, :],
                            antimask,
                            params,
                        )
                        if isLearn:
                            mismatch_now = profileTriNucMismatches(
                                readSet,
                                ref_np[start_ind:end_ind],
                                antimask,
                                params,
                            )
                            for trinuc in mismatch_dict.keys():
                                for nnn in range(4):
                                    mismatch_dict[trinuc][nnn] += mismatch_now[trinuc][
                                        nnn
                                    ]
                        ref_int = ref_np[start_ind:end_ind]
                        refs_ind = np.nonzero(
                            np.logical_and(
                                F1R2_ARLR
                                <= -params["pcutoff"] + 2 * log10(params["mutRate"]),
                                F2R1_ARLR
                                <= -params["pcutoff"] + 2 * log10(params["mutRate"]),
                            )
                        )[0].tolist()
                        muts_ind = np.nonzero(
                            np.logical_and(
                                F1R2_ARLR >= params["pcutoff"],
                                F2R1_ARLR >= params["pcutoff"],
                            )
                        )[0].tolist()
                        pass_bool = np.full(F1R2_ARLR.size, False, dtype=bool)
                        pass_bool[refs_ind] = True
                        pass_bool[muts_ind] = True
                        pos = [
                            mut_ind + start_ind + reference_mat_start
                            for mut_ind in muts_ind
                        ]
                        # muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                        mut_positions = [
                            mut_ind + start_ind + reference_mat_start + 1
                            for mut_ind in muts_ind
                        ]
                        """
                        if len(mut_positions) > 1:
                            break_flag = 0
                            for nn in range(1, len(mut_positions)):
                                if mut_positions[nn] - mut_positions[nn - 1] != 1:
                                    break_flag = 1
                            if break_flag:
                                continue
                        """
                        NMs = [seq.get_tag("NM") for seq in readSet]
                        averageNM = sum(NMs) / len(NMs)
                        if averageNM - len(mut_positions) > 1:
                            continue
                        if len(mut_positions) > params["maxMnv"]:
                            continue
                        for nn in range(len(mut_positions)):
                            mut_chrom = reference_mat_chrom
                            mut_pos = mut_positions[nn]
                            mut_ref = num2base[ref_int[muts_ind[nn]]]
                            mut_alt = num2base[alt_int[muts_ind[nn]]]
                            """
                            # if muts_dict.get(
                            # "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                            # ):
                            # continue
                            ta, tr, tdp = extractDepthSnv(
                                tumorBam, mut_chrom, mut_pos, mut_ref, mut_alt, params
                            )
                            # if ta + tr != tdp:
                            # continue
                            if ta == 0:
                                continue
                            if ta > params["maxAltCount"]:
                                continue
                            if normalBam:
                                na, nr, ndp = extractDepthSnv(
                                    normalBam,
                                    mut_chrom,
                                    mut_pos,
                                    mut_ref,
                                    mut_alt,
                                    params,
                                )
                                if na > 0:
                                    continue
                                if ndp < params["minNdepth"]:
                                    continue
                            else:
                                na, nr, ndp = (0, 0, 0)
                                if ta / tdp > params["maxAF"]:
                                    continue
                            
                            """
                            if readSet[0].template_length > 0:
                                readPos5p = min(
                                    muts_ind[nn] + 1,
                                    abs(readSet[0].template_length) - muts_ind[nn],
                                )
                                readPos3p = min(
                                    abs(readSet[0].reference_length) - muts_ind[nn],
                                    muts_ind[nn] + 1,
                                )
                            else:
                                readPos5p = min(
                                    readSet[0].reference_length - muts_ind[nn],
                                    abs(readSet[0].template_length)
                                    - readSet[0].reference_length
                                    + muts_ind[nn],
                                )
                                readPos3p = min(
                                    abs(readSet[0].reference_length) - muts_ind[nn],
                                    muts_ind[nn] + 1,
                                )
                            FPs.append(readPos5p)
                            RPs.append(readPos3p)
                            mut = {
                                "chrom": mut_chrom,
                                "pos": mut_pos,
                                "ref": mut_ref,
                                "alt": mut_alt,
                                "infos": {
                                    "F1R2": F1R2,
                                    "F2R1": F2R1,
                                    "TG": F1R2_ARLR[muts_ind[nn]],
                                    "BG": F2R1_ARLR[muts_ind[nn]],
                                    "TC": ",".join(
                                        [
                                            str(_)
                                            for _ in F1R2_count[
                                                :, muts_ind[nn]
                                            ].tolist()
                                        ]
                                    ),
                                    "BC": ",".join(
                                        [
                                            str(_)
                                            for _ in F2R1_count[
                                                :, muts_ind[nn]
                                            ].tolist()
                                        ]
                                    ),
                                },
                                "formats": ["AC", "RC", "DP"],
                                # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_dict[
                                "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                            ] = 0
                            muts.append(mut)
                        if isLearn:
                            continue
                        coverage[start_ind:end_ind][pass_bool] += 1
                        duplex_read_num_dict[duplex_no][1] += np.count_nonzero(
                            pass_bool
                        )
                        duplex_read_num_dict[duplex_no][0] += 1
                        duplex_count += 1
            """
            Calling block ends
            """
            currentReadDict = {
                bc: {"seqs": [rec], "F1R2": 0, "F2R1": 0, "names": {rec.query_name: 0}}
            }
            if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                currentReadDict[bc1 + "+" + bc2]["F1R2"] += 1
            else:
                currentReadDict[bc1 + "+" + bc2]["F2R1"] += 1
            currentStart = start
    """
    Calling block starts
    """
    for key in currentReadDict.keys():
        readSet = currentReadDict[key]["seqs"]
        mean_mapq = sum([seq.mapping_quality for seq in readSet]) / len(readSet)
        if mean_mapq < params["mapq"]:
            continue
        meanASXS = sum(
            [seq.get_tag("AS") - seq.get_tag("XS") for seq in readSet]
        ) / len(readSet)
        if meanASXS < 50:
            continue
        setBc = key.split("+")
        setBc1 = setBc[0]
        setBc2 = setBc[1]
        F2R1 = currentReadDict[key]["F2R1"]
        F1R2 = currentReadDict[key]["F1R2"]
        duplex_no = f"{min([F1R2,F2R1])}+{max([F1R2,F2R1])}"
        if duplex_read_num_dict.get(duplex_no) is None:
            duplex_read_num_dict[duplex_no] = [0, 0]
        unique_read_num += 1
        if setBc1 != setBc2 and F2R1 >= 1 and F1R2 >= 1:
            rs_reference_end = max([r.reference_end for r in readSet])
            chromNow = readSet[0].reference_name
            if chromNow != reference_mat_chrom or rs_reference_end >= reference_mat_end:
                ### Output coverage
                if "coverage" in locals():
                    non_zero_positions = np.nonzero(coverage + coverage_indel)
                    for pos in non_zero_positions[0].tolist():
                        locus_bed.write(
                            (
                                "\t".join(
                                    [
                                        reference_mat_chrom,
                                        str(pos + reference_mat_start),
                                        str(pos + 1 + reference_mat_start),
                                        str(coverage[pos]),
                                        str(coverage_indel[pos]),
                                    ]
                                )
                                + "\n"
                            ).encode("utf-8")
                        )
                    total_coverage += np.sum(coverage)
                    total_coverage_indel += np.sum(coverage_indel)
                # if chromNow != reference_mat_chrom:
                reference_mat_chrom = chromNow
                # current_reference = str(fasta[reference_mat_chrom].seq)
                reference_mat_start = readSet[0].reference_start
                try:
                    region_end = region[2]
                except:
                    region_end = 10e10
                contig_end = tumorBam.get_reference_length(chromNow)
                reference_mat_end = min(
                    readSet[0].reference_start + 100000,
                    max(region_end, max([seq.reference_end for seq in readSet])),
                    contig_end,
                )
                current_fasta = fasta.fetch(
                    reference_mat_chrom, reference_mat_start, reference_mat_end
                )
                (
                    prior_mat,
                    snp_mask,
                    indel_mask,
                    noise_mask,
                    n_cov_mask,
                    ref_np,
                ) = prepare_reference_mats(
                    reference_mat_chrom,
                    reference_mat_start,
                    reference_mat_end,
                    current_fasta,
                    germline,
                    noise,
                    indel_bed,
                    normalBam,
                    tumorBam,
                    params,
                )
                coverage = np.zeros(100000, dtype=int)
                coverage_indel = np.zeros(100000, dtype=int)

            ### Record read names to check if mate has been processed
            processed_flag = 0
            for seq in readSet:
                if seq.query_name in processed_read_names:
                    processed_read_names.remove(seq.query_name)
                    processed_flag = 1
                    break
            if processed_flag == 0:
                processed_read_names.add(readSet[0].query_name)

            start_ind = readSet[0].reference_start - reference_mat_start
            reference_length_min = min([read.reference_length for read in readSet])
            end_ind = (
                readSet[0].reference_start + reference_length_min - reference_mat_start
            )

            reference_length_max = max([read.reference_length for read in readSet])
            end_ind_max = (
                readSet[0].reference_start + reference_length_max - reference_mat_start
            )
            masks = np.zeros([4, end_ind - start_ind], dtype=bool)
            masks[0, :] = snp_mask[start_ind:end_ind]
            masks[1, :] = noise_mask[start_ind:end_ind]
            masks[2, :] = n_cov_mask[start_ind:end_ind]
            left, right = determineTrimLength(
                readSet[0], params=params, processed_flag=processed_flag
            )
            masks[3, :left] = True
            masks[3, -right:] = True
            antimask = np.all(~masks, axis=0)
            ### If the whole reads are masked:
            if not np.any(antimask):
                continue
            indel_bool = [
                ("I" in seq.cigarstring or "D" in seq.cigarstring) for seq in readSet
            ]
            if any(indel_bool):
                masks_indel = np.zeros([4, end_ind_max - start_ind], dtype=bool)
                masks_indel[0, :] = indel_mask[start_ind:end_ind_max]
                masks_indel[1, :] = noise_mask[start_ind:end_ind_max]
                masks_indel[2, :] = n_cov_mask[start_ind:end_ind_max]
                left, right = determineTrimLength(
                    readSet[0], params=params, processed_flag=processed_flag
                )
                masks_indel[3, :left] = True
                masks_indel[3, -right:] = True
                antimask_indel = np.all(~masks_indel, axis=0)
                (
                    F1R2_ARLR,
                    F2R1_ARLR,
                    indels,
                    F1R2_ref_count,
                    F1R2_alt_count,
                    F2R1_ref_count,
                    F2R1_alt_count,
                ) = genotypeDSIndel(readSet, tumorBam, antimask_indel, params)
                DCS = np.logical_and(
                    F1R2_ARLR >= params["pcutoff"],
                    F2R1_ARLR >= params["pcutoff"],
                )
                pass_inds = np.nonzero(DCS)[0].tolist()
                indels_pass = [indels[_] for _ in pass_inds]
                coverage_indel[start_ind:end_ind_max][antimask_indel] += 1
                if len(indels_pass) == 1:
                    indel = indels_pass[0]
                    indel_chrom = chromNow
                    indel_pos = int(indel.split(":")[0])
                    indel_size = int(indel.split(":")[1])
                    NMs = [seq.get_tag("NM") for seq in readSet]
                    averageNM = sum(NMs) / len(NMs)
                    if averageNM - abs(indel_size) > 1:
                        continue
                    if indel_size < 0:
                        indel_ref = str(
                            current_fasta[
                                indel_pos
                                - reference_mat_start : indel_pos
                                - reference_mat_start
                                - indel_size
                                + 1
                            ]
                        ).upper()
                        indel_alt = current_fasta[
                            indel_pos - reference_mat_start
                        ].upper()
                    else:
                        indel_ref = current_fasta[
                            indel_pos - reference_mat_start
                        ].upper()
                        indel_alt = indel_ref + indel.split(":")[2]
                    indel_str = (
                        str(indel_chrom)
                        + ":"
                        + str(indel_pos)
                        + str(indel_ref)
                        + ":"
                        + str(indel_alt)
                    )
                    # if indel_dict.get(indel_str):
                    # continue
                    """
                    ta, tr, tdp = extractDepthIndel(
                        tumorBam,
                        indel_chrom,
                        indel_pos + 1,
                        indel_ref,
                        indel_alt,
                        params,
                    )
                    if ta + tr != tdp:
                        continue
                    if ta == 0:
                        continue
                    if ta > params["maxAltCount"]:
                        continue
                    if IndelFilterByWindows(tumorBam,indel_chrom,indel_pos+1,5,params):
                        continue
                    if normalBam:
                        na, nr, ndp = extractDepthIndel(
                            normalBam,
                            indel_chrom,
                            indel_pos + 1,
                            indel_ref,
                            indel_alt,
                            params,
                        )
                        if ndp - nr > 0:
                            continue
                        if ndp < params["minNdepth"]:
                            continue
                        if IndelFilterByWindows(normalBam,indel_chrom,indel_pos+1,5,params):
                            continue
                    else:
                        na, nr, ndp = (0, 0, 0)
                        if ta / tdp > params["maxAF"]:
                            continue
                    """
                    indel_rec = {
                        "chrom": chromNow,
                        "pos": indel_pos + 1,
                        "ref": indel_ref,
                        "alt": indel_alt,
                        "infos": {
                            "F1R2": int(
                                F1R2_alt_count[pass_inds[0]]
                                + F1R2_ref_count[pass_inds[0]]
                            ),
                            "F2R1": int(
                                F2R1_alt_count[pass_inds[0]]
                                + F2R1_ref_count[pass_inds[0]]
                            ),
                            "TG": F1R2_ARLR[pass_inds[0]],
                            "BG": F2R1_ARLR[pass_inds[0]],
                            "TC": ",".join(
                                [
                                    str(F1R2_alt_count[pass_inds[0]]),
                                    str(F1R2_ref_count[pass_inds[0]]),
                                ]
                            ),
                            "BC": ",".join(
                                [
                                    str(F2R1_alt_count[pass_inds[0]]),
                                    str(F2R1_ref_count[pass_inds[0]]),
                                ]
                            ),
                        },
                        "formats": ["AC", "RC", "DP"],
                        # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_indels.append(indel_rec)
                    indel_dict[indel_str] = 1
            else:
                if isLearn:
                    if sum([int(seq.get_tag("NM")) for seq in readSet]) == 0:
                        continue
                ### Calculate genotype probability
                (
                    F1R2_ARLR,
                    F2R1_ARLR,
                    ref_int,
                    alt_int,
                    antimask,
                    F1R2_count,
                    F2R1_count,
                ) = genotypeDSSnv(
                    readSet,
                    ref_np[start_ind:end_ind],
                    prior_mat[start_ind:end_ind, :],
                    antimask,
                    params,
                )
                if isLearn:
                    mismatch_now = profileTriNucMismatches(
                        readSet,
                        ref_np[start_ind:end_ind],
                        antimask,
                        params,
                    )
                    for trinuc in mismatch_dict.keys():
                        for nnn in range(4):
                            mismatch_dict[trinuc][nnn] += mismatch_now[trinuc][nnn]
                ref_int = ref_np[start_ind:end_ind]
                refs_ind = np.nonzero(
                    np.logical_and(
                        F1R2_ARLR <= -params["pcutoff"] + 2 * log10(params["mutRate"]),
                        F2R1_ARLR <= -params["pcutoff"] + 2 * log10(params["mutRate"]),
                    )
                )[0].tolist()
                muts_ind = np.nonzero(
                    np.logical_and(
                        F1R2_ARLR >= params["pcutoff"],
                        F2R1_ARLR >= params["pcutoff"],
                    )
                )[0].tolist()
                pass_bool = np.full(F1R2_ARLR.size, False, dtype=bool)
                pass_bool[refs_ind] = True
                pass_bool[muts_ind] = True
                pos = [
                    mut_ind + start_ind + reference_mat_start for mut_ind in muts_ind
                ]
                # muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                mut_positions = [
                    mut_ind + start_ind + reference_mat_start + 1
                    for mut_ind in muts_ind
                ]
                """
                if len(mut_positions) > 1:
                    break_flag = 0
                    for nn in range(1, len(mut_positions)):
                        if mut_positions[nn] - mut_positions[nn - 1] != 1:
                            break_flag = 1
                    if break_flag:
                        continue
                """
                NMs = [seq.get_tag("NM") for seq in readSet]
                averageNM = sum(NMs) / len(NMs)
                if averageNM - len(mut_positions) > 1:
                    continue
                if len(mut_positions) > params["maxMnv"]:
                    continue
                for nn in range(len(mut_positions)):
                    mut_chrom = reference_mat_chrom
                    mut_pos = mut_positions[nn]
                    mut_ref = num2base[ref_int[muts_ind[nn]]]
                    mut_alt = num2base[alt_int[muts_ind[nn]]]
                    """
                    # if muts_dict.get(
                    # "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                    # ):
                    # continue
                    ta, tr, tdp = extractDepthSnv(
                        tumorBam, mut_chrom, mut_pos, mut_ref, mut_alt, params
                    )
                    # if ta + tr != tdp:
                    # continue
                    if ta == 0:
                        continue
                    if ta > params["maxAltCount"]:
                        continue
                    if normalBam:
                        na, nr, ndp = extractDepthSnv(
                            normalBam,
                            mut_chrom,
                            mut_pos,
                            mut_ref,
                            mut_alt,
                            params,
                        )
                        if na > 0:
                            continue
                        if ndp < params["minNdepth"]:
                            continue
                    else:
                        na, nr, ndp = (0, 0, 0)
                        if ta / tdp > params["maxAF"]:
                            continue
                    """

                    if readSet[0].template_length > 0:
                        readPos5p = min(
                            muts_ind[nn] + 1,
                            abs(readSet[0].template_length) - muts_ind[nn],
                        )
                        readPos3p = min(
                            abs(readSet[0].reference_length) - muts_ind[nn],
                            muts_ind[nn] + 1,
                        )
                    else:
                        readPos5p = min(
                            readSet[0].reference_length - muts_ind[nn],
                            abs(readSet[0].template_length)
                            - readSet[0].reference_length
                            + muts_ind[nn],
                        )
                        readPos3p = min(
                            abs(readSet[0].reference_length) - muts_ind[nn],
                            muts_ind[nn] + 1,
                        )
                    FPs.append(readPos5p)
                    RPs.append(readPos3p)
                    mut = {
                        "chrom": mut_chrom,
                        "pos": mut_pos,
                        "ref": mut_ref,
                        "alt": mut_alt,
                        "infos": {
                            "F1R2": F1R2,
                            "F2R1": F2R1,
                            "TG": F1R2_ARLR[muts_ind[nn]],
                            "BG": F2R1_ARLR[muts_ind[nn]],
                            "TC": ",".join(
                                [str(_) for _ in F1R2_count[:, muts_ind[nn]].tolist()]
                            ),
                            "BC": ",".join(
                                [str(_) for _ in F2R1_count[:, muts_ind[nn]].tolist()]
                            ),
                        },
                        "formats": ["AC", "RC", "DP"],
                        # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_dict["_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])] = 0
                    muts.append(mut)
                if isLearn:
                    continue
                coverage[start_ind:end_ind][pass_bool] += 1
                duplex_read_num_dict[duplex_no][1] += np.count_nonzero(pass_bool)
                duplex_read_num_dict[duplex_no][0] += 1
                duplex_count += 1
    """
    Calling block ends
    """
    mut_dict = dict()
    mut_pass_filter = []
    for mut in muts:
        chrom = mut["chrom"]
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]

        if not mut_dict.get(":".join([chrom, str(pos), ref, alt])):
            ta, tr, tdp = extractDepthSnv(tumorBam, chrom, pos, ref, alt, params)
            window_filter = False
            if IndelFilterByWindows(tumorBam, chrom, pos, 3, params):
                window_filter = True
            if normalBam:
                na, nr, ndp = extractDepthSnv(normalBam, chrom, pos, ref, alt, params)
                if IndelFilterByWindows(normalBam, chrom, pos, 3, params):
                    window_filter = True
            else:
                na, nr, ndp = (0, 0, 0)
            mut_dict[":".join([chrom, str(pos), ref, alt])] = (
                ta,
                tr,
                tdp,
                na,
                nr,
                ndp,
                window_filter,
            )
        else:
            ta, tr, tdp, na, nr, ndp, window_filter = mut_dict[
                ":".join([chrom, str(pos), ref, alt])
            ]
        if window_filter:
            continue
        if ta > params["maxAltCount"]:
            continue
        if ta == 0:
            continue
        if normalBam:
            if na > 0:
                continue
            if ndp < params["minNdepth"]:
                continue
        else:
            if ta / tdp > params["maxAF"]:
                continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        mut_pass_filter.append(mut)

    muts_indels_dict = dict()
    muts_indels_pass_filter = []
    for mut in muts_indels:
        chrom = mut["chrom"]
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]

        if not muts_indels_dict.get(":".join([chrom, str(pos), ref, alt])):
            ta, tr, tdp = extractDepthIndel(tumorBam, chrom, pos, ref, alt, params)
            window_filter = False
            if IndelFilterByWindows(tumorBam, chrom, pos, 3, params):
                window_filter = True
            if normalBam:
                na, nr, ndp = extractDepthIndel(normalBam, chrom, pos, ref, alt, params)
                if IndelFilterByWindows(normalBam, chrom, pos, 3, params):
                    window_filter = True
            else:
                na, nr, ndp = (0, 0, 0)
            muts_indels_dict[":".join([chrom, str(pos), ref, alt])] = (
                ta,
                tr,
                tdp,
                na,
                nr,
                ndp,
                window_filter,
            )
        else:
            ta, tr, tdp, na, nr, ndp, window_filter = muts_indels_dict[
                ":".join([chrom, str(pos), ref, alt])
            ]
        if window_filter:
            continue
        if ta > params["maxAltCount"]:
            continue
        if ta == 0:
            continue
        if normalBam:
            if ndp - nr > 0:
                continue
            if ndp < params["minNdepth"]:
                continue
        else:
            if ta / tdp > params["maxAF"]:
                continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        muts_indels_pass_filter.append(mut)
    
    if isLearn:
        return mismatch_dict, FPs, RPs
    
    if "coverage" in locals():
        non_zero_positions = np.nonzero(coverage + coverage_indel)
        for pos in non_zero_positions[0].tolist():
            locus_bed.write(
                (
                    "\t".join(
                        [
                            reference_mat_chrom,
                            str(pos + reference_mat_start),
                            str(pos + 1 + reference_mat_start),
                            str(coverage[pos]),
                            str(coverage_indel[pos]),
                        ]
                    )
                    + "\n"
                ).encode("utf-8")
            )
        total_coverage += np.sum(coverage)
        total_coverage_indel += np.sum(coverage_indel)
        
    print(
        f"Process {processNo} finished. Time: {(time.time()-starttime)/60: .2f} minutes"
    )

    return (
        mut_pass_filter,
        total_coverage,
        recCount,
        duplex_count,
        duplex_read_num_dict,
        muts_indels_pass_filter,
        total_coverage_indel,
        unique_read_num,
        pass_read_num,
        FPs,
        RPs,
    )
