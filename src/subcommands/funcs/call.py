#!/usr/bin/env python3

import gzip

import time
import math

import numpy as np
import pandas as pd
from pysam import FastaFile as FASTA
from pysam import TabixFile as BED
from pysam import VariantFile as VCF
import pysam
import re
import h5py

from .depth import (
    extractDepthRegion,
    extractDepthSnv,
    extractDepthIndel,
    detectOverlapDiscord,
)
from .prob import genotypeDSSnv, genotypeDSIndel
from .learn import profileTriNucMismatches
from .misc import getAlignmentObject as BAM
from .indels import findIndels

# from . misc import IndelFilterByWindows


def prepare_reference_mats(
    chrom,
    start,
    end,
    reference_int,
    trinuc_int,
    germline_bed,
    noise_bed,
    indel_bed,
    include_bed,
    feature_beds,
    nbams,
    tbam,
    params,
):
    ### Define and Initialize
    af_miss = params["mutRate"]
    af_cutoff = params["germline_cutoff"]
    m = end - start
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    trinuc2num = params["trinuc2num_dict"]

    snp_mask = np.full(m, False, dtype=bool)
    indel_mask = np.full(m, False, dtype=bool)
    noise_mask = np.full(m, False, dtype=bool)
    n_cov_mask = np.full(m, False, dtype=bool)

    ### Initialize prior mat as no germline resource
    noise_mask[reference_int == 4] = True
    noise_mask[trinuc_int == 96] = True
    reference_int[reference_int == 4] = 0
    prior_mat = np.full([m, 4], af_miss, dtype=float)
    prior_mat[np.ogrid[:m], reference_int] = 1 - 3 * af_miss

    ### Adjust by germline
    if germline_bed != None:
        for rec in germline_bed.fetch(chrom, start, end):
            ind = rec.pos - 1 - start
            ref = rec.ref
            try:
                afs = rec.info["AF"]
            except:
                afs = [1 for _ in range(len(rec.alts))]
            if len(ref) == 1:
                for ii, alt in enumerate(rec.alts):
                    if len(alt) == 1:
                        # prior_mat[ind, base2num[alt]] = afs[ii]
                        # has_snp = True
                        if afs[ii] >= af_cutoff:
                            snp_mask[ind] = True
                    elif afs[ii] >= af_cutoff:
                        indel_mask[ind] = True
            if len(ref) != 1:
                for ii, alt in enumerate(rec.alts):
                    if len(alt) == len(ref):
                        diff = np.array([a != b for a, b in zip(list(ref), list(alt))])
                        if ind >= 0:
                            if afs[ii] >= af_cutoff:
                                snp_mask[ind : ind + len(alt)][
                                    diff[: (snp_mask.size - ind)]
                                ] = True
                        else:
                            if afs[ii] >= af_cutoff:
                                snp_mask[: min(snp_mask.size - ind, diff.size) + ind][
                                    diff[-ind : min(snp_mask.size - ind, diff.size)]
                                ] = True
                    elif afs[ii] >= af_cutoff:
                        indel_mask[max(ind, 0) : ind + len(ref)] = True
    ### Prepare noise mask
    if noise_bed != None:
        for rec in noise_bed.fetch(chrom, start, end, parser=pysam.asBed()):
            interval_start = max(rec.start, start)
            interval_end = min(rec.end, end)
            interval_len = interval_end - interval_start
            interval_start_ind = interval_start - start
            noise_mask[interval_start_ind : interval_start_ind + interval_len] = True
    if indel_bed != None:
        for rec in indel_bed.fetch(chrom, start, end, parser=pysam.asBed()):
            interval_start = max(rec.start, start)
            interval_end = min(rec.end, end)
            interval_len = interval_end - interval_start
            interval_start_ind = interval_start - start
            indel_mask[interval_start_ind : interval_start_ind + interval_len] = True
    if include_bed != None:
        include_arr = np.zeros(end - start, dtype=bool)
        for rec in include_bed.fetch(chrom, start, end, parser=pysam.asBed()):
            interval_start = max(rec.start, start)
            interval_end = min(rec.end, end)
            interval_len = interval_end - interval_start
            interval_start_ind = interval_start - start
            include_arr[interval_start_ind : interval_start_ind + interval_len] = True
        include_mask = ~include_arr
    else:
        include_mask = np.zeros(end - start, dtype=bool)
    if feature_beds != None:
        feature_mat = np.zeros(
            (len(feature_beds), end - start),
        )
        for nn, rec in range(
            feature_files.fetch(chrom, start, end, parser=pysam.asBed())
        ):
            interval_start = max(rec.start, start)
            interval_end = min(rec.end, end)
            interval_len = interval_end - interval_start
            interval_start_ind = interval_start - start
            feature_mat[
                nn, interval_start_ind : interval_start_ind + interval_len
            ] = True
    else:
        feature_mat = None

    ### Prepare normal coverage mask
    if not params["isLearn"]:
        if nbams:
            depth = np.zeros(end - start)
            for nbam in nbams:
                depth_now, indel_mask_out = extractDepthRegion(
                    nbam, chrom, start, end, params
                )
                indel_mask[indel_mask_out] = True
                depth += depth_now
            n_cov_mask = depth < params["minNdepth"]

        if params["maxAF"] < 1:
            depth, indel_mask_out = extractDepthRegion(tbam, chrom, start, end, params)
            ma = params["maxAF"]
            min_depth = math.ceil(1 / ma)
            n_cov_mask = depth < min_depth
            indel_mask[indel_mask_out] = True
    return (
        prior_mat,
        snp_mask,
        indel_mask,
        noise_mask,
        n_cov_mask,
        include_mask,
        feature_mat,
    )  # , reference_int, trinuc_int


def determineTrimLength(seq, params, processed_flag):
    if seq.template_length > 0 and not processed_flag:
        overlap = 0  # Never mask overlap of forward read
        left = params["trim5"]
        right_frag = params["trim5"] - min(
            params["trim5"], abs(seq.template_length) - seq.reference_length
        )
        right_read = params["trim3"]
        right = max(right_frag, right_read)
    else:
        ### Mask overlap of reverse read
        if processed_flag:
            mate_cigar = seq.get_tag("MC")
            cigar_m = re.findall(r"(\d+)M", mate_cigar)
            cigar_d = re.findall(r"(\d+)D", mate_cigar)
            mate_reference_length = sum([int(_) for _ in cigar_m]) + sum(
                [int(_) for _ in cigar_d]
            )
            overlap = max(
                0,
                seq.reference_length + mate_reference_length - abs(seq.template_length),
            )
        else:
            overlap = 0
        right_frag = params["trim5"]
        right_read = params["trim3"]
        right = max(right_frag, right_read)
        left_frag = params["trim5"] - min(
            params["trim5"], abs(seq.template_length) - seq.reference_length
        )
        left = max(left_frag, overlap, params["trim3"])
    return left, right


def nums2str(nums, num2base="ATCG"):
    bases = [num2base[_] for _ in nums]
    return "".join(bases)


def bamIterateMultipleRegion(bam, regions):
    bamObject = BAM(bam, "rb")
    for region in regions:
        for rec in bamObject.fetch(*region):
            if len(region) >= 2:
                if rec.reference_start < region[1]:
                    continue
            yield rec, region


def callBam(params, processNo):
    # Get parameters
    bam = params["tumorBam"]
    nbams = params["normalBams"]
    regions = params["regions"]
    if params["germline"]:
        germline = VCF(params["germline"], "rb")
    else:
        germline = None
    all_chroms = [_[0] for _ in regions]
    start_time = time.time()
    ##for ch in all_chroms:
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
    if params["region_file"]:
        include_bed = BED(params["region_file"])
    else:
        include_bed = None
    if params["feature_files"]:
        feature_beds = [BED(_) for _ in params["feature_files"]]
    else:
        feature_beds = None
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    muts = []
    muts_dict = dict()
    muts_indels = []
    duplex_read_num_dict = dict()
    duplex_read_num_dict_trinuc = dict()
    if feature_beds:
        duplex_read_num_dict_trinuc_features = [dict() for _ in feature_beds]
    else:
        duplex_read_num_dict_trinuc_features = None
    unique_read_num = 0
    pass_read_num = 0
    FPs = []
    RPs = []
    indel_dict = dict()
    mismatch_mat = np.zeros([64, 4])
    indelerr_mat = np.zeros([40, 11])
    mismatch_dmg_mat = np.zeros([64, 4])
    indel_dmg_mat = np.zeros([40, 11])
    trinuc2num = dict()
    num2trinuc = list()
    trinuc_order = 0
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                trinuc = minus_base + ref_base + plus_base
                trinuc2num[trinuc] = trinuc_order
                num2trinuc.append(trinuc)
                trinuc_order += 1
    for plus_base in ["T", "A", "G", "C"]:
        for ref_base in ["G", "A"]:
            for minus_base in ["T", "A", "G", "C"]:
                trinuc = minus_base + ref_base + plus_base
                trinuc2num[trinuc] = trinuc_order
                num2trinuc.append(trinuc)
                trinuc_order += 1

    trinuc_convert_np = np.zeros([64, 4], dtype=np.uint8)
    for trinuc in trinuc2num.keys():
        row = np.zeros(4)
        row_num = trinuc2num[trinuc]
        row[0] = trinuc2num[trinuc[0] + "A" + trinuc[2]]
        row[1] = trinuc2num[trinuc[0] + "T" + trinuc[2]]
        row[2] = trinuc2num[trinuc[0] + "C" + trinuc[2]]
        row[3] = trinuc2num[trinuc[0] + "G" + trinuc[2]]
        trinuc_convert_np[row_num, :] = row
    params["trinuc_convert"] = trinuc_convert_np
    params["trinuc2num_dict"] = trinuc2num
    params["num2trinuc_list"] = num2trinuc

    ### Load amp error matrix
    if not params["amperr_file"]:
        prob_amp = params["amperr"]
        row1 = np.array([[prob_amp / 3, prob_amp / 3, 1 - prob_amp, prob_amp / 3]])
        row5 = np.array([[prob_amp / 3, 1 - prob_amp, prob_amp / 3, prob_amp / 3]])
        first_32rows = np.tile(
            np.repeat(np.concatenate([row1, row5], axis=0), 4, axis=0), (4, 1)
        )
        second_32rows = first_32rows[:, np.array([1, 0, 3, 2])]
        ampmat = np.vstack((first_32rows, second_32rows))
    else:
        ampmat = pd.read_csv(params["amperr_file"], sep="\t", index_col=0).to_numpy()
    ampmat = ampmat / ampmat.sum(axis=1, keepdims=True)
    # ampmat_avg_error = (1 - ampmat.max(axis=1,keepdims=True))/3
    ampmat_min_error = ampmat.min(axis=1, keepdims=True)
    ampmat = np.concatenate([ampmat, ampmat_min_error], axis=1)
    params["ampmat"] = ampmat

    ampmat_rev = np.zeros([64, 4])
    for trinuc in trinuc2num.keys():
        refbase = trinuc[1]
        for nn, altbase in enumerate(["A", "T", "C", "G"]):
            ampmat_rev[trinuc2num[trinuc], nn] = ampmat[
                trinuc2num[trinuc[0] + altbase + trinuc[2]], base2num[refbase]
            ]
    # ampmat_rev_avg_error = (1 - ampmat_rev.max(axis=1,keepdims=True))/3
    ampmat_rev_min_error = ampmat_rev.min(axis=1, keepdims=True)
    ampmat_rev = np.concatenate([ampmat_rev, ampmat_rev_min_error], axis=1)
    params["ampmat_rev"] = ampmat_rev

    if params["amperri_file"]:
        ampmat_indel = np.loadtxt(params["amperri_file"], delimiter="\t")
        ampmat_indel = ampmat_indel / np.sum(ampmat_indel, axis=1, keepdims=True)
        for nn in range(1, 40):
            current_row = ampmat_indel[nn, :]
            current_row[current_row == 0] = ampmat_indel[nn - 1, :][current_row == 0]
            ampmat_indel[nn, :] = current_row
        params["ampmat_indel"] = ampmat_indel
    else:
        params["ampmat_indel"] = np.ones([40, 11]) * params["amperri"]
        ampmat_indel = params["ampmat_indel"]

    ampmat_indel_rev = np.fliplr(ampmat_indel)
    params["ampmat_indel_rev"] = ampmat_indel_rev

    # params["ampmat_indel_mean"] = np.mean(ampmat_indel,axis=1)
    # params["ampmat_indel_rev_mean"] = np.mean(ampmat_indel_rev,axis=1)

    ### Load damage matrix
    if not params["dmgerr_file"]:
        prob_dmg = params["dmgerr"]
        row1 = np.array([[prob_dmg / 3, prob_dmg / 3, 1 - prob_dmg, prob_dmg / 3]])
        row5 = np.array([[prob_dmg / 3, 1 - prob_dmg, prob_dmg / 3, prob_dmg / 3]])
        first_32rows = np.tile(
            np.repeat(np.concatenate([row1, row5], axis=0), 4, axis=0), (4, 1)
        )
        second_32rows = first_32rows[:, np.array([1, 0, 3, 2])]
        dmgmat = np.vstack((first_32rows, second_32rows))
    else:
        dmgmat = pd.read_csv(params["dmgerr_file"], sep="\t", index_col=0).to_numpy()
        dmgmat += 1
    dmgmat = dmgmat / dmgmat.sum(axis=1, keepdims=True)
    dmgmat_min_error = dmgmat.min(axis=1, keepdims=True)
    dmgmat = np.concatenate([dmgmat, dmgmat_min_error], axis=1)
    # dmgmat += 1e-9

    params["dmgmat_top"] = dmgmat
    params["trinuc2num_dict"] = trinuc2num

    dmgmat_rev = np.zeros([64, 4])
    for trinuc in trinuc2num.keys():
        refbase = trinuc[1]
        for nn, altbase in enumerate(["A", "T", "C", "G"]):
            dmgmat_rev[trinuc2num[trinuc], nn] = dmgmat[
                trinuc2num[trinuc[0] + altbase + trinuc[2]], base2num[refbase]
            ]
    dmgmat_rev_min_error = dmgmat_rev.min(axis=1, keepdims=True)
    dmgmat_rev = np.concatenate([dmgmat_rev, dmgmat_rev_min_error], axis=1)
    # dmgmat_rev += 1e-9
    params["dmgmat_rev_top"] = dmgmat_rev

    dmgmat_b = np.vstack((dmgmat[32:64, [1, 0, 3, 2]], dmgmat[:32, [1, 0, 3, 2]]))
    dmgmat_b_min_error = dmgmat_b.min(axis=1, keepdims=True)
    dmgmat_rev_b = np.vstack(
        (dmgmat_rev[32:64, [1, 0, 3, 2]], dmgmat_rev[:32, [1, 0, 3, 2]])
    )
    dmgmat_rev_b_min_error = dmgmat_rev_b.min(axis=1, keepdims=True)
    dmgmat_b = np.concatenate([dmgmat_b, dmgmat_b_min_error], axis=1)
    dmgmat_rev_b = np.concatenate([dmgmat_rev_b, dmgmat_rev_b_min_error], axis=1)
    # dmgmat_b += 1e-9
    # dmgmat_rev_b += 1e-9
    params["dmgmat_bot"] = dmgmat_b
    params["dmgmat_rev_bot"] = dmgmat_rev_b

    if params["dmgerri_file"]:
        dmgmat_indel = np.loadtxt(params["dmgerri_file"], delimiter="\t")
        dmgmat_indel[0, :] += 1
        dmgmat_indel[:, 5] += 1
        dmgmat_indel = dmgmat_indel / np.sum(dmgmat_indel, axis=1, keepdims=True)
        for nn in range(1, 40):
            current_row = dmgmat_indel[nn, :]
            current_row[current_row == 0] = dmgmat_indel[nn - 1, :][current_row == 0]
            dmgmat_indel[nn, :] = current_row
        # dmgmat_indel += 1e-9
        params["dmgmat_indel"] = dmgmat_indel
    else:
        params["dmgmat_indel"] = np.ones([40, 11]) * params["dmgerri"]
        dmgmat_indel = params["dmgmat_indel"]

    dmgmat_indel_rev = np.fliplr(dmgmat_indel)
    # dmgmat_indel_rev += 1e-9
    params["dmgmat_indel_rev"] = dmgmat_indel_rev
    params["dmgmat_indel_mean"] = np.mean(dmgmat_indel, axis=1)
    params["dmgmat_indel_rev_mean"] = np.mean(dmgmat_indel_rev, axis=1)

    # Initialize
    total_coverage = 0
    total_coverage_indel = 0
    starttime = time.time()
    tumorBam = BAM(bam, "rb")
    if nbams:
        normalBams = list()
        for nbam in nbams:
            normalBams.append(BAM(nbam, "rb"))
    else:
        normalBams = None
    currentStart = -1
    currentReadDict = {}
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
    if params["maxNM"]:
        print(
            f"Process {str(processNo)}: screening for highly-damaged or misaligned reads"
        )
        read_blacklist = set()
        rec_num = 0
        for rec, region in bamIterateMultipleRegion(bam, regions):
            rec_num += 1
            if rec.query_name in read_blacklist or rec.is_unmapped:
                continue
            id_length = 0
            for cigar in rec.cigartuples:
                if cigar[0] == 1 or cigar[0] == 2:
                    id_length += cigar[1]
            NM_no_id = rec.get_tag("NM") - id_length
            if NM_no_id >= params["maxNM"]:
                read_blacklist.add(rec.query_name)
        currentTime = (time.time() - starttime) / 60
        starttime = time.time()
        if rec_num > 0:
            percent_blocked = len(read_blacklist) * 2 / rec_num * 100
        else:
            percent_blocked = 0
        print(
            f"Process {str(processNo)}: finished screening highly damaged reads in {currentTime: .2f} minutes. Blacklisted {len(read_blacklist)} ({percent_blocked: .2f}%) possible highly damaged read and started variant calling."
        )
    for rec, region in bamIterateMultipleRegion(bam, regions):
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
        if "I" in rec.cigarstring or "D" in rec.cigarstring:
            if len(findIndels(rec)) >= 2:
                continue
        pass_read_num += 1
        start = rec.reference_start
        bc = rec.query_name.split("_")[1]
        bcsplit = bc.split("+")
        bc1 = bcsplit[0]
        bc2 = bcsplit[1]
        if (rec.is_read1 and rec.is_forward) or (rec.is_read2 and rec.is_reverse):
            label = bc1 + "+" + bc2  # + ":" + "f1r2"
            # label = bc1 + "+" + bc2 + ":" + "f2r1"
        else:
            label = bc2 + "+" + bc1  # + ":" + "f1r2"
        chrom = tumorBam.get_reference_name(rec.reference_id)
        if currentStart == -1:
            currentStart = start
        if start == currentStart:
            if currentReadDict.get(label):
                if currentReadDict[label]["names"].get(rec.query_name):
                    if rec.is_read2:
                        continue
                    else:
                        ind = currentReadDict[label]["names"].get(rec.query_name)
                        currentReadDict[label]["seqs"][ind] = rec
                else:
                    currentReadDict[label]["seqs"].append(rec)
                    currentReadDict[label]["names"][rec.query_name] = (
                        len(currentReadDict[label]["seqs"]) - 1
                    )
                    if (rec.is_forward and rec.is_read1) or (
                        rec.is_reverse and rec.is_read2
                    ):
                        currentReadDict[label]["F1R2"] += 1
                    else:
                        currentReadDict[label]["F2R1"] += 1
            else:
                currentReadDict.update(
                    {
                        label: {
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
                    currentReadDict[label]["F1R2"] += 1
                else:
                    currentReadDict[label]["F2R1"] += 1

            """
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
            """
        else:
            """
            Calling block starts
            """
            for key in currentReadDict.keys():
                readSet = currentReadDict[key]["seqs"]
                all_dup = True
                for _ in readSet:
                    if not _.is_duplicate:
                        all_dup = False
                        break
                if all_dup:
                    continue
                if params["maxNM"]:
                    blacklist_num = 0
                    for seq in readSet:
                        if seq.query_name in read_blacklist:
                            blacklist_num += 1
                    if blacklist_num / len(readSet) >= 0.5:
                        continue
                mean_mapq = sum([seq.mapping_quality for seq in readSet]) / len(readSet)
                if mean_mapq < params["mapq"]:
                    continue
                meanASXS = sum(
                    [seq.get_tag("AS") - seq.get_tag("XS") for seq in readSet]
                ) / len(readSet)
                if meanASXS < 50:
                    continue
                setBc = key.split(":")[0].split("+")
                setBc1 = setBc[0]
                setBc2 = setBc[1]
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]
                duplex_no = f"{F1R2}+{F2R1}"
                if duplex_read_num_dict.get(duplex_no) is None:
                    duplex_read_num_dict[duplex_no] = [0, 0]
                    duplex_read_num_dict_trinuc[duplex_no] = np.zeros(96, dtype=int)
                if feature_beds:
                    for nn in range(len(duplex_read_num_dict_trinuc_features)):
                        if duplex_read_num_dict[nn].get(duplex_no) is None:
                            duplex_read_num_dict_trinuc_features[nn][
                                duplex_no
                            ] = np.zeros(96, dtype=int)
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
                            readSet[0].reference_start + 1000000,
                            max(
                                region_end, max([seq.reference_end for seq in readSet])
                            ),
                            contig_end,
                        )
                        #
                        ref_h5 = h5py.File(params["reference"] + ".ref.h5", "r")
                        tn_h5 = h5py.File(params["reference"] + ".tn.h5", "r")
                        hp_h5 = h5py.File(params["reference"] + ".hp.h5", "r")
                        ref_np = ref_h5[reference_mat_chrom][
                            reference_mat_start:reference_mat_end
                        ]
                        trinuc_np = tn_h5[reference_mat_chrom][
                            reference_mat_start:reference_mat_end
                        ]
                        hp_np = hp_h5[reference_mat_chrom][
                            :, reference_mat_start:reference_mat_end
                        ]
                        (
                            prior_mat,
                            snp_mask,
                            indel_mask,
                            noise_mask,
                            n_cov_mask,
                            include_mask,
                            feature_mat
                            # ref_np,
                            # trinuc_np
                        ) = prepare_reference_mats(
                            reference_mat_chrom,
                            reference_mat_start,
                            reference_mat_end,
                            # current_fasta,
                            ref_np,
                            trinuc_np,
                            germline,
                            noise,
                            indel_bed,
                            include_bed,
                            feature_beds,
                            nbams,
                            bam,
                            params,
                        )
                        # print(ref_np,reference_mat_start)
                        coverage = np.zeros(1000000, dtype=int)
                        coverage_indel = np.zeros(1000000, dtype=int)
                        if feature_beds:
                            trinuc_features = [
                                np.zeros((96, 1000000), dtype=int) for _ in feature_beds
                            ]
                        else:
                            trinuc_features = list()

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
                    reference_length_max = max(
                        [read.reference_length for read in readSet]
                    )
                    end_ind = (
                        readSet[0].reference_start
                        + reference_length_max
                        - reference_mat_start
                    )

                    end_ind_max = (
                        readSet[0].reference_start
                        + reference_length_max
                        - reference_mat_start
                    )
                    masks = np.zeros([5, end_ind - start_ind], dtype=bool)
                    masks[0, :] = snp_mask[start_ind:end_ind]
                    masks[1, :] = noise_mask[start_ind:end_ind]
                    masks[2, :] = n_cov_mask[start_ind:end_ind]
                    masks[3, :] = include_mask[start_ind:end_ind]
                    left, right = determineTrimLength(
                        readSet[0], params=params, processed_flag=processed_flag
                    )
                    masks[4, :left] = True
                    masks[4, -right:] = True
                    antimask = np.all(~masks, axis=0)
                    antimask[trinuc_np[start_ind:end_ind] > 64] = False
                    ### If the whole reads are masked:
                    if not np.any(antimask):
                        continue
                    indel_bool = [
                        ("I" in seq.cigarstring or "D" in seq.cigarstring)
                        for seq in readSet
                    ]
                    # if any(indel_bool):
                    if not isLearn:
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
                            LR,
                            indels,
                            hps,
                            F1R2_ref_count,
                            F1R2_alt_count,
                            F2R1_ref_count,
                            F2R1_alt_count,
                        ) = genotypeDSIndel(
                            readSet,
                            tumorBam,
                            antimask_indel,
                            hp_np[0, start_ind:end_ind],
                            params,
                        )
                        """
                        DCS = np.logical_and(
                            F1R2_BLR >= params["pcutoff"],
                            F2R1_BLR >= params["pcutoff"],
                        )
                        """
                        pass_inds = np.nonzero(LR <= params["pcutoff"])[0].tolist()
                        indels_pass = [indels[_] for _ in pass_inds]
                        coverage_indel[start_ind:end_ind_max][antimask_indel] += 1
                        for nn in range(len(indels_pass)):
                            indel = indels_pass[nn]
                            indel_chrom = chromNow
                            indel_pos = int(indel.split(":")[0])
                            indel_size = int(indel.split(":")[1])
                            NMs = [seq.get_tag("NM") for seq in readSet]
                            if indel_size < 0:
                                indel_ref = nums2str(
                                    ref_np[
                                        indel_pos
                                        - reference_mat_start : indel_pos
                                        - reference_mat_start
                                        - indel_size
                                        + 1
                                    ]
                                ).upper()
                                indel_alt = nums2str(
                                    ref_np[[indel_pos - reference_mat_start]]
                                ).upper()
                            else:
                                indel_ref = nums2str(
                                    ref_np[[indel_pos - reference_mat_start]]
                                ).upper()
                                indel_alt = indel_ref + indel.split(":")[2]
                            indel_str = (
                                str(indel_chrom)
                                + ":"
                                + str(indel_pos)
                                + str(indel_ref)
                                + ":"
                                + str(indel_alt)
                            )
                            readPos = indel_pos - readSet[0].reference_start
                            if not processed_flag:
                                readPos5p = min(
                                    readPos + 1,
                                    abs(readSet[0].template_length) - readPos,
                                )
                                readPos3p = min(
                                    abs(readSet[0].reference_length) - readPos,
                                    readPos + 1,
                                )
                            else:
                                readPos5p = min(
                                    readSet[0].reference_length - readPos,
                                    abs(readSet[0].template_length)
                                    - readSet[0].reference_length
                                    + readPos,
                                )
                                readPos3p = min(
                                    abs(readSet[0].reference_length) - readPos,
                                    readPos + 1,
                                )
                            if (
                                F1R2_alt_count[pass_inds[nn]]
                                + F1R2_ref_count[pass_inds[nn]]
                                == 0
                            ):
                                continue
                            if (
                                F2R1_alt_count[pass_inds[nn]]
                                + F2R1_ref_count[pass_inds[nn]]
                                == 0
                            ):
                                continue
                            indel_rec = {
                                "chrom": chromNow,
                                "pos": indel_pos + 1,
                                "ref": indel_ref,
                                "alt": indel_alt,
                                "infos": {
                                    "F1R2": int(
                                        F1R2_alt_count[pass_inds[nn]]
                                        + F1R2_ref_count[pass_inds[nn]]
                                    ),
                                    "F2R1": int(
                                        F2R1_alt_count[pass_inds[nn]]
                                        + F2R1_ref_count[pass_inds[nn]]
                                    ),
                                    # "LR": LR[pass_inds[0]],
                                    "LR": LR[pass_inds[0]],
                                    # "BLR": F2R1_LR[pass_inds[0]],
                                    "TC": ",".join(
                                        [
                                            str(F1R2_alt_count[pass_inds[nn]]),
                                            str(F1R2_ref_count[pass_inds[nn]]),
                                        ]
                                    ),
                                    "BC": ",".join(
                                        [
                                            str(F2R1_alt_count[pass_inds[nn]]),
                                            str(F2R1_ref_count[pass_inds[nn]]),
                                        ]
                                    ),
                                    "TAG1": setBc[0],
                                    "TAG2": setBc[1],
                                    "SP": currentStart,
                                    "DF": readPos5p,
                                    "DR": readPos3p,
                                    "TN": ".",
                                    "HP": hps[pass_inds[nn]],
                                },
                                "formats": ["AC", "RC", "DP"],
                                # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_indels.append(indel_rec)
                            indel_dict[indel_str] = 1
                    # else:
                    ### Calculate genotype probability
                    if not any(indel_bool) or isLearn:
                        # if 1:
                        if isLearn:
                            (
                                mismatch_now,
                                indelerr_now,
                                mismatch_dmg_now,
                                indel_dmg_now,
                            ) = profileTriNucMismatches(
                                readSet,
                                ref_np[start_ind:end_ind],
                                trinuc_np[start_ind:end_ind],
                                hp_np[0, start_ind:end_ind],
                                np.copy(antimask),
                                params,
                            )
                            mismatch_mat += mismatch_now
                            indelerr_mat += indelerr_now
                            mismatch_dmg_mat += mismatch_dmg_now
                            indel_dmg_mat += indel_dmg_now
                            continue
                        (
                            LR,
                            b1_int,
                            antimask,
                            F1R2_count,
                            F2R1_count,
                        ) = genotypeDSSnv(
                            readSet,
                            ref_np[start_ind:end_ind],
                            trinuc_np[start_ind:end_ind],
                            prior_mat[start_ind:end_ind, :],
                            np.copy(antimask),
                            params,
                        )
                        ref_int = ref_np[start_ind:end_ind]
                        refs_ind = np.nonzero(
                            np.logical_and(
                                LR
                                <= params["pcutoff"],  # - np.log10(params["mutRate"]),
                                b1_int == ref_int,
                            )
                        )[0].tolist()
                        muts_ind = np.nonzero(
                            np.logical_and(LR <= params["pcutoff"], b1_int != ref_int)
                        )[0].tolist()
                        alt_int = b1_int
                        pass_bool = np.full(LR.size, False, dtype=bool)
                        pass_bool[refs_ind] = True
                        pass_bool[muts_ind] = True
                        # if F1R2 >=2 and F2R1 >=2:
                        # print("pass_bool",pass_bool,"antimask",antimask,"LR",LR,F1R2_count,F2R1_count)
                        pos = [
                            mut_ind + start_ind + reference_mat_start
                            for mut_ind in muts_ind
                        ]
                        # muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                        mut_positions = [
                            mut_ind + start_ind + reference_mat_start + 1
                            for mut_ind in muts_ind
                        ]
                        NMs = [seq.get_tag("NM") for seq in readSet]
                        averageNM = sum(NMs) / len(NMs)
                        # if averageNM - len(mut_positions) > 1:
                        # continue
                        # if len(mut_positions) > params["maxMnv"]:
                        # continue
                        for nn in range(len(mut_positions)):
                            mut_chrom = reference_mat_chrom
                            mut_pos = mut_positions[nn]
                            mut_ref = num2base[ref_int[muts_ind[nn]]]
                            mut_alt = num2base[alt_int[muts_ind[nn]]]
                            mut_trinuc = num2trinuc[
                                trinuc_np[start_ind:end_ind][muts_ind[nn]]
                            ]
                            if not processed_flag:
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
                            if F1R2_count[:, muts_ind[nn]].sum() == 0:
                                continue
                            if F2R1_count[:, muts_ind[nn]].sum() == 0:
                                continue
                            mut = {
                                "chrom": mut_chrom,
                                "pos": mut_pos,
                                "ref": mut_ref,
                                "alt": mut_alt,
                                "infos": {
                                    "F1R2": F1R2,
                                    "F2R1": F2R1,
                                    "LR": LR[muts_ind[nn]],
                                    # "BLR": F2R1_LR[muts_ind[nn]],
                                    # "LR": LR[muts_ind[nn]],
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
                                    "TAG1": setBc[0],
                                    "TAG2": setBc[1],
                                    "SP": currentStart,
                                    "DF": readPos5p,
                                    "DR": readPos3p,
                                    "TN": mut_trinuc,
                                    "HP": "0",
                                },
                                "formats": ["AC", "RC", "DP"],
                                # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_dict[
                                "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                            ] = 0
                            muts.append(mut)
                        """
                        if isLearn:
                            continue
                        """
                        coverage[start_ind:end_ind][pass_bool] += 1
                        duplex_read_num_dict[duplex_no][1] += np.count_nonzero(
                            pass_bool
                        )
                        trinuc_pass = trinuc_np[start_ind:end_ind][pass_bool]
                        duplex_read_num_dict_trinuc[duplex_no] += np.bincount(
                            trinuc_pass, minlength=96
                        ).astype(int)
                        if feature_beds:
                            for nn in range(feature_beds):
                                trinuc_pass_feature = trinuc_np[start_ind:end_ind][
                                    pass_bool[feature_mat[nn, :]]
                                ]
                                duplex_read_num_dict_trinuc[nn][
                                    duplex_no
                                ] += np.bincount(
                                    trinuc_pass_feature, minlength=96
                                ).astype(
                                    int
                                )
                        duplex_read_num_dict[duplex_no][0] += 1
                        duplex_count += 1
            """
            Calling block ends
            """
            currentReadDict = {
                label: {
                    "seqs": [rec],
                    "F1R2": 0,
                    "F2R1": 0,
                    "names": {rec.query_name: 0},
                }
            }
            if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
                currentReadDict[label]["F1R2"] += 1
            else:
                currentReadDict[label]["F2R1"] += 1
            currentStart = start
    """
    Calling block starts
    """
    for key in currentReadDict.keys():
        readSet = currentReadDict[key]["seqs"]
        all_dup = True
        for _ in readSet:
            if not _.is_duplicate:
                all_dup = False
                break
        if all_dup:
            continue
        if params["maxNM"]:
            blacklist_num = 0
            for seq in readSet:
                if seq.query_name in read_blacklist:
                    blacklist_num += 1
            if blacklist_num / len(readSet) >= 0.5:
                continue
        mean_mapq = sum([seq.mapping_quality for seq in readSet]) / len(readSet)
        if mean_mapq < params["mapq"]:
            continue
        meanASXS = sum(
            [seq.get_tag("AS") - seq.get_tag("XS") for seq in readSet]
        ) / len(readSet)
        if meanASXS < 50:
            continue
        setBc = key.split(":")[0].split("+")
        setBc1 = setBc[0]
        setBc2 = setBc[1]
        F2R1 = currentReadDict[key]["F2R1"]
        F1R2 = currentReadDict[key]["F1R2"]
        duplex_no = f"{F1R2}+{F2R1}"
        if duplex_read_num_dict.get(duplex_no) is None:
            duplex_read_num_dict[duplex_no] = [0, 0]
            duplex_read_num_dict_trinuc[duplex_no] = np.zeros(96, dtype=int)
        if feature_beds:
            for nn in range(len(duplex_read_num_dict_trinuc_features)):
                if duplex_read_num_dict[nn].get(duplex_no) is None:
                    duplex_read_num_dict_trinuc_features[nn][duplex_no] = np.zeros(
                        96, dtype=int
                    )
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
                    readSet[0].reference_start + 1000000,
                    max(region_end, max([seq.reference_end for seq in readSet])),
                    contig_end,
                )
                #
                ref_h5 = h5py.File(params["reference"] + ".ref.h5", "r")
                tn_h5 = h5py.File(params["reference"] + ".tn.h5", "r")
                hp_h5 = h5py.File(params["reference"] + ".hp.h5", "r")
                ref_np = ref_h5[reference_mat_chrom][
                    reference_mat_start:reference_mat_end
                ]
                trinuc_np = tn_h5[reference_mat_chrom][
                    reference_mat_start:reference_mat_end
                ]
                hp_np = hp_h5[reference_mat_chrom][
                    :, reference_mat_start:reference_mat_end
                ]
                (
                    prior_mat,
                    snp_mask,
                    indel_mask,
                    noise_mask,
                    n_cov_mask,
                    include_mask,
                    feature_mat
                    # ref_np,
                    # trinuc_np
                ) = prepare_reference_mats(
                    reference_mat_chrom,
                    reference_mat_start,
                    reference_mat_end,
                    # current_fasta,
                    ref_np,
                    trinuc_np,
                    germline,
                    noise,
                    indel_bed,
                    include_bed,
                    feature_beds,
                    nbams,
                    bam,
                    params,
                )
                # print(ref_np,reference_mat_start)
                coverage = np.zeros(1000000, dtype=int)
                coverage_indel = np.zeros(1000000, dtype=int)
                if feature_beds:
                    trinuc_features = [
                        np.zeros((96, 1000000), dtype=int) for _ in feature_beds
                    ]
                else:
                    trinuc_features = list()

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
            reference_length_max = max([read.reference_length for read in readSet])
            end_ind = (
                readSet[0].reference_start + reference_length_max - reference_mat_start
            )

            end_ind_max = (
                readSet[0].reference_start + reference_length_max - reference_mat_start
            )
            masks = np.zeros([5, end_ind - start_ind], dtype=bool)
            masks[0, :] = snp_mask[start_ind:end_ind]
            masks[1, :] = noise_mask[start_ind:end_ind]
            masks[2, :] = n_cov_mask[start_ind:end_ind]
            masks[3, :] = include_mask[start_ind:end_ind]
            left, right = determineTrimLength(
                readSet[0], params=params, processed_flag=processed_flag
            )
            masks[4, :left] = True
            masks[4, -right:] = True
            antimask = np.all(~masks, axis=0)
            antimask[trinuc_np[start_ind:end_ind] > 64] = False
            ### If the whole reads are masked:
            if not np.any(antimask):
                continue
            indel_bool = [
                ("I" in seq.cigarstring or "D" in seq.cigarstring) for seq in readSet
            ]
            # if any(indel_bool):
            if not isLearn:
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
                    LR,
                    indels,
                    hps,
                    F1R2_ref_count,
                    F1R2_alt_count,
                    F2R1_ref_count,
                    F2R1_alt_count,
                ) = genotypeDSIndel(
                    readSet,
                    tumorBam,
                    antimask_indel,
                    hp_np[0, start_ind:end_ind],
                    params,
                )
                """
                DCS = np.logical_and(
                    F1R2_BLR >= params["pcutoff"],
                    F2R1_BLR >= params["pcutoff"],
                )
                """
                pass_inds = np.nonzero(LR <= params["pcutoff"])[0].tolist()
                indels_pass = [indels[_] for _ in pass_inds]
                coverage_indel[start_ind:end_ind_max][antimask_indel] += 1
                for nn in range(len(indels_pass)):
                    indel = indels_pass[nn]
                    indel_chrom = chromNow
                    indel_pos = int(indel.split(":")[0])
                    indel_size = int(indel.split(":")[1])
                    NMs = [seq.get_tag("NM") for seq in readSet]
                    if indel_size < 0:
                        indel_ref = nums2str(
                            ref_np[
                                indel_pos
                                - reference_mat_start : indel_pos
                                - reference_mat_start
                                - indel_size
                                + 1
                            ]
                        ).upper()
                        indel_alt = nums2str(
                            ref_np[[indel_pos - reference_mat_start]]
                        ).upper()
                    else:
                        indel_ref = nums2str(
                            ref_np[[indel_pos - reference_mat_start]]
                        ).upper()
                        indel_alt = indel_ref + indel.split(":")[2]
                    indel_str = (
                        str(indel_chrom)
                        + ":"
                        + str(indel_pos)
                        + str(indel_ref)
                        + ":"
                        + str(indel_alt)
                    )
                    readPos = indel_pos - readSet[0].reference_start
                    if not processed_flag:
                        readPos5p = min(
                            readPos + 1,
                            abs(readSet[0].template_length) - readPos,
                        )
                        readPos3p = min(
                            abs(readSet[0].reference_length) - readPos,
                            readPos + 1,
                        )
                    else:
                        readPos5p = min(
                            readSet[0].reference_length - readPos,
                            abs(readSet[0].template_length)
                            - readSet[0].reference_length
                            + readPos,
                        )
                        readPos3p = min(
                            abs(readSet[0].reference_length) - readPos,
                            readPos + 1,
                        )
                    if (
                        F1R2_alt_count[pass_inds[nn]] + F1R2_ref_count[pass_inds[nn]]
                        == 0
                    ):
                        continue
                    if (
                        F2R1_alt_count[pass_inds[nn]] + F2R1_ref_count[pass_inds[nn]]
                        == 0
                    ):
                        continue
                    indel_rec = {
                        "chrom": chromNow,
                        "pos": indel_pos + 1,
                        "ref": indel_ref,
                        "alt": indel_alt,
                        "infos": {
                            "F1R2": int(
                                F1R2_alt_count[pass_inds[nn]]
                                + F1R2_ref_count[pass_inds[nn]]
                            ),
                            "F2R1": int(
                                F2R1_alt_count[pass_inds[nn]]
                                + F2R1_ref_count[pass_inds[nn]]
                            ),
                            # "LR": LR[pass_inds[0]],
                            "LR": LR[pass_inds[0]],
                            # "BLR": F2R1_LR[pass_inds[0]],
                            "TC": ",".join(
                                [
                                    str(F1R2_alt_count[pass_inds[nn]]),
                                    str(F1R2_ref_count[pass_inds[nn]]),
                                ]
                            ),
                            "BC": ",".join(
                                [
                                    str(F2R1_alt_count[pass_inds[nn]]),
                                    str(F2R1_ref_count[pass_inds[nn]]),
                                ]
                            ),
                            "TAG1": setBc[0],
                            "TAG2": setBc[1],
                            "SP": currentStart,
                            "DF": readPos5p,
                            "DR": readPos3p,
                            "TN": ".",
                            "HP": hps[pass_inds[nn]],
                        },
                        "formats": ["AC", "RC", "DP"],
                        # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_indels.append(indel_rec)
                    indel_dict[indel_str] = 1
            # else:
            ### Calculate genotype probability
            # if 1:
            if not any(indel_bool) or isLearn:  # or len(indels_pass) == 0:
                if isLearn:
                    (
                        mismatch_now,
                        indelerr_now,
                        mismatch_dmg_now,
                        indel_dmg_now,
                    ) = profileTriNucMismatches(
                        readSet,
                        ref_np[start_ind:end_ind],
                        trinuc_np[start_ind:end_ind],
                        hp_np[0, start_ind:end_ind],
                        np.copy(antimask),
                        params,
                    )
                    mismatch_mat += mismatch_now
                    indelerr_mat += indelerr_now
                    mismatch_dmg_mat += mismatch_dmg_now
                    indel_dmg_mat += indel_dmg_now
                    continue
                (
                    LR,
                    b1_int,
                    antimask,
                    F1R2_count,
                    F2R1_count,
                ) = genotypeDSSnv(
                    readSet,
                    ref_np[start_ind:end_ind],
                    trinuc_np[start_ind:end_ind],
                    prior_mat[start_ind:end_ind, :],
                    np.copy(antimask),
                    params,
                )
                ref_int = ref_np[start_ind:end_ind]
                refs_ind = np.nonzero(
                    np.logical_and(
                        LR <= params["pcutoff"],  # - np.log10(params["mutRate"]),
                        b1_int == ref_int,
                    )
                )[0].tolist()
                muts_ind = np.nonzero(
                    np.logical_and(LR <= params["pcutoff"], b1_int != ref_int)
                )[0].tolist()
                alt_int = b1_int
                pass_bool = np.full(LR.size, False, dtype=bool)
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
                NMs = [seq.get_tag("NM") for seq in readSet]
                averageNM = sum(NMs) / len(NMs)
                # if averageNM - len(mut_positions) > 1:
                # continue
                # if len(mut_positions) > params["maxMnv"]:
                # continue
                for nn in range(len(mut_positions)):
                    mut_chrom = reference_mat_chrom
                    mut_pos = mut_positions[nn]
                    mut_ref = num2base[ref_int[muts_ind[nn]]]
                    mut_alt = num2base[alt_int[muts_ind[nn]]]
                    mut_trinuc = num2trinuc[trinuc_np[start_ind:end_ind][muts_ind[nn]]]
                    if not processed_flag:
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
                    if F1R2_count[:, muts_ind[nn]].sum() == 0:
                        continue
                    if F2R1_count[:, muts_ind[nn]].sum() == 0:
                        continue
                    mut = {
                        "chrom": mut_chrom,
                        "pos": mut_pos,
                        "ref": mut_ref,
                        "alt": mut_alt,
                        "infos": {
                            "F1R2": F1R2,
                            "F2R1": F2R1,
                            "LR": LR[muts_ind[nn]],
                            # "BLR": F2R1_LR[muts_ind[nn]],
                            # "LR": LR[muts_ind[nn]],
                            "TC": ",".join(
                                [str(_) for _ in F1R2_count[:, muts_ind[nn]].tolist()]
                            ),
                            "BC": ",".join(
                                [str(_) for _ in F2R1_count[:, muts_ind[nn]].tolist()]
                            ),
                            "TAG1": setBc[0],
                            "TAG2": setBc[1],
                            "SP": currentStart,
                            "DF": readPos5p,
                            "DR": readPos3p,
                            "TN": mut_trinuc,
                            "HP": "0",
                        },
                        "formats": ["AC", "RC", "DP"],
                        # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_dict["_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])] = 0
                    muts.append(mut)
                """
                if isLearn:
                    continue
                """
                coverage[start_ind:end_ind][pass_bool] += 1
                duplex_read_num_dict[duplex_no][1] += np.count_nonzero(pass_bool)
                trinuc_pass = trinuc_np[start_ind:end_ind][pass_bool]
                duplex_read_num_dict_trinuc[duplex_no] += np.bincount(
                    trinuc_pass, minlength=96
                ).astype(int)
                if feature_beds:
                    for nn in range(feature_beds):
                        trinuc_pass_feature = trinuc_np[start_ind:end_ind][
                            pass_bool[feature_mat[nn, :]]
                        ]
                        duplex_read_num_dict_trinuc_features[nn][
                            duplex_no
                        ] += np.bincount(trinuc_pass_feature, minlength=96).astype(int)
                duplex_read_num_dict[duplex_no][0] += 1
                duplex_count += 1
    """
    Calling block ends
    """
    if isLearn:
        return mismatch_mat, indelerr_mat, mismatch_dmg_mat, indel_dmg_mat
    mut_dict = dict()
    mut_pass_filter = []
    for mut in muts:
        chrom = mut["chrom"]
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]
        bc1 = mut["infos"]["TAG1"]
        bc2 = mut["infos"]["TAG2"]
        readStartCoord = mut["infos"]["SP"]

        if not mut_dict.get(":".join([chrom, str(pos), ref, alt])):
            ta, tr, ti, tdp = extractDepthSnv(
                tumorBam, chrom, pos, ref, alt, params, minbq=params["minBq"]
            )
            # overlap_error = detectOverlapDiscord(
            # tumorBam, chrom, pos, ref, alt, params, bc1, bc2, readStartCoord
            # )
            # window_filter = False
            # if IndelFilterByWindows(tumorBam, chrom, pos, 3, params):
            # window_filter = True
            if normalBams:
                na = 0
                nr = 0
                ni = 0
                ndp = 0
                for normalBam in normalBams:
                    na_now, nr_now, ni_now, ndp_now = extractDepthSnv(
                        normalBam, chrom, pos, ref, alt, params, minbq=params["minBq"]
                    )
                    na += na_now
                    nr += nr_now
                    ni += ni_now
                    ndp += ndp_now
                # if IndelFilterByWindows(normalBam, chrom, pos, 3, params):
                # window_filter = True
            else:
                na, nr, ni, ndp = (0, 0, 0, 0)
            mut_dict[":".join([chrom, str(pos), ref, alt])] = (
                ta,
                tr,
                ti,
                tdp,
                # overlap_error,
                na,
                nr,
                ni,
                ndp,
                # window_filter,
            )
        else:
            ta, tr, ti, tdp, na, nr, ni, ndp = mut_dict[
                ":".join([chrom, str(pos), ref, alt])
            ]
        # if window_filter:
        # continue
        # if ta > params["maxAltCount"]:
        # continue
        if ta == 0:
            continue
        if ta / tdp > params["maxAF"]:
            continue
        if ti >= 1:
            continue
        # if overlap_error:
        # continue
        if normalBams:
            if na > 0:
                continue
            if ndp < params["minNdepth"]:
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
            ta, tr, ti, tdp = extractDepthIndel(tumorBam, chrom, pos, ref, alt, params)
            # window_filter = False
            # if IndelFilterByWindows(tumorBam, chrom, pos, 3, params):
            # window_filter = True
            if normalBams:
                na = 0
                nr = 0
                ni = 0
                ndp = 0
                for normalBam in normalBams:
                    na_now, nr_now, ni_now, ndp_now = extractDepthIndel(
                        normalBam, chrom, pos, ref, alt, params, minbq=params["minBq"]
                    )
                    na += na_now
                    nr += nr_now
                    ni += ni_now
                    ndp += ndp_now
                # if IndelFilterByWindows(normalBam, chrom, pos, 3, params):
                # window_filter = True
            else:
                na, nr, ndp = (0, 0, 0)
            muts_indels_dict[":".join([chrom, str(pos), ref, alt])] = (
                ta,
                tr,
                tdp,
                na,
                nr,
                ndp,
                # window_filter,
            )
        else:
            ta, tr, tdp, na, nr, ndp = muts_indels_dict[
                ":".join([chrom, str(pos), ref, alt])
            ]
        # if window_filter:
        # continue
        # if ta > params["maxAltCount"]:
        # continue
        if ta == 0:
            continue
        if ti > 0:
            continue
        if ta / tdp > params["maxAF"]:
            continue
        if normalBams:
            if na > 0:
                continue
            if ndp < params["minNdepth"]:
                continue
            if ni > 0:
                continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        muts_indels_pass_filter.append(mut)

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
        f"Process {processNo} finished in {(time.time()-starttime)/60: .2f} minutes and processed {recCount} reads"
    )

    for duplex_no in duplex_read_num_dict_trinuc.keys():
        trinuc_profile = (
            duplex_read_num_dict_trinuc[duplex_no][:32]
            + duplex_read_num_dict_trinuc[duplex_no][32:64]
        )
        duplex_read_num_dict_trinuc[duplex_no] = trinuc_profile

    if feature_beds:
        for nn in range(feature_beds):
            for duplex_no in duplex_read_num_dict_trinuc.keys():
                trinuc_profile = (
                    duplex_read_num_dict_trinuc[nn][duplex_no][:32]
                    + duplex_read_num_dict_trinuc[nn][duplex_no][32:64]
                )
                duplex_read_num_dict_trinuc[nn][duplex_no] = trinuc_profile

    return (
        mut_pass_filter,
        total_coverage,
        recCount,
        duplex_count,
        duplex_read_num_dict,
        duplex_read_num_dict_trinuc,
        muts_indels_pass_filter,
        total_coverage_indel,
        unique_read_num,
        pass_read_num,
        FPs,
        RPs,
    )
