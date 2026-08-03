#!/usr/bin/env python3

from Bio import bgzf

import time
import math
import numpy as np
import pandas as pd
from pysam import FastaFile as FASTA
from pysam import TabixFile as BED
from pysam import VariantFile as VCF
import pysam
import os
import re
import h5py
import errno
import copy

from .depth import (
    extractDepthRegion,
    extractDepthBatchSnv,
    extractDepthBatchIndel,
    extractDepthBatchDbs,
    prepareAlignMask,
    detectOverlapDiscord,
)
from .prob import (
    genotypeDSSnv,
    genotypeDSIndel,
    indelErrorProbs,
    calculateSSPosterior,
    calculateDSPosterior,
)
from .learn import profileTriNucMismatches
from .misc import getAlignmentObject as BAM
from .misc import build_trinuc64_order
from .misc import load_repeat_context
from .misc import log_progress
from .misc import build_dbs_raw144_labels, build_dbs_raw144_index_grid
from .indels import findIndels

# [4,4,4,4] (ref1_num, ref2_num, alt1_num, alt2_num) -> row index into the
# 144-raw-class DBS opportunity space, all in the same "ATCG" 0-3 encoding
# ref_int/alt_int already use below (base2num/num2base) -- built once at
# import time (cheap: 144 valid entries out of 256) rather than per-call,
# and reused by every duplex-group DBS-opportunity accumulation below.
_DBS_RAW144_IDX_GRID = build_dbs_raw144_index_grid()
_, DBS_RAW144_LABELS = build_dbs_raw144_labels()

# from .prob_old import genotypeDSSnv,genotypeDSIndel# as oldDSSnv

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
    nm_mask = np.full(m, False, dtype=bool)

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
        for n in noise_bed:
            for rec in n.fetch(chrom, start, end, parser=pysam.asBed()):
                interval_start = max(rec.start, start)
                interval_end = min(rec.end, end)
                interval_len = interval_end - interval_start
                interval_start_ind = interval_start - start
                noise_mask[
                    interval_start_ind : interval_start_ind + interval_len
                ] = True
    if indel_bed != None:
        for rec in indel_bed.fetch(chrom, start, end, parser=pysam.asBed()):
            interval_start = max(rec.start, start)
            interval_end = min(rec.end, end)
            interval_len = interval_end - interval_start
            interval_start_ind = interval_start - start
            indel_mask[interval_start_ind : interval_start_ind + interval_len] = True
    if include_bed != None:
        include_arr = np.zeros(end - start, dtype=bool)
        if chrom in include_bed.contigs:
            for rec in include_bed.fetch(chrom, start, end, parser=pysam.asBed()):
                interval_start = max(rec.start, start)
                interval_end = min(rec.end, end)
                interval_len = interval_end - interval_start
                interval_start_ind = interval_start - start
                include_arr[
                    interval_start_ind : interval_start_ind + interval_len
                ] = True
        include_mask = ~include_arr
    else:
        include_mask = np.zeros(end - start, dtype=bool)

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

        nm_mask = prepareAlignMask(tbam, chrom, start, end, params)
        # nm_mask = nm_avg >= params["maxNM"]/2

    return (
        prior_mat,
        snp_mask,
        indel_mask,
        noise_mask,
        n_cov_mask,
        include_mask,
        nm_mask,
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
                seq.reference_length
                + mate_reference_length
                - abs(seq.template_length)
                - params["trim3"],
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


def get_bed_file_for_position(
    pos,
    chrom,
    regions_start_chrom,
    regions_start_pos,
    regions_end_chrom,
    regions_end_pos,
    locus_bed,
    locus_bed_prev,
    locus_bed_next,
):
    """
    Determine which bed file to write to based on position relative to region boundaries
    Only compare positions within the same chromosome
    """
    if chrom == regions_start_chrom and pos < regions_start_pos:
        return locus_bed_prev
    elif chrom == regions_end_chrom and pos > regions_end_pos:
        return locus_bed_next
    else:
        return locus_bed


def bamIterateMultipleRegion(bam, regions, ref):  # , regionFile):
    bamObject = BAM(bam, "rb", ref)
    # if not regionFile:
    for region in regions:
        for rec in bamObject.fetch(*region):
            if len(region) >= 2:
                if rec.reference_start < region[1]:
                    continue
            yield rec, region
    """
    else:
        for region in regions:
            chrom = region[0]
            if chrom in regionFile.contigs:
                for interval in regionFile.fetch(*region):
                    for rec in bamObject.fetch(interval.contig, interval.start, interval.end):
                        if len(region) >= 2:
                            if rec.reference_start < region[1]:
                                continue
                    yield rec, region
            else: continue
    """


def regularizeErrorMat(mat, minerr):
    """Replace invalid cells — exact zero, or NaN from a 0/0 row-sum
    division upstream (a trinuc/indel context with zero observed counts)
    — with the smallest valid value found elsewhere in the matrix, so one
    unobserved context can't leave a zero or NaN in the error model.
    Falls back to `minerr` only if the matrix has no valid cells at all."""
    invalid = (mat == 0) | np.isnan(mat)
    valid = mat[~invalid]
    mat_min = valid.min() if valid.size > 0 else minerr
    mat[invalid] = mat_min
    return mat


def _detect_dbs_pairs(
    mut_chrom,
    mut_positions,
    muts_ind,
    ref_int,
    alt_int,
    pass_bool,
    LR_pass_bool,
    flt_rs,
    F1R2,
    F2R1,
    setBc,
    currentStart,
    template_length,
    num2base,
):
    """Detect DBS (dinucleotide substitution) events within one duplex
    read group's SNV candidate list: two positions changing
    SIMULTANEOUSLY in the SAME physical molecule, not two independent
    SNVs that happen to sit next to each other in the genome. Because
    mut_positions/muts_ind/ref_int/alt_int all come from THIS one read
    group's own consensus calling, a genomically-adjacent pair found here
    is, by construction, supported by the same underlying reads — unlike
    scanning the already-written SNV VCF after the fact, which has no way
    to tell same-molecule adjacency from coincidence.

    Only pairs where BOTH constituent SNVs individually reach full PASS
    status are called as a DBS (recomputed here from pass_bool/
    LR_pass_bool/flt_rs — the same three inputs the existing per-position
    loop above already uses to decide flt — rather than threading a new
    list through that loop, so this can be added without touching any of
    its existing lines). Weaker DBS-specific filter tiers (mirroring
    "masked"/"underpowered" for SNVs) aren't attempted in this first cut.

    Returns a list of DBS mut dicts (ref/alt are 2-character
    dinucleotides), in the same shape as the SNV mut dict built above,
    minus the SNV-only INFO fields (CS/LR/LM/TC/BC/TN/HP/STR) that don't
    apply to a 2-base event.
    """

    def _flt(i):
        if pass_bool[i] and LR_pass_bool[i]:
            return flt_rs
        if pass_bool[i]:
            return "underpowered"
        return "masked"

    dbs_muts = []
    for k in range(len(mut_positions) - 1):
        if mut_positions[k + 1] != mut_positions[k] + 1:
            continue
        i0, i1 = muts_ind[k], muts_ind[k + 1]
        if _flt(i0) != "PASS" or _flt(i1) != "PASS":
            continue
        ref_dinuc = num2base[ref_int[i0]] + num2base[ref_int[i1]]
        alt_dinuc = num2base[alt_int[i0]] + num2base[alt_int[i1]]
        dbs_muts.append(
            {
                "chrom": mut_chrom,
                "pos": mut_positions[k],
                "ref": ref_dinuc,
                "alt": alt_dinuc,
                "filter": "PASS",
                "infos": {
                    "F1R2": F1R2,
                    "F2R1": F2R1,
                    "TAG1": setBc[0],
                    "TAG2": setBc[1],
                    "SP": currentStart,
                    "TL": template_length,
                },
                "formats": ["AC", "RC", "DP"],
                # No independent raw-BAM depth re-verification for DBS
                # yet (mirroring extractDepthSnv/extractDepthBatchSnv for
                # SNVs) — both constituent bases already passed full
                # duplex-consensus CS/LR calling, which is a strong
                # signal on its own, but there's currently no equivalent
                # of the SNV/indel maxAF / normal-VAF post-processing
                # check for DBS. Placeholder AC=1/RC=0/DP=1 (tumor),
                # 0/0/0 (normal) so downstream VCF writing has values to
                # write; a real depth-extraction pass would replace this.
                "samples": [[1, 0, 1], [0, 0, 0]],
            }
        )
    return dbs_muts


def _compute_dbs_opportunity(cov_mat, ref_int, pass_bool):
    """Per-duplex-group DBS opportunity contribution from one family's
    window: for every pair of directly ADJACENT in-window positions
    (i, i+1) that both individually pass calling (pass_bool[i] and
    pass_bool[i+1]) and have a valid A/T/C/G reference base at both,
    opportunity is coverage(alt1 at i) * coverage(alt2 at i+1) summed over
    all 9 non-ref (alt1, alt2) combinations -- the same product-of-
    coverages approximation Estimate.py's calculate_dbs_opportunity uses
    genome-wide from the coverage bed, just computed here per-family (so
    it can be attributed to that family's own duplex_no) from the same
    cov_mat/ref_int/pass_bool the trinuc accumulator right after this
    function's call site already uses -- window-local arrays, aligned the
    same way.

    cov_mat  : (window_len, 4) L-weighted per-alt-base coverage.
    ref_int  : (window_len,) reference base, 0-3=A/T/C/G, 4=N/invalid.
    pass_bool: (window_len,) bool, calling-eligible positions.
    Returns a (144,) raw-DBS-class contribution vector (funcs/misc.py's
    build_dbs_raw144_labels order), zero wherever no adjacent pass_bool
    pair with valid reference bases exists.
    """
    contribution = np.zeros(144)
    if cov_mat.shape[0] < 2:
        return contribution
    pair_ok = pass_bool[:-1] & pass_bool[1:]
    if not np.any(pair_ok):
        return contribution
    r1 = ref_int[:-1]
    r2 = ref_int[1:]
    valid = pair_ok & (r1 <= 3) & (r2 <= 3)
    if not np.any(valid):
        return contribution
    c1 = cov_mat[:-1][valid]  # (k, 4)
    c2 = cov_mat[1:][valid]  # (k, 4)
    pair_contrib = c1[:, :, None] * c2[:, None, :]  # (k, 4, 4)
    raw_idx = _DBS_RAW144_IDX_GRID[r1[valid], r2[valid]]  # (k, 4, 4)
    flat_idx = raw_idx.reshape(-1)
    keep = flat_idx >= 0
    contribution += np.bincount(
        flat_idx[keep], weights=pair_contrib.reshape(-1)[keep], minlength=144
    )
    return contribution


def callBam(params, processNo):
    # Get parameters
    bam = params["tumorBam"]
    nbams = params["normalBams"]
    regions = params["regions"]
    if len(regions[0]) == 1:
        regions_start_chrom = regions[0][0]
        regions_start_pos = 0
    else:
        regions_start_chrom = regions[0][0]
        original_start_pos = regions[0][1]

        # If starting position is 0, skip modification and keep as 0
        if original_start_pos == 0:
            regions_start_pos = 0
        else:
            # Determine regions_start_pos using the specified approach:
            # 1) Find locus 1bp before the starting position
            target_position = original_start_pos - 1

            # 2) Fetch reads that contain that position
            # 3) Find the last position in those reads
            # 4) Define regions_start_pos as position 1bp after that last position
            try:
                bamObject = BAM(bam, "rb", params.get("reference"))
                max_reference_end = 0

                # Fetch reads overlapping the target position
                for read in bamObject.fetch(
                    regions_start_chrom, target_position, target_position + 1
                ):
                    if not read.is_unmapped and read.reference_end is not None:
                        max_reference_end = max(max_reference_end, read.reference_end)

                # Set regions_start_pos as 1bp after the last position in those reads
                if max_reference_end > 0:
                    regions_start_pos = max_reference_end + 1
                else:
                    # Fallback to original position if no reads found
                    regions_start_pos = original_start_pos

                bamObject.close()
            except Exception:
                # Fallback to original position if BAM access fails
                regions_start_pos = original_start_pos
    if len(regions[-1]) <= 2:
        regions_end_chrom = regions[-1][0]
        regions_end_pos = 10e10
    else:
        regions_end_chrom = regions[-1][0]
        regions_end_pos = regions[-1][2]
    if params["germline"]:
        germline = VCF(params["germline"], mode="rb")
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
    tmp_dir = params["tmp_dir"]
    sample_name = os.path.basename(params["output"])
    output = os.path.join(tmp_dir, sample_name + "_" + str(nn))
    sample_dir = tmp_dir
    if not os.path.exists(sample_dir):
        try:
            os.makedirs(sample_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    if params["noise"]:
        noise = list()
        for bb in params["noise"]:
            noise.append(BED(bb))
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

    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    muts = []
    muts_dict = dict()
    muts_indels = []
    muts_dbs = []
    duplex_read_num_dict = dict()
    duplex_read_num_dict_trinuc = dict()
    duplex_read_num_dict_indel = dict()
    duplex_read_num_dict_dbs = dict()
    # unmasked_duplex_read_num_dict_trinuc = dict()
    unique_read_num = 0
    pass_read_num = 0
    FPs = []
    RPs = []
    indel_dict = dict()
    mismatch_mat = np.zeros([64, 4])
    indelerr_mat = np.zeros([23, 23])
    mismatch_dmg_mat = np.zeros([64, 4])
    indel_dmg_mat = np.zeros([23, 23])
    trinuc2num, num2trinuc = build_trinuc64_order()

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
    if isLearn:
        ampmat = np.zeros([64, 4])
    else:
        ampmat = (
            pd.read_csv(params["amperr_file"], sep="\t", index_col=0)
            .to_numpy()
            .astype(float)
        )
    # ampmat += 0.5
    # A trinuc context with zero observed counts across all 4 alt bases
    # has row sum 0; dividing those rows would produce NaN (0/0) instead
    # of leaving them as exact zero, which regularizeErrorMat below can
    # then correctly patch with the matrix's smallest valid value.
    ampmat_row_sum = ampmat.sum(axis=1, keepdims=True)
    ampmat_row_nonzero = (ampmat_row_sum != 0).flatten()
    ampmat_normalized = np.zeros_like(ampmat)
    ampmat_normalized[ampmat_row_nonzero, :] = (
        ampmat[ampmat_row_nonzero, :] / ampmat_row_sum[ampmat_row_nonzero, :]
    )
    ampmat = ampmat_normalized
    # ampmat_avg_error = (1 - ampmat.max(axis=1,keepdims=True))/3
    ampmat_min_error = ampmat.min(axis=1, keepdims=True)
    ampmat = np.concatenate([ampmat, ampmat_min_error], axis=1)
    ampmat = regularizeErrorMat(ampmat, 1e-6)
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

    if isLearn:
        ampmat_indel = np.zeros([23, 23])
    else:
        ampmat_indel = np.loadtxt(params["amperri_file"], delimiter="\t", dtype=float)
        # ampmat_indel += 0.5
        ampmats = list()
        ampmats.append(ampmat_indel[:, 0:11])
        ampmats.append(ampmat_indel[:, 11:14])
        ampmats.append(ampmat_indel[:, 14:17])
        ampmats.append(ampmat_indel[:, 17:20])
        ampmats.append(ampmat_indel[:, 20:23])
        ampmats_new = list()
        ampmats_new_rev = list()
        for mm, mat in enumerate(ampmats):
            row_non_zero = np.sum(mat, axis=1, keepdims=False) != 0
            mat_new = np.zeros(mat.shape)
            mat_new[row_non_zero, :] = mat[row_non_zero, :] / np.sum(
                mat[row_non_zero, :], axis=1, keepdims=True
            )
            for nn in range(1, 20):
                current_row = mat_new[nn, :]
                smaller_entries = current_row <= mat_new[nn - 1, :]
                current_row[smaller_entries] = mat_new[nn - 1, :][smaller_entries]
                mat_new[nn, :] = current_row
            if mm == 0:
                str_mat = mat[20:23, 0:11]
                str_mat_total = np.vstack(
                    [str_mat.sum(axis=0) / str_mat.sum(axis=0).sum()] * 3
                )
                mat_new[20:23, 0:11][str_mat == 0] = str_mat_total[str_mat == 0]
            ampmats_new.append(mat_new.copy())
            ampmats_new_rev.append(np.fliplr(mat_new.copy()))
        ampmat_indel = np.hstack(ampmats_new)
        ampmat_indel = regularizeErrorMat(ampmat_indel, 1e-6)
        params["ampmat_indel"] = ampmat_indel

        ampmat_indel_rev = np.hstack(ampmats_new_rev)
        ampmat_indel_rev = regularizeErrorMat(ampmat_indel_rev, 1e-6)
        params["ampmat_indel_rev"] = ampmat_indel_rev

    # params["ampmat_indel_mean"] = np.mean(ampmat_indel,axis=1)
    # params["ampmat_indel_rev_mean"] = np.mean(ampmat_indel_rev,axis=1)

    ### Load damage matrix
    if isLearn:
        dmgmat = np.zeros([64, 4])
    else:
        dmgmat = (
            pd.read_csv(params["dmgerr_file"], sep="\t", index_col=0)
            .to_numpy()
            .astype(float)
        )
        # dmgmat += 1

    # dmgmat += 0.5
    # Same zero-row guard as ampmat above: a trinuc context with zero
    # observed counts must stay exact zero here, not become NaN via 0/0,
    # so regularizeErrorMat below can patch it correctly.
    dmgmat_row_sum = dmgmat.sum(axis=1, keepdims=True)
    dmgmat_row_nonzero = (dmgmat_row_sum != 0).flatten()
    dmgmat_normalized = np.zeros_like(dmgmat)
    dmgmat_normalized[dmgmat_row_nonzero, :] = (
        dmgmat[dmgmat_row_nonzero, :] / dmgmat_row_sum[dmgmat_row_nonzero, :]
    )
    dmgmat = dmgmat_normalized
    # dmgmat_ref_error = 1 -  dmgmat.max(axis=1, keepdims=True)
    dmgmat_ref_error = dmgmat.min(axis=1, keepdims=True)
    dmgmat = np.concatenate([dmgmat, dmgmat_ref_error], axis=1)
    dmgmat = regularizeErrorMat(dmgmat, 1e-8)

    params["dmgmat_top"] = dmgmat
    params["trinuc2num_dict"] = trinuc2num

    dmgmat_rev = np.zeros([64, 4])
    dmgmat_rev_ref_error = np.zeros(64)
    for trinuc in trinuc2num.keys():
        refbase = trinuc[1]
        for nn, altbase in enumerate(["A", "T", "C", "G"]):
            dmgmat_rev[trinuc2num[trinuc], nn] = dmgmat[
                trinuc2num[trinuc[0] + altbase + trinuc[2]], base2num[refbase]
            ]
    # dmgmat_rev_ref_error = 1 -  dmgmat_rev.max(axis=1, keepdims=True)
    dmgmat_rev_ref_error = dmgmat_rev.min(axis=1, keepdims=True)
    dmgmat_rev = np.concatenate([dmgmat_rev, dmgmat_rev_ref_error], axis=1)
    params["dmgmat_rev_top"] = dmgmat_rev

    dmgmat_b = np.vstack((dmgmat[32:64, [1, 0, 3, 2]], dmgmat[:32, [1, 0, 3, 2]]))
    # dmgmat_b_ref_error = 1 - dmgmat_b.max(axis=1, keepdims=True)
    dmgmat_b_ref_error = dmgmat_b.min(axis=1, keepdims=True)
    dmgmat_rev_b = np.vstack(
        (dmgmat_rev[32:64, [1, 0, 3, 2]], dmgmat_rev[:32, [1, 0, 3, 2]])
    )
    dmgmat_rev_b_ref_error = dmgmat_rev_b.min(axis=1, keepdims=True)
    dmgmat_b = np.concatenate([dmgmat_b, dmgmat_b_ref_error], axis=1)
    dmgmat_rev_b = np.concatenate([dmgmat_rev_b, dmgmat_rev_b_ref_error], axis=1)
    params["dmgmat_bot"] = dmgmat_b
    params["dmgmat_rev_bot"] = dmgmat_rev_b

    if isLearn:
        dmgmat_indel = np.zeros([21, 23])
    else:
        # dmgmat_indel = np.loadtxt(params["dmgerri_file"], delimiter="\t")
        dmgmat_indel = np.loadtxt(params["dmgerri_file"], delimiter="\t", dtype=float)
        # dmgmat_indel+= 0.5
        dmgmats = list()
        dmgmats.append(dmgmat_indel[:, 0:11])
        dmgmats.append(dmgmat_indel[:, 11:14])
        dmgmats.append(dmgmat_indel[:, 14:17])
        dmgmats.append(dmgmat_indel[:, 17:20])
        dmgmats.append(dmgmat_indel[:, 20:23])
        dmgmats_new = list()
        dmgmats_new_rev = list()
        for mm, mat in enumerate(dmgmats):
            row_non_zero = np.sum(mat, axis=1, keepdims=False) != 0
            mat_new = np.zeros(mat.shape)
            mat_new[row_non_zero, :] = mat[row_non_zero, :] / np.sum(
                mat[row_non_zero, :], axis=1, keepdims=True
            )
            for nn in range(1, 20):
                current_row = mat_new[nn, :]
                smaller_entries = current_row <= mat_new[nn - 1, :]
                current_row[smaller_entries] = mat_new[nn - 1, :][smaller_entries]
                mat_new[nn, :] = current_row
            if mm == 0:
                str_mat = mat[20:23, 0:11]
                str_mat_total = np.vstack(
                    [str_mat.sum(axis=0) / str_mat.sum(axis=0).sum()] * 3
                )
                mat_new[20:23, 0:11][str_mat == 0] = str_mat_total[str_mat == 0]
            dmgmats_new.append(mat_new.copy())
            dmgmats_new_rev.append(np.fliplr(mat_new.copy()))
        dmgmat_indel = np.hstack(dmgmats_new)
        dmgmat_indel = regularizeErrorMat(dmgmat_indel, 1e-8)
        params["dmgmat_indel"] = dmgmat_indel
        dmgmat_indel_rev = np.hstack(dmgmats_new_rev)
        dmgmat_indel_rev = regularizeErrorMat(dmgmat_indel_rev, 1e-8)
        params["dmgmat_indel_rev"] = dmgmat_indel_rev
    # Initialize

    total_coverage = np.zeros(4)
    total_coverage_indel_cat = np.zeros(14)
    total_unmasked_coverage = np.zeros(4)
    total_unmasked_coverage_indel_cat = np.zeros(14)
    starttime = time.time()
    tumorBam = BAM(bam, "rb", params.get("reference"))
    if nbams:
        normalBams = list()
        for nbam in nbams:
            normalBams.append(BAM(nbam, "rb", params.get("reference")))
    else:
        normalBams = None
    currentStart = -1
    # currentReadDict = {}
    recCount = 0
    currentCheckPoint = 1000000
    lastTime = 0
    duplex_count = 0
    reference_mat_chrom = "anyChrom"
    reference_mat_start = 0
    locus_bed = bgzf.open(output + "_coverage.bed.gz", "wt")
    locus_bed_prev = bgzf.open(output + "_coverage_prev_region.tmp.bed.gz", "wt")
    locus_bed_next = bgzf.open(output + "_coverage_next_region.tmp.bed.gz", "wt")
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
        for rec, region in bamIterateMultipleRegion(
            bam, regions, params.get("reference")
        ):
            rec_num += 1
            if rec.query_name in read_blacklist or rec.is_unmapped:
                continue
            id_length = 0
            id_num = 0
            for cigar in rec.cigartuples:
                if cigar[0] == 1 or cigar[0] == 2:
                    id_length += cigar[1]
                    id_num += 1
            NM_no_id = rec.get_tag("NM") - id_length + id_num
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
    # Build 10x10x64x4 LR lookup matrix L for coverage/power calculation.
    # Axis 0: top strand (F1R2) read count 0-9
    # Axis 1: bottom strand (F2R1) read count 0-9
    # Axis 2: trinuc context index 0-63
    # Axis 3: converted (alt) base index 0-3 (A/T/C/G)
    L = np.zeros([10, 10, 64, 4])

    # Sample base quality distribution from the BAM regions
    all_quals = []
    reads_sampled = 0
    max_qual_reads = 1000
    for region in regions:
        for read in tumorBam.fetch(*region):
            if not read.is_unmapped and read.query_alignment_qualities is not None:
                quals = np.array(read.query_alignment_qualities, dtype=float)
                valid_quals = quals[quals > params["minBq"]]
                all_quals.extend(valid_quals.tolist())
                reads_sampled += 1
                if reads_sampled >= max_qual_reads:
                    break
        if reads_sampled >= max_qual_reads:
            break
    if not all_quals:
        all_quals = [30]
    all_quals = np.array(all_quals, dtype=float)

    prob_amp_mat_L = params["ampmat"]
    prob_amp_mat_rev_L = params["ampmat_rev"]
    prob_dmg_t_L = params["dmgmat_top"]
    prob_dmg_rev_t_L = params["dmgmat_rev_top"]
    prob_dmg_b_L = params["dmgmat_bot"]
    prob_dmg_rev_b_L = params["dmgmat_rev_bot"]
    trinuc_conv_np_L = params["trinuc_convert"]
    ln10_L = np.log(10)
    N_SIM = 100
    rng_L = np.random.default_rng()

    for t in range(64):
        ref_base_idx_L = base2num[num2trinuc[t][1]]
        for b in range(4):
            if b == ref_base_idx_L:
                # A base "changing" to itself isn't a mutation opportunity;
                # leave L[:, :, t, ref_base_idx_L] at its zero-initialized
                # value instead of simulating a self-to-self detection power.
                continue
            tc = int(trinuc_conv_np_L[t, b])
            P_arr = np.full(N_SIM, prob_amp_mat_L[tc, ref_base_idx_L])
            P_rev_arr = np.full(N_SIM, prob_amp_mat_rev_L[tc, ref_base_idx_L])
            Pt_arr = np.full(N_SIM, prob_dmg_t_L[tc, ref_base_idx_L])
            Prev_t_arr = np.full(N_SIM, prob_dmg_rev_t_L[tc, ref_base_idx_L])
            Pb_arr = np.full(N_SIM, prob_dmg_b_L[tc, ref_base_idx_L])
            Prev_b_arr = np.full(N_SIM, prob_dmg_rev_b_L[tc, ref_base_idx_L])

            for i in range(10):
                for j in range(10):
                    if i == 0 and j == 0:
                        continue

                    if i > 0:
                        q_top = rng_L.choice(all_quals, size=(i, N_SIM), replace=True)
                        F1R2_Pseq = -q_top / 10 * ln10_L
                        F1R2_bin = np.ones([i, N_SIM], dtype=bool)
                    else:
                        F1R2_Pseq = np.zeros([0, N_SIM])
                        F1R2_bin = np.zeros([0, N_SIM], dtype=bool)

                    if j > 0:
                        q_bot = rng_L.choice(all_quals, size=(j, N_SIM), replace=True)
                        F2R1_Pseq = -q_bot / 10 * ln10_L
                        F2R1_bin = np.ones([j, N_SIM], dtype=bool)
                    else:
                        F2R1_Pseq = np.zeros([0, N_SIM])
                        F2R1_bin = np.zeros([0, N_SIM], dtype=bool)

                    F1R2_b1, F1R2_b2 = calculateSSPosterior(
                        P_arr, P_rev_arr, F1R2_bin, F1R2_Pseq
                    )
                    F2R1_b1, F2R1_b2 = calculateSSPosterior(
                        P_arr, P_rev_arr, F2R1_bin, F2R1_Pseq
                    )
                    LL_B1, LL_B2 = calculateDSPosterior(
                        Pt_arr,
                        Prev_t_arr,
                        Pb_arr,
                        Prev_b_arr,
                        F1R2_b1,
                        F2R1_b1,
                        F1R2_b2,
                        F2R1_b2,
                    )
                    L[i, j, t, b] = ((LL_B1 - LL_B2) >= pcut).mean()

    # Build indel detection-power tables, analogous to L above but for indel
    # calling. genotypeDSIndel classifies each candidate indel by repeat
    # context (hps: homopolymer run length 0-20, strs: STR length bin 0-3)
    # and indel length (idLen), then selects Pamp/Pdmg via indelErrorProbs.
    # These tables simulate that same selection across every (top count,
    # bottom count) combination and precompute the fraction of simulations
    # that would clear indel calling's tiered pcutoffi threshold, so
    # per-position indel coverage can be looked up by depth+context instead
    # of re-simulated. Only built outside learn mode: ampmat_indel/
    # dmgmat_indel (and pcutoffi-driven calling) aren't used in learn mode.
    if not isLearn:
        prob_amp_indel_L = params["ampmat_indel"]
        prob_amp_indel_rev_L = params["ampmat_indel_rev"]
        prob_dmg_indel_L = params["dmgmat_indel"]
        prob_dmg_indel_rev_L = params["dmgmat_indel_rev"]
        pcutoffi_L = params["pcutoffi"]

        def indelTierThreshold(hps_c, strs_c):
            # Mirrors the OR'd tiering used to build LR_pass_bool for real
            # indel calls: hps<=5 -> tier0, 5<hps<=10 -> tier1, hps>10 ->
            # tier2, and strs>0 additionally ORs in tier2. OR'ing two
            # `LR >= x` comparisons on the same LR is equivalent to
            # comparing against the smaller threshold.
            if hps_c <= 5:
                threshold = pcutoffi_L[0]
            elif hps_c <= 10:
                threshold = pcutoffi_L[1]
            else:
                threshold = pcutoffi_L[2]
            if strs_c > 0:
                threshold = min(threshold, pcutoffi_L[2])
            return threshold

        def simulateIndelPowerGrid(
            Pamp_c,
            Pamp_rev_c,
            Pdmg_c,
            Pdmg_rev_c,
            Pdmg_bot_c,
            Pdmg_rev_bot_c,
            threshold,
        ):
            grid = np.zeros([10, 10])
            P_arr = np.full(N_SIM, Pamp_c)
            P_rev_arr = np.full(N_SIM, Pamp_rev_c)
            Pt_arr = np.full(N_SIM, Pdmg_c)
            Prev_t_arr = np.full(N_SIM, Pdmg_rev_c)
            Pb_arr = np.full(N_SIM, Pdmg_bot_c)
            Prev_b_arr = np.full(N_SIM, Pdmg_rev_bot_c)
            for i in range(10):
                for j in range(10):
                    if i == 0 and j == 0:
                        continue
                    if i > 0:
                        q_top = rng_L.choice(all_quals, size=(i, N_SIM), replace=True)
                        F1R2_Pseq = -q_top / 10 * ln10_L
                        F1R2_bin = np.ones([i, N_SIM], dtype=bool)
                    else:
                        F1R2_Pseq = np.zeros([0, N_SIM])
                        F1R2_bin = np.zeros([0, N_SIM], dtype=bool)
                    if j > 0:
                        q_bot = rng_L.choice(all_quals, size=(j, N_SIM), replace=True)
                        F2R1_Pseq = -q_bot / 10 * ln10_L
                        F2R1_bin = np.ones([j, N_SIM], dtype=bool)
                    else:
                        F2R1_Pseq = np.zeros([0, N_SIM])
                        F2R1_bin = np.zeros([0, N_SIM], dtype=bool)
                    F1R2_b1, F1R2_b2 = calculateSSPosterior(
                        P_arr, P_rev_arr, F1R2_bin, F1R2_Pseq
                    )
                    F2R1_b1, F2R1_b2 = calculateSSPosterior(
                        P_arr, P_rev_arr, F2R1_bin, F2R1_Pseq
                    )
                    LL_B1, LL_B2 = calculateDSPosterior(
                        Pt_arr,
                        Prev_t_arr,
                        Pb_arr,
                        Prev_b_arr,
                        F1R2_b1,
                        F2R1_b1,
                        F1R2_b2,
                        F2R1_b2,
                    )
                    grid[i, j] = ((LL_B1 - LL_B2) >= threshold).mean()
            return grid

        # 1bp indels (deletion/insertion of a homopolymer's repeat unit):
        # axes are top count, bottom count, hps (0-20), ref_allele (0-3),
        # sign (0=deletion idLen=-1, 1=insertion idLen=+1).
        L_indel_1bp = np.zeros([10, 10, 21, 4, 2])
        for hps_c in range(21):
            threshold = indelTierThreshold(hps_c, 0)
            for ref_allele_c in range(4):
                for sign_idx, idLen_c in enumerate((-1, 1)):
                    (
                        Pamp_c,
                        Pamp_rev_c,
                        Pdmg_c,
                        Pdmg_rev_c,
                        Pdmg_bot_c,
                        Pdmg_rev_bot_c,
                    ) = indelErrorProbs(
                        hps_c,
                        0,
                        idLen_c,
                        ref_allele_c,
                        prob_amp_indel_L,
                        prob_amp_indel_rev_L,
                        prob_dmg_indel_L,
                        prob_dmg_indel_rev_L,
                    )
                    L_indel_1bp[
                        :, :, hps_c, ref_allele_c, sign_idx
                    ] = simulateIndelPowerGrid(
                        Pamp_c,
                        Pamp_rev_c,
                        Pdmg_c,
                        Pdmg_rev_c,
                        Pdmg_bot_c,
                        Pdmg_rev_bot_c,
                        threshold,
                    )

        # Length-bin indels (raw length 2/3/4/5+, either direction): axes are
        # top count, bottom count, row (0-22, same hps-1/19+strs convention
        # as ampmat_indel), idLen+5 (0-10, only 8 of 11 slots populated).
        # Rows 0-19 are HP-context (hps 1-20); rows 20-22 are STR-context
        # (strs 1-3, hps forced to 1, matching genotypeDSIndel/indelErrorProbs).
        L_indel_len = np.zeros([10, 10, 23, 11])
        len_idLens = (-5, -4, -3, -2, 2, 3, 4, 5)
        for hps_c in range(1, 21):
            row = hps_c - 1
            threshold = indelTierThreshold(hps_c, 0)
            for idLen_c in len_idLens:
                (
                    Pamp_c,
                    Pamp_rev_c,
                    Pdmg_c,
                    Pdmg_rev_c,
                    Pdmg_bot_c,
                    Pdmg_rev_bot_c,
                ) = indelErrorProbs(
                    hps_c,
                    0,
                    idLen_c,
                    0,
                    prob_amp_indel_L,
                    prob_amp_indel_rev_L,
                    prob_dmg_indel_L,
                    prob_dmg_indel_rev_L,
                )
                L_indel_len[:, :, row, idLen_c + 5] = simulateIndelPowerGrid(
                    Pamp_c,
                    Pamp_rev_c,
                    Pdmg_c,
                    Pdmg_rev_c,
                    Pdmg_bot_c,
                    Pdmg_rev_bot_c,
                    threshold,
                )
        for strs_c in range(1, 4):
            row = 19 + strs_c
            threshold = indelTierThreshold(1, strs_c)
            for idLen_c in len_idLens:
                (
                    Pamp_c,
                    Pamp_rev_c,
                    Pdmg_c,
                    Pdmg_rev_c,
                    Pdmg_bot_c,
                    Pdmg_rev_bot_c,
                ) = indelErrorProbs(
                    1,
                    strs_c,
                    idLen_c,
                    0,
                    prob_amp_indel_L,
                    prob_amp_indel_rev_L,
                    prob_dmg_indel_L,
                    prob_dmg_indel_rev_L,
                )
                L_indel_len[:, :, row, idLen_c + 5] = simulateIndelPowerGrid(
                    Pamp_c,
                    Pamp_rev_c,
                    Pdmg_c,
                    Pdmg_rev_c,
                    Pdmg_bot_c,
                    Pdmg_rev_bot_c,
                    threshold,
                )

    retain_base = 5
    currentReadDictList = [
        dict() for _ in range(retain_base)
    ]  # Adjustable parameter pending
    for rec, region in bamIterateMultipleRegion(bam, regions, params.get("reference")):
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
            or not rec.is_proper_pair
            or rec.is_qcfail
        ):
            continue
        if rec.has_tag("DT"):
            if rec.get_tag("DT") == "SQ":
                continue
        if rec.has_tag("dt"):
            if rec.get_tag("dt") == "SQ":
                continue
        # If 5 prime is soft clipped
        if (rec.is_forward and rec.cigartuples[0][0] != 0) or (
            rec.is_reverse and rec.cigartuples[-1][0] != 0
        ):
            continue
        pass_read_num += 1
        start = rec.reference_start
        bc = rec.query_name.split("_")[-1]
        bcsplit = bc.split("+")
        bc1 = bcsplit[0]
        bc2 = bcsplit[1]
        if (rec.is_read1 and rec.is_forward) or (rec.is_read2 and rec.is_reverse):
            label = bc1 + "+" + bc2 + "+" + str(rec.template_length)
        else:
            label = bc2 + "+" + bc1 + "+" + str(rec.template_length)
        chrom = tumorBam.get_reference_name(rec.reference_id)
        if currentStart == -1:
            currentStart = start
        if start == currentStart:
            has_same_label_flag = False
            for rb in range(retain_base):
                if currentReadDictList[rb].get(label):
                    currentReadDictList[rb][label]["seqs"].append(rec)
                    currentReadDictList[rb][label]["names"][rec.query_name] = (
                        len(currentReadDictList[rb][label]["seqs"]) - 1
                    )
                    if (rec.is_forward and rec.is_read1) or (
                        rec.is_reverse and rec.is_read2
                    ):
                        currentReadDictList[rb][label]["F1R2"] += 1
                    else:
                        currentReadDictList[rb][label]["F2R1"] += 1
                    has_same_label_flag = True
                    break

            if not has_same_label_flag:
                # else:
                currentReadDictList[-1].update(
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
                    currentReadDictList[retain_base - 1][label]["F1R2"] += 1
                else:
                    currentReadDictList[retain_base - 1][label]["F2R1"] += 1
            # print(currentStart,start)
        else:
            # print(currentReadDictList)
            """
            Calling block starts
            """
            # print(currentReadDictList):
            currentReadDict = dict()
            for _ in range(min(start - currentStart, retain_base)):
                currentReadDict |= currentReadDictList.pop(0)
                currentReadDictList.append(dict())
            for key in currentReadDict.keys():
                flt_rs = "PASS"
                readSet = currentReadDict[key]["seqs"]
                all_dup = True
                for _ in readSet:
                    if not _.is_duplicate:
                        all_dup = False
                        break
                if all_dup:
                    continue
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]

                mean_mapq = sum([seq.mapping_quality for seq in readSet]) / len(readSet)
                if mean_mapq < params["mapq"]:
                    if params["rescue"]:
                        flt_rs = "low_mapq"
                    else:
                        continue
                meanASXS = sum(
                    [seq.get_tag("AS") - seq.get_tag("XS") for seq in readSet]
                ) / len(readSet)
                if meanASXS < params["minMeanASXS"]:
                    if params["rescue"]:
                        flt_rs = "low_ASXS"
                    else:
                        continue
                setBc = key.split(":")[0].split("+")
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]
                duplex_no = f"{F1R2}+{F2R1}"
                if duplex_read_num_dict.get(duplex_no) is None:
                    duplex_read_num_dict[duplex_no] = [0, 0, 0]
                    duplex_read_num_dict_trinuc[duplex_no] = np.zeros((64, 4))
                    duplex_read_num_dict_indel[duplex_no] = np.zeros(100)
                    duplex_read_num_dict_dbs[duplex_no] = np.zeros(144)
                    # unmasked_duplex_read_num_dict_trinuc[duplex_no] = np.zeros(96, dtype=int)
                duplex_read_num_dict[duplex_no][2] += 1
                unique_read_num += 1
                if F2R1 >= 1 and F1R2 >= 1:
                    if params["maxNM"]:
                        f1r2_blacklist_num = 0
                        f2r1_blacklist_num = 0
                        f1r2_count = 0
                        f2r1_count = 0
                        for _ in readSet:
                            if (_.is_forward and _.is_read1) or (
                                _.is_reverse and _.is_read2
                            ):
                                f1r2_count += 1
                                if _.query_name in read_blacklist:
                                    f1r2_blacklist_num += 1
                            else:
                                f2r1_count += 1
                                if _.query_name in read_blacklist:
                                    f2r1_blacklist_num += 1
                        if (
                            (
                                (f1r2_blacklist_num + f2r1_blacklist_num) / len(readSet)
                                >= 0.5
                            )
                            or f1r2_blacklist_num == f1r2_count
                            or f2r1_blacklist_num == f2r1_count
                        ):
                            if params["rescue"]:
                                flt_rs = "high_nm"
                            else:
                                continue
                    rs_reference_end = max([r.reference_end for r in readSet])
                    rs_reference_start = min([r.reference_start for r in readSet])

                    chromNow = readSet[0].reference_name
                    if (
                        chromNow != reference_mat_chrom
                        or rs_reference_end > reference_mat_end
                    ):
                        ### Output coverage
                        if "coverage" in locals():
                            if "coverage_leftover" in locals():
                                coverage[
                                    0 : coverage_leftover.shape[0]
                                ] += coverage_leftover
                                coverage_indel_cat[
                                    0 : coverage_leftover.shape[0]
                                ] += coverage_indel_cat_leftover
                                unmasked_coverage[
                                    0 : coverage_leftover.shape[0]
                                ] += unmasked_coverage_leftover
                                unmasked_coverage_indel_cat[
                                    0 : coverage_leftover.shape[0]
                                ] += unmasked_coverage_indel_cat_leftover
                                unmasked_coverage_leftover = np.zeros((1, 4))
                                unmasked_coverage_indel_cat_leftover = np.zeros((1, 14))
                                coverage_leftover = np.zeros((1, 4))
                                coverage_indel_cat_leftover = np.zeros((1, 14))
                            if chromNow == reference_mat_chrom:
                                coverage_leftover = copy.deepcopy(
                                    coverage[
                                        (rs_reference_start - reference_mat_start) : (
                                            reference_mat_end - reference_mat_start
                                        )
                                    ]
                                )
                                coverage_indel_cat_leftover = copy.deepcopy(
                                    coverage_indel_cat[
                                        (rs_reference_start - reference_mat_start) : (
                                            reference_mat_end - reference_mat_start
                                        )
                                    ]
                                )
                                unmasked_coverage_leftover = copy.deepcopy(
                                    unmasked_coverage[
                                        (rs_reference_start - reference_mat_start) : (
                                            reference_mat_end - reference_mat_start
                                        )
                                    ]
                                )
                                unmasked_coverage_indel_cat_leftover = copy.deepcopy(
                                    unmasked_coverage_indel_cat[
                                        (rs_reference_start - reference_mat_start) : (
                                            reference_mat_end - reference_mat_start
                                        )
                                    ]
                                )
                                non_zero_positions = np.nonzero(
                                    coverage[
                                        0 : (rs_reference_start - reference_mat_start)
                                    ].sum(axis=1)
                                    + coverage_indel_cat[
                                        0 : (rs_reference_start - reference_mat_start)
                                    ].sum(axis=1)
                                )
                            else:
                                non_zero_positions = np.nonzero(
                                    coverage.sum(axis=1)
                                    + coverage_indel_cat.sum(axis=1)
                                )
                            for pos in non_zero_positions[0].tolist():
                                current_pos = pos + reference_mat_start
                                bed_file = get_bed_file_for_position(
                                    current_pos,
                                    reference_mat_chrom,
                                    regions_start_chrom,
                                    regions_start_pos,
                                    regions_end_chrom,
                                    regions_end_pos,
                                    locus_bed,
                                    locus_bed_prev,
                                    locus_bed_next,
                                )
                                bed_file.write(
                                    "\t".join(
                                        [
                                            reference_mat_chrom,
                                            str(current_pos),
                                            str(current_pos + 1),
                                            "\t".join(str(v) for v in coverage[pos]),
                                            "\t".join(
                                                str(v) for v in coverage_indel_cat[pos]
                                            ),
                                        ]
                                    )
                                    + "\n"
                                )
                                total_coverage += coverage[pos]
                                total_coverage_indel_cat += coverage_indel_cat[pos]
                                total_unmasked_coverage += unmasked_coverage[pos]
                                total_unmasked_coverage_indel_cat += (
                                    unmasked_coverage_indel_cat[pos]
                                )
                        # if chromNow != reference_mat_chrom:
                        reference_mat_chrom = chromNow
                        # current_reference = str(fasta[reference_mat_chrom].seq)
                        reference_mat_start = rs_reference_start
                        try:
                            region_end = region[2]
                        except:
                            region_end = 10e10
                        contig_end = tumorBam.get_reference_length(chromNow)
                        reference_mat_end = min(
                            rs_reference_start + 1000000,
                            max(
                                region_end, max([seq.reference_end for seq in readSet])
                            ),
                            contig_end,
                        )
                        #
                        ref_h5 = h5py.File(params["reference"] + ".ref.h5", "r")
                        tn_h5 = h5py.File(params["reference"] + ".tn.h5", "r")
                        hp_h5 = h5py.File(params["reference"] + ".hp.h5", "r")
                        str_h5 = h5py.File(params["reference"] + ".str.h5", "r")
                        ref_np = ref_h5[reference_mat_chrom][
                            reference_mat_start:reference_mat_end
                        ]
                        trinuc_np = tn_h5[reference_mat_chrom][
                            reference_mat_start:reference_mat_end
                        ]
                        # (3, window) — unit_len, cut (start-of-run), and
                        # repeat_count, merged from the self-derived
                        # homopolymer index and the BED-derived STR index
                        # (see funcs/misc.py's load_repeat_context).
                        hp_np = load_repeat_context(
                            hp_h5[reference_mat_chrom],
                            str_h5[reference_mat_chrom],
                            reference_mat_start,
                            reference_mat_end,
                        )
                        # Raw, unmerged slices for the learn-phase error-rate
                        # estimation (profileTriNucMismatches), which needs
                        # the self-derived homopolymer run and the BED-derived
                        # STR annotation as independent sources rather than
                        # load_repeat_context's STR-priority merge (that merge
                        # is for indel *calling*'s classification, where a
                        # position can only be one category at a time).
                        hp_raw_np = hp_h5[reference_mat_chrom][
                            :, reference_mat_start:reference_mat_end
                        ]
                        str_raw_np = str_h5[reference_mat_chrom][
                            :, reference_mat_start:reference_mat_end
                        ]
                        (
                            prior_mat,
                            snp_mask,
                            indel_mask,
                            noise_mask,
                            n_cov_mask,
                            include_mask,
                            nm_mask
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
                            nbams,
                            bam,
                            params,
                        )
                        # print(ref_np,reference_mat_start)
                        coverage = np.zeros((1000000, 4))
                        coverage_indel_cat = np.zeros((1000000, 14))
                        unmasked_coverage = np.zeros((1000000, 4))
                        unmasked_coverage_indel_cat = np.zeros((1000000, 14))
                    ### Record read names to check if mate has been processed
                    processed_flag = 0
                    for seq in readSet:
                        if seq.query_name in processed_read_names:
                            processed_read_names.remove(seq.query_name)
                            processed_flag = 1
                            break
                    if processed_flag == 0 and flt_rs == "PASS":
                        processed_read_names.add(readSet[0].query_name)
                    start_ind = rs_reference_start - reference_mat_start
                    # reference_length_min = min(
                    # [read.reference_length for read in readSet]
                    # )
                    # reference_length_max = max(
                    # [read.reference_length for read in readSet]
                    # )
                    # end_ind = (
                    # rs_reference_start
                    # + reference_length_max
                    # - reference_mat_start
                    # )
                    end_ind = rs_reference_end - reference_mat_start

                    end_ind_max = end_ind
                    ref_lens = [_.reference_length for _ in readSet]
                    max_ref_len_abs = max([abs(_) for _ in ref_lens])
                    max_ref_num = ref_lens.index(max_ref_len_abs)
                    max_ref_len = ref_lens[max_ref_num]

                    masks = np.zeros([6, end_ind - start_ind], dtype=bool)
                    masks[0, :] = snp_mask[start_ind:end_ind]
                    masks[1, :] = noise_mask[start_ind:end_ind]
                    masks[2, :] = n_cov_mask[start_ind:end_ind]
                    masks[3, :] = include_mask[start_ind:end_ind]
                    masks[4, :] = nm_mask[start_ind:end_ind]
                    left, right = determineTrimLength(
                        readSet[max_ref_num],
                        params=params,
                        processed_flag=processed_flag,
                    )
                    masks[5, :left] = True
                    masks[5, -right:] = True
                    antimask = np.all(~masks, axis=0)
                    antimask[trinuc_np[start_ind:end_ind] > 64] = False
                    antimask[ref_np[start_ind:end_ind] == 4] = False
                    # Create unmasked version that only excludes trinuc > 64
                    unmasked_antimask = np.all(~masks[2:, :], axis=0)
                    unmasked_antimask[trinuc_np[start_ind:end_ind] > 64] = False
                    unmasked_antimask[ref_np[start_ind:end_ind] == 4] = False
                    learn_antimask = np.all(~masks[3:, :], axis=0)
                    learn_antimask[trinuc_np[start_ind:end_ind] > 64] = False
                    learn_antimask[ref_np[start_ind:end_ind] == 4] = False
                    ### If the whole reads are masked:
                    if not np.any(unmasked_antimask):
                        continue
                    indel_bool = [
                        ("I" in seq.cigarstring or "D" in seq.cigarstring)
                        for seq in readSet
                    ]
                    # if any(indel_bool):
                    if not isLearn:
                        masks_indel = np.zeros([6, end_ind_max - start_ind], dtype=bool)
                        masks_indel[0, :] = indel_mask[start_ind:end_ind_max]
                        masks_indel[1, :] = noise_mask[start_ind:end_ind_max]
                        masks_indel[2, :] = n_cov_mask[start_ind:end_ind_max]
                        left, right = determineTrimLength(
                            readSet[max_ref_num],
                            params=params,
                            processed_flag=processed_flag,
                        )
                        masks_indel[3, :left] = True
                        masks_indel[3, -right:] = True
                        masks_indel[4, :] = include_mask[start_ind:end_ind_max]
                        masks_indel[5, :] = nm_mask[start_ind:end_ind_max]
                        antimask_indel = np.all(~masks_indel, axis=0)
                        unmasked_antimask_indel = np.all(~masks_indel[2:, :], axis=0)
                        (
                            CS,
                            LR_raw,
                            LR_max,
                            indels,
                            hps,
                            strs,
                            F1R2_ref_count,
                            F1R2_alt_count,
                            F2R1_ref_count,
                            F2R1_alt_count,
                        ) = genotypeDSIndel(
                            readSet,
                            rs_reference_start,
                            rs_reference_end,
                            ref_np[start_ind:end_ind],
                            unmasked_antimask_indel,
                            hp_raw_np[:, start_ind:end_ind],
                            str_raw_np[:, start_ind:end_ind],
                            params,
                        )
                        # print(LR_max)
                        # pass_inds = np.nonzero(LR <= params["pcutoffi"])[0].tolist()

                        # pass_inds = np.nonzero(LR_raw >= params["pcutoffi"])[0].tolist()
                        # pass_inds = np.nonzero(LR_raw >= params["pcutoffi"])[0].tolist()
                        pass_inds = np.nonzero(CS >= params["cscutoffi"])[0].tolist()
                        # LR_pass_bool = (LR_raw >= params["pcutoffi"])
                        if len(LR_raw) > 0:
                            LR_pass_bool = np.zeros(len(LR_raw), dtype=bool)
                            LR_pass_bool[
                                np.logical_and(
                                    LR_raw >= params["pcutoffi"][0], hps <= 5
                                )
                            ] = True
                            LR_pass_bool[
                                np.all(
                                    np.vstack(
                                        [
                                            LR_raw >= params["pcutoffi"][1],
                                            hps <= 10,
                                            hps > 5,
                                        ]
                                    ),
                                    axis=0,
                                )
                            ] = True
                            LR_pass_bool[
                                np.logical_and(
                                    LR_raw >= params["pcutoffi"][2], hps > 10
                                )
                            ] = True
                            LR_pass_bool[
                                np.logical_and(
                                    LR_raw >= params["pcutoffi"][2], strs > 0
                                )
                            ] = True
                        indels_pass = [indels[_] for _ in pass_inds]
                        for nn in range(len(indels_pass)):
                            indel = indels_pass[nn]
                            indel_chrom = chromNow
                            indel_pos = int(indel.split(":")[0])
                            indel_size = int(indel.split(":")[1])
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
                            readPos = indel_pos - rs_reference_start
                            if readSet[0].template_length > 0:
                                readPos5p = min(
                                    readPos + 1,
                                    abs(readSet[0].template_length) - readPos,
                                )
                                readPos3p = min(
                                    abs(max_ref_len) - readPos,
                                    readPos + 1,
                                )
                            else:
                                readPos5p = min(
                                    max_ref_len - readPos,
                                    abs(readSet[0].template_length)
                                    - max_ref_len
                                    + readPos
                                    + 1,
                                )
                                readPos3p = min(
                                    abs(max_ref_len) - readPos,
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
                            if indel_size > 0:
                                offset = 0
                            else:
                                offset = -indel_size
                            unmasked_flag = antimask_indel[
                                indel_pos
                                - reference_mat_start
                                - start_ind : indel_pos
                                - reference_mat_start
                                - start_ind
                                + offset
                                + 1
                            ].all()
                            if unmasked_flag and LR_pass_bool[nn]:
                                flt = flt_rs
                            elif unmasked_flag:
                                flt = "underpowered"
                            else:
                                flt = "masked"

                            indel_rec = {
                                "chrom": chromNow,
                                "pos": indel_pos + 1,
                                "ref": indel_ref,
                                "alt": indel_alt,
                                "filter": flt,
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
                                    "CS": CS[pass_inds[nn]],
                                    "LR": LR_raw[pass_inds[nn]],
                                    "LM": LR_max[pass_inds[nn]],
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
                                    "TL": readSet[0].template_length,
                                    "STR": strs[nn],
                                },
                                "formats": ["AC", "RC", "DP"],
                                # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_indels.append(indel_rec)
                            indel_dict[indel_str] = 1
                    # else:
                    ### Calculate genotype probability
                    # if not any(indel_bool) or isLearn:
                    if 1:
                        if isLearn and F1R2 > 2 and F2R1 > 2:
                            (
                                mismatch_now,
                                indelerr_now,
                                mismatch_dmg_now,
                                indel_dmg_now,
                            ) = profileTriNucMismatches(
                                readSet,
                                rs_reference_start,
                                ref_np[start_ind:end_ind],
                                trinuc_np[start_ind:end_ind],
                                hp_raw_np[:, start_ind:end_ind],
                                str_raw_np[:, start_ind:end_ind],
                                np.copy(learn_antimask),
                                params,
                            )
                            mismatch_mat += mismatch_now
                            indelerr_mat += indelerr_now
                            mismatch_dmg_mat += mismatch_dmg_now
                            indel_dmg_mat += indel_dmg_now
                        if isLearn:
                            continue
                        (
                            cov_mat,
                            CS_mut,
                            LR_raw_mut,
                            LR_max_mut,
                            mut_mask,
                            b1_int,
                            unmasked_antimask,
                            F1R2_count,
                            F2R1_count,
                        ) = genotypeDSSnv(
                            readSet,
                            rs_reference_start,
                            ref_np[start_ind:end_ind],
                            trinuc_np[start_ind:end_ind],
                            prior_mat[start_ind:end_ind, :],
                            np.copy(unmasked_antimask),
                            params,
                            L,
                        )
                        # Per-position indel category coverage (del_unit,
                        # del_2, del_3, del_4, del_5plus, ins_unit, ins_2,
                        # ins_3, ins_4, ins_5plus): look up
                        # detection power from the L_indel_1bp/L_indel_len
                        # tables built above, using this window's repeat
                        # context (hp_np) and the same per-base read-family
                        # counts (F1R2_count/F2R1_count) SBS coverage uses
                        # as its depth proxy.
                        if not isLearn:
                            unit_len_arr = hp_np[0, start_ind:end_ind].astype(int)
                            repeat_count_arr = hp_np[2, start_ind:end_ind].astype(int)
                            ref_allele_arr = ref_np[start_ind:end_ind]
                            ref_allele_safe = np.where(
                                ref_allele_arr > 3, 0, ref_allele_arr
                            )
                            is_hp = unit_len_arr <= 1
                            hps_for_row = np.where(
                                is_hp, np.minimum(repeat_count_arr, 20), 1
                            )
                            total_len_arr = unit_len_arr * repeat_count_arr
                            strs_for_row = np.zeros_like(total_len_arr)
                            strs_for_row[
                                np.logical_and(~is_hp, total_len_arr >= 10)
                            ] = 1
                            strs_for_row[
                                np.logical_and(~is_hp, total_len_arr >= 25)
                            ] = 2
                            strs_for_row[
                                np.logical_and(~is_hp, total_len_arr >= 40)
                            ] = 3
                            row_len = np.where(
                                strs_for_row > 0, 19 + strs_for_row, hps_for_row - 1
                            )
                            unit_len_clamped = np.clip(unit_len_arr, 2, 5)
                            n_top_indel = np.minimum(F1R2_count.sum(axis=0), 9).astype(
                                int
                            )
                            n_bot_indel = np.minimum(F2R1_count.sum(axis=0), 9).astype(
                                int
                            )
                            # "Last cut"/"first cut" mask: repeat-based
                            # classes (ins_unit/del_unit here, and the
                            # homopolymer/STR raw classes below) can't
                            # reliably confirm a run's true length once this
                            # read family no longer extends past that run's
                            # start (hp_np row-1 "cut", from
                            # load_repeat_context), so everything from the
                            # last cut within this family's span onward is
                            # excluded — and symmetrically, the entire run
                            # touching the family's own LEFT edge is
                            # excluded too: whether an indel occurred within
                            # a run right at the start of what this family's
                            # reads cover is exactly as undecidable as one
                            # at the end, since the family has no visible
                            # flanking context on that side either way. This
                            # is the whole run containing window position 0
                            # — not just its first base — so if the window
                            # itself starts exactly on a cut (that run's
                            # start is confirmed), the excluded span still
                            # runs to the *next* cut, not just past the
                            # first one. Does not apply to del_2..del_5plus
                            # / microhomology, which use a fixed context
                            # regardless of any actual run.
                            cut_arr = hp_np[1, start_ind:end_ind]
                            last_cut_valid = np.ones(end_ind - start_ind, dtype=bool)
                            cut_positions = np.nonzero(cut_arr)[0]
                            if cut_positions.size > 0:
                                last_cut_valid[cut_positions[-1] :] = False
                                if cut_positions[0] == 0:
                                    first_run_end = (
                                        cut_positions[1]
                                        if cut_positions.size > 1
                                        else end_ind - start_ind
                                    )
                                else:
                                    first_run_end = cut_positions[0]
                                last_cut_valid[:first_run_end] = False
                            # Reference base immediately following each
                            # position, needed for the "Insertion A/T/C/G"
                            # columns below and the matching fine-grained
                            # 1bpins{base} columns (96-99): rep0 for
                            # inserting base N requires only that this next
                            # base differs from N — nothing to do with the
                            # position's own reference base.
                            window_len = end_ind - start_ind
                            next_ref_arr = np.full(window_len, -1, dtype=int)
                            avail = min(end_ind + 1, ref_np.shape[0]) - (start_ind + 1)
                            if avail > 0:
                                next_ref_arr[:avail] = ref_np[
                                    start_ind + 1 : start_ind + 1 + avail
                                ]
                            next_ref_valid = (next_ref_arr >= 0) & (next_ref_arr <= 3)

                            cov_mat_indel = np.zeros([window_len, 14])
                            cov_mat_indel[:, 0] = (
                                np.where(
                                    is_hp,
                                    L_indel_1bp[
                                        n_top_indel,
                                        n_bot_indel,
                                        hps_for_row,
                                        ref_allele_safe,
                                        0,
                                    ],
                                    L_indel_len[
                                        n_top_indel,
                                        n_bot_indel,
                                        row_len,
                                        -unit_len_clamped + 5,
                                    ],
                                )
                                * last_cut_valid
                            )
                            # Insertion of Repeat Unit: same shape as column
                            # 0 (deletion of repeat unit) but the insertion
                            # direction (last L_indel_1bp index 1, and
                            # +unit_len_clamped+5 instead of
                            # -unit_len_clamped+5 into L_indel_len).
                            cov_mat_indel[:, 1] = (
                                np.where(
                                    is_hp,
                                    L_indel_1bp[
                                        n_top_indel,
                                        n_bot_indel,
                                        hps_for_row,
                                        ref_allele_safe,
                                        1,
                                    ],
                                    L_indel_len[
                                        n_top_indel,
                                        n_bot_indel,
                                        row_len,
                                        unit_len_clamped + 5,
                                    ],
                                )
                                * last_cut_valid
                            )
                            # del_2..del_5plus: fixed homopolymer-length-1
                            # context (microhomology-style, row 0 = hps=1),
                            # not the position's actual repeat context.
                            cov_mat_indel[:, 2] = L_indel_len[
                                n_top_indel, n_bot_indel, 0, 3
                            ]
                            cov_mat_indel[:, 3] = L_indel_len[
                                n_top_indel, n_bot_indel, 0, 2
                            ]
                            cov_mat_indel[:, 4] = L_indel_len[
                                n_top_indel, n_bot_indel, 0, 1
                            ]
                            cov_mat_indel[:, 5] = L_indel_len[
                                n_top_indel, n_bot_indel, 0, 0
                            ]
                            # Insertion A/T/C/G: rep0 opportunity to insert
                            # exactly that base — fixed hps=1 context, using
                            # THAT base (not the position's own reference
                            # base) as the L_indel_1bp ref_allele index —
                            # zeroed wherever the next reference base already
                            # equals it (that's a repeat extension, not a
                            # novel insertion).
                            for b in range(4):
                                cov_mat_indel[:, 6 + b] = (
                                    L_indel_1bp[n_top_indel, n_bot_indel, 1, b, 1]
                                    * next_ref_valid
                                    * (next_ref_arr != b)
                                    * antimask_indel
                                )
                            # Insertion Length 2..5+: power to see a novel
                            # U-bp unit inserted where no repeat exists yet
                            # (rep0 in classify_indel_record's is_hp/indel_len
                            # >=2 branch, misc.py) — fixed hps=1 context (row
                            # 0), same as del_2..del_5plus, but only where a
                            # "no repeat here yet" call is even possible
                            # (is_hp) and the read family can still confirm
                            # that (last_cut_valid). Without this, novel
                            # multi-bp insertions that DO get classified and
                            # counted in the fine-grained per-duplex-group
                            # totals had no per-locus coverage anywhere.
                            ins_len_valid = is_hp & last_cut_valid & antimask_indel
                            cov_mat_indel[:, 10] = (
                                L_indel_len[n_top_indel, n_bot_indel, 0, 7]
                                * ins_len_valid
                            )
                            cov_mat_indel[:, 11] = (
                                L_indel_len[n_top_indel, n_bot_indel, 0, 8]
                                * ins_len_valid
                            )
                            cov_mat_indel[:, 12] = (
                                L_indel_len[n_top_indel, n_bot_indel, 0, 9]
                                * ins_len_valid
                            )
                            cov_mat_indel[:, 13] = (
                                L_indel_len[n_top_indel, n_bot_indel, 0, 10]
                                * ins_len_valid
                            )
                            # Fine-grained 100-class raw indel classification
                            # (mirrors the SBS 192-class raw trinuc scheme):
                            # base-specific homopolymer del (24) + ins (20,
                            # no rep0 bin — see below), unit-size/repeat-
                            # count-specific STR del/ins (24+24),
                            # microhomology deletion length (4). See
                            # funcs/misc.py's build_indel100_labels for
                            # the exact column layout.
                            indel100 = np.zeros(100)

                            def _accumulate_indel100(col_idx, weight, valid):
                                if np.any(valid):
                                    indel100[:] += np.bincount(
                                        col_idx[valid].astype(int),
                                        weights=weight[valid],
                                        minlength=100,
                                    )

                            to_cgta = np.array([3, 2, 0, 1])  # base2num->CGTA order
                            base_cgta = to_cgta[ref_allele_safe]

                            # HP and STR opportunity below are credited
                            # INDEPENDENTLY, from the raw (unmerged)
                            # hp.h5/str.h5 arrays (hp_raw_np/str_raw_np) —
                            # not from is_hp/unit_len_arr/cut_arr above,
                            # which come from load_repeat_context's
                            # STR-priority merge and silently drop a
                            # position's homopolymer opportunity wherever
                            # it also sits inside an annotated STR tract
                            # (e.g. the "AA" in a (AAT)n repeat). That merge
                            # is correct for classifying an *observed*
                            # mutation event (necessarily one category or
                            # the other) but wrong here: such a position is
                            # a real candidate for both a 1bp homopolymer
                            # slip and a larger STR-unit slip. Mirrors
                            # funcs/misc.py's
                            # indel100_reference_bucket_indices, which
                            # already double-counts these positions into
                            # both buckets on the genome-composition side —
                            # see [[feedback-hp-str-independence]].
                            def _run_boundary_valid(cut_bool_arr):
                                valid = np.ones(cut_bool_arr.shape[0], dtype=bool)
                                cut_positions = np.nonzero(cut_bool_arr)[0]
                                if cut_positions.size > 0:
                                    valid[cut_positions[-1] :] = False
                                    if cut_positions[0] == 0:
                                        first_run_end = (
                                            cut_positions[1]
                                            if cut_positions.size > 1
                                            else cut_bool_arr.shape[0]
                                        )
                                    else:
                                        first_run_end = cut_positions[0]
                                    valid[:first_run_end] = False
                                return valid

                            hp_run_arr = np.minimum(
                                hp_raw_np[0, start_ind:end_ind].astype(int), 20
                            )
                            hp_cut_bool = hp_raw_np[1, start_ind:end_ind].astype(bool)
                            hp_repeat_valid = (
                                _run_boundary_valid(hp_cut_bool) & antimask_indel
                            )

                            str_unit_raw = str_raw_np[0, start_ind:end_ind].astype(int)
                            str_repeat_raw = str_raw_np[1, start_ind:end_ind].astype(
                                int
                            )
                            str_cut_bool = str_raw_np[2, start_ind:end_ind].astype(bool)
                            is_real_str = str_unit_raw >= 2
                            str_repeat_valid = (
                                _run_boundary_valid(str_cut_bool) & antimask_indel
                            )
                            total_len_str = str_unit_raw * str_repeat_raw
                            strs_for_row_str = np.zeros_like(total_len_str)
                            strs_for_row_str[is_real_str & (total_len_str >= 10)] = 1
                            strs_for_row_str[is_real_str & (total_len_str >= 25)] = 2
                            strs_for_row_str[is_real_str & (total_len_str >= 40)] = 3
                            row_len_str = np.where(
                                strs_for_row_str > 0, 19 + strs_for_row_str, 0
                            )
                            unit_len_clamped_str = np.clip(str_unit_raw, 2, 5)
                            unit_bucket = np.clip(str_unit_raw, 2, 5) - 2
                            count_bucket_del = np.clip(str_repeat_raw, 1, 6) - 1

                            # 1-4: homopolymer deletion (cols 0-23) — always
                            # from the position's own raw hp run/cut,
                            # regardless of any STR annotation there too.
                            hp_del_power = L_indel_1bp[
                                n_top_indel,
                                n_bot_indel,
                                hp_run_arr,
                                ref_allele_safe,
                                0,
                            ]
                            del_bucket = np.clip(hp_run_arr, 1, 6) - 1
                            _accumulate_indel100(
                                base_cgta * 6 + del_bucket,
                                hp_del_power,
                                hp_repeat_valid & hp_cut_bool,
                            )

                            # 5-8: homopolymer insertion (cols 24-43) —
                            # length bins 1..5+ only; rep0 (hp_run_arr==0,
                            # never actually occurs) is deliberately not
                            # credited here, since that value is always
                            # discarded and replaced by the exact
                            # 1bpins{base} next-base computation below (see
                            # Estimate.py's
                            # override_inshp0_with_next_base_opportunity).
                            hp_ins_power = L_indel_1bp[
                                n_top_indel,
                                n_bot_indel,
                                hp_run_arr,
                                ref_allele_safe,
                                1,
                            ]
                            ins_bucket = np.clip(hp_run_arr, 1, 5) - 1
                            _accumulate_indel100(
                                24 + base_cgta * 5 + ins_bucket,
                                hp_ins_power,
                                (hp_run_arr >= 1) & hp_repeat_valid & hp_cut_bool,
                            )

                            # 9: STR deletion (cols 44-67) — real part, from
                            # the position's own raw STR annotation,
                            # independent of any homopolymer overlap.
                            str_del_power = L_indel_len[
                                n_top_indel,
                                n_bot_indel,
                                row_len_str,
                                -unit_len_clamped_str + 5,
                            ]
                            _accumulate_indel100(
                                44 + unit_bucket * 6 + count_bucket_del,
                                str_del_power,
                                is_real_str & str_repeat_valid & str_cut_bool,
                            )
                            # Flat rep1 opportunity at positions with no
                            # real STR annotation (whether or not they're
                            # also a homopolymer run — that's credited
                            # separately above): every such position is a
                            # candidate for an arbitrary, non-repeating
                            # U-bp deletion — a fixed hps=1/strs=0 power
                            # lookup (row 0), not an attempt to discover a
                            # true higher count that might coincidentally
                            # exist there (that would need rescanning the
                            # reference; see funcs/misc.py's
                            # indel100_reference_bucket_indices for the
                            # matching scan-free credit on the genome-wide
                            # side). Gated by hp_repeat_valid, not
                            # str_repeat_valid: there's no real STR tract
                            # here to have its own boundary, so this reuses
                            # the run-visibility this position already has
                            # under its own (possibly trivial, length-1)
                            # homopolymer context.
                            for u_idx, U in enumerate([2, 3, 4, 5]):
                                hyp_del_power = L_indel_len[
                                    n_top_indel, n_bot_indel, 0, -U + 5
                                ]
                                _accumulate_indel100(
                                    np.full(end_ind - start_ind, 44 + u_idx * 6),
                                    hyp_del_power,
                                    (~is_real_str) & hp_repeat_valid,
                                )

                            # 10: STR insertion (cols 68-91) — real part
                            str_ins_power = L_indel_len[
                                n_top_indel,
                                n_bot_indel,
                                row_len_str,
                                unit_len_clamped_str + 5,
                            ]
                            count_bucket_ins = np.clip(str_repeat_raw, 0, 5)
                            _accumulate_indel100(
                                68 + unit_bucket * 6 + count_bucket_ins,
                                str_ins_power,
                                is_real_str & str_repeat_valid & str_cut_bool,
                            )
                            is_rep1 = str_repeat_raw == 1
                            _accumulate_indel100(
                                68 + unit_bucket * 6 + 0,
                                str_ins_power,
                                is_real_str & is_rep1 & str_repeat_valid & str_cut_bool,
                            )
                            # Flat rep0/rep1 opportunity at positions with
                            # no real STR annotation ("what if a U-bp unit
                            # were inserted here"): fixed hps=1/strs=0
                            # context (row 0), same value added to both the
                            # "0" and "1" repeat-count buckets for every
                            # unit size — scan-free, matching the deletion
                            # case above, and gated the same way (see
                            # comment there for why hp_repeat_valid, not
                            # str_repeat_valid).
                            for u_idx, U in enumerate([2, 3, 4, 5]):
                                hyp_power = L_indel_len[
                                    n_top_indel, n_bot_indel, 0, U + 5
                                ]
                                base_col = 68 + u_idx * 6
                                _accumulate_indel100(
                                    np.full(end_ind - start_ind, base_col),
                                    hyp_power,
                                    (~is_real_str) & hp_repeat_valid,
                                )
                                _accumulate_indel100(
                                    np.full(end_ind - start_ind, base_col + 1),
                                    hyp_power,
                                    (~is_real_str) & hp_repeat_valid,
                                )

                            # 11: microhomology deletion (cols 92-95) —
                            # fixed context, every position, no last-cut
                            # restriction; reuse cov_mat_indel[:,2:6].
                            for k in range(4):
                                _accumulate_indel100(
                                    np.full(end_ind - start_ind, 92 + k),
                                    cov_mat_indel[:, 2 + k],
                                    antimask_indel,
                                )

                            # 12: 1bpins{base} rep0 opportunity (cols
                            # 96-99) — identical per-position values to
                            # cov_mat_indel[:,6:10] (Insertion A/T/C/G), so
                            # just reuse them rather than recomputing; the
                            # antimask_indel/next-base gating is already
                            # baked into those columns.
                            for b in range(4):
                                _accumulate_indel100(
                                    np.full(window_len, 96 + b),
                                    cov_mat_indel[:, 6 + b],
                                    np.ones(window_len, dtype=bool),
                                )
                        """
                        (
                            LR_old,
                            b1_int_old,
                            unmasked_antimask_old,
                            F1R2_count_old,
                            F2R1_count_old,
                            prob1_old,prob2_old,prob3_old,prob4_old
                        ) = oldDSSnv(
                            readSet,
                            ref_np[start_ind:end_ind],
                            trinuc_np[start_ind:end_ind],
                            prior_mat[start_ind:end_ind, :],
                            np.copy(unmasked_antimask),
                            params,
                        )
                        """
                        # prob1_diff = prob1-prob1_old
                        # notice1 = np.abs(prob1_diff) >=1
                        # prob3_diff = prob3-prob3_old
                        # notice2 = np.abs(prob3_diff) >=1
                        # if np.any(notice1):
                        # print("pr1",prob1[notice1],prob1_old[notice1],F1R2_count[:,unmasked_antimask_old][:,notice1])

                        ref_int = ref_np[start_ind:end_ind]
                        n_win = ref_int.size
                        mut_pos_in_win = np.nonzero(mut_mask)[0]
                        if CS_mut.size > 0:
                            muts_ind_compressed = np.nonzero(
                                CS_mut >= params["cscutoff"]
                            )[0]
                            muts_ind = mut_pos_in_win[muts_ind_compressed].tolist()
                        else:
                            muts_ind_compressed = np.zeros(0, dtype=int)
                            muts_ind = []
                        refs_ind = np.nonzero(
                            np.logical_and(unmasked_antimask, b1_int == ref_int)
                        )[0].tolist()
                        LR_pass_bool = np.zeros(n_win, dtype=bool)
                        LR_pass_bool[mut_mask] = LR_raw_mut >= params["pcutoff"]
                        alt_int = b1_int
                        unmasked_pass_bool = np.full(n_win, False, dtype=bool)
                        unmasked_pass_bool[refs_ind] = True
                        unmasked_pass_bool[muts_ind] = True
                        pass_bool = np.copy(unmasked_pass_bool)
                        pass_bool[~antimask] = False
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
                        for nn in range(len(mut_positions)):
                            mut_chrom = reference_mat_chrom
                            mut_pos = mut_positions[nn]
                            mut_ref = num2base[ref_int[muts_ind[nn]]]
                            mut_alt = num2base[alt_int[muts_ind[nn]]]
                            mut_trinuc = num2trinuc[
                                trinuc_np[start_ind:end_ind][muts_ind[nn]]
                            ]
                            if readSet[0].template_length > 0:
                                readPos5p = min(
                                    muts_ind[nn] + 1,
                                    abs(readSet[0].template_length) - muts_ind[nn],
                                )
                                readPos3p = min(
                                    abs(max_ref_len) - muts_ind[nn],
                                    muts_ind[nn] + 1,
                                )
                            else:
                                readPos5p = min(
                                    max_ref_len - muts_ind[nn],
                                    abs(readSet[0].template_length)
                                    - max_ref_len
                                    + muts_ind[nn]
                                    + 1,
                                )
                                readPos3p = min(
                                    abs(max_ref_len) - muts_ind[nn],
                                    muts_ind[nn] + 1,
                                )
                            FPs.append(readPos5p)
                            RPs.append(readPos3p)
                            if F1R2_count[:, muts_ind[nn]].sum() == 0:
                                continue
                            if F2R1_count[:, muts_ind[nn]].sum() == 0:
                                continue
                            if pass_bool[muts_ind[nn]] and LR_pass_bool[muts_ind[nn]]:
                                flt = flt_rs
                            elif pass_bool[muts_ind[nn]]:
                                flt = "underpowered"
                            else:
                                flt = "masked"
                            mut = {
                                "chrom": mut_chrom,
                                "pos": mut_pos,
                                "ref": mut_ref,
                                "alt": mut_alt,
                                "filter": flt,
                                "infos": {
                                    "F1R2": F1R2,
                                    "F2R1": F2R1,
                                    "CS": CS_mut[muts_ind_compressed[nn]],
                                    "LR": LR_raw_mut[muts_ind_compressed[nn]],
                                    "LM": LR_max_mut[muts_ind_compressed[nn]],
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
                                    "TL": readSet[0].template_length,
                                    "STR": 0,
                                },
                                "formats": ["AC", "RC", "DP"],
                                # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                            }
                            muts_dict[
                                "_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])
                            ] = 0
                            muts.append(mut)
                        muts_dbs.extend(
                            _detect_dbs_pairs(
                                reference_mat_chrom,
                                mut_positions,
                                muts_ind,
                                ref_int,
                                alt_int,
                                pass_bool,
                                LR_pass_bool,
                                flt_rs,
                                F1R2,
                                F2R1,
                                setBc,
                                currentStart,
                                readSet[0].template_length,
                                num2base,
                            )
                        )
                        """
                        if isLearn:
                            continue
                        """
                        if flt_rs == "PASS":
                            coverage[start_ind:end_ind][pass_bool] += cov_mat[pass_bool]
                            duplex_read_num_dict[duplex_no][1] += np.count_nonzero(
                                pass_bool
                            )
                            trinuc_pass = trinuc_np[start_ind:end_ind][pass_bool]
                            cov_pass = cov_mat[pass_bool]
                            for b in range(4):
                                duplex_read_num_dict_trinuc[duplex_no][
                                    :, b
                                ] += np.bincount(
                                    trinuc_pass, weights=cov_pass[:, b], minlength=64
                                )
                            duplex_read_num_dict_dbs[
                                duplex_no
                            ] += _compute_dbs_opportunity(cov_mat, ref_int, pass_bool)

                            # Update unmasked coverage and trinuc counts (includes all passing sites)

                            unmasked_coverage[start_ind:end_ind][
                                unmasked_pass_bool
                            ] += cov_mat[unmasked_pass_bool]
                            if pass_bool.any():
                                duplex_read_num_dict[duplex_no][0] += 1
                                duplex_count += 1
                            if not isLearn:
                                coverage_indel_cat[start_ind:end_ind][
                                    antimask_indel
                                ] += cov_mat_indel[antimask_indel]
                                unmasked_coverage_indel_cat[start_ind:end_ind][
                                    unmasked_antimask_indel
                                ] += cov_mat_indel[unmasked_antimask_indel]
                                duplex_read_num_dict_indel[duplex_no] += indel100
            """
            Calling block ends
            """
            has_same_label_flag = False
            for rb in range(retain_base):
                if currentReadDictList[rb].get(label):
                    """
                    if currentReadDictList[rb][label]["names"].get(rec.query_name):
                        if rec.is_read2:
                            continue
                        else:
                            ind = currentReadDictList[rb][label]["names"].get(rec.query_name)
                            currentReadDictList[rb][label]["seqs"][ind] = rec
                    else:
                    """
                    currentReadDictList[rb][label]["seqs"].append(rec)
                    currentReadDictList[rb][label]["names"][rec.query_name] = (
                        len(currentReadDictList[rb][label]["seqs"]) - 1
                    )
                    if (rec.is_forward and rec.is_read1) or (
                        rec.is_reverse and rec.is_read2
                    ):
                        currentReadDictList[rb][label]["F1R2"] += 1
                    else:
                        currentReadDictList[rb][label]["F2R1"] += 1
                    has_same_label_flag = True
                    break

            if not has_same_label_flag:
                # else:
                currentReadDictList[retain_base - 1].update(
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
                    currentReadDictList[retain_base - 1][label]["F1R2"] += 1
                else:
                    currentReadDictList[retain_base - 1][label]["F2R1"] += 1
            currentStart = start
    """
    Calling block starts
    """
    # print(currentReadDictList):
    currentReadDict = dict()
    for _ in range(len(currentReadDictList)):
        currentReadDict |= currentReadDictList.pop(0)
        currentReadDictList.append(dict())

    for key in currentReadDict.keys():
        flt_rs = "PASS"
        readSet = currentReadDict[key]["seqs"]
        all_dup = True
        for _ in readSet:
            if not _.is_duplicate:
                all_dup = False
                break
        if all_dup:
            continue
        F2R1 = currentReadDict[key]["F2R1"]
        F1R2 = currentReadDict[key]["F1R2"]

        mean_mapq = sum([seq.mapping_quality for seq in readSet]) / len(readSet)
        if mean_mapq < params["mapq"]:
            if params["rescue"]:
                flt_rs = "low_mapq"
            else:
                continue
        meanASXS = sum(
            [seq.get_tag("AS") - seq.get_tag("XS") for seq in readSet]
        ) / len(readSet)
        if meanASXS < params["minMeanASXS"]:
            if params["rescue"]:
                flt_rs = "low_ASXS"
            else:
                continue
        setBc = key.split(":")[0].split("+")
        F2R1 = currentReadDict[key]["F2R1"]
        F1R2 = currentReadDict[key]["F1R2"]
        duplex_no = f"{F1R2}+{F2R1}"
        if duplex_read_num_dict.get(duplex_no) is None:
            duplex_read_num_dict[duplex_no] = [0, 0, 0]
            duplex_read_num_dict_trinuc[duplex_no] = np.zeros((64, 4))
            duplex_read_num_dict_indel[duplex_no] = np.zeros(100)
            duplex_read_num_dict_dbs[duplex_no] = np.zeros(144)
            # unmasked_duplex_read_num_dict_trinuc[duplex_no] = np.zeros(96, dtype=int)
        duplex_read_num_dict[duplex_no][2] += 1
        unique_read_num += 1
        if F2R1 >= 1 and F1R2 >= 1:
            if params["maxNM"]:
                f1r2_blacklist_num = 0
                f2r1_blacklist_num = 0
                f1r2_count = 0
                f2r1_count = 0
                for _ in readSet:
                    if (_.is_forward and _.is_read1) or (_.is_reverse and _.is_read2):
                        f1r2_count += 1
                        if _.query_name in read_blacklist:
                            f1r2_blacklist_num += 1
                    else:
                        f2r1_count += 1
                        if _.query_name in read_blacklist:
                            f2r1_blacklist_num += 1
                if (
                    ((f1r2_blacklist_num + f2r1_blacklist_num) / len(readSet) >= 0.5)
                    or f1r2_blacklist_num == f1r2_count
                    or f2r1_blacklist_num == f2r1_count
                ):
                    if params["rescue"]:
                        flt_rs = "high_nm"
                    else:
                        continue
            rs_reference_end = max([r.reference_end for r in readSet])
            rs_reference_start = min([r.reference_start for r in readSet])

            chromNow = readSet[0].reference_name
            if chromNow != reference_mat_chrom or rs_reference_end > reference_mat_end:
                ### Output coverage
                if "coverage" in locals():
                    if "coverage_leftover" in locals():
                        coverage[0 : coverage_leftover.shape[0]] += coverage_leftover
                        coverage_indel_cat[
                            0 : coverage_leftover.shape[0]
                        ] += coverage_indel_cat_leftover
                        unmasked_coverage[
                            0 : coverage_leftover.shape[0]
                        ] += unmasked_coverage_leftover
                        unmasked_coverage_indel_cat[
                            0 : coverage_leftover.shape[0]
                        ] += unmasked_coverage_indel_cat_leftover
                        unmasked_coverage_leftover = np.zeros((1, 4))
                        unmasked_coverage_indel_cat_leftover = np.zeros((1, 14))
                        coverage_leftover = np.zeros((1, 4))
                        coverage_indel_cat_leftover = np.zeros((1, 14))
                    if chromNow == reference_mat_chrom:
                        coverage_leftover = copy.deepcopy(
                            coverage[
                                (rs_reference_start - reference_mat_start) : (
                                    reference_mat_end - reference_mat_start
                                )
                            ]
                        )
                        coverage_indel_cat_leftover = copy.deepcopy(
                            coverage_indel_cat[
                                (rs_reference_start - reference_mat_start) : (
                                    reference_mat_end - reference_mat_start
                                )
                            ]
                        )
                        unmasked_coverage_leftover = copy.deepcopy(
                            unmasked_coverage[
                                (rs_reference_start - reference_mat_start) : (
                                    reference_mat_end - reference_mat_start
                                )
                            ]
                        )
                        unmasked_coverage_indel_cat_leftover = copy.deepcopy(
                            unmasked_coverage_indel_cat[
                                (rs_reference_start - reference_mat_start) : (
                                    reference_mat_end - reference_mat_start
                                )
                            ]
                        )
                        non_zero_positions = np.nonzero(
                            coverage[
                                0 : (rs_reference_start - reference_mat_start)
                            ].sum(axis=1)
                            + coverage_indel_cat[
                                0 : (rs_reference_start - reference_mat_start)
                            ].sum(axis=1)
                        )
                    else:
                        non_zero_positions = np.nonzero(
                            coverage.sum(axis=1) + coverage_indel_cat.sum(axis=1)
                        )
                    for pos in non_zero_positions[0].tolist():
                        current_pos = pos + reference_mat_start
                        bed_file = get_bed_file_for_position(
                            current_pos,
                            reference_mat_chrom,
                            regions_start_chrom,
                            regions_start_pos,
                            regions_end_chrom,
                            regions_end_pos,
                            locus_bed,
                            locus_bed_prev,
                            locus_bed_next,
                        )
                        bed_file.write(
                            "\t".join(
                                [
                                    reference_mat_chrom,
                                    str(current_pos),
                                    str(current_pos + 1),
                                    "\t".join(str(v) for v in coverage[pos]),
                                    "\t".join(str(v) for v in coverage_indel_cat[pos]),
                                ]
                            )
                            + "\n"
                        )
                        total_coverage += coverage[pos]
                        total_coverage_indel_cat += coverage_indel_cat[pos]
                        total_unmasked_coverage += unmasked_coverage[pos]
                        total_unmasked_coverage_indel_cat += (
                            unmasked_coverage_indel_cat[pos]
                        )
                # if chromNow != reference_mat_chrom:
                reference_mat_chrom = chromNow
                # current_reference = str(fasta[reference_mat_chrom].seq)
                reference_mat_start = rs_reference_start
                try:
                    region_end = region[2]
                except:
                    region_end = 10e10
                contig_end = tumorBam.get_reference_length(chromNow)
                reference_mat_end = min(
                    rs_reference_start + 1000000,
                    max(region_end, max([seq.reference_end for seq in readSet])),
                    contig_end,
                )
                #
                ref_h5 = h5py.File(params["reference"] + ".ref.h5", "r")
                tn_h5 = h5py.File(params["reference"] + ".tn.h5", "r")
                hp_h5 = h5py.File(params["reference"] + ".hp.h5", "r")
                str_h5 = h5py.File(params["reference"] + ".str.h5", "r")
                ref_np = ref_h5[reference_mat_chrom][
                    reference_mat_start:reference_mat_end
                ]
                trinuc_np = tn_h5[reference_mat_chrom][
                    reference_mat_start:reference_mat_end
                ]
                # (3, window) — unit_len, cut (start-of-run), and
                # repeat_count, merged from the self-derived homopolymer
                # index and the BED-derived STR index (see
                # funcs/misc.py's load_repeat_context).
                hp_np = load_repeat_context(
                    hp_h5[reference_mat_chrom],
                    str_h5[reference_mat_chrom],
                    reference_mat_start,
                    reference_mat_end,
                )
                # Raw, unmerged slices for the learn-phase error-rate
                # estimation (profileTriNucMismatches), which needs the
                # self-derived homopolymer run and the BED-derived STR
                # annotation as independent sources rather than
                # load_repeat_context's STR-priority merge (that merge is
                # for indel *calling*'s classification, where a position
                # can only be one category at a time).
                hp_raw_np = hp_h5[reference_mat_chrom][
                    :, reference_mat_start:reference_mat_end
                ]
                str_raw_np = str_h5[reference_mat_chrom][
                    :, reference_mat_start:reference_mat_end
                ]
                (
                    prior_mat,
                    snp_mask,
                    indel_mask,
                    noise_mask,
                    n_cov_mask,
                    include_mask,
                    nm_mask
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
                    nbams,
                    bam,
                    params,
                )
                # print(ref_np,reference_mat_start)
                coverage = np.zeros((1000000, 4))
                coverage_indel_cat = np.zeros((1000000, 14))
                unmasked_coverage = np.zeros((1000000, 4))
                unmasked_coverage_indel_cat = np.zeros((1000000, 14))
            ### Record read names to check if mate has been processed
            processed_flag = 0
            for seq in readSet:
                if seq.query_name in processed_read_names:
                    processed_read_names.remove(seq.query_name)
                    processed_flag = 1
                    break
            if processed_flag == 0 and flt_rs == "PASS":
                processed_read_names.add(readSet[0].query_name)
            start_ind = rs_reference_start - reference_mat_start
            end_ind = rs_reference_end - reference_mat_start

            end_ind_max = end_ind
            ref_lens = [_.reference_length for _ in readSet]
            max_ref_len_abs = max([abs(_) for _ in ref_lens])
            max_ref_num = ref_lens.index(max_ref_len_abs)
            max_ref_len = ref_lens[max_ref_num]

            masks = np.zeros([6, end_ind - start_ind], dtype=bool)
            masks[0, :] = snp_mask[start_ind:end_ind]
            masks[1, :] = noise_mask[start_ind:end_ind]
            masks[2, :] = n_cov_mask[start_ind:end_ind]
            masks[3, :] = include_mask[start_ind:end_ind]
            masks[4, :] = nm_mask[start_ind:end_ind]
            left, right = determineTrimLength(
                readSet[max_ref_num], params=params, processed_flag=processed_flag
            )
            masks[5, :left] = True
            masks[5, -right:] = True
            antimask = np.all(~masks, axis=0)
            antimask[trinuc_np[start_ind:end_ind] > 64] = False
            antimask[ref_np[start_ind:end_ind] == 4] = False
            # Create unmasked version that only excludes trinuc > 64
            unmasked_antimask = np.all(~masks[2:, :], axis=0)
            unmasked_antimask[trinuc_np[start_ind:end_ind] > 64] = False
            unmasked_antimask[ref_np[start_ind:end_ind] == 4] = False
            learn_antimask = np.all(~masks[3:, :], axis=0)
            learn_antimask[trinuc_np[start_ind:end_ind] > 64] = False
            learn_antimask[ref_np[start_ind:end_ind] == 4] = False
            ### If the whole reads are masked:
            if not np.any(unmasked_antimask):
                continue
            indel_bool = [
                ("I" in seq.cigarstring or "D" in seq.cigarstring) for seq in readSet
            ]
            # if any(indel_bool):
            if not isLearn:
                masks_indel = np.zeros([6, end_ind_max - start_ind], dtype=bool)
                masks_indel[0, :] = indel_mask[start_ind:end_ind_max]
                masks_indel[1, :] = noise_mask[start_ind:end_ind_max]
                masks_indel[2, :] = n_cov_mask[start_ind:end_ind_max]
                left, right = determineTrimLength(
                    readSet[max_ref_num], params=params, processed_flag=processed_flag
                )
                masks_indel[3, :left] = True
                masks_indel[3, -right:] = True
                masks_indel[4, :] = include_mask[start_ind:end_ind_max]
                masks_indel[5, :] = nm_mask[start_ind:end_ind_max]
                antimask_indel = np.all(~masks_indel, axis=0)
                unmasked_antimask_indel = np.all(~masks_indel[2:, :], axis=0)
                (
                    CS,
                    LR_raw,
                    LR_max,
                    indels,
                    hps,
                    strs,
                    F1R2_ref_count,
                    F1R2_alt_count,
                    F2R1_ref_count,
                    F2R1_alt_count,
                ) = genotypeDSIndel(
                    readSet,
                    rs_reference_start,
                    rs_reference_end,
                    ref_np[start_ind:end_ind],
                    unmasked_antimask_indel,
                    hp_raw_np[:, start_ind:end_ind],
                    str_raw_np[:, start_ind:end_ind],
                    params,
                )
                # print(LR_max)
                # pass_inds = np.nonzero(LR <= params["pcutoffi"])[0].tolist()

                # pass_inds = np.nonzero(LR_raw >= params["pcutoffi"])[0].tolist()
                # pass_inds = np.nonzero(LR_raw >= params["pcutoffi"])[0].tolist()
                pass_inds = np.nonzero(CS >= params["cscutoffi"])[0].tolist()
                # LR_pass_bool = (LR_raw >= params["pcutoffi"])
                if len(LR_raw) > 0:
                    LR_pass_bool = np.zeros(len(LR_raw), dtype=bool)
                    LR_pass_bool[
                        np.logical_and(LR_raw >= params["pcutoffi"][0], hps <= 5)
                    ] = True
                    LR_pass_bool[
                        np.all(
                            np.vstack(
                                [LR_raw >= params["pcutoffi"][1], hps <= 10, hps > 5]
                            ),
                            axis=0,
                        )
                    ] = True
                    LR_pass_bool[
                        np.logical_and(LR_raw >= params["pcutoffi"][2], hps > 10)
                    ] = True
                    LR_pass_bool[
                        np.logical_and(LR_raw >= params["pcutoffi"][2], strs > 0)
                    ] = True
                indels_pass = [indels[_] for _ in pass_inds]
                for nn in range(len(indels_pass)):
                    indel = indels_pass[nn]
                    indel_chrom = chromNow
                    indel_pos = int(indel.split(":")[0])
                    indel_size = int(indel.split(":")[1])
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
                    readPos = indel_pos - rs_reference_start
                    if readSet[0].template_length > 0:
                        readPos5p = min(
                            readPos + 1,
                            abs(readSet[0].template_length) - readPos,
                        )
                        readPos3p = min(
                            abs(max_ref_len) - readPos,
                            readPos + 1,
                        )
                    else:
                        readPos5p = min(
                            max_ref_len - readPos,
                            abs(readSet[0].template_length) - max_ref_len + readPos + 1,
                        )
                        readPos3p = min(
                            abs(max_ref_len) - readPos,
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
                    if indel_size > 0:
                        offset = 0
                    else:
                        offset = -indel_size
                    unmasked_flag = antimask_indel[
                        indel_pos
                        - reference_mat_start
                        - start_ind : indel_pos
                        - reference_mat_start
                        - start_ind
                        + offset
                        + 1
                    ].all()
                    if unmasked_flag and LR_pass_bool[nn]:
                        flt = flt_rs
                    elif unmasked_flag:
                        flt = "underpowered"
                    else:
                        flt = "masked"

                    indel_rec = {
                        "chrom": chromNow,
                        "pos": indel_pos + 1,
                        "ref": indel_ref,
                        "alt": indel_alt,
                        "filter": flt,
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
                            "CS": CS[pass_inds[nn]],
                            "LR": LR_raw[pass_inds[nn]],
                            "LM": LR_max[pass_inds[nn]],
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
                            "TL": readSet[0].template_length,
                            "STR": strs[nn],
                        },
                        "formats": ["AC", "RC", "DP"],
                        # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_indels.append(indel_rec)
                    indel_dict[indel_str] = 1
            # else:
            ### Calculate genotype probability
            # if not any(indel_bool) or isLearn:
            if 1:
                if isLearn and F1R2 > 2 and F2R1 > 2:
                    (
                        mismatch_now,
                        indelerr_now,
                        mismatch_dmg_now,
                        indel_dmg_now,
                    ) = profileTriNucMismatches(
                        readSet,
                        rs_reference_start,
                        ref_np[start_ind:end_ind],
                        trinuc_np[start_ind:end_ind],
                        hp_raw_np[:, start_ind:end_ind],
                        str_raw_np[:, start_ind:end_ind],
                        np.copy(learn_antimask),
                        params,
                    )
                    mismatch_mat += mismatch_now
                    indelerr_mat += indelerr_now
                    mismatch_dmg_mat += mismatch_dmg_now
                    indel_dmg_mat += indel_dmg_now
                if isLearn:
                    continue
                (
                    cov_mat,
                    CS_mut,
                    LR_raw_mut,
                    LR_max_mut,
                    mut_mask,
                    b1_int,
                    unmasked_antimask,
                    F1R2_count,
                    F2R1_count,
                ) = genotypeDSSnv(
                    readSet,
                    rs_reference_start,
                    ref_np[start_ind:end_ind],
                    trinuc_np[start_ind:end_ind],
                    prior_mat[start_ind:end_ind, :],
                    np.copy(unmasked_antimask),
                    params,
                    L,
                )
                # Per-position indel category coverage (del_unit, del_2,
                # del_3, del_4, del_5plus, ins_unit, ins_2, ins_3, ins_4,
                # ins_5plus): look up detection power from the
                # L_indel_1bp/L_indel_len tables built above, using this
                # window's repeat context (hp_np) and the same per-base
                # read-family counts (F1R2_count/F2R1_count) SBS coverage
                # uses as its depth proxy.
                if not isLearn:
                    unit_len_arr = hp_np[0, start_ind:end_ind].astype(int)
                    repeat_count_arr = hp_np[2, start_ind:end_ind].astype(int)
                    ref_allele_arr = ref_np[start_ind:end_ind]
                    ref_allele_safe = np.where(ref_allele_arr > 3, 0, ref_allele_arr)
                    is_hp = unit_len_arr <= 1
                    hps_for_row = np.where(is_hp, np.minimum(repeat_count_arr, 20), 1)
                    total_len_arr = unit_len_arr * repeat_count_arr
                    strs_for_row = np.zeros_like(total_len_arr)
                    strs_for_row[np.logical_and(~is_hp, total_len_arr >= 10)] = 1
                    strs_for_row[np.logical_and(~is_hp, total_len_arr >= 25)] = 2
                    strs_for_row[np.logical_and(~is_hp, total_len_arr >= 40)] = 3
                    row_len = np.where(
                        strs_for_row > 0, 19 + strs_for_row, hps_for_row - 1
                    )
                    unit_len_clamped = np.clip(unit_len_arr, 2, 5)
                    n_top_indel = np.minimum(F1R2_count.sum(axis=0), 9).astype(int)
                    n_bot_indel = np.minimum(F2R1_count.sum(axis=0), 9).astype(int)
                    # "Last cut"/"first cut" mask: repeat-based classes
                    # (ins_unit/del_unit here, and the homopolymer/STR raw
                    # classes below) can't reliably confirm a run's true
                    # length once this read family no longer extends past
                    # that run's start (hp_np row-1 "cut", from
                    # load_repeat_context), so everything from the last cut
                    # within this family's span onward is excluded — and
                    # symmetrically, the entire run touching the family's
                    # own LEFT edge is excluded too: whether an indel
                    # occurred within a run right at the start of what this
                    # family's reads cover is exactly as undecidable as one
                    # at the end, since the family has no visible flanking
                    # context on that side either way. This is the whole
                    # run containing window position 0 — not just its first
                    # base — so if the window itself starts exactly on a
                    # cut (that run's start is confirmed), the excluded
                    # span still runs to the *next* cut, not just past the
                    # first one. Does not apply to del_2..del_5plus /
                    # microhomology, which use a fixed context regardless
                    # of any actual run.
                    cut_arr = hp_np[1, start_ind:end_ind]
                    last_cut_valid = np.ones(end_ind - start_ind, dtype=bool)
                    cut_positions = np.nonzero(cut_arr)[0]
                    if cut_positions.size > 0:
                        last_cut_valid[cut_positions[-1] :] = False
                        if cut_positions[0] == 0:
                            first_run_end = (
                                cut_positions[1]
                                if cut_positions.size > 1
                                else end_ind - start_ind
                            )
                        else:
                            first_run_end = cut_positions[0]
                        last_cut_valid[:first_run_end] = False
                    # Reference base immediately following each position,
                    # needed for the "Insertion A/T/C/G" columns below and
                    # the matching fine-grained 1bpins{base} columns
                    # (96-99): rep0 for inserting base N requires only
                    # that this next base differs from N — nothing to do
                    # with the position's own reference base.
                    window_len = end_ind - start_ind
                    next_ref_arr = np.full(window_len, -1, dtype=int)
                    avail = min(end_ind + 1, ref_np.shape[0]) - (start_ind + 1)
                    if avail > 0:
                        next_ref_arr[:avail] = ref_np[
                            start_ind + 1 : start_ind + 1 + avail
                        ]
                    next_ref_valid = (next_ref_arr >= 0) & (next_ref_arr <= 3)

                    cov_mat_indel = np.zeros([window_len, 14])
                    cov_mat_indel[:, 0] = (
                        np.where(
                            is_hp,
                            L_indel_1bp[
                                n_top_indel,
                                n_bot_indel,
                                hps_for_row,
                                ref_allele_safe,
                                0,
                            ],
                            L_indel_len[
                                n_top_indel,
                                n_bot_indel,
                                row_len,
                                -unit_len_clamped + 5,
                            ],
                        )
                        * last_cut_valid
                    )
                    # Insertion of Repeat Unit: same shape as column 0
                    # (deletion of repeat unit) but the insertion direction
                    # (last L_indel_1bp index 1, and +unit_len_clamped+5
                    # instead of -unit_len_clamped+5 into L_indel_len).
                    cov_mat_indel[:, 1] = (
                        np.where(
                            is_hp,
                            L_indel_1bp[
                                n_top_indel,
                                n_bot_indel,
                                hps_for_row,
                                ref_allele_safe,
                                1,
                            ],
                            L_indel_len[
                                n_top_indel,
                                n_bot_indel,
                                row_len,
                                unit_len_clamped + 5,
                            ],
                        )
                        * last_cut_valid
                    )
                    # del_2..del_5plus: fixed homopolymer-length-1 context
                    # (microhomology-style, row 0 = hps=1), not the
                    # position's actual repeat context.
                    cov_mat_indel[:, 2] = L_indel_len[n_top_indel, n_bot_indel, 0, 3]
                    cov_mat_indel[:, 3] = L_indel_len[n_top_indel, n_bot_indel, 0, 2]
                    cov_mat_indel[:, 4] = L_indel_len[n_top_indel, n_bot_indel, 0, 1]
                    cov_mat_indel[:, 5] = L_indel_len[n_top_indel, n_bot_indel, 0, 0]
                    # Insertion A/T/C/G: rep0 opportunity to insert exactly
                    # that base — fixed hps=1 context, using THAT base (not
                    # the position's own reference base) as the
                    # L_indel_1bp ref_allele index — zeroed wherever the
                    # next reference base already equals it (that's a
                    # repeat extension, not a novel insertion).
                    for b in range(4):
                        cov_mat_indel[:, 6 + b] = (
                            L_indel_1bp[n_top_indel, n_bot_indel, 1, b, 1]
                            * next_ref_valid
                            * (next_ref_arr != b)
                            * antimask_indel
                        )
                    # Insertion Length 2..5+: power to see a novel U-bp unit
                    # inserted where no repeat exists yet (rep0 in
                    # classify_indel_record's is_hp/indel_len>=2 branch,
                    # misc.py) — fixed hps=1 context (row 0), same as
                    # del_2..del_5plus, but only where a "no repeat here yet"
                    # call is even possible (is_hp) and the read family can
                    # still confirm that (last_cut_valid).
                    ins_len_valid = is_hp & last_cut_valid & antimask_indel
                    cov_mat_indel[:, 10] = (
                        L_indel_len[n_top_indel, n_bot_indel, 0, 7] * ins_len_valid
                    )
                    cov_mat_indel[:, 11] = (
                        L_indel_len[n_top_indel, n_bot_indel, 0, 8] * ins_len_valid
                    )
                    cov_mat_indel[:, 12] = (
                        L_indel_len[n_top_indel, n_bot_indel, 0, 9] * ins_len_valid
                    )
                    cov_mat_indel[:, 13] = (
                        L_indel_len[n_top_indel, n_bot_indel, 0, 10] * ins_len_valid
                    )
                    # Fine-grained 100-class raw indel classification
                    # (mirrors the SBS 192-class raw trinuc scheme):
                    # base-specific homopolymer del (24) + ins (20, no
                    # rep0 bin — see below), unit-size/repeat-count-
                    # specific STR del/ins (24+24), microhomology deletion
                    # length (4). See funcs/misc.py's build_indel100_labels
                    # for the exact column layout.
                    indel100 = np.zeros(100)

                    def _accumulate_indel100(col_idx, weight, valid):
                        if np.any(valid):
                            indel100[:] += np.bincount(
                                col_idx[valid].astype(int),
                                weights=weight[valid],
                                minlength=100,
                            )

                    to_cgta = np.array([3, 2, 0, 1])  # base2num->CGTA order
                    base_cgta = to_cgta[ref_allele_safe]

                    # HP and STR opportunity below are credited
                    # INDEPENDENTLY, from the raw (unmerged) hp.h5/str.h5
                    # arrays (hp_raw_np/str_raw_np) — not from
                    # is_hp/unit_len_arr/cut_arr above, which come from
                    # load_repeat_context's STR-priority merge and
                    # silently drop a position's homopolymer opportunity
                    # wherever it also sits inside an annotated STR tract
                    # (e.g. the "AA" in a (AAT)n repeat). That merge is
                    # correct for classifying an *observed* mutation event
                    # (necessarily one category or the other) but wrong
                    # here: such a position is a real candidate for both a
                    # 1bp homopolymer slip and a larger STR-unit slip.
                    # Mirrors funcs/misc.py's
                    # indel100_reference_bucket_indices, which already
                    # double-counts these positions into both buckets on
                    # the genome-composition side — see
                    # [[feedback-hp-str-independence]].
                    def _run_boundary_valid(cut_bool_arr):
                        valid = np.ones(cut_bool_arr.shape[0], dtype=bool)
                        cut_positions = np.nonzero(cut_bool_arr)[0]
                        if cut_positions.size > 0:
                            valid[cut_positions[-1] :] = False
                            if cut_positions[0] == 0:
                                first_run_end = (
                                    cut_positions[1]
                                    if cut_positions.size > 1
                                    else cut_bool_arr.shape[0]
                                )
                            else:
                                first_run_end = cut_positions[0]
                            valid[:first_run_end] = False
                        return valid

                    hp_run_arr = np.minimum(
                        hp_raw_np[0, start_ind:end_ind].astype(int), 20
                    )
                    hp_cut_bool = hp_raw_np[1, start_ind:end_ind].astype(bool)
                    hp_repeat_valid = _run_boundary_valid(hp_cut_bool) & antimask_indel

                    str_unit_raw = str_raw_np[0, start_ind:end_ind].astype(int)
                    str_repeat_raw = str_raw_np[1, start_ind:end_ind].astype(int)
                    str_cut_bool = str_raw_np[2, start_ind:end_ind].astype(bool)
                    is_real_str = str_unit_raw >= 2
                    str_repeat_valid = (
                        _run_boundary_valid(str_cut_bool) & antimask_indel
                    )
                    total_len_str = str_unit_raw * str_repeat_raw
                    strs_for_row_str = np.zeros_like(total_len_str)
                    strs_for_row_str[is_real_str & (total_len_str >= 10)] = 1
                    strs_for_row_str[is_real_str & (total_len_str >= 25)] = 2
                    strs_for_row_str[is_real_str & (total_len_str >= 40)] = 3
                    row_len_str = np.where(
                        strs_for_row_str > 0, 19 + strs_for_row_str, 0
                    )
                    unit_len_clamped_str = np.clip(str_unit_raw, 2, 5)
                    unit_bucket = np.clip(str_unit_raw, 2, 5) - 2
                    count_bucket_del = np.clip(str_repeat_raw, 1, 6) - 1

                    # 1-4: homopolymer deletion (cols 0-23) — always from
                    # the position's own raw hp run/cut, regardless of any
                    # STR annotation there too.
                    hp_del_power = L_indel_1bp[
                        n_top_indel, n_bot_indel, hp_run_arr, ref_allele_safe, 0
                    ]
                    del_bucket = np.clip(hp_run_arr, 1, 6) - 1
                    _accumulate_indel100(
                        base_cgta * 6 + del_bucket,
                        hp_del_power,
                        hp_repeat_valid & hp_cut_bool,
                    )

                    # 5-8: homopolymer insertion (cols 24-43) — length
                    # bins 1..5+ only; rep0 (hp_run_arr==0, never actually
                    # occurs) is deliberately not credited here, since that
                    # value is always discarded and replaced by the exact
                    # 1bpins{base} next-base computation below (see
                    # Estimate.py's
                    # override_inshp0_with_next_base_opportunity).
                    hp_ins_power = L_indel_1bp[
                        n_top_indel, n_bot_indel, hp_run_arr, ref_allele_safe, 1
                    ]
                    ins_bucket = np.clip(hp_run_arr, 1, 5) - 1
                    _accumulate_indel100(
                        24 + base_cgta * 5 + ins_bucket,
                        hp_ins_power,
                        (hp_run_arr >= 1) & hp_repeat_valid & hp_cut_bool,
                    )

                    # 9: STR deletion (cols 44-67) — real part, from the
                    # position's own raw STR annotation, independent of
                    # any homopolymer overlap.
                    str_del_power = L_indel_len[
                        n_top_indel, n_bot_indel, row_len_str, -unit_len_clamped_str + 5
                    ]
                    _accumulate_indel100(
                        44 + unit_bucket * 6 + count_bucket_del,
                        str_del_power,
                        is_real_str & str_repeat_valid & str_cut_bool,
                    )
                    # Flat rep1 opportunity at positions with no real STR
                    # annotation (whether or not they're also a
                    # homopolymer run — that's credited separately above):
                    # every such position is a candidate for an arbitrary,
                    # non-repeating U-bp deletion — a fixed hps=1/strs=0
                    # power lookup (row 0), not an attempt to discover a
                    # true higher count that might coincidentally exist
                    # there (that would need rescanning the reference; see
                    # funcs/misc.py's indel100_reference_bucket_indices for
                    # the matching scan-free credit on the genome-wide
                    # side). Gated by hp_repeat_valid, not
                    # str_repeat_valid: there's no real STR tract here to
                    # have its own boundary, so this reuses the
                    # run-visibility this position already has under its
                    # own (possibly trivial, length-1) homopolymer context.
                    for u_idx, U in enumerate([2, 3, 4, 5]):
                        hyp_del_power = L_indel_len[n_top_indel, n_bot_indel, 0, -U + 5]
                        _accumulate_indel100(
                            np.full(end_ind - start_ind, 44 + u_idx * 6),
                            hyp_del_power,
                            (~is_real_str) & hp_repeat_valid,
                        )

                    # 10: STR insertion (cols 68-91) — real part
                    str_ins_power = L_indel_len[
                        n_top_indel, n_bot_indel, row_len_str, unit_len_clamped_str + 5
                    ]
                    count_bucket_ins = np.clip(str_repeat_raw, 0, 5)
                    _accumulate_indel100(
                        68 + unit_bucket * 6 + count_bucket_ins,
                        str_ins_power,
                        is_real_str & str_repeat_valid & str_cut_bool,
                    )
                    is_rep1 = str_repeat_raw == 1
                    _accumulate_indel100(
                        68 + unit_bucket * 6 + 0,
                        str_ins_power,
                        is_real_str & is_rep1 & str_repeat_valid & str_cut_bool,
                    )
                    # Flat rep0/rep1 opportunity at positions with no real
                    # STR annotation ("what if a U-bp unit were inserted
                    # here"): fixed hps=1/strs=0 context (row 0), same
                    # value added to both the "0" and "1" repeat-count
                    # buckets for every unit size — scan-free, matching the
                    # deletion case above, and gated the same way (see
                    # comment there for why hp_repeat_valid, not
                    # str_repeat_valid).
                    for u_idx, U in enumerate([2, 3, 4, 5]):
                        hyp_power = L_indel_len[n_top_indel, n_bot_indel, 0, U + 5]
                        base_col = 68 + u_idx * 6
                        _accumulate_indel100(
                            np.full(end_ind - start_ind, base_col),
                            hyp_power,
                            (~is_real_str) & hp_repeat_valid,
                        )
                        _accumulate_indel100(
                            np.full(end_ind - start_ind, base_col + 1),
                            hyp_power,
                            (~is_real_str) & hp_repeat_valid,
                        )

                    # 11: microhomology deletion (cols 92-95) — fixed
                    # context, every position, no last-cut restriction;
                    # reuse cov_mat_indel[:,2:6].
                    for k in range(4):
                        _accumulate_indel100(
                            np.full(end_ind - start_ind, 92 + k),
                            cov_mat_indel[:, 2 + k],
                            antimask_indel,
                        )

                    # 12: 1bpins{base} rep0 opportunity (cols 96-99) —
                    # identical per-position values to
                    # cov_mat_indel[:,6:10] (Insertion A/T/C/G), so just
                    # reuse them rather than recomputing; the
                    # antimask_indel/next-base gating is already baked
                    # into those columns.
                    for b in range(4):
                        _accumulate_indel100(
                            np.full(window_len, 96 + b),
                            cov_mat_indel[:, 6 + b],
                            np.ones(window_len, dtype=bool),
                        )
                """
                (
                    LR_old,
                    b1_int_old,
                    unmasked_antimask_old,
                    F1R2_count_old,
                    F2R1_count_old,
                    prob1_old,prob2_old,prob3_old,prob4_old
                ) = oldDSSnv(
                    readSet,
                    ref_np[start_ind:end_ind],
                    trinuc_np[start_ind:end_ind],
                    prior_mat[start_ind:end_ind, :],
                    np.copy(unmasked_antimask),
                    params,
                )
                """
                # prob1_diff = prob1-prob1_old
                # notice1 = np.abs(prob1_diff) >=1
                # prob3_diff = prob3-prob3_old
                # notice2 = np.abs(prob3_diff) >=1
                # if np.any(notice1):
                # print("pr1",prob1[notice1],prob1_old[notice1],F1R2_count[:,unmasked_antimask_old][:,notice1])

                ref_int = ref_np[start_ind:end_ind]
                n_win = ref_int.size
                mut_pos_in_win = np.nonzero(mut_mask)[0]
                if CS_mut.size > 0:
                    muts_ind_compressed = np.nonzero(CS_mut >= params["cscutoff"])[0]
                    muts_ind = mut_pos_in_win[muts_ind_compressed].tolist()
                else:
                    muts_ind_compressed = np.zeros(0, dtype=int)
                    muts_ind = []
                refs_ind = np.nonzero(
                    np.logical_and(unmasked_antimask, b1_int == ref_int)
                )[0].tolist()
                LR_pass_bool = np.zeros(n_win, dtype=bool)
                LR_pass_bool[mut_mask] = LR_raw_mut >= params["pcutoff"]
                alt_int = b1_int
                unmasked_pass_bool = np.full(n_win, False, dtype=bool)
                unmasked_pass_bool[refs_ind] = True
                unmasked_pass_bool[muts_ind] = True
                pass_bool = np.copy(unmasked_pass_bool)
                pass_bool[~antimask] = False
                pos = [
                    mut_ind + start_ind + reference_mat_start for mut_ind in muts_ind
                ]
                # muts_ind = np.nonzero(np.logical_and(mut_bool,pass_bool))[0].tolist()
                mut_positions = [
                    mut_ind + start_ind + reference_mat_start + 1
                    for mut_ind in muts_ind
                ]
                NMs = [seq.get_tag("NM") for seq in readSet]
                for nn in range(len(mut_positions)):
                    mut_chrom = reference_mat_chrom
                    mut_pos = mut_positions[nn]
                    mut_ref = num2base[ref_int[muts_ind[nn]]]
                    mut_alt = num2base[alt_int[muts_ind[nn]]]
                    mut_trinuc = num2trinuc[trinuc_np[start_ind:end_ind][muts_ind[nn]]]
                    if readSet[0].template_length > 0:
                        readPos5p = min(
                            muts_ind[nn] + 1,
                            abs(readSet[0].template_length) - muts_ind[nn],
                        )
                        readPos3p = min(
                            abs(max_ref_len) - muts_ind[nn],
                            muts_ind[nn] + 1,
                        )
                    else:
                        readPos5p = min(
                            max_ref_len - muts_ind[nn],
                            abs(readSet[0].template_length)
                            - max_ref_len
                            + muts_ind[nn]
                            + 1,
                        )
                        readPos3p = min(
                            abs(max_ref_len) - muts_ind[nn],
                            muts_ind[nn] + 1,
                        )
                    FPs.append(readPos5p)
                    RPs.append(readPos3p)
                    if F1R2_count[:, muts_ind[nn]].sum() == 0:
                        continue
                    if F2R1_count[:, muts_ind[nn]].sum() == 0:
                        continue
                    if pass_bool[muts_ind[nn]] and LR_pass_bool[muts_ind[nn]]:
                        flt = flt_rs
                    elif pass_bool[muts_ind[nn]]:
                        flt = "underpowered"
                    else:
                        flt = "masked"
                    mut = {
                        "chrom": mut_chrom,
                        "pos": mut_pos,
                        "ref": mut_ref,
                        "alt": mut_alt,
                        "filter": flt,
                        "infos": {
                            "F1R2": F1R2,
                            "F2R1": F2R1,
                            "CS": CS_mut[muts_ind_compressed[nn]],
                            "LR": LR_raw_mut[muts_ind_compressed[nn]],
                            "LM": LR_max_mut[muts_ind_compressed[nn]],
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
                            "TL": readSet[0].template_length,
                            "STR": 0,
                        },
                        "formats": ["AC", "RC", "DP"],
                        # "samples": [[ta, tr, tdp], [na, nr, ndp]],
                    }
                    muts_dict["_".join([mut_chrom, str(mut_pos), mut_ref, mut_alt])] = 0
                    muts.append(mut)
                muts_dbs.extend(
                    _detect_dbs_pairs(
                        reference_mat_chrom,
                        mut_positions,
                        muts_ind,
                        ref_int,
                        alt_int,
                        pass_bool,
                        LR_pass_bool,
                        flt_rs,
                        F1R2,
                        F2R1,
                        setBc,
                        currentStart,
                        readSet[0].template_length,
                        num2base,
                    )
                )
                """
                if isLearn:
                    continue
                """
                if flt_rs == "PASS":
                    coverage[start_ind:end_ind][pass_bool] += cov_mat[pass_bool]
                    duplex_read_num_dict[duplex_no][1] += np.count_nonzero(pass_bool)
                    trinuc_pass = trinuc_np[start_ind:end_ind][pass_bool]
                    cov_pass = cov_mat[pass_bool]
                    for b in range(4):
                        duplex_read_num_dict_trinuc[duplex_no][:, b] += np.bincount(
                            trinuc_pass, weights=cov_pass[:, b], minlength=64
                        )
                    duplex_read_num_dict_dbs[duplex_no] += _compute_dbs_opportunity(
                        cov_mat, ref_int, pass_bool
                    )

                    # Update unmasked coverage and trinuc counts (includes all passing sites)

                    unmasked_coverage[start_ind:end_ind][unmasked_pass_bool] += cov_mat[
                        unmasked_pass_bool
                    ]
                    if pass_bool.any():
                        duplex_read_num_dict[duplex_no][0] += 1
                        duplex_count += 1
                    if not isLearn:
                        coverage_indel_cat[start_ind:end_ind][
                            antimask_indel
                        ] += cov_mat_indel[antimask_indel]
                        unmasked_coverage_indel_cat[start_ind:end_ind][
                            unmasked_antimask_indel
                        ] += cov_mat_indel[unmasked_antimask_indel]
                        duplex_read_num_dict_indel[duplex_no] += indel100
    """
    Calling block ends
    """
    if isLearn:
        return mismatch_mat, indelerr_mat, mismatch_dmg_mat, indel_dmg_mat

    ### SNV depth extraction, batched: gather every unique PASS/masked/
    ### underpowered candidate first, then resolve the tumor BAM and each
    ### normal BAM with a handful of sequential pileup() scans
    ### (extractDepthBatchSnv) instead of one random-access pileup() call
    ### per candidate.
    mut_dict = dict()
    snv_candidate_keys = []
    for mut in muts:
        if mut["filter"] not in ("PASS", "masked", "underpowered"):
            continue
        key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
        if key not in mut_dict:
            mut_dict[key] = None
            snv_candidate_keys.append(key)

    if snv_candidate_keys:
        tumor_snv_depths = extractDepthBatchSnv(
            tumorBam, snv_candidate_keys, params, minbq=params["minBq"]
        )
        normal_snv_depths = (
            [
                extractDepthBatchSnv(
                    normalBam, snv_candidate_keys, params, minbq=params["minBq"]
                )
                for normalBam in normalBams
            ]
            if normalBams
            else []
        )
        for key in snv_candidate_keys:
            ta, tr, ti, tdp = tumor_snv_depths.get(key, (0, 0, 0, 0))
            if normalBams:
                na = nr = ni = ndp = 0
                for normal_depths in normal_snv_depths:
                    na_now, nr_now, ni_now, ndp_now = normal_depths.get(
                        key, (0, 0, 0, 0)
                    )
                    na += na_now
                    nr += nr_now
                    ni += ni_now
                    ndp += ndp_now
            else:
                na, nr, ni, ndp = (0, 0, 0, 0)
            mut_dict[key] = (ta, tr, ti, tdp, na, nr, ni, ndp)

    mut_pass_filter = []
    n_muts = len(muts)
    snv_filter_progress = {"start": time.time(), "last": time.time()}
    for i, mut in enumerate(muts):
        chrom = mut["chrom"]
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]
        flt = mut["filter"]
        log_progress(
            snv_filter_progress,
            f"Process {processNo} filtering SNV mutations",
            i,
            n_muts,
            extra=f"{chrom}:{pos}",
        )
        if flt == "PASS" or flt == "masked" or flt == "underpowered":
            ta, tr, ti, tdp, na, nr, ni, ndp = mut_dict[(chrom, pos, ref, alt)]
        else:
            ta, tr, ti, tdp, na, nr, ni, ndp = (0, 0, 0, 0, 0, 0, 0, 0)
        # if window_filter:
        # continue
        # if ta > params["maxAltCount"]:
        # continue
        if flt == "PASS" or flt == "masked" or flt == "underpowered":
            if ta == 0:
                continue
            if ta / tdp > params["maxAF"]:
                continue
            # if ti >= 1:
            # continue
            if normalBams:
                if ndp < params["minNdepth"]:
                    continue
                if na / ndp > params["normalVAF"]:
                    continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        mut_pass_filter.append(mut)

    ### Indel depth extraction, batched the same way.
    muts_indels_dict = dict()
    indel_candidate_keys = []
    for mut in muts_indels:
        if mut["filter"] not in ("PASS", "masked", "underpowered"):
            continue
        key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
        if key not in muts_indels_dict:
            muts_indels_dict[key] = None
            indel_candidate_keys.append(key)

    if indel_candidate_keys:
        tumor_indel_depths = extractDepthBatchIndel(
            tumorBam, indel_candidate_keys, params, minbq=params["minBq"]
        )
        normal_indel_depths = (
            [
                extractDepthBatchIndel(
                    normalBam, indel_candidate_keys, params, minbq=params["minBq"]
                )
                for normalBam in normalBams
            ]
            if normalBams
            else []
        )
        for key in indel_candidate_keys:
            ta, tr, ti, tdp = tumor_indel_depths.get(key, (0, 0, 0, 0))
            if normalBams:
                na = nr = ni = ndp = 0
                for normal_depths in normal_indel_depths:
                    na_now, nr_now, ni_now, ndp_now = normal_depths.get(
                        key, (0, 0, 0, 0)
                    )
                    na += na_now
                    nr += nr_now
                    ni += ni_now
                    ndp += ndp_now
            else:
                na, nr, ni, ndp = (0, 0, 0, 0)
            # na/nr/ndp match the pre-batching dict layout; ni (summed
            # normal-BAM otherIndelCount) is additionally kept so the
            # "if ni > 0" filter below always sees the real value instead
            # of a stale ni left over from whatever mutation the
            # single-pass loop happened to process previously (see
            # commit message for details).
            muts_indels_dict[key] = (ta, tr, tdp, na, nr, ndp, ni)

    muts_indels_pass_filter = []
    n_muts_indels = len(muts_indels)
    indel_filter_progress = {"start": time.time(), "last": time.time()}
    for i, mut in enumerate(muts_indels):
        chrom = mut["chrom"]
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]
        flt = mut["filter"]
        log_progress(
            indel_filter_progress,
            f"Process {processNo} filtering indel mutations",
            i,
            n_muts_indels,
            extra=f"{chrom}:{pos}",
        )
        if flt == "PASS" or flt == "masked" or flt == "underpowered":
            ta, tr, tdp, na, nr, ndp, ni = muts_indels_dict[(chrom, pos, ref, alt)]
        else:
            ta, tr, ti, tdp, na, nr, ni, ndp = (0, 0, 0, 0, 0, 0, 0, 0)
        if flt == "PASS" or flt == "masked" or flt == "underpowered":
            if ta == 0:
                continue
            # if ti > 0:
            # continue
            if ta / tdp > params["maxAF"]:
                continue
            if normalBams:
                if ndp < params["minNdepth"]:
                    continue
                if na / ndp > params["normalVAF"]:
                    continue
                if ni > 0:
                    continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        muts_indels_pass_filter.append(mut)

    ### DBS depth extraction, batched the same way, and filtered under the
    ### same maxAF/minNdepth/normalVAF strategy as SNVs/indels above.
    ### _detect_dbs_pairs only ever emits filter=="PASS" candidates (both
    ### constituent bases already independently passed full duplex-consensus
    ### calling), so — unlike the SNV/indel loops — there's no
    ### masked/underpowered filter state to preserve here: a candidate either
    ### clears this gate and is kept, or it's dropped.
    muts_dbs_dict = dict()
    dbs_candidate_keys = []
    for mut in muts_dbs:
        key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
        if key not in muts_dbs_dict:
            muts_dbs_dict[key] = None
            dbs_candidate_keys.append(key)

    if dbs_candidate_keys:
        tumor_dbs_depths = extractDepthBatchDbs(
            tumorBam, dbs_candidate_keys, params, minbq=params["minBq"]
        )
        normal_dbs_depths = (
            [
                extractDepthBatchDbs(
                    normalBam, dbs_candidate_keys, params, minbq=params["minBq"]
                )
                for normalBam in normalBams
            ]
            if normalBams
            else []
        )
        for key in dbs_candidate_keys:
            ta, tr, ti, tdp = tumor_dbs_depths.get(key, (0, 0, 0, 0))
            if normalBams:
                na = nr = ni = ndp = 0
                for normal_depths in normal_dbs_depths:
                    na_now, nr_now, ni_now, ndp_now = normal_depths.get(
                        key, (0, 0, 0, 0)
                    )
                    na += na_now
                    nr += nr_now
                    ni += ni_now
                    ndp += ndp_now
            else:
                na, nr, ni, ndp = (0, 0, 0, 0)
            muts_dbs_dict[key] = (ta, tr, ti, tdp, na, nr, ni, ndp)

    muts_dbs_pass_filter = []
    for mut in muts_dbs:
        chrom = mut["chrom"]
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]
        ta, tr, ti, tdp, na, nr, ni, ndp = muts_dbs_dict[(chrom, pos, ref, alt)]
        if ta == 0:
            continue
        if ta / tdp > params["maxAF"]:
            continue
        if normalBams:
            if ndp < params["minNdepth"]:
                continue
            if na / ndp > params["normalVAF"]:
                continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        muts_dbs_pass_filter.append(mut)

    if "coverage" in locals():
        if "coverage_leftover" in locals():
            coverage[0 : coverage_leftover.shape[0]] += coverage_leftover
            coverage_indel_cat[
                0 : coverage_leftover.shape[0]
            ] += coverage_indel_cat_leftover
            unmasked_coverage[
                0 : coverage_leftover.shape[0]
            ] += unmasked_coverage_leftover
            unmasked_coverage_indel_cat[
                0 : coverage_leftover.shape[0]
            ] += unmasked_coverage_indel_cat_leftover
        non_zero_positions = np.nonzero(
            coverage.sum(axis=1) + coverage_indel_cat.sum(axis=1)
        )
        for pos in non_zero_positions[0].tolist():
            current_pos = pos + reference_mat_start
            bed_file = get_bed_file_for_position(
                current_pos,
                reference_mat_chrom,
                regions_start_chrom,
                regions_start_pos,
                regions_end_chrom,
                regions_end_pos,
                locus_bed,
                locus_bed_prev,
                locus_bed_next,
            )
            bed_file.write(
                "\t".join(
                    [
                        reference_mat_chrom,
                        str(current_pos),
                        str(current_pos + 1),
                        "\t".join(str(v) for v in coverage[pos]),
                        "\t".join(str(v) for v in coverage_indel_cat[pos]),
                    ]
                )
                + "\n"
            )
            total_coverage += coverage[pos]
            total_coverage_indel_cat += coverage_indel_cat[pos]
            total_unmasked_coverage += unmasked_coverage[pos]
            total_unmasked_coverage_indel_cat += unmasked_coverage_indel_cat[pos]

    # Close the bed files
    locus_bed.close()
    locus_bed_prev.close()
    locus_bed_next.close()

    print(
        f"Process {processNo} finished in {(time.time()-starttime)/60: .2f} minutes and processed {recCount} reads"
    )

    # duplex_read_num_dict_trinuc[duplex_no] stays a raw (64, 4) matrix here
    # (64 trinuc contexts, un-folded by reverse complement); Caller.py's
    # writer selects the 192 non-self (context, alt) cells for output, and
    # reverse-complement pairs are only combined at estimation time.

    return (
        mut_pass_filter,
        total_coverage.sum(),
        recCount,
        duplex_count,
        duplex_read_num_dict,
        duplex_read_num_dict_trinuc,
        muts_indels_pass_filter,
        unique_read_num,
        pass_read_num,
        FPs,
        RPs,
        total_unmasked_coverage.sum(),
        total_coverage_indel_cat,
        total_unmasked_coverage_indel_cat,
        duplex_read_num_dict_indel,
        muts_dbs_pass_filter,
        duplex_read_num_dict_dbs,
    )
