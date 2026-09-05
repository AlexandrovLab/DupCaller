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
import h5py
import errno
import copy
from scipy.optimize import brentq

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
    indelMaxLR,
    calculateSSPosterior,
    calculateDSPosterior,
)
from .learn import profileTriNucMismatches, NUM_BQ
from .misc import getAlignmentObject as BAM
from .misc import build_trinuc64_order
from .misc import load_repeat_context
from .misc import log_progress
from .misc import get_duplex_barcode
from .misc import (
    determineTrimLength,
    nums2str,
    get_bed_file_for_position,
    _filter_reason_label,
    bamIterateMultipleRegion,
    _detect_dbs_pairs,
    _compute_dbs_opportunity,
    _accumulate_depth_matrix,
    threshold_rng,
    simulate_power_grid,
    load_error_matrices,
    _compute_read_label,
    _place_read_in_dict,
    _index_rugged_mates,
    _drain_rugged_pool,
)
from .indels import findIndels


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
    af_cutoff = params["germline_cutoff"]
    m = end - start
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    trinuc2num = params["trinuc2num_dict"]

    snp_mask = np.full(m, False, dtype=bool)
    indel_mask = np.full(m, False, dtype=bool)
    noise_mask = np.full(m, False, dtype=bool)
    n_cov_mask = np.full(m, False, dtype=bool)
    nm_mask = np.full(m, False, dtype=bool)

    noise_mask[reference_int == 4] = True
    noise_mask[trinuc_int == 96] = True
    reference_int[reference_int == 4] = 0

    ### Adjust by germline
    if germline_bed != None:
        for rec in germline_bed.fetch(chrom, start, end):
            ind = rec.pos - 1 - start
            ref = rec.ref
            try:
                afs = rec.info["AF"]
            except KeyError:
                afs = [1 for _ in range(len(rec.alts))]
            if len(ref) == 1:
                for ii, alt in enumerate(rec.alts):
                    if len(alt) == 1:
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
        snp_mask,
        indel_mask,
        noise_mask,
        n_cov_mask,
        include_mask,
        nm_mask,
    )  # , reference_int, trinuc_int


def _process_duplex_family(
    F1R2,
    F2R1,
    FPs,
    L,
    L_indel_1bp,
    L_indel_len,
    RPs,
    chromNow,
    coverage,
    coverage_indel_cat,
    currentStart,
    depth_by_dinuc,
    depth_by_hpstr,
    depth_by_trinuc,
    duplex_no,
    duplex_read_num_dict,
    duplex_read_num_dict_dbs,
    duplex_read_num_dict_indel,
    duplex_read_num_dict_trinuc,
    flt_rs,
    hp_alt_mat,
    hp_dmg_mat,
    hp_np,
    hp_raw_np,
    include_mask,
    indel_dict,
    indel_mask,
    isLearn,
    mismatch_dmg_mat,
    mismatch_mat,
    muts,
    muts_dbs,
    muts_dict,
    muts_indels,
    n_cov_mask,
    nm_mask,
    noise_mask,
    num2base,
    num2trinuc,
    params,
    processed_read_names,
    readSet,
    ref_np,
    reference_mat_chrom,
    reference_mat_start,
    rs_reference_end,
    rs_reference_start,
    sbs_alt_bq_hist_mat,
    setBc,
    snp_mask,
    str_alt_mat,
    str_dmg_mat,
    str_raw_np,
    trinuc_np,
    unmasked_coverage,
    unmasked_coverage_indel_cat,
):
    """Genotype and record one candidate duplex family, mutating the
    caller's per-process accumulators (coverage, indel100, VCF-record
    lists, error-learning matrices, depth-composition tables) in place.

    Extracted from callBam's two near-identical "process a completed
    duplex family" sites (the mid-stream trigger and the end-of-region
    flush) -- both call this with identical arguments, verified via AST
    analysis to read exactly the same 50 names from outside their own
    scope and touch no window-management state (reference_mat_*/coverage
    allocation/leftover-carry, still duplicated inline at each call site
    since that part is a locals()-based state machine spanning multiple
    trigger firings, not safely extractable without first replacing that
    state machine).

    All array/dict/list arguments (coverage, muts, duplex_read_num_dict*,
    depth_by_*, indel_dict, muts_dict, processed_read_names, etc.) are
    the caller's own persistent objects, mutated here in place -- nothing
    needs to flow back out through a return value except whether this
    family counts toward the caller's duplex_count tally, since a plain
    scalar int (unlike these containers) can't be mutated through a
    function call.

    Returns True if this family produced at least one passing (pass_bool)
    genotype call -- the caller does `if _process_duplex_family(...):
    duplex_count += 1` -- False for every early-exit path (masked/learn-
    mode/no-pass-call), mirroring the original inline code's `continue`
    at those same points.
    """
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
    # Rescue eligibility -- see the matching block earlier in this
    # function for the full explanation.
    rescue_antimask = ~masks[3, :]
    rescue_antimask[trinuc_np[start_ind:end_ind] > 64] = False
    rescue_antimask[ref_np[start_ind:end_ind] == 4] = False
    learn_antimask = np.all(~masks[3:, :], axis=0)
    learn_antimask[trinuc_np[start_ind:end_ind] > 64] = False
    learn_antimask[ref_np[start_ind:end_ind] == 4] = False
    ### If the whole reads are masked: under --rescue, only include_mask is
    ### truly non-rescuable (rescue_antimask == ~include_mask, matching the
    ### per-site rescue-eligibility checks below), so gate on that wider
    ### definition of "has an analyzable position" instead of
    ### unmasked_antimask (which also requires clearing n_cov/nm/trim).
    ### Otherwise a family blocked only by those masks across its whole
    ### span would be dropped here before ever reaching the per-site
    ### logic that's supposed to report it as a filtered-but-rescued call.
    if params["rescue"]:
        if not np.any(rescue_antimask):
            return False  # no site in this family survives even rescue
    elif not np.any(unmasked_antimask):
        return False  # this family never reaches pass_bool -- doesn't count
    indel_bool = [("I" in seq.cigarstring or "D" in seq.cigarstring) for seq in readSet]
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
        # Ignores noise_mask only -- see the matching block
        # earlier in this function for the full explanation.
        unmasked_antimask_indel = np.all(~masks_indel[[0, 2, 3, 4, 5], :], axis=0)
        # Same rescue-eligibility split as the SNV side above.
        rescue_antimask_indel = ~masks_indel[4, :]
        (
            LR_raw,
            LR_max,
            indels,
            hps,
            strs,
            F1R2_ref_count,
            F1R2_alt_count,
            F2R1_ref_count,
            F2R1_alt_count,
            hp_match,
        ) = genotypeDSIndel(
            readSet,
            rs_reference_start,
            rs_reference_end,
            ref_np[start_ind:end_ind],
            # Wide (rescue) scope -- see the matching call site
            # earlier in this function for the full explanation.
            rescue_antimask_indel,
            hp_raw_np[:, start_ind:end_ind],
            str_raw_np[:, start_ind:end_ind],
            params,
        )
        # No CS-based prefilter -- every genotyped candidate indel
        # is emitted, tagged "masked" or "PASS" by the
        # antimask/LR_pass_bool checks below.
        pass_inds = list(range(len(indels)))
        if len(LR_raw) > 0:
            # ts/ti removed: any candidate with a positive LR is
            # treated as a candidate mutation, capped per-call by
            # that call's own theoretical LR ceiling (LR_max/"LM",
            # already context-specific -- see indelMaxLR).
            LR_pass_bool = LR_raw > np.minimum(params["pcutoffi"], LR_max)
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
                indel_alt = nums2str(ref_np[[indel_pos - reference_mat_start]]).upper()
            else:
                indel_ref = nums2str(ref_np[[indel_pos - reference_mat_start]]).upper()
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
            # simplex branch: this local-coverage re-check only makes sense
            # for a family that is genuinely dual-strand overall (catches
            # e.g. soft-clip-induced dropout at this specific position); for
            # a pure simplex family the absent side's count is trivially 0
            # everywhere, so skip the re-check rather than reject every
            # candidate outright.
            if F1R2 > 0 and F2R1 > 0:
                if F1R2_alt_count[pass_inds[nn]] + F1R2_ref_count[pass_inds[nn]] == 0:
                    continue
                if F2R1_alt_count[pass_inds[nn]] + F2R1_ref_count[pass_inds[nn]] == 0:
                    continue
            if indel_size > 0:
                offset = 0
            else:
                offset = -indel_size
            indel_slice_start = indel_pos - reference_mat_start - start_ind
            indel_slice = slice(indel_slice_start, indel_slice_start + offset + 1)
            # See the matching block earlier in this function for
            # the full explanation.
            if not LR_pass_bool[nn]:
                continue
            if antimask_indel[indel_slice].all():
                flt = flt_rs
            elif unmasked_antimask_indel[indel_slice].all():
                flt = "masked"
            elif params["rescue"] and rescue_antimask_indel[indel_slice].all():
                flt = _filter_reason_label(
                    masks_indel[2, indel_slice],
                    masks_indel[5, indel_slice],
                    masks_indel[3, indel_slice],
                    extra_flag=masks_indel[0, indel_slice],
                    extra_label="indel_mask",
                )
            else:
                continue

            indel_rec = {
                "chrom": chromNow,
                "pos": indel_pos + 1,
                "ref": indel_ref,
                "alt": indel_alt,
                "filter": flt,
                "infos": {
                    # Read-bundle (whole duplex family) totals, matching the
                    # documented F1R2/F2R1 INFO contract and the SBS side's
                    # own convention (call.py's SBS mut dict, above) -- NOT
                    # the locus-specific alt+ref count at this indel's own
                    # position (that's TC/BC, below). The two counts can
                    # differ when a read doesn't span this exact position or
                    # fails the post-anchor base-quality check in
                    # getIndelArr, and Estimate.py's duplex-group lookup
                    # keys strictly off F1R2+F2R1, so using the locus count
                    # here could point at a duplex group that was never
                    # recorded as having any coverage.
                    "F1R2": int(F1R2),
                    "F2R1": int(F2R1),
                    # "LR": LR[pass_inds[0]],
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
                    "HM": int(hp_match[pass_inds[nn]]),
                    "SNPM": 0,
                    "NOISEM": int(np.any(masks_indel[1, indel_slice])),
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
                hp_alt_now,
                str_alt_now,
                mismatch_dmg_now,
                hp_dmg_now,
                str_dmg_now,
                sbs_alt_bq_hist_now,
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
            hp_alt_mat += hp_alt_now
            str_alt_mat += str_alt_now
            mismatch_dmg_mat += mismatch_dmg_now
            hp_dmg_mat += hp_dmg_now
            str_dmg_mat += str_dmg_now
            sbs_alt_bq_hist_mat += sbs_alt_bq_hist_now
        if isLearn:
            return False  # learn mode: no calling, doesn't count toward duplex_count
        (
            cov_mat,
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
            np.copy(unmasked_antimask),
            # Wide (rescue) scope -- see the matching call site
            # earlier in this function for the full explanation.
            np.copy(rescue_antimask),
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
            hps_for_row = np.where(is_hp, np.minimum(repeat_count_arr, 10), 1)
            total_len_arr = unit_len_arr * repeat_count_arr
            # STR-length bin for L_indel_len's row axis (1="2-9"
            # through 4="40+") -- always >=1 wherever ~is_hp
            # (unit_len>=2 there implies a real annotated repeat by
            # construction, so there's no more "not really a
            # repeat" fallback to an hp-length row; row 0 of
            # L_indel_len is reserved for genuinely unannotated
            # positions, which this branch never represents since
            # is_hp already covers those via L_indel_1bp below).
            strs_for_row = np.where(~is_hp, 1, 0)
            strs_for_row[np.logical_and(~is_hp, total_len_arr >= 10)] = 2
            strs_for_row[np.logical_and(~is_hp, total_len_arr >= 25)] = 3
            strs_for_row[np.logical_and(~is_hp, total_len_arr >= 40)] = 4
            row_len = strs_for_row
            # depth_by_hpstr bucketing: base-specific hp length
            # (1-10, capped at 10) x A/T/C/G, plus 4 STR-length
            # buckets -- a different (finer, base-split) encoding
            # from row_len above, which stays base-agnostic and
            # indexes L_indel_len's power-table rows (left alone
            # here, still STR-priority-merged, since it feeds
            # cov_mat_indel below).
            #
            # Unlike row_len, the two depth_by_hpstr buckets below
            # are credited INDEPENDENTLY from the raw (unmerged)
            # hp_raw_np/str_raw_np arrays, not from
            # is_hp/unit_len_arr above -- see the matching call
            # site earlier in this function for the full
            # explanation (mirrors the indel100 accumulation
            # further down).
            hp_run_capped10 = np.minimum(
                hp_raw_np[0, start_ind:end_ind].astype(int), 10
            )
            hp_bucket_row = ref_allele_safe * 10 + hp_run_capped10 - 1
            str_unit_row = str_raw_np[0, start_ind:end_ind].astype(int)
            str_repeat_row = str_raw_np[1, start_ind:end_ind].astype(int)
            is_real_str_row = str_unit_row >= 2
            total_len_row = str_unit_row * str_repeat_row
            str_bin_row = np.where(is_real_str_row, 1, 0)
            str_bin_row[is_real_str_row & (total_len_row >= 10)] = 2
            str_bin_row[is_real_str_row & (total_len_row >= 25)] = 3
            str_bin_row[is_real_str_row & (total_len_row >= 40)] = 4
            str_bucket_row = 39 + str_bin_row
            unit_len_clamped = np.clip(unit_len_arr, 2, 5)
            n_top_indel = np.minimum(F1R2_count.sum(axis=0), 9).astype(int)
            n_bot_indel = np.minimum(F2R1_count.sum(axis=0), 9).astype(int)

            # Run-boundary visibility -- see the matching call site
            # earlier in this function for the full explanation of
            # why this moved up here (shared with cov_mat_indel and
            # the indel100 block below, instead of two
            # independently-approximated gates).
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

            hp_cut_bool = hp_raw_np[1, start_ind:end_ind].astype(bool)
            hp_repeat_valid = _run_boundary_valid(hp_cut_bool) & antimask_indel
            str_cut_bool = str_raw_np[2, start_ind:end_ind].astype(bool)
            str_repeat_valid = _run_boundary_valid(str_cut_bool) & antimask_indel
            unit_len_clamped_str = np.clip(str_unit_row, 2, 5)
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
                next_ref_arr[:avail] = ref_np[start_ind + 1 : start_ind + 1 + avail]
            next_ref_valid = (next_ref_arr >= 0) & (next_ref_arr <= 3)

            cov_mat_indel = np.zeros([window_len, 16])
            # Deletion/Insertion of Repeat Unit -- see the matching
            # call site earlier in this function for the full
            # explanation (0 outside a real STR, not the HP/1bp
            # fallback; that's columns 14/15's job now).
            cov_mat_indel[:, 0] = (
                np.where(
                    is_real_str_row,
                    L_indel_len[
                        n_top_indel,
                        n_bot_indel,
                        str_bin_row,
                        -unit_len_clamped_str + 5,
                    ],
                    0,
                )
                * str_repeat_valid
                * str_cut_bool
            )
            cov_mat_indel[:, 1] = (
                np.where(
                    is_real_str_row,
                    L_indel_len[
                        n_top_indel,
                        n_bot_indel,
                        str_bin_row,
                        unit_len_clamped_str + 5,
                    ],
                    0,
                )
                * str_repeat_valid
                * str_cut_bool
            )
            # del_2..del_5plus: fixed homopolymer-length-1 context
            # (microhomology-style, row 0 = hps=1), not the
            # position's actual repeat context.
            cov_mat_indel[:, 2] = L_indel_len[n_top_indel, n_bot_indel, 0, 3]
            cov_mat_indel[:, 3] = L_indel_len[n_top_indel, n_bot_indel, 0, 2]
            cov_mat_indel[:, 4] = L_indel_len[n_top_indel, n_bot_indel, 0, 1]
            cov_mat_indel[:, 5] = L_indel_len[n_top_indel, n_bot_indel, 0, 0]
            # Insertion A/T/C/G: rep0 opportunity to insert exactly
            # that base where it's genuinely novel (next reference
            # base != that base -- a real repeat extension is
            # credited under hp.txt/L_indel_1bp instead, not here)
            # -- the row-0/mismatch context, which has no base
            # axis of its own, hence the same L_indel_len[...,0,6]
            # value for all 4 columns, just gated by whichever
            # base's own next-ref match excludes it. Deliberately
            # not STR-gated -- see the matching call site earlier
            # in this function.
            for b in range(4):
                cov_mat_indel[:, 6 + b] = (
                    L_indel_len[n_top_indel, n_bot_indel, 0, 6]
                    * next_ref_valid
                    * (next_ref_arr != b)
                    * antimask_indel
                )
            # Insertion Length 2..5+: the flat, scan-free "arbitrary
            # novel U-bp insertion" opportunity -- gated only by
            # antimask_indel, not hp_repeat_valid/str_repeat_valid
            # or ~is_real_str -- see the matching call site earlier
            # in this function for the full explanation.
            cov_mat_indel[:, 10] = (
                L_indel_len[n_top_indel, n_bot_indel, 0, 7] * antimask_indel
            )
            cov_mat_indel[:, 11] = (
                L_indel_len[n_top_indel, n_bot_indel, 0, 8] * antimask_indel
            )
            cov_mat_indel[:, 12] = (
                L_indel_len[n_top_indel, n_bot_indel, 0, 9] * antimask_indel
            )
            cov_mat_indel[:, 13] = (
                L_indel_len[n_top_indel, n_bot_indel, 0, 10] * antimask_indel
            )
            # Deletion/Insertion of Homopolymer -- see the matching
            # call site earlier in this function for the full
            # explanation.
            cov_mat_indel[:, 14] = (
                L_indel_1bp[
                    n_top_indel, n_bot_indel, hp_run_capped10, ref_allele_safe, 0
                ]
                * hp_repeat_valid
                * hp_cut_bool
            )
            cov_mat_indel[:, 15] = (
                L_indel_1bp[
                    n_top_indel, n_bot_indel, hp_run_capped10, ref_allele_safe, 1
                ]
                * hp_repeat_valid
                * hp_cut_bool
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
            # hp_run_capped10, hp_cut_bool, hp_repeat_valid,
            # str_unit_row, str_repeat_row, str_cut_bool,
            # is_real_str_row, str_repeat_valid, str_bin_row,
            # unit_len_clamped_str are all computed once, further up
            # in this function -- see the matching call site earlier
            # for the full explanation.
            unit_bucket = np.clip(str_unit_row, 2, 5) - 2
            count_bucket_del = np.clip(str_repeat_row, 1, 6) - 1

            # 1-4: homopolymer deletion (cols 0-23) — always from
            # the position's own raw hp run/cut, regardless of any
            # STR annotation there too.
            hp_del_power = L_indel_1bp[
                n_top_indel, n_bot_indel, hp_run_capped10, ref_allele_safe, 0
            ]
            del_bucket = np.clip(hp_run_capped10, 1, 6) - 1
            _accumulate_indel100(
                base_cgta * 6 + del_bucket,
                hp_del_power,
                hp_repeat_valid & hp_cut_bool,
            )

            # 5-8: homopolymer insertion (cols 24-43) — length
            # bins 1..5+ only; rep0 (hp_run_capped10==0, never
            # actually occurs) is deliberately not credited here,
            # since that value is always discarded and replaced by
            # the exact 1bpins{base} next-base computation below
            # (see Estimate.py's
            # override_inshp0_with_next_base_opportunity).
            hp_ins_power = L_indel_1bp[
                n_top_indel, n_bot_indel, hp_run_capped10, ref_allele_safe, 1
            ]
            ins_bucket = np.clip(hp_run_capped10, 1, 5) - 1
            _accumulate_indel100(
                24 + base_cgta * 5 + ins_bucket,
                hp_ins_power,
                (hp_run_capped10 >= 1) & hp_repeat_valid & hp_cut_bool,
            )

            # 9: STR deletion (cols 44-67) — real part, from the
            # position's own raw STR annotation, independent of
            # any homopolymer overlap. Same value/gate as
            # cov_mat_indel[:,0] above.
            str_del_power = L_indel_len[
                n_top_indel, n_bot_indel, str_bin_row, -unit_len_clamped_str + 5
            ]
            _accumulate_indel100(
                44 + unit_bucket * 6 + count_bucket_del,
                str_del_power,
                is_real_str_row & str_repeat_valid & str_cut_bool,
            )
            # Flat rep1 opportunity, for an arbitrary, non-repeating
            # U-bp deletion — a fixed hps=1/strs=0 power lookup (row
            # 0), not an attempt to discover a true higher count
            # that might coincidentally exist there (that would
            # need rescanning the reference; see funcs/misc.py's
            # indel100_reference_bucket_indices for the matching
            # scan-free credit on the genome-wide side). Gated only
            # by antimask_indel (mapped, unmasked), not
            # hp_repeat_valid/str_repeat_valid -- see the matching
            # call site earlier in this function for the full
            # explanation. Deliberately NOT gated by
            # ~is_real_str_row either (a hypothetical, non-matching
            # U-bp deletion is a real, independent opportunity even
            # at a position that also has a real, different-unit
            # STR annotation).
            for u_idx, U in enumerate([2, 3, 4, 5]):
                hyp_del_power = L_indel_len[n_top_indel, n_bot_indel, 0, -U + 5]
                _accumulate_indel100(
                    np.full(end_ind - start_ind, 44 + u_idx * 6),
                    hyp_del_power,
                    antimask_indel,
                )

            # 10: STR insertion (cols 68-91) — real part
            str_ins_power = L_indel_len[
                n_top_indel, n_bot_indel, str_bin_row, unit_len_clamped_str + 5
            ]
            count_bucket_ins = np.clip(str_repeat_row, 0, 5)
            _accumulate_indel100(
                68 + unit_bucket * 6 + count_bucket_ins,
                str_ins_power,
                is_real_str_row & str_repeat_valid & str_cut_bool,
            )
            is_rep1 = str_repeat_row == 1
            _accumulate_indel100(
                68 + unit_bucket * 6 + 0,
                str_ins_power,
                is_real_str_row & is_rep1 & str_repeat_valid & str_cut_bool,
            )
            # Flat rep0/rep1 opportunity ("what if a U-bp unit were
            # inserted here"): fixed hps=1/strs=0 context (row 0),
            # same value added to both the "0" and "1" repeat-count
            # buckets for every unit size — scan-free, matching the
            # deletion case above, and gated only by antimask_indel
            # for the same reason (see the deletion loop above).
            for u_idx, U in enumerate([2, 3, 4, 5]):
                hyp_power = L_indel_len[n_top_indel, n_bot_indel, 0, U + 5]
                base_col = 68 + u_idx * 6
                _accumulate_indel100(
                    np.full(end_ind - start_ind, base_col),
                    hyp_power,
                    antimask_indel,
                )
                _accumulate_indel100(
                    np.full(end_ind - start_ind, base_col + 1),
                    hyp_power,
                    antimask_indel,
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
        ref_int = ref_np[start_ind:end_ind]
        n_win = ref_int.size
        # No CS-based prefilter -- every genotyped candidate
        # mismatch is emitted, tagged "masked" or "PASS" by the
        # antimask/LR_pass_bool checks below.
        mut_pos_in_win = np.nonzero(mut_mask)[0]
        muts_ind_compressed = np.arange(mut_pos_in_win.size)
        muts_ind = mut_pos_in_win.tolist()
        refs_ind = np.nonzero(np.logical_and(unmasked_antimask, b1_int == ref_int))[
            0
        ].tolist()
        LR_pass_bool = np.zeros(n_win, dtype=bool)
        mut_thresholds = params["pcutoff_sbs"][
            trinuc_np[start_ind:end_ind][mut_mask],
            b1_int[mut_mask],
        ]
        # ts removed: any candidate with a positive LR is treated
        # as a candidate mutation.
        LR_pass_bool[mut_mask] = LR_raw_mut > mut_thresholds
        alt_int = b1_int
        unmasked_pass_bool = np.full(n_win, False, dtype=bool)
        unmasked_pass_bool[refs_ind] = True
        unmasked_pass_bool[muts_ind] = True
        pass_bool = np.copy(unmasked_pass_bool)
        pass_bool[~antimask] = False
        pos = [mut_ind + start_ind + reference_mat_start for mut_ind in muts_ind]
        mut_positions = [
            mut_ind + start_ind + reference_mat_start + 1 for mut_ind in muts_ind
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
                    abs(readSet[0].template_length) - max_ref_len + muts_ind[nn] + 1,
                )
                readPos3p = min(
                    abs(max_ref_len) - muts_ind[nn],
                    muts_ind[nn] + 1,
                )
            FPs.append(readPos5p)
            RPs.append(readPos3p)
            # simplex branch: see matching comment at the indel re-check above.
            if F1R2 > 0 and F2R1 > 0:
                if F1R2_count[:, muts_ind[nn]].sum() == 0:
                    continue
                if F2R1_count[:, muts_ind[nn]].sum() == 0:
                    continue
            # See the matching block earlier in this function for
            # the full explanation.
            if not LR_pass_bool[muts_ind[nn]]:
                continue
            if pass_bool[muts_ind[nn]]:
                flt = flt_rs
            elif masks[0, muts_ind[nn]] or masks[1, muts_ind[nn]]:
                # snp_mask/noise_mask specifically (masks rows 0/1, same as
                # SNPM/NOISEM below) -- was `unmasked_pass_bool[muts_ind[nn]]`,
                # which is unconditionally True for every muts_ind entry, so it
                # mislabeled any ncov/nm/trim/include block as "masked" too and
                # made the rescue branch below unreachable for those cases.
                flt = "masked"
            elif params["rescue"] and rescue_antimask[muts_ind[nn]]:
                flt = _filter_reason_label(
                    masks[2, muts_ind[nn]],
                    masks[4, muts_ind[nn]],
                    masks[5, muts_ind[nn]],
                )
            else:
                continue
            mut = {
                "chrom": mut_chrom,
                "pos": mut_pos,
                "ref": mut_ref,
                "alt": mut_alt,
                "filter": flt,
                "infos": {
                    "F1R2": F1R2,
                    "F2R1": F2R1,
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
                    "SNPM": int(masks[0, muts_ind[nn]]),
                    "NOISEM": int(masks[1, muts_ind[nn]]),
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
            if not isLearn:
                coverage_indel_cat[start_ind:end_ind][antimask_indel] += cov_mat_indel[
                    antimask_indel
                ]
                unmasked_coverage_indel_cat[start_ind:end_ind][
                    unmasked_antimask_indel
                ] += cov_mat_indel[unmasked_antimask_indel]
                duplex_read_num_dict_indel[duplex_no] += indel100
                # Read-family depth-composition totals (see the
                # matching block earlier in this function for the
                # full explanation).
                trinuc_ctx = trinuc_np[start_ind:end_ind]
                _accumulate_depth_matrix(
                    depth_by_trinuc,
                    n_top_indel,
                    n_bot_indel,
                    trinuc_ctx,
                    unmasked_pass_bool & (trinuc_ctx < 64),
                    64,
                )
                _accumulate_depth_matrix(
                    depth_by_hpstr,
                    n_top_indel,
                    n_bot_indel,
                    hp_bucket_row,
                    unmasked_pass_bool,
                    44,
                )
                _accumulate_depth_matrix(
                    depth_by_hpstr,
                    n_top_indel,
                    n_bot_indel,
                    str_bucket_row,
                    unmasked_pass_bool & is_real_str_row,
                    44,
                )
                dinuc_ctx = ref_allele_safe * 4 + np.where(
                    next_ref_arr >= 0, next_ref_arr, 0
                )
                # A dinuc spans two positions, so a mask at either
                # one has to exclude it: the dinuc "at" a masked
                # position (it's the first base) and the dinuc
                # "before" it (it's the second base) both need to
                # be dropped, not just the former.
                next_unmasked_pass_bool = np.zeros_like(unmasked_pass_bool)
                next_unmasked_pass_bool[:-1] = unmasked_pass_bool[1:]
                _accumulate_depth_matrix(
                    depth_by_dinuc,
                    n_top_indel,
                    n_bot_indel,
                    dinuc_ctx,
                    unmasked_pass_bool
                    & next_unmasked_pass_bool
                    & (ref_allele_arr <= 3)
                    & next_ref_valid,
                    16,
                )

    return bool(pass_bool.any())


def _collect_call_barcode(call_barcodes, key, mut):
    """Record one more founding duplex family's (TAG1, TAG2) barcode pair
    as supporting this candidate -- shared by the SNV/indel/DBS depth-
    extraction loops below, each of which needs this same accounting so
    extractDepthBatchSnv/Indel/Dbs's call_barcodes arg can exempt a
    founding-family read from the minBq filter regardless of base
    quality."""
    call_barcodes.setdefault(key, set()).add(
        (mut["infos"]["TAG1"], mut["infos"]["TAG2"])
    )


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
    hp_alt_mat = np.zeros([10, 12])
    str_alt_mat = np.zeros([5, 11])
    mismatch_dmg_mat = np.zeros([64, 4])
    hp_dmg_mat = np.zeros([10, 12])
    str_dmg_mat = np.zeros([5, 11])
    # Amp-error BQ histogram accumulated across families (SBS only) --
    # see funcs/learn.py's profileTriNucMismatches for the (channel, BQ)
    # shape convention and NUM_BQ, and estimate_sbs_srd_rates for how it
    # feeds the SRD rate matrix.
    sbs_alt_bq_hist_mat = np.zeros([64, 4, NUM_BQ])
    # Loads/regularizes ampmat/ampmat_rev/dmgmat_top/rev_top/bot/rev_bot/
    # ampmat_indel/ampmat_indel_rev/dmgmat_indel/dmgmat_indel_rev and sets
    # trinuc_convert/trinuc2num_dict/num2trinuc_list on params -- see
    # load_error_matrices's docstring for why this is a standalone,
    # importable function rather than inlined here.
    load_error_matrices(params)
    trinuc2num = params["trinuc2num_dict"]
    num2trinuc = params["num2trinuc_list"]
    # Initialize

    total_coverage = np.zeros(4)
    total_coverage_indel_cat = np.zeros(16)
    total_unmasked_coverage = np.zeros(4)
    total_unmasked_coverage_indel_cat = np.zeros(16)
    # Read-family depth-composition totals: depth_mat[top_count, bot_count,
    # category] = number of passing loci with that exact (F1R2, F2R1) count
    # composition (each capped at 9) in that context bucket. Unlike the old
    # per-locus dict keyed by (chrom, pos), these are pure aggregate counts
    # -- no locus identity is kept, so there is nothing to carry across
    # window boundaries and no cross-process boundary merging to do; see
    # _accumulate_depth_matrix.
    depth_by_trinuc = np.zeros((10, 10, 64), dtype=np.int64)  # 64 trinuc contexts
    depth_by_hpstr = np.zeros(
        (10, 10, 44), dtype=np.int64
    )  # 40 hp (A/T/C/G x len 1-10, capped) + 4 str buckets
    depth_by_dinuc = np.zeros((10, 10, 16), dtype=np.int64)  # 16 dinuc contexts
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
    # Highest position already written to locus_bed for reference_mat_chrom
    # (reset on chromosome change, below). A backward window re-init (see
    # the rs_reference_start < reference_mat_start branch further down)
    # can revisit positions whose coverage was already flushed under the
    # window active before the re-init -- writing them again would corrupt
    # the bed file's required sort order and double-count their coverage,
    # so any flush only ever emits strictly-increasing positions and
    # silently drops the (already-accounted-for) rest.
    max_flushed_pos = -1
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
            bam, regions, params.get("reference"), params.get("region_file")
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
    # Only ever populated (via simulate_power_grid, below) in coverage_only
    # mode -- round-1/genotyping mode leaves this all-zero, which is fine
    # since genotypeDSSnv treats it exactly like L=None (cov_mat stays 0,
    # see prob.py) and round 1 doesn't write coverage.bed at all.
    coverage_only = params.get("coverage_only", False)
    L = np.zeros([10, 10, 64, 4])

    # Sample base quality distribution from the BAM regions -- only needed
    # to feed simulate_power_grid below, so skip entirely outside
    # coverage_only mode.
    all_quals = []
    if coverage_only:
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
    N_SIM = 100
    # Each simulate_power_grid call below derives its own RNG fresh from
    # (params["seed"], threshold) via threshold_rng instead of sharing one
    # live generator across every context/threshold this call processes --
    # see threshold_rng's docstring for why a shared, order-dependent RNG
    # broke reproducibility across different -p thread counts.
    base_seed_L = params["seed"]

    # Per-(trinuc context, alt base) SBS calling threshold: max(pcutoff,
    # maxLR), where maxLR = log10(1-Pdmg_t) + log10(1-Pdmg_b)
    # - log10(Pdmg_rev_t) - log10(Pdmg_rev_b) is the theoretical ceiling of
    # LR_masked for that SBS type (same formula as genotypeDSSnv's
    # LR_max/"LM" INFO field, here evaluated once per context instead of
    # once per called family, since Pdmg_t/Pdmg_rev_t/Pdmg_b/Pdmg_rev_b
    # only depend on context, never on read depth). Contexts t and
    # (t+32)%64 with alt bases b and b^1 are always exact reverse
    # complements (see build_trinuc64_order), and dmgmat_bot/
    # dmgmat_rev_bot are themselves built from dmgmat_top/dmgmat_rev_top
    # via that same reverse-complement mapping, so maxLR at an RC pair is
    # always identical -- only the canonical (t<32) half needs computing;
    # the other half is filled in by mirroring.
    pcutoff_sbs = np.full((64, 4), pcut, dtype=float)
    for t_ctx in range(32):
        ref_base_idx_ctx = base2num[num2trinuc[t_ctx][1]]
        for b_ctx in range(4):
            if b_ctx == ref_base_idx_ctx:
                continue
            tc_ctx = int(trinuc_conv_np_L[t_ctx, b_ctx])
            Pdmg_t_ctx = prob_dmg_t_L[tc_ctx, ref_base_idx_ctx]
            Pdmg_rev_t_ctx = prob_dmg_rev_t_L[tc_ctx, ref_base_idx_ctx]
            Pdmg_b_ctx = prob_dmg_b_L[tc_ctx, ref_base_idx_ctx]
            Pdmg_rev_b_ctx = prob_dmg_rev_b_L[tc_ctx, ref_base_idx_ctx]
            if Pdmg_t_ctx == 0:
                Pdmg_t_ctx = 1e-9
            if Pdmg_rev_t_ctx == 0:
                Pdmg_rev_t_ctx = 1e-9
            if Pdmg_b_ctx == 0:
                Pdmg_b_ctx = 1e-9
            if Pdmg_rev_b_ctx == 0:
                Pdmg_rev_b_ctx = 1e-9
            maxLR_ctx = (
                np.log10(1 - Pdmg_t_ctx)
                + np.log10(1 - Pdmg_b_ctx)
                - np.log10(Pdmg_rev_t_ctx)
                - np.log10(Pdmg_rev_b_ctx)
            )
            chosen_threshold = min(pcut, maxLR_ctx)
            pcutoff_sbs[t_ctx, b_ctx] = chosen_threshold
            pcutoff_sbs[t_ctx + 32, b_ctx ^ 1] = chosen_threshold
    # Per-channel FDR-driven threshold overrides (Caller.py's post-hoc
    # per-channel refinement), keyed by raw (trinuc, alt) index pairs --
    # applied after the maxLR ceiling above so a rerun restricted to a
    # handful of adjusted channels only touches L at those (t, b) cells,
    # leaving every other cell's power table exactly as the full run built
    # it. Absent (the normal, non-rerun path) for every other caller.
    for (t_ov, b_ov), new_threshold in params.get("pcutoff_sbs_override", {}).items():
        pcutoff_sbs[t_ov, b_ov] = new_threshold
    params["pcutoff_sbs"] = pcutoff_sbs

    if coverage_only:
        for t in range(64):
            ref_base_idx_L = base2num[num2trinuc[t][1]]
            for b in range(4):
                if b == ref_base_idx_L:
                    # A base "changing" to itself isn't a mutation
                    # opportunity; leave L[:, :, t, ref_base_idx_L] at its
                    # zero-initialized value instead of simulating a
                    # self-to-self detection power.
                    continue
                tc = int(trinuc_conv_np_L[t, b])
                L[:, :, t, b] = simulate_power_grid(
                    prob_amp_mat_L[tc, ref_base_idx_L],
                    prob_amp_mat_rev_L[tc, ref_base_idx_L],
                    prob_dmg_t_L[tc, ref_base_idx_L],
                    prob_dmg_rev_t_L[tc, ref_base_idx_L],
                    prob_dmg_b_L[tc, ref_base_idx_L],
                    prob_dmg_rev_b_L[tc, ref_base_idx_L],
                    pcutoff_sbs[t, b],
                    all_quals,
                    threshold_rng(base_seed_L, pcutoff_sbs[t, b]),
                    N_SIM,
                )

    # Build indel detection-power tables, analogous to L above but for indel
    # calling. genotypeDSIndel classifies each candidate indel by repeat
    # context (hps: homopolymer run length 0-10, capped; strs: STR length
    # bin 0-3)
    # and indel length (idLen), then selects Pamp/Pdmg via indelErrorProbs.
    # These tables simulate that same selection across every (top count,
    # bottom count) combination and precompute the fraction of simulations
    # that would clear indel calling's per-context pcutoffi threshold, so
    # per-position indel coverage can be looked up by depth+context instead
    # of re-simulated. Only built outside learn mode: ampmat_indel/
    # dmgmat_indel (and pcutoffi-driven calling) aren't used in learn mode.
    # Defaults for isLearn mode, where neither table is built below --
    # kept as None (rather than omitted) so callBam's return signature is
    # identical in both modes.
    L_indel_1bp = None
    L_indel_len = None
    if not isLearn:
        prob_amp_hp_L = params["ampmat_hp"]
        prob_dmg_hp_L = params["dmgmat_hp"]
        prob_amp_str_L = params["ampmat_str"]
        prob_dmg_str_L = params["dmgmat_str"]
        pcutoffi_L = params["pcutoffi"]

        def indelContextThreshold(Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c):
            # Uniform default cutoff, capped per-context by the theoretical
            # ceiling of masked LR for that context (same formula as
            # genotypeDSSnv's LR_max/"LM" field -- see indelMaxLR), mirroring
            # pcutoff_sbs's min(pcut, maxLR_ctx) treatment above.
            return min(
                pcutoffi_L,
                float(indelMaxLR(Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c)),
            )

        def simulateIndelPowerGrid(
            Pamp_c,
            Pamp_rev_c,
            Pdmg_c,
            Pdmg_rev_c,
            Pdmg_bot_c,
            Pdmg_rev_bot_c,
            threshold,
        ):
            return simulate_power_grid(
                Pamp_c,
                Pamp_rev_c,
                Pdmg_c,
                Pdmg_rev_c,
                Pdmg_bot_c,
                Pdmg_rev_bot_c,
                threshold,
                all_quals,
                threshold_rng(base_seed_L, threshold),
                N_SIM,
            )

        # 1bp indels (deletion/insertion of a homopolymer's repeat unit):
        # axes are top count, bottom count, hps (0-10, capped), ref_allele
        # (0-3), sign (0=deletion idLen=-1, 1=insertion idLen=+1).
        # inserted_base is passed equal to ref_allele_c throughout -- this
        # grid represents the true same-base homopolymer-extension
        # opportunity at a position with this (hps, ref_allele) context
        # (what genotypeDSIndel/depth_by_hpstr's base axis actually means),
        # not an arbitrary mismatched-base scenario. A mismatched 1bp
        # insertion (str.txt row 0) has no per-(hps, base) coverage grid of
        # its own -- known gap, see funcs/call.py's indelErrorProbs
        # docstring and Caller.py's channel job-building for the same
        # scoping note; such candidates still get scored correctly by
        # genotypeDSIndel, they just don't get an FDR-refined threshold.
        L_indel_1bp = np.zeros([10, 10, 11, 4, 2])
        indel_1bp_threshold_override = params.get("indel_1bp_threshold_override", {})
        if coverage_only:
            for hps_c in range(11):
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
                            ref_allele_c,
                            prob_amp_hp_L,
                            prob_dmg_hp_L,
                            prob_amp_str_L,
                            prob_dmg_str_L,
                        )
                        # Per-(hp length, sign, base pool) FDR-driven
                        # override, same idea as pcutoff_sbs_override
                        # above; falls back to the per-context threshold
                        # when absent. Pool ("C"=C/G, "T"=A/T) matches
                        # Caller.py's raw_lr_hp accumulation -- base2num
                        # order A,T,C,G, so ref_allele_c in (2,3) is "C".
                        pool_c = "C" if ref_allele_c in (2, 3) else "T"
                        threshold = indel_1bp_threshold_override.get(
                            (hps_c, sign_idx, pool_c),
                            indelContextThreshold(
                                Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c
                            ),
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

        # Length-bin indels (raw length 2/3/4/5+, either direction): axes
        # are top count, bottom count, STR-length-bin row (0="0-1"/not a
        # real repeat through 4="40+"), idLen+5. Always STR-context -- a
        # length>=2 indel outside an annotated repeat is row 0, so this is
        # a single loop over all 5 str.txt rows.
        L_indel_len = np.zeros([10, 10, 5, 11])
        len_idLens = (-5, -4, -3, -2, 2, 3, 4, 5)
        indel_len_threshold_override = params.get("indel_len_threshold_override", {})
        if coverage_only:
            for strs_c in range(5):
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
                        0,
                        prob_amp_hp_L,
                        prob_dmg_hp_L,
                        prob_amp_str_L,
                        prob_dmg_str_L,
                    )
                    threshold = indel_len_threshold_override.get(
                        (strs_c, idLen_c),
                        indelContextThreshold(
                            Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c
                        ),
                    )
                    L_indel_len[:, :, strs_c, idLen_c + 5] = simulateIndelPowerGrid(
                        Pamp_c,
                        Pamp_rev_c,
                        Pdmg_c,
                        Pdmg_rev_c,
                        Pdmg_bot_c,
                        Pdmg_rev_bot_c,
                        threshold,
                    )
            # Row 0 ("not a real repeat") at idLen=+1 specifically: the
            # mismatched-insertion background scenario (str.txt row 0's
            # only real column), needed below for the "Insertion A/T/C/G"
            # coverage columns -- a *novel* insertion (next reference base
            # != the inserted base, i.e. not a repeat extension; see that
            # code's own comment), which is exactly the row-0/mismatch
            # context, not the hp.txt/matching one L_indel_1bp represents.
            # Not part of len_idLens above since every other +-1bp
            # scenario is HP-context, handled by L_indel_1bp instead.
            (
                Pamp_c,
                Pamp_rev_c,
                Pdmg_c,
                Pdmg_rev_c,
                Pdmg_bot_c,
                Pdmg_rev_bot_c,
            ) = indelErrorProbs(
                1,
                0,
                1,
                0,
                1,
                prob_amp_hp_L,
                prob_dmg_hp_L,
                prob_amp_str_L,
                prob_dmg_str_L,
            )
            threshold = indel_len_threshold_override.get(
                (0, 1),
                indelContextThreshold(Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c),
            )
            L_indel_len[:, :, 0, 6] = simulateIndelPowerGrid(
                Pamp_c,
                Pamp_rev_c,
                Pdmg_c,
                Pdmg_rev_c,
                Pdmg_bot_c,
                Pdmg_rev_bot_c,
                threshold,
            )

    # Families are keyed by label and require an exact reference_start match.
    currentReadDict = {}
    # rugged_reads_index/rugged_reads_pool: reroute a downstream mate whose
    # own alignment start disagrees with the rest of its anchor family's
    # next_reference_start into that family's majority-position group.
    rugged_reads_index = {}
    rugged_reads_pool = {}
    # Chromosome of the entries currently held in rugged_reads_index/
    # rugged_reads_pool, so they can be dropped in bulk once we move past
    # it (see the purge below) instead of leaking for the rest of the
    # process's lifetime.
    rugged_pool_chrom = None
    for rec, region in bamIterateMultipleRegion(
        bam, regions, params.get("reference"), params.get("region_file")
    ):
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
        if rec.template_length < 0 and rec.query_name in rugged_reads_index:
            rugged_reads_pool[rugged_reads_index.pop(rec.query_name)].append(rec)
            continue
        start = rec.reference_start
        label = _compute_read_label(rec, params)
        chrom = tumorBam.get_reference_name(rec.reference_id)
        if chrom != rugged_pool_chrom:
            # BAM is coordinate-sorted, so once we've moved past a
            # chromosome any redirect still sitting in rugged_reads_pool/
            # rugged_reads_index for it (its downstream mate was filtered
            # out or never turned up in this process's region) can never
            # be matched -- drop it now rather than hold onto it, and
            # every prior chromosome's leftovers, for the rest of the
            # process.
            rugged_reads_pool.clear()
            rugged_reads_index.clear()
            rugged_pool_chrom = chrom
        if currentStart == -1:
            currentStart = start
            _drain_rugged_pool(
                chrom, currentStart, rugged_reads_pool, params, currentReadDict
            )
        if start == currentStart:
            _place_read_in_dict(rec, label, currentReadDict)
            # print(currentStart,start)
        else:
            """
            Calling block starts
            """
            _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool)
            # _drain_rugged_pool (above, at the previous batch transition) can
            # merge a read into a family at this batch's currentStart even
            # though the read's own true reference_start is a few bp earlier
            # (that's the whole point -- reconciling a minority mate whose
            # alignment disagrees slightly with the rest of its family). So
            # the first key that ends up establishing/extending the
            # reference window below isn't guaranteed to have the earliest
            # start of any key still to come in this same batch, and a
            # later key's start_ind can go negative. Floor the window at the
            # batch's true minimum read start to rule that out up front.
            batch_min_start = min(
                r.reference_start
                for entry in currentReadDict.values()
                for r in entry["seqs"]
            )
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
                # simplex branch: round-0 error learning (isLearn=True) keeps
                # the strict duplex requirement unchanged (profileTriNucMismatches'
                # amp/damage decomposition is fundamentally a strand-discordance
                # comparison and has no meaning for a single-strand family);
                # only the calling rounds (isLearn=False) admit simplex families.
                family_has_evidence = (
                    (F2R1 >= 1 and F1R2 >= 1) if isLearn else (F2R1 >= 1 or F1R2 >= 1)
                )
                if family_has_evidence:
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
                    contig_len_now = tumorBam.get_reference_length(chromNow)
                    if rs_reference_start >= contig_len_now:
                        # A read set can't be placed in any window bounded by
                        # its own chromosome's length -- skip it rather than
                        # building a zero/negative-length reference window
                        # that crashes the downstream mask slicing.
                        print(
                            f"WARNING: skipping read set at {chromNow}:"
                            f"{rs_reference_start}-{rs_reference_end}, which "
                            f"starts at or beyond the contig length "
                            f"({contig_len_now})"
                        )
                        continue
                    if (
                        chromNow != reference_mat_chrom
                        or rs_reference_end > reference_mat_end
                        or rs_reference_start < reference_mat_start
                    ):
                        ### Output coverage
                        # Original per-locus multi-column coverage write-out
                        # (per-subprocess only -- Caller.py does not merge
                        # these across process boundaries yet).
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
                                unmasked_coverage_indel_cat_leftover = np.zeros((1, 16))
                                coverage_leftover = np.zeros((1, 4))
                                coverage_indel_cat_leftover = np.zeros((1, 16))
                            if (
                                chromNow == reference_mat_chrom
                                and rs_reference_start >= reference_mat_start
                            ):
                                # Forward progress within the same
                                # chromosome: the new window starts inside
                                # the old window's range, so the old
                                # window's tail is still relevant and gets
                                # carried forward as "leftover" rather than
                                # flushed immediately. A *backward* trigger
                                # (rs_reference_start < reference_mat_start)
                                # can't reuse this -- the new window's
                                # positions aren't a re-based slice of the
                                # old coverage array -- so it falls through
                                # to the full-flush branch below instead,
                                # same as a chromosome change.
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
                                if current_pos <= max_flushed_pos:
                                    # Already flushed under the window
                                    # active before a backward re-init (see
                                    # the max_flushed_pos comment above) --
                                    # skip to avoid an out-of-order bed
                                    # line and double-counting.
                                    continue
                                max_flushed_pos = current_pos
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
                        if chromNow != reference_mat_chrom:
                            max_flushed_pos = -1
                        reference_mat_chrom = chromNow
                        # current_reference = str(fasta[reference_mat_chrom].seq)
                        # batch_min_start <= rs_reference_start always (it's
                        # the min over every key in this batch, of which
                        # this key's readSet is one) -- use it so a later
                        # key in this same batch can never start before the
                        # window.
                        reference_mat_start = batch_min_start
                        try:
                            region_end = region[2]
                        except IndexError:
                            region_end = 10e10
                        contig_end = contig_len_now
                        # The +1000000 cap must be measured from
                        # reference_mat_start (batch_min_start), not from
                        # this key's own rs_reference_start -- coverage/
                        # coverage_indel_cat/etc. below are fixed at
                        # 1,000,000 rows, so basing the cap on the
                        # unfloored rs_reference_start could let the
                        # window (reference_mat_end - reference_mat_start)
                        # exceed 1,000,000 whenever batch_min_start pulled
                        # the start back, silently truncating those
                        # arrays relative to ref_np/trinuc_np/masks (which
                        # aren't capped at 1,000,000).
                        reference_mat_end = min(
                            reference_mat_start + 1000000,
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
                        coverage_indel_cat = np.zeros((1000000, 16))
                        unmasked_coverage = np.zeros((1000000, 4))
                        unmasked_coverage_indel_cat = np.zeros((1000000, 16))
                    if _process_duplex_family(
                        F1R2,
                        F2R1,
                        FPs,
                        L,
                        L_indel_1bp,
                        L_indel_len,
                        RPs,
                        chromNow,
                        coverage,
                        coverage_indel_cat,
                        currentStart,
                        depth_by_dinuc,
                        depth_by_hpstr,
                        depth_by_trinuc,
                        duplex_no,
                        duplex_read_num_dict,
                        duplex_read_num_dict_dbs,
                        duplex_read_num_dict_indel,
                        duplex_read_num_dict_trinuc,
                        flt_rs,
                        hp_alt_mat,
                        hp_dmg_mat,
                        hp_np,
                        hp_raw_np,
                        include_mask,
                        indel_dict,
                        indel_mask,
                        isLearn,
                        mismatch_dmg_mat,
                        mismatch_mat,
                        muts,
                        muts_dbs,
                        muts_dict,
                        muts_indels,
                        n_cov_mask,
                        nm_mask,
                        noise_mask,
                        num2base,
                        num2trinuc,
                        params,
                        processed_read_names,
                        readSet,
                        ref_np,
                        reference_mat_chrom,
                        reference_mat_start,
                        rs_reference_end,
                        rs_reference_start,
                        sbs_alt_bq_hist_mat,
                        setBc,
                        snp_mask,
                        str_alt_mat,
                        str_dmg_mat,
                        str_raw_np,
                        trinuc_np,
                        unmasked_coverage,
                        unmasked_coverage_indel_cat,
                    ):
                        duplex_count += 1
            """
            Calling block ends
            """
            currentReadDict = {}
            _place_read_in_dict(rec, label, currentReadDict)
            currentStart = start
            _drain_rugged_pool(
                chrom, currentStart, rugged_reads_pool, params, currentReadDict
            )
    """
    Calling block starts
    """
    _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool)
    # See the matching comment on the equivalent block earlier in this
    # function: floor the window at this final batch's true minimum read
    # start too, for the same reason. currentReadDict can be empty here --
    # unlike the equivalent mid-stream flush, this final flush runs
    # unconditionally even when this process's assigned region (e.g. a
    # --rescue run with a sparse -R bed covering only a few chromosomes)
    # never accumulated a single read -- so guard the min() rather than
    # crash on an empty iterable; batch_min_start is only read inside the
    # loop below, which is a no-op when currentReadDict is empty.
    batch_min_start = (
        min(
            r.reference_start
            for entry in currentReadDict.values()
            for r in entry["seqs"]
        )
        if currentReadDict
        else None
    )

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
        # simplex branch: see matching comment at the mid-batch flush site above.
        family_has_evidence = (
            (F2R1 >= 1 and F1R2 >= 1) if isLearn else (F2R1 >= 1 or F1R2 >= 1)
        )
        if family_has_evidence:
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
            contig_len_now = tumorBam.get_reference_length(chromNow)
            if rs_reference_start >= contig_len_now:
                # See the matching guard earlier in this function -- skip a
                # read set that starts at or beyond its own chromosome's
                # length instead of building a degenerate reference window.
                print(
                    f"WARNING: skipping read set at {chromNow}:"
                    f"{rs_reference_start}-{rs_reference_end}, which "
                    f"starts at or beyond the contig length "
                    f"({contig_len_now})"
                )
                continue
            if (
                chromNow != reference_mat_chrom
                or rs_reference_end > reference_mat_end
                or rs_reference_start < reference_mat_start
            ):
                ### Output coverage
                # Original per-locus multi-column coverage write-out
                # (per-subprocess only -- see the matching block earlier in
                # this function).
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
                        unmasked_coverage_indel_cat_leftover = np.zeros((1, 16))
                        coverage_leftover = np.zeros((1, 4))
                        coverage_indel_cat_leftover = np.zeros((1, 16))
                    if (
                        chromNow == reference_mat_chrom
                        and rs_reference_start >= reference_mat_start
                    ):
                        # See the matching comment on the equivalent block
                        # earlier in this function.
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
                        if current_pos <= max_flushed_pos:
                            # See the matching comment on the equivalent
                            # block earlier in this function.
                            continue
                        max_flushed_pos = current_pos
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
                if chromNow != reference_mat_chrom:
                    max_flushed_pos = -1
                reference_mat_chrom = chromNow
                # current_reference = str(fasta[reference_mat_chrom].seq)
                reference_mat_start = batch_min_start
                try:
                    region_end = region[2]
                except IndexError:
                    region_end = 10e10
                contig_end = contig_len_now
                # See the matching comment on the equivalent block earlier
                # in this function.
                reference_mat_end = min(
                    reference_mat_start + 1000000,
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
                coverage_indel_cat = np.zeros((1000000, 16))
                unmasked_coverage = np.zeros((1000000, 4))
                unmasked_coverage_indel_cat = np.zeros((1000000, 16))
            if _process_duplex_family(
                F1R2,
                F2R1,
                FPs,
                L,
                L_indel_1bp,
                L_indel_len,
                RPs,
                chromNow,
                coverage,
                coverage_indel_cat,
                currentStart,
                depth_by_dinuc,
                depth_by_hpstr,
                depth_by_trinuc,
                duplex_no,
                duplex_read_num_dict,
                duplex_read_num_dict_dbs,
                duplex_read_num_dict_indel,
                duplex_read_num_dict_trinuc,
                flt_rs,
                hp_alt_mat,
                hp_dmg_mat,
                hp_np,
                hp_raw_np,
                include_mask,
                indel_dict,
                indel_mask,
                isLearn,
                mismatch_dmg_mat,
                mismatch_mat,
                muts,
                muts_dbs,
                muts_dict,
                muts_indels,
                n_cov_mask,
                nm_mask,
                noise_mask,
                num2base,
                num2trinuc,
                params,
                processed_read_names,
                readSet,
                ref_np,
                reference_mat_chrom,
                reference_mat_start,
                rs_reference_end,
                rs_reference_start,
                sbs_alt_bq_hist_mat,
                setBc,
                snp_mask,
                str_alt_mat,
                str_dmg_mat,
                str_raw_np,
                trinuc_np,
                unmasked_coverage,
                unmasked_coverage_indel_cat,
            ):
                duplex_count += 1
    """
    Calling block ends
    """
    if isLearn:
        return (
            mismatch_mat,
            hp_alt_mat,
            str_alt_mat,
            mismatch_dmg_mat,
            hp_dmg_mat,
            str_dmg_mat,
            sbs_alt_bq_hist_mat,
        )

    ### SNV depth extraction, batched: gather every eligible candidate
    ### first, then resolve the tumor BAM and each normal BAM with a
    ### handful of sequential pileup() scans (extractDepthBatchSnv) instead
    ### of one random-access pileup() call per candidate.
    ###
    ### Eligibility differs by round: round 1 (not coverage_only) only
    ### depth-extracts true PASS candidates -- masked ones (and any
    ### rescue-reason label) get no depth-extraction here, matching this
    ### round's default-threshold classification. Round 2 (coverage_only)
    ### instead only depth-extracts the deferred set Caller.py identified in
    ### its re-filter step: candidates masked in round 1 whose raw LR now
    ### clears the channel's final refined threshold. PASS candidates from
    ### round 1 (even ones later downgraded to "underpowered") already have
    ### real depth from round 1 and are never re-extracted here.
    deferred_depth_keys = params.get("deferred_depth_keys", frozenset())
    mut_dict = dict()
    snv_candidate_keys = []
    # (TAG1, TAG2) of every duplex family that itself supports each
    # candidate -- tumor depth extraction below counts a primary read
    # toward a candidate's depth regardless of base quality when its own
    # duplex barcode pair matches one of these (either orientation), since
    # it's one of the founding reads of the call being verified. A
    # candidate can be supported by more than one duplex family, so this
    # collects every one seen among this candidate's own eligible mut
    # records, not just the first.
    call_barcodes = {}
    for mut in muts:
        key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
        if coverage_only:
            if key not in deferred_depth_keys:
                continue
        elif mut["filter"] != "PASS":
            continue
        _collect_call_barcode(call_barcodes, key, mut)
        if key not in mut_dict:
            mut_dict[key] = None
            snv_candidate_keys.append(key)

    if snv_candidate_keys:
        _snv_depth_t0 = time.time()
        tumor_snv_depths = extractDepthBatchSnv(
            tumorBam,
            snv_candidate_keys,
            params,
            minbq=params["minBq"],
            call_barcodes=call_barcodes,
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
        print(
            f"Process {processNo}: SNV depth extraction for "
            f"{len(snv_candidate_keys)} candidates took "
            f"{(time.time() - _snv_depth_t0) / 60:.2f} minutes"
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
        key = (chrom, pos, ref, alt)
        eligible = key in deferred_depth_keys if coverage_only else flt == "PASS"
        if eligible:
            ta, tr, ti, tdp, na, nr, ni, ndp = mut_dict[key]
        else:
            ta, tr, ti, tdp, na, nr, ni, ndp = (0, 0, 0, 0, 0, 0, 0, 0)
        # if window_filter:
        # continue
        # if ta > params["maxAltCount"]:
        # continue
        # Real depth was extracted (eligible), but it fails a post-hoc
        # sanity check PASS candidates are also held to. Rather than
        # vanishing silently here, this is always tagged with which check
        # it failed and kept (with real AC/RC/DP) -- Caller.py prunes these
        # reject-reason records out after round 2's merge, unless --rescue
        # is on (same "junk mutations that may carry real biological
        # signal" opt-in as every other rescue reason). Not dropped here
        # directly: for round 2 (coverage_only) this record is a
        # region-local copy that only feeds Caller.py's merge-by-position,
        # never the final mutsAll/indelsAll directly, so a `continue` here
        # would be invisible to the real keep/drop decision.
        if eligible:
            if ta == 0:
                mut["filter"] = "no_good_alt_read"
            elif ta / tdp > params["maxAF"]:
                mut["filter"] = "duplex_vaf"
            # if ti >= 1:
            # continue
            elif normalBams:
                if ndp < params["minNdepth"]:
                    mut["filter"] = "n_cov_mask"
                elif na / ndp > params["normalVAF"]:
                    mut["filter"] = "normal_vaf"
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        mut_pass_filter.append(mut)

    ### Indel depth extraction, batched the same way, with the same
    ### round-1-PASS-only / round-2-deferred-set eligibility split as SNVs
    ### above.
    muts_indels_dict = dict()
    indel_candidate_keys = []
    # Same founding-family barcode collection as the SNV loop above.
    indel_call_barcodes = {}
    for mut in muts_indels:
        key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
        if coverage_only:
            if key not in deferred_depth_keys:
                continue
        elif mut["filter"] != "PASS":
            continue
        _collect_call_barcode(indel_call_barcodes, key, mut)
        if key not in muts_indels_dict:
            muts_indels_dict[key] = None
            indel_candidate_keys.append(key)

    if indel_candidate_keys:
        _indel_depth_t0 = time.time()
        tumor_indel_depths = extractDepthBatchIndel(
            tumorBam,
            indel_candidate_keys,
            params,
            minbq=params["minBq"],
            call_barcodes=indel_call_barcodes,
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
        print(
            f"Process {processNo}: indel depth extraction for "
            f"{len(indel_candidate_keys)} candidates took "
            f"{(time.time() - _indel_depth_t0) / 60:.2f} minutes"
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
            # normal-BAM otherIndelCount) is additionally kept per key so
            # the "if ni > 0" filter below always sees this mutation's own
            # value, not a stale value from whichever mutation a
            # single-pass loop last happened to process.
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
        key = (chrom, pos, ref, alt)
        eligible = key in deferred_depth_keys if coverage_only else flt == "PASS"
        if eligible:
            ta, tr, tdp, na, nr, ndp, ni = muts_indels_dict[key]
        else:
            ta, tr, ti, tdp, na, nr, ni, ndp = (0, 0, 0, 0, 0, 0, 0, 0)
        # Same real-depth-but-rejected reporting as the SNV loop above --
        # see its comment (tagged here, pruned centrally by Caller.py after
        # round 2's merge). ni>0 (other-indel evidence in normal) is left
        # as an unconditional drop, not one of the labeled reasons -- never
        # reported regardless of --rescue.
        if eligible:
            if ta == 0:
                mut["filter"] = "no_good_alt_read"
            # if ti > 0:
            # continue
            elif ta / tdp > params["maxAF"]:
                mut["filter"] = "duplex_vaf"
            elif normalBams:
                if ndp < params["minNdepth"]:
                    mut["filter"] = "n_cov_mask"
                elif na / ndp > params["normalVAF"]:
                    mut["filter"] = "normal_vaf"
                elif ni > 0:
                    continue
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        muts_indels_pass_filter.append(mut)

    ### DBS depth extraction, batched the same way, and filtered under the
    ### same maxAF/minNdepth/normalVAF strategy as SNVs/indels above.
    ### _detect_dbs_pairs only ever emits filter=="PASS" candidates (both
    ### constituent bases already independently passed full duplex-consensus
    ### calling), so — unlike the SNV/indel loops — there's no masked filter
    ### state to preserve here (a candidate either clears this gate and is
    ### kept, or it's dropped) and therefore no round-2 deferred set either
    ### -- DBS depth-extraction only runs in round 1.
    muts_dbs_pass_filter = []
    if not coverage_only:
        muts_dbs_dict = dict()
        dbs_candidate_keys = []
        # Same founding-family barcode collection as the SNV loop above.
        dbs_call_barcodes = {}
        for mut in muts_dbs:
            key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
            _collect_call_barcode(dbs_call_barcodes, key, mut)
            if key not in muts_dbs_dict:
                muts_dbs_dict[key] = None
                dbs_candidate_keys.append(key)

        if dbs_candidate_keys:
            _dbs_depth_t0 = time.time()
            tumor_dbs_depths = extractDepthBatchDbs(
                tumorBam,
                dbs_candidate_keys,
                params,
                minbq=params["minBq"],
                call_barcodes=dbs_call_barcodes,
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
            print(
                f"Process {processNo}: DBS depth extraction for "
                f"{len(dbs_candidate_keys)} candidates took "
                f"{(time.time() - _dbs_depth_t0) / 60:.2f} minutes"
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

        for mut in muts_dbs:
            chrom = mut["chrom"]
            pos = mut["pos"]
            ref = mut["ref"]
            alt = mut["alt"]
            ta, tr, ti, tdp, na, nr, ni, ndp = muts_dbs_dict[(chrom, pos, ref, alt)]
            if tdp > 0 and ta / tdp > params["maxAF"]:
                continue
            if normalBams:
                if ndp < params["minNdepth"]:
                    continue
                if na / ndp > params["normalVAF"]:
                    continue
            mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
            muts_dbs_pass_filter.append(mut)

    # Original per-locus multi-column coverage write-out: final,
    # unconditional flush of whatever remains buffered (per-subprocess
    # only -- see the matching blocks earlier in this function).
    _final_flush_t0 = time.time()
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
            if current_pos <= max_flushed_pos:
                # See the matching comment on the equivalent blocks
                # earlier in this function.
                continue
            max_flushed_pos = current_pos
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
        f"Process {processNo}: final coverage flush took "
        f"{(time.time() - _final_flush_t0) / 60:.2f} minutes"
    )
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
        # Read-family depth-composition totals by trinuc/hp/str context
        # (see _accumulate_depth_matrix above) and the detection-power
        # lookup tables (L for SBS, L_indel_1bp/L_indel_len for indels;
        # None in isLearn mode). Caller.py uses these two ways: (1) at the
        # initial mutation-rate estimate, combined with PASS calls' raw
        # (non-log) likelihood ratios; (2) for the per-channel FDR
        # threshold refinement, re-simulating just that channel's (10,10)
        # power grid at a new, stricter LR threshold (via
        # load_error_matrices + simulate_power_grid) and recombining with
        # these same depth counts for a new effective coverage, without
        # rescanning the bam.
        depth_by_trinuc,
        depth_by_hpstr,
        L,
        L_indel_1bp,
        L_indel_len,
    )
