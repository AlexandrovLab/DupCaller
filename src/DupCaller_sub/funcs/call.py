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
from array import array
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
    indel_strand_evidence,
    calculateSSPosterior,
    calculateDSPosterior,
)
from .learn import profileTriNucMismatches, NUM_BQ
from .misc import getAlignmentObject as BAM
from .misc import build_trinuc64_order
from .misc import load_repeat_context
from .misc import log_progress
from .misc import get_duplex_barcode
from .misc import family_fragment_starts
from .misc import (
    DEPTH_HPSTR_NCAT,
    DEPTH_HPSTR_STR0_ANY,
    STR_EEFF_LENGTHS,
    depth_hpstr_str_bucket,
    depth_hpstr_str0_bucket,
    indel_eeff_site_masks,
    indel_mu_channel,
    str_tract_valid,
    determineTrimLength,
    trim_mask,
    nums2str,
    get_bed_file_for_position,
    _filter_reason_label,
    bamIterateMultipleRegion,
    bamIterateMultipleRegionWithOverflow,
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
from .indels import STR_NCOL, findIndels, indel_mask_span
from .power_tally import (
    PowerTally,
    contract_power_tally,
    indel_keys,
    sbs_power_index,
)


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
            if rec.alts is None:
                # ALT="." (non-variant record): nothing to mask.
                continue
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
                depth_now = extractDepthRegion(nbam, chrom, start, end, params)
                depth += depth_now
            n_cov_mask = depth < params["minNdepth"]

        if params["maxAF"] < 1:
            depth = extractDepthRegion(tbam, chrom, start, end, params)
            ma = params["maxAF"]
            min_depth = math.ceil(1 / ma)
            # OR, not overwrite: a position failing minNdepth in the
            # normal must stay masked when --maxAF also applies.
            n_cov_mask = n_cov_mask | (depth < min_depth)

        nm_mask = prepareAlignMask(tbam, chrom, start, end, params)

    return (
        snp_mask,
        indel_mask,
        noise_mask,
        n_cov_mask,
        include_mask,
        nm_mask,
    )


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
    neg_lr_sbs=None,
    neg_lr_indel=None,
    power_tally=None,
    coverage_totals=None,
):
    """Genotype and record one duplex family, mutating callBam's
    accumulators (coverage arrays, VCF-record lists, learning matrices,
    depth-composition tables, ...) in place.

    Window management (reference_mat_*, coverage allocation/carry) stays
    in callBam. Returns True if any site of the family is in pass_bool
    (reference or alt), i.e. whether it counts toward duplex_count; False
    on every early exit (fully masked, learn mode).
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
    masks[5, :] = trim_mask(masks.shape[1], left, right)
    antimask = np.all(~masks, axis=0)
    antimask[trinuc_np[start_ind:end_ind] > 64] = False
    antimask[ref_np[start_ind:end_ind] == 4] = False
    # "Unmasked": ignores snp_mask/noise_mask (rows 0-1) only.
    unmasked_antimask = np.all(~masks[2:, :], axis=0)
    unmasked_antimask[trinuc_np[start_ind:end_ind] > 64] = False
    unmasked_antimask[ref_np[start_ind:end_ind] == 4] = False
    # Rescue eligibility: under --rescue only include_mask is non-rescuable.
    rescue_antimask = ~masks[3, :]
    rescue_antimask[trinuc_np[start_ind:end_ind] > 64] = False
    rescue_antimask[ref_np[start_ind:end_ind] == 4] = False
    learn_antimask = np.all(~masks[3:, :], axis=0)
    learn_antimask[trinuc_np[start_ind:end_ind] > 64] = False
    learn_antimask[ref_np[start_ind:end_ind] == 4] = False
    ### Skip a fully masked family. Under --rescue gate on rescue_antimask,
    ### so a family blocked only by n_cov/nm/trim still reaches the per-site
    ### logic that reports it as a rescued call.
    if params["rescue"]:
        if not np.any(rescue_antimask):
            return False
    elif not np.any(unmasked_antimask):
        return False
    indel_bool = [("I" in seq.cigarstring or "D" in seq.cigarstring) for seq in readSet]
    if not isLearn:
        masks_indel = np.zeros([6, end_ind_max - start_ind], dtype=bool)
        masks_indel[0, :] = indel_mask[start_ind:end_ind_max]
        masks_indel[1, :] = noise_mask[start_ind:end_ind_max]
        masks_indel[2, :] = n_cov_mask[start_ind:end_ind_max]
        left, right = determineTrimLength(
            readSet[max_ref_num], params=params, processed_flag=processed_flag
        )
        masks_indel[3, :] = trim_mask(masks_indel.shape[1], left, right)
        masks_indel[4, :] = include_mask[start_ind:end_ind_max]
        masks_indel[5, :] = nm_mask[start_ind:end_ind_max]
        antimask_indel = np.all(~masks_indel, axis=0)
        # Ignores noise_mask only.
        unmasked_antimask_indel = np.all(~masks_indel[[0, 2, 3, 4, 5], :], axis=0)
        # Same rescue-eligibility split as the SNV side (include_mask only).
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
            rev_rate,
            str_units,
        ) = genotypeDSIndel(
            readSet,
            rs_reference_start,
            rs_reference_end,
            ref_np[start_ind:end_ind],
            # Genotype over the wide (rescue) scope; masks are applied
            # per call below.
            rescue_antimask_indel,
            hp_raw_np[:, start_ind:end_ind],
            str_raw_np[:, start_ind:end_ind],
            params,
        )
        # Every genotyped candidate is considered; the mask/LR checks
        # below decide its filter label.
        pass_inds = list(range(len(indels)))
        if len(LR_raw) > 0:
            # Threshold capped by the call's own LR ceiling (LM, see
            # indelMaxLR). >= matches simulate_power_grid and Caller.py.
            LR_pass_bool = LR_raw >= np.minimum(params["pcutoffi"], LR_max)
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
            if F1R2_alt_count[pass_inds[nn]] + F1R2_ref_count[pass_inds[nn]] == 0:
                continue
            if F2R1_alt_count[pass_inds[nn]] + F2R1_ref_count[pass_inds[nn]] == 0:
                continue
            # Anchor plus the deleted bases / the insertion's context base,
            # as in learning and the Eeff site masks.
            indel_slice = slice(
                *indel_mask_span(
                    indel_pos - reference_mat_start - start_ind, indel_size
                )
            )
            # LR < 0: only where every other filter passes, its log10 LR
            # goes to neg_lr_indel[channel] as a mu-solve input; no record
            # (it skips depth extraction and is never written to a VCF).
            if LR_raw[pass_inds[nn]] < 0:
                if (
                    neg_lr_indel is None
                    or params.get("coverage_only", False)
                    or flt_rs != "PASS"
                    or not antimask_indel[indel_slice].all()
                ):
                    continue
                channel = indel_mu_channel(
                    indel_size,
                    int(hps[pass_inds[nn]]),
                    int(strs[nn]),
                    int(hp_match[pass_inds[nn]]),
                    indel_ref[1] if indel_size < 0 else indel_alt[1],
                    int(str_units[pass_inds[nn]]),
                )
                if channel is not None:
                    neg_lr_indel.setdefault(channel, array("f")).append(
                        LR_raw[pass_inds[nn]]
                    )
                continue
            elif not LR_pass_bool[nn]:
                continue
            elif antimask_indel[indel_slice].all():
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
                    # Whole-family read totals (as for SBS), not the
                    # locus-specific counts in TC/BC: Estimate.py's
                    # duplex-group lookup keys off F1R2+F2R1.
                    "F1R2": int(F1R2),
                    "F2R1": int(F2R1),
                    "LR": LR_raw[pass_inds[nn]],
                    "LM": LR_max[pass_inds[nn]],
                    # Filled in by Caller.py for PASS calls from _strand
                    # (per strand: informative reads, REF reads, reversion
                    # rate); _strand is not written to the VCF.
                    "MSP": ".",
                    "_strand": (
                        (
                            indel_strand_evidence(
                                F1R2_alt_count[pass_inds[nn]],
                                F1R2_ref_count[pass_inds[nn]],
                                rev_rate[pass_inds[nn]],
                            ),
                            indel_strand_evidence(
                                F2R1_alt_count[pass_inds[nn]],
                                F2R1_ref_count[pass_inds[nn]],
                                rev_rate[pass_inds[nn]],
                            ),
                        )
                        if flt == "PASS" and not params.get("coverage_only", False)
                        else None
                    ),
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
                    "FS": family_fragment_starts(readSet),
                    "STR": strs[nn],
                    "HM": int(hp_match[pass_inds[nn]]),
                    # Repeat unit of the STR tract the indel slips (0: not
                    # a slip); routes the call to its STR channel.
                    "SU": int(str_units[pass_inds[nn]]),
                    "SNPM": 0,
                    "NOISEM": int(np.any(masks_indel[1, indel_slice])),
                },
                "formats": ["AC", "RC", "DP"],
            }
            muts_indels.append(indel_rec)
            indel_dict[indel_str] = 1
    ### SNV genotyping (and error learning in learn mode)
    if 1:
        if isLearn:
            # Per-strand/per-family gating (srdMinRead/ssmMinRead) is
            # applied inside profileTriNucMismatches, not here.
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
            return False
        (
            cov_mat,
            LR_raw_mut,
            LR_max_mut,
            mut_mask,
            b1_int,
            unmasked_antimask,
            F1R2_count,
            F2R1_count,
            strand_ev_F1R2,
            strand_ev_F2R1,
            neg_pos,
            neg_alt,
            neg_LR,
        ) = genotypeDSSnv(
            readSet,
            rs_reference_start,
            ref_np[start_ind:end_ind],
            trinuc_np[start_ind:end_ind],
            np.copy(unmasked_antimask),
            np.copy(rescue_antimask),
            params,
            L,
        )
        if power_tally is not None:
            # The L index cov_mat would have been looked up at.
            sbs_pidx = sbs_power_index(
                F1R2_count,
                F2R1_count,
                trinuc_np[start_ind:end_ind],
                unmasked_antimask,
            )
        # Per-position indel opportunity: detection power looked up in
        # L_indel_1bp/L_indel_len by repeat context and per-strand
        # family read counts (capped at 9).
        if not isLearn:
            unit_len_arr = hp_np[0, start_ind:end_ind].astype(int)
            repeat_count_arr = hp_np[2, start_ind:end_ind].astype(int)
            ref_allele_arr = ref_np[start_ind:end_ind]
            ref_allele_safe = np.where(ref_allele_arr > 3, 0, ref_allele_arr)
            is_hp = unit_len_arr <= 1
            hps_for_row = np.where(is_hp, np.minimum(repeat_count_arr, 10), 1)
            total_len_arr = unit_len_arr * repeat_count_arr
            # STR-length bin (1="2-9" .. 4="40+"), >=1 wherever ~is_hp.
            strs_for_row = np.where(~is_hp, 1, 0)
            strs_for_row[np.logical_and(~is_hp, total_len_arr >= 10)] = 2
            strs_for_row[np.logical_and(~is_hp, total_len_arr >= 25)] = 3
            strs_for_row[np.logical_and(~is_hp, total_len_arr >= 40)] = 4
            row_len = strs_for_row
            # Everything below uses the raw, unmerged hp_raw_np/str_raw_np,
            # so HP and STR opportunity are credited independently (a
            # position inside an STR tract still has HP opportunity).
            # depth_by_hpstr hp bucket: base (A/T/C/G) x run length 1-10.
            hp_run_capped10 = np.minimum(
                hp_raw_np[0, start_ind:end_ind].astype(int), 10
            )
            hp_bucket_row = ref_allele_safe * 10 + hp_run_capped10 - 1
            # L_indel_1bp's run axis keeps exactly-10 and 11+ apart: a
            # deletion from an 11+ run reverts as a 10+ run.
            hp_run_power = np.minimum(hp_raw_np[0, start_ind:end_ind].astype(int), 11)
            str_unit_row = str_raw_np[0, start_ind:end_ind].astype(int)
            str_repeat_row = str_raw_np[1, start_ind:end_ind].astype(int)
            is_real_str_row = str_unit_row >= 2
            total_len_row = str_unit_row * str_repeat_row
            str_bin_row = np.where(is_real_str_row, 1, 0)
            str_bin_row[is_real_str_row & (total_len_row >= 10)] = 2
            str_bin_row[is_real_str_row & (total_len_row >= 25)] = 3
            str_bin_row[is_real_str_row & (total_len_row >= 40)] = 4
            unit_len_clamped = np.clip(unit_len_arr, 2, 5)
            n_top_indel = np.minimum(F1R2_count.sum(axis=0), 9).astype(int)
            n_bot_indel = np.minimum(F2R1_count.sum(axis=0), 9).astype(int)

            # Run-boundary visibility: positions in the first and last run
            # of the family's span are invalid, since the reads may not
            # cover those runs in full. Shared by cov_mat_indel and indel100.
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
            # str.h5 flags only tract starts, so validity comes from each
            # tract's own end, not from the next cut (str_tract_valid).
            str_repeat_valid = (
                str_tract_valid(str_cut_bool, total_len_row) & antimask_indel
            )
            unit_len_clamped_str = np.clip(str_unit_row, 2, 5)
            # Next reference base: a rep0 (novel) insertion of base N
            # requires only that it differs from N.
            window_len = end_ind - start_ind
            next_ref_arr = np.full(window_len, -1, dtype=int)
            avail = min(end_ind + 1, ref_np.shape[0]) - (start_ind + 1)
            if avail > 0:
                next_ref_arr[:avail] = ref_np[start_ind + 1 : start_ind + 1 + avail]
            next_ref_valid = (next_ref_arr >= 0) & (next_ref_arr <= 3)

            cov_mat_indel = np.zeros([window_len, 16])
            # Cols 0/1: deletion/insertion of an STR repeat unit, at tract
            # starts only; 0 outside a real STR (HP is cols 14/15).
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
            # del_2..del_5plus: fixed row-0 context, not the position's
            # actual repeat context.
            cov_mat_indel[:, 2] = L_indel_len[n_top_indel, n_bot_indel, 0, 3]
            cov_mat_indel[:, 3] = L_indel_len[n_top_indel, n_bot_indel, 0, 2]
            cov_mat_indel[:, 4] = L_indel_len[n_top_indel, n_bot_indel, 0, 1]
            cov_mat_indel[:, 5] = L_indel_len[n_top_indel, n_bot_indel, 0, 0]
            # Insertion A/T/C/G: novel 1bp insertion (next ref base !=
            # inserted base; repeat extensions are HP opportunity). The
            # row-0 context has no base axis, so all 4 columns share
            # L_indel_len[..., 0, 6]. Not STR-gated.
            for b in range(4):
                cov_mat_indel[:, 6 + b] = (
                    L_indel_len[n_top_indel, n_bot_indel, 0, 6]
                    * next_ref_valid
                    * (next_ref_arr != b)
                    * antimask_indel
                )
            # Insertion length 2..5+: flat "novel U-bp insertion"
            # opportunity, gated only by antimask_indel (a hypothetical
            # insertion does not depend on the local repeat context).
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
            # Cols 14/15: deletion/insertion of a homopolymer base, at
            # run starts only.
            cov_mat_indel[:, 14] = (
                L_indel_1bp[n_top_indel, n_bot_indel, hp_run_power, ref_allele_safe, 0]
                * hp_repeat_valid
                * hp_cut_bool
            )
            cov_mat_indel[:, 15] = (
                L_indel_1bp[n_top_indel, n_bot_indel, hp_run_power, ref_allele_safe, 1]
                * hp_repeat_valid
                * hp_cut_bool
            )
            # 100-class indel opportunity: HP del (24) + ins (20, no rep0),
            # STR del/ins by unit x count (24+24), MH deletion length (4),
            # 1bpins{base} (4). Layout: misc.py's build_indel100_labels.
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

            # HP and STR opportunity are credited independently from the
            # raw hp/str arrays, not load_repeat_context's STR-priority
            # merge: a homopolymer inside an STR tract (the "AA" in
            # (AAT)n) can slip either way. That merge is only right for
            # classifying an observed event. Mirrors misc.py's
            # indel100_reference_bucket_indices.
            unit_bucket = np.clip(str_unit_row, 2, 5) - 2
            count_bucket_del = np.clip(str_repeat_row, 1, 6) - 1

            # Homopolymer deletion (cols 0-23).
            hp_del_power = L_indel_1bp[
                n_top_indel, n_bot_indel, hp_run_power, ref_allele_safe, 0
            ]
            del_bucket = np.clip(hp_run_capped10, 1, 6) - 1
            _accumulate_indel100(
                base_cgta * 6 + del_bucket,
                hp_del_power,
                hp_repeat_valid & hp_cut_bool,
            )

            # Homopolymer insertion (cols 24-43), length bins 1..5+ only;
            # rep0 comes from the 1bpins{base} columns instead (Estimate.py's
            # override_inshp0_with_next_base_opportunity).
            hp_ins_power = L_indel_1bp[
                n_top_indel, n_bot_indel, hp_run_power, ref_allele_safe, 1
            ]
            ins_bucket = np.clip(hp_run_capped10, 1, 5) - 1
            _accumulate_indel100(
                24 + base_cgta * 5 + ins_bucket,
                hp_ins_power,
                (hp_run_capped10 >= 1) & hp_repeat_valid & hp_cut_bool,
            )

            # STR deletion (cols 44-67), annotated tracts: same value/gate
            # as cov_mat_indel[:, 0].
            str_del_power = L_indel_len[
                n_top_indel, n_bot_indel, str_bin_row, -unit_len_clamped_str + 5
            ]
            _accumulate_indel100(
                44 + unit_bucket * 6 + count_bucket_del,
                str_del_power,
                is_real_str_row & str_repeat_valid & str_cut_bool,
            )
            # Flat rep1 opportunity for an arbitrary non-repeating U-bp
            # deletion: fixed row-0 power, scan-free (as on the genome-wide
            # side), gated only by antimask_indel -- also inside STRs, where
            # a different-unit deletion is an independent opportunity.
            for u_idx, U in enumerate([2, 3, 4, 5]):
                hyp_del_power = L_indel_len[n_top_indel, n_bot_indel, 0, -U + 5]
                _accumulate_indel100(
                    np.full(end_ind - start_ind, 44 + u_idx * 6),
                    hyp_del_power,
                    antimask_indel,
                )

            # STR insertion (cols 68-91), annotated tracts.
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
            # Flat rep0/rep1 opportunity for a novel U-bp insertion: row-0
            # power added to both the "0" and "1" count buckets, gated
            # only by antimask_indel (as for deletions above).
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

            # Microhomology deletion (cols 92-95): reuse cov_mat_indel[:, 2:6].
            for k in range(4):
                _accumulate_indel100(
                    np.full(end_ind - start_ind, 92 + k),
                    cov_mat_indel[:, 2 + k],
                    antimask_indel,
                )

            # 1bpins{base} rep0 (cols 96-99): reuse cov_mat_indel[:, 6:10],
            # already gated by antimask_indel and next base.
            for b in range(4):
                _accumulate_indel100(
                    np.full(window_len, 96 + b),
                    cov_mat_indel[:, 6 + b],
                    np.ones(window_len, dtype=bool),
                )
            if power_tally is not None:
                indel_tally_keys = indel_keys(
                    n_top_indel,
                    n_bot_indel,
                    hp_run_power,
                    ref_allele_safe,
                    hp_run_capped10,
                    hp_repeat_valid,
                    hp_cut_bool,
                    is_real_str_row,
                    str_repeat_valid,
                    str_cut_bool,
                    str_bin_row,
                    str_unit_row,
                    str_repeat_row,
                    antimask_indel,
                    next_ref_arr,
                    next_ref_valid,
                )
        ref_int = ref_np[start_ind:end_ind]
        n_win = ref_int.size
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
        # >=, matching simulate_power_grid and Caller.py's final check.
        LR_pass_bool[mut_mask] = LR_raw_mut >= mut_thresholds
        alt_int = b1_int
        unmasked_pass_bool = np.full(n_win, False, dtype=bool)
        unmasked_pass_bool[refs_ind] = True
        unmasked_pass_bool[muts_ind] = True
        pass_bool = np.copy(unmasked_pass_bool)
        pass_bool[~antimask] = False
        # LR < 0 sites that pass every other filter (unmasked, family PASS,
        # reads on both strands) go to neg_lr_sbs[trinuc * 4 + base] for
        # the per-channel mu solve. Round 1 only.
        if (
            neg_lr_sbs is not None
            and neg_pos.size
            and flt_rs == "PASS"
            and not params.get("coverage_only", False)
        ):
            neg_trinuc = trinuc_np[start_ind:end_ind][neg_pos]
            neg_ok = (
                pass_bool[neg_pos]
                & (F1R2_count[:, neg_pos].sum(axis=0) > 0)
                & (F2R1_count[:, neg_pos].sum(axis=0) > 0)
                & (neg_trinuc < 64)
            )
            for k, lr in zip(
                (neg_trinuc * 4 + neg_alt)[neg_ok].tolist(), neg_LR[neg_ok].tolist()
            ):
                neg_lr_sbs[k].append(lr)
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
            if F1R2_count[:, muts_ind[nn]].sum() == 0:
                continue
            if F2R1_count[:, muts_ind[nn]].sum() == 0:
                continue
            if not LR_pass_bool[muts_ind[nn]]:
                continue
            if pass_bool[muts_ind[nn]]:
                flt = flt_rs
            elif masks[0, muts_ind[nn]] or masks[1, muts_ind[nn]]:
                # "masked" = snp_mask/noise_mask only (as SNPM/NOISEM);
                # other mask reasons go to the rescue branch.
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
                    # Filled in by Caller.py for PASS calls from _strand
                    # (not written to the VCF).
                    "MSP": ".",
                    "_strand": (
                        (
                            strand_ev_F1R2[muts_ind_compressed[nn]],
                            strand_ev_F2R1[muts_ind_compressed[nn]],
                        )
                        if flt == "PASS" and not params.get("coverage_only", False)
                        else None
                    ),
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
                    "FS": family_fragment_starts(readSet),
                    "STR": 0,
                    "SNPM": int(masks[0, muts_ind[nn]]),
                    "NOISEM": int(masks[1, muts_ind[nn]]),
                },
                "formats": ["AC", "RC", "DP"],
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
                readSet,
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
            if power_tally is not None:
                power_tally.add_sbs(
                    duplex_no, sbs_pidx, pass_bool, unmasked_pass_bool, ref_int
                )

            unmasked_coverage[start_ind:end_ind][unmasked_pass_bool] += cov_mat[
                unmasked_pass_bool
            ]
            if coverage_totals is not None:
                coverage_totals[0] += cov_mat[pass_bool].sum()
                coverage_totals[1] += cov_mat[unmasked_pass_bool].sum()
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
                if power_tally is not None:
                    power_tally.add_indel(duplex_no, indel_tally_keys)
                # Read-family depth-composition totals (Eeff). They cover
                # the same sites as the mu solve: SBS uses pass_bool; indels
                # use the indel masks and one site per run/tract
                # (indel_eeff_site_masks), with depth taken at the anchor
                # (the base before the context base), which every read that
                # can show the event covers -- a deletion read has no base
                # at a deleted context base.
                trinuc_ctx = trinuc_np[start_ind:end_ind]
                n_top_anchor = np.zeros_like(n_top_indel)
                n_top_anchor[1:] = n_top_indel[:-1]
                n_bot_anchor = np.zeros_like(n_bot_indel)
                n_bot_anchor[1:] = n_bot_indel[:-1]
                (
                    hp_eeff_bool,
                    str_eeff_by_unit,
                    str0_any_eeff_bool,
                    str0_eeff_by_len,
                ) = indel_eeff_site_masks(
                    antimask_indel,
                    ref_int,
                    hp_cut_bool,
                    hp_repeat_valid,
                    str_unit_row,
                    total_len_row,
                    str_cut_bool,
                    str_repeat_valid,
                )
                _accumulate_depth_matrix(
                    depth_by_trinuc,
                    n_top_indel,
                    n_bot_indel,
                    trinuc_ctx,
                    pass_bool & (trinuc_ctx < 64),
                    64,
                )
                _accumulate_depth_matrix(
                    depth_by_hpstr,
                    n_top_anchor,
                    n_bot_anchor,
                    hp_bucket_row,
                    hp_eeff_bool,
                    DEPTH_HPSTR_NCAT,
                )
                _accumulate_depth_matrix(
                    depth_by_hpstr,
                    n_top_anchor,
                    n_bot_anchor,
                    np.full(window_len, DEPTH_HPSTR_STR0_ANY),
                    str0_any_eeff_bool,
                    DEPTH_HPSTR_NCAT,
                )
                for k in STR_EEFF_LENGTHS:
                    for signed_k in (-k, k):
                        _accumulate_depth_matrix(
                            depth_by_hpstr,
                            n_top_anchor,
                            n_bot_anchor,
                            depth_hpstr_str_bucket(
                                np.maximum(str_bin_row, 1), signed_k
                            ),
                            str_eeff_by_unit[signed_k],
                            DEPTH_HPSTR_NCAT,
                        )
                        _accumulate_depth_matrix(
                            depth_by_hpstr,
                            n_top_anchor,
                            n_bot_anchor,
                            np.full(window_len, depth_hpstr_str0_bucket(signed_k)),
                            str0_eeff_by_len[signed_k],
                            DEPTH_HPSTR_NCAT,
                        )
                dinuc_ctx = ref_allele_safe * 4 + np.where(
                    next_ref_arr >= 0, next_ref_arr, 0
                )
                # A dinuc is excluded if either of its two positions is masked.
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
    """Add a founding family's (TAG1, TAG2, |TL|, fragment start) ids to
    call_barcodes[key], so extractDepthBatch*'s call_barcodes can exempt
    its reads from the minBq and --mapq filters."""
    infos = mut["infos"]
    ids = call_barcodes.setdefault(key, set())
    for fs in infos["FS"]:
        ids.add((infos["TAG1"], infos["TAG2"], abs(infos["TL"]), fs))


COVERAGE_WINDOW_ROWS = 1000000


def _cut_coverage_leftover(arrays, old_start, new_start, same_chrom):
    """Split a window's per-locus coverage arrays at the next window's start.

    arrays are row-aligned (row i = genomic old_start + i). Returns
    (flush_rows, leftovers, leftover_gstart): rows [0, flush_rows) lie
    before new_start and are final; leftovers are copies of the rows from
    genomic leftover_gstart onward (trailing all-zero rows trimmed), to be
    added back into the next window with _carry_in_coverage. On a
    chromosome change everything is flushed and leftovers is None. When
    new_start < old_start nothing is flushed and the whole window is
    carried.
    """
    rows = arrays[0].shape[0]
    if not same_chrom:
        return rows, None, None
    cut = min(max(new_start - old_start, 0), rows)
    tail = np.zeros(rows - cut)
    for arr in arrays:
        tail += arr[cut:].sum(axis=1)
    nz = np.flatnonzero(tail)
    if nz.size == 0:
        return cut, None, None
    n = int(nz[-1]) + 1
    return cut, [arr[cut : cut + n].copy() for arr in arrays], old_start + cut


def _carry_in_coverage(arrays, leftovers, leftover_gstart, window_start):
    """Add carried rows back at their genomic position in the current window."""
    off = leftover_gstart - window_start
    for arr, rows in zip(arrays, leftovers):
        arr[off : off + rows.shape[0]] += rows


def _coverage_window_rows(leftovers, leftover_gstart, window_start):
    """Rows for a new window's coverage arrays: the usual 1 Mb, or more if
    the carried rows (which start at or after window_start) reach past it."""
    if leftovers is None:
        return COVERAGE_WINDOW_ROWS
    return max(
        COVERAGE_WINDOW_ROWS,
        leftover_gstart - window_start + leftovers[0].shape[0],
    )


def _early_deferred_depth(records, tumorBam, normalBams, params, extract, indel):
    """--skipCoveragePass: depth-extract, in round 1, every masked candidate
    that could join Caller.py's deferred set (raw LR >= --minLR, the floor
    of every final threshold), as round 2 would, and store it on the
    record as "deferred_depth": a list of (LR cut, (samples, filter)).
    Records round 2 would drop (indel ni > 0) get none.

    Founding families (exempt from minBq/mapq) depend on the final SNV
    threshold, unknown here, so SNV tumor depth is extracted once per
    distinct family LR, with the families at or above it; Caller.py picks
    the entry matching the final threshold. Indels: one entry, cut -inf."""
    floor = params.get("minLR", 0.0)
    keys = {}
    for mut in records:
        infos = mut["infos"]
        if (
            mut["filter"] == "masked"
            and infos["LR"] >= floor
            and (infos.get("NOISEM") or (not indel and infos.get("SNPM")))
        ):
            keys[(mut["chrom"], mut["pos"], mut["ref"], mut["alt"])] = None
    if not keys:
        return
    at_key = {}
    for mut in records:
        key = (mut["chrom"], mut["pos"], mut["ref"], mut["alt"])
        if key in keys:
            at_key.setdefault(key, []).append(mut)
    cuts = {
        key: (
            [-np.inf]
            if indel
            else sorted(
                {m["infos"]["LR"] for m in muts if m["infos"]["LR"] >= floor},
                reverse=True,
            )
        )
        for key, muts in at_key.items()
    }
    tumor_levels = []
    for level in range(max(len(c) for c in cuts.values())):
        level_keys = [k for k, c in cuts.items() if len(c) > level]
        barcodes = {}
        for key in level_keys:
            for mut in at_key[key]:
                if mut["infos"]["LR"] >= cuts[key][level]:
                    _collect_call_barcode(barcodes, key, mut)
        tumor_levels.append(
            extract(
                tumorBam,
                level_keys,
                params,
                minbq=params["minBq"],
                call_barcodes=barcodes,
            )
        )
    key_list = list(keys)
    normals = (
        [extract(nb, key_list, params, minbq=params["minBq"]) for nb in normalBams]
        if normalBams
        else []
    )
    for key, muts in at_key.items():
        na = nr = ni = ndp = 0
        for depths in normals:
            na_now, nr_now, ni_now, ndp_now = depths.get(key, (0, 0, 0, 0))
            na += na_now
            nr += nr_now
            ni += ni_now
            ndp += ndp_now
        for mut in muts:
            entries = []
            for level, cut in enumerate(cuts[key]):
                ta, tr, ti, tdp = tumor_levels[level].get(key, (0, 0, 0, 0))
                flt = mut["filter"]
                if ta == 0:
                    flt = "no_good_alt_read"
                elif ta / tdp > params["maxAF"]:
                    flt = "duplex_vaf"
                elif normalBams:
                    if ndp < params["minNdepth"]:
                        flt = "n_cov_mask"
                    elif na / ndp > params["normalVAF"]:
                        flt = "normal_vaf"
                    elif indel and ni > 0:
                        break  # round 2 drops the record
                entries.append((cut, ([[ta, tr, tdp], [na, nr, ndp]], flt)))
            else:
                mut["deferred_depth"] = entries


def sample_power_quals(tumorBam, regions, min_bq, max_reads=1000):
    """Base qualities (>= min_bq) of the first max_reads mapped reads in
    regions, sampled by simulate_power_grid; [30] if none."""
    all_quals = []
    reads_sampled = 0
    for region in regions:
        for read in tumorBam.fetch(*region):
            if not read.is_unmapped and read.query_alignment_qualities is not None:
                quals = np.array(read.query_alignment_qualities, dtype=float)
                valid_quals = quals[quals >= min_bq]
                all_quals.extend(valid_quals.tolist())
                reads_sampled += 1
                if reads_sampled >= max_reads:
                    break
        if reads_sampled >= max_reads:
            break
    if not all_quals:
        all_quals = [30]
    return np.array(all_quals, dtype=float)


def build_power_tables(params, all_quals, build=True):
    """Per-(trinuc, alt) SBS thresholds pcutoff_sbs (also set in
    params["pcutoff_sbs"]) and detection-power tables L (SBS), L_indel_1bp,
    L_indel_len (indels; None in learn mode), honoring per-channel
    overrides. build=False leaves the power tables all-zero (round 1).
    """
    pcut = params["pcutoff"]
    isLearn = params.get("isLearn", False)
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2trinuc = params["num2trinuc_list"]
    # SBS power table L[top count 0-9, bottom count 0-9, trinuc 0-63,
    # alt base A/T/C/G]. All-zero unless build; genotypeDSSnv then treats it
    # like L=None (cov_mat stays 0).
    L = np.zeros([10, 10, 64, 4])

    prob_amp_mat_L = params["ampmat"]
    prob_amp_mat_rev_L = params["ampmat_rev"]
    prob_dmg_t_L = params["dmgmat_top"]
    prob_dmg_rev_t_L = params["dmgmat_rev_top"]
    prob_dmg_b_L = params["dmgmat_bot"]
    prob_dmg_rev_b_L = params["dmgmat_rev_bot"]
    trinuc_conv_np_L = params["trinuc_convert"]
    N_SIM = 100
    # Each simulate_power_grid call derives its own RNG from (seed,
    # threshold) via threshold_rng, so results don't depend on -p.
    base_seed_L = params["seed"]

    # Per-(trinuc, alt) SBS threshold: min(pcutoff, maxLR), where
    # maxLR = log10(1-Pdmg_t) + log10(1-Pdmg_b) - log10(Pdmg_rev_t)
    # - log10(Pdmg_rev_b) is the LR ceiling for that context (genotypeDSSnv's
    # LM). Context t with alt b and context t+32 with alt b^1 are reverse
    # complements with identical damage rates, so only t < 32 is computed
    # and mirrored.
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
    # Per-channel FDR-refined thresholds from Caller.py, keyed by (trinuc,
    # alt); applied after (and not capped by) the maxLR ceiling.
    for (t_ov, b_ov), new_threshold in params.get("pcutoff_sbs_override", {}).items():
        pcutoff_sbs[t_ov, b_ov] = new_threshold
    params["pcutoff_sbs"] = pcutoff_sbs

    if build:
        for t in range(64):
            ref_base_idx_L = base2num[num2trinuc[t][1]]
            for b in range(4):
                if b == ref_base_idx_L:
                    continue  # ref -> ref is no opportunity; stays 0
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

    # Indel power tables: fraction of simulated families, per (top count,
    # bottom count) and repeat context (hps, strs, idLen -> indelErrorProbs,
    # as in genotypeDSIndel), that clear the context's threshold. None in
    # learn mode.
    L_indel_1bp = None
    L_indel_len = None
    if not isLearn:
        prob_amp_hp_L = params["ampmat_hp"]
        prob_dmg_hp_L = params["dmgmat_hp"]
        prob_amp_str_L = params["ampmat_str"]
        prob_dmg_str_L = params["dmgmat_str"]
        pcutoffi_L = params["pcutoffi"]
        min_lr_L = params.get("minLR", 0.0)

        def indelContextThreshold(Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c):
            # min(pcutoffi, context LR ceiling (indelMaxLR)), as for SBS,
            # floored at --minLR like Caller.py's cutoff for contexts
            # without a channel threshold.
            return max(
                min_lr_L,
                min(
                    pcutoffi_L,
                    float(indelMaxLR(Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c)),
                ),
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

        # 1bp HP indels: axes top count, bottom count, hps (0-11, 11 = 11+;
        # a deletion from an 11+ run reverts as 10+, from exactly 10 as 9 --
        # see _indel_error_cells), ref_allele (0-3), sign (0 = del, 1 = ins).
        # inserted_base == ref_allele: same-base run extension only. A
        # mismatched 1bp insertion has no per-(hps, base) grid (see
        # prob.py's indelErrorProbs); it is scored but not FDR-refined.
        L_indel_1bp = np.zeros([10, 10, 12, 4, 2])
        indel_1bp_threshold_override = params.get("indel_1bp_threshold_override", {})
        if build:
            for hps_c in range(12):
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
                        # Per-(hp length, sign, base pool) override; pool
                        # "C" = C/G (ref_allele 2/3 in A,T,C,G order), "T" =
                        # A/T, as in Caller.py's raw_lr_hp.
                        pool_c = "C" if ref_allele_c in (2, 3) else "T"
                        # Channels stop at HP10 (10+).
                        threshold = indel_1bp_threshold_override.get(
                            (min(hps_c, 10), sign_idx, pool_c),
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

        # >=2bp indels: axes top count, bottom count, STR bin (0 = not a
        # repeat .. 4 = "40+"), idLen+5. Rows 1-4: |idLen| is the repeat
        # unit (one-unit slip, as the STR channels); row 0: indel length.
        L_indel_len = np.zeros([10, 10, 5, 11])
        len_idLens = (-5, -4, -3, -2, 2, 3, 4, 5)
        indel_len_threshold_override = params.get("indel_len_threshold_override", {})
        if build:
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
                        unit_len=abs(idLen_c),
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
            # Row 0, idLen=+1: novel (mismatched) 1bp insertion, used by
            # the "Insertion A/T/C/G" coverage columns. Other +-1bp cases
            # are HP context (L_indel_1bp).
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
    return pcutoff_sbs, L, L_indel_1bp, L_indel_len


def contract_worker_tally(params, summary):
    """--skipCoveragePass: contract one round-1 worker's power tallies with
    the power tables round 2 would have built (same thresholds, seed and
    quality sample). Returns contract_power_tally's result."""
    _, L, L_indel_1bp, L_indel_len = build_power_tables(params.copy(), summary["quals"])
    return contract_power_tally(summary, L, L_indel_1bp, L_indel_len)


def callBam(params, processNo):
    bam = params["tumorBam"]
    nbams = params["normalBams"]
    regions = params["regions"]
    # Raw start of this worker's first region (None if no preceding worker),
    # used to keep family processing exactly-once across workers. Distinct
    # from regions_start_pos, which routes coverage.bed output.
    prev_boundary_chrom = None
    prev_boundary_pos = None
    if len(regions[0]) > 1 and regions[0][1] != 0:
        prev_boundary_chrom = regions[0][0]
        prev_boundary_pos = regions[0][1]
    coverage_start = params.get("coverage_start")
    if coverage_start is not None and coverage_start[0] == regions[0][0]:
        regions_start_chrom, regions_start_pos = coverage_start
    elif len(regions[0]) == 1:
        regions_start_chrom = regions[0][0]
        regions_start_pos = 0
    else:
        regions_start_chrom = regions[0][0]
        original_start_pos = regions[0][1]

        if original_start_pos == 0:
            regions_start_pos = 0
        else:
            # Coverage starts at max(reference_end) + 1 over the reads that
            # overlap the base just before the region (falls back to the
            # region start if there are none or the BAM can't be read).
            target_position = original_start_pos - 1

            try:
                bamObject = BAM(bam, "rb", params.get("reference"))
                max_reference_end = 0

                for read in bamObject.fetch(
                    regions_start_chrom, target_position, target_position + 1
                ):
                    if not read.is_unmapped and read.reference_end is not None:
                        max_reference_end = max(max_reference_end, read.reference_end)

                if max_reference_end > 0:
                    regions_start_pos = max_reference_end + 1
                else:
                    regions_start_pos = original_start_pos

                bamObject.close()
            except Exception:
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
    unique_read_num = 0
    pass_read_num = 0
    FPs = []
    RPs = []
    indel_dict = dict()
    # float32 log10 LRs of LR < 0 sites for the mu solve: SBS by
    # trinuc * 4 + base, indels by indel_mu_channel name.
    neg_lr_sbs = [array("f") for _ in range(256)]
    neg_lr_indel = {}
    mismatch_mat = np.zeros([64, 4])
    hp_alt_mat = np.zeros([10, 12])
    str_alt_mat = np.zeros([5, STR_NCOL])
    mismatch_dmg_mat = np.zeros([64, 4])
    hp_dmg_mat = np.zeros([10, 12])
    str_dmg_mat = np.zeros([5, STR_NCOL])
    # SBS amp-error BQ histogram (learn mode; see learn.py's
    # profileTriNucMismatches).
    sbs_alt_bq_hist_mat = np.zeros([64, 4, NUM_BQ])
    # Loads error matrices and trinuc lookups into params.
    load_error_matrices(params)
    trinuc2num = params["trinuc2num_dict"]
    num2trinuc = params["num2trinuc_list"]

    # SNV coverage totals (masked, unmasked), summed per family rather than
    # at the coverage.bed flush, which skips noise-masked-only positions.
    coverage_totals = np.zeros(2)
    total_coverage_indel_cat = np.zeros(16)
    total_unmasked_coverage_indel_cat = np.zeros(16)
    # Read-family depth composition: depth_mat[top_count, bot_count,
    # category] = number of passing loci with that (F1R2, F2R1) count
    # (each capped at 9) in that context. Pure aggregates, so nothing needs
    # carrying across windows or merging across processes.
    depth_by_trinuc = np.zeros((10, 10, 64), dtype=np.int64)  # 64 trinuc contexts
    # 40 hp (A/T/C/G x len 1-10, capped) + 16 STR (bin x length) + 5 STR0
    # buckets, see DEPTH_HPSTR_NCAT.
    depth_by_hpstr = np.zeros((10, 10, DEPTH_HPSTR_NCAT), dtype=np.int64)
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
    recCount = 0
    currentCheckPoint = 1000000
    lastTime = 0
    duplex_count = 0
    reference_mat_chrom = "anyChrom"
    reference_mat_start = 0
    # Highest position already written to the coverage beds on
    # reference_mat_chrom. A backward window re-init can revisit flushed
    # positions; flushes emit only strictly increasing positions to keep the
    # bed sorted and avoid double counting.
    max_flushed_pos = -1
    # Coverage rows carried from one window to the next (see
    # _cut_coverage_leftover) and the genomic position of their first row.
    coverage_leftovers = None
    coverage_leftover_gstart = None
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
        # Same overflow-aware iterator as the main pass, so overflow reads
        # are screened too.
        for rec, region in bamIterateMultipleRegionWithOverflow(
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
    coverage_only = params.get("coverage_only", False)
    # --skipCoveragePass: tally coverage opportunities by power-table index
    # (PowerTally); Caller.py contracts them with power tables built at the
    # final thresholds.
    power_tally = PowerTally() if params.get("power_tally", False) else None
    all_quals = (
        sample_power_quals(tumorBam, regions, params["minBq"])
        if coverage_only or power_tally is not None
        else None
    )
    pcutoff_sbs, L, L_indel_1bp, L_indel_len = build_power_tables(
        params, all_quals, build=coverage_only
    )
    # Families are keyed by label and require an exact reference_start match.
    currentReadDict = {}
    # Reroute a downstream mate whose alignment start disagrees with its
    # family's next_reference_start into that family's group (see
    # _index_rugged_mates).
    rugged_reads_index = {}
    rugged_reads_pool = {}
    # Chromosome of the pending rugged entries, purged on chromosome change.
    rugged_pool_chrom = None
    for rec, region in bamIterateMultipleRegionWithOverflow(
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
        # Skip reads whose 5' end is clipped
        if (rec.is_forward and rec.cigartuples[0][0] != 0) or (
            rec.is_reverse and rec.cigartuples[-1][0] != 0
        ):
            continue
        if (
            prev_boundary_chrom is not None
            and rec.reference_name == prev_boundary_chrom
            and rec.next_reference_start < prev_boundary_pos
        ):
            # Mate starts in the previous worker's region, which claims
            # the family via its overflow pass; skip to count it once.
            continue
        pass_read_num += 1
        if rec.template_length < 0 and rec.query_name in rugged_reads_index:
            rugged_reads_pool[rugged_reads_index.pop(rec.query_name)].append(rec)
            continue
        start = rec.reference_start
        label = _compute_read_label(rec, params)
        chrom = tumorBam.get_reference_name(rec.reference_id)
        if chrom != rugged_pool_chrom:
            # Coordinate-sorted BAM: redirects left from a previous
            # chromosome can never be matched.
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
        else:
            """
            Calling block starts
            """
            _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool)
            # Rerouted rugged mates can start a few bp before currentStart,
            # so floor the window at the batch's minimum read start or a
            # later key's start_ind could go negative. currentReadDict can
            # be empty if its only read was rerouted.
            batch_min_start = min(
                (
                    r.reference_start
                    for entry in currentReadDict.values()
                    for r in entry["seqs"]
                ),
                default=currentStart,
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
                    contig_len_now = tumorBam.get_reference_length(chromNow)
                    if rs_reference_start >= contig_len_now:
                        # Would give an empty reference window.
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
                        ### Flush the finished part of the window's coverage
                        if "coverage" in locals():
                            if coverage_leftovers is not None:
                                _carry_in_coverage(
                                    (
                                        coverage,
                                        coverage_indel_cat,
                                        unmasked_coverage,
                                        unmasked_coverage_indel_cat,
                                    ),
                                    coverage_leftovers,
                                    coverage_leftover_gstart,
                                    reference_mat_start,
                                )
                            # Split at the next window's start (batch_min_start):
                            # rows before it are final, the rest is carried over.
                            (
                                flush_rows,
                                coverage_leftovers,
                                coverage_leftover_gstart,
                            ) = _cut_coverage_leftover(
                                (
                                    coverage,
                                    coverage_indel_cat,
                                    unmasked_coverage,
                                    unmasked_coverage_indel_cat,
                                ),
                                reference_mat_start,
                                batch_min_start,
                                chromNow == reference_mat_chrom,
                            )
                            non_zero_positions = np.nonzero(
                                coverage[0:flush_rows].sum(axis=1)
                                + coverage_indel_cat[0:flush_rows].sum(axis=1)
                            )
                            for pos in non_zero_positions[0].tolist():
                                current_pos = pos + reference_mat_start
                                if current_pos <= max_flushed_pos:
                                    continue  # see max_flushed_pos
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
                                total_coverage_indel_cat += coverage_indel_cat[pos]
                                total_unmasked_coverage_indel_cat += (
                                    unmasked_coverage_indel_cat[pos]
                                )
                        if chromNow != reference_mat_chrom:
                            max_flushed_pos = -1
                        reference_mat_chrom = chromNow
                        # batch_min_start <= every key's start in this batch.
                        reference_mat_start = batch_min_start
                        try:
                            region_end = region[2]
                        except IndexError:
                            region_end = 10e10
                        contig_end = contig_len_now
                        # Cap measured from reference_mat_start so the window
                        # never exceeds the coverage arrays' 1 Mb rows.
                        reference_mat_end = min(
                            reference_mat_start + 1000000,
                            max(
                                region_end, max([seq.reference_end for seq in readSet])
                            ),
                            contig_end,
                        )
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
                        # (3, window): unit_len, cut (run start), repeat_count;
                        # STR-priority merge of hp.h5 and str.h5.
                        hp_np = load_repeat_context(
                            hp_h5[reference_mat_chrom],
                            str_h5[reference_mat_chrom],
                            reference_mat_start,
                            reference_mat_end,
                        )
                        # Raw, unmerged slices, for learning and opportunity
                        # where HP and STR must be treated independently.
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
                            nm_mask,
                        ) = prepare_reference_mats(
                            reference_mat_chrom,
                            reference_mat_start,
                            reference_mat_end,
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
                        cov_rows = _coverage_window_rows(
                            coverage_leftovers,
                            coverage_leftover_gstart,
                            reference_mat_start,
                        )
                        coverage = np.zeros((cov_rows, 4))
                        coverage_indel_cat = np.zeros((cov_rows, 16))
                        unmasked_coverage = np.zeros((cov_rows, 4))
                        unmasked_coverage_indel_cat = np.zeros((cov_rows, 16))
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
                        neg_lr_sbs=neg_lr_sbs,
                        neg_lr_indel=neg_lr_indel,
                        power_tally=power_tally,
                        coverage_totals=coverage_totals,
                    ):
                        duplex_count += 1
            """
            Calling block ends
            """
            currentReadDict = {}
            currentStart = start
            if rec.template_length < 0 and rec.query_name in rugged_reads_index:
                # Re-check: this read was tested before the flushed batch
                # was indexed by _index_rugged_mates.
                rugged_reads_pool[rugged_reads_index.pop(rec.query_name)].append(rec)
            else:
                _place_read_in_dict(rec, label, currentReadDict)
            _drain_rugged_pool(
                chrom, currentStart, rugged_reads_pool, params, currentReadDict
            )
    """
    Calling block starts
    """
    _index_rugged_mates(currentReadDict, rugged_reads_index, rugged_reads_pool)
    # Final batch: same window floor as above. currentReadDict is empty if
    # the worker's regions had no reads.
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
            contig_len_now = tumorBam.get_reference_length(chromNow)
            if rs_reference_start >= contig_len_now:
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
                ### Flush the finished part of the window's coverage
                if "coverage" in locals():
                    if coverage_leftovers is not None:
                        _carry_in_coverage(
                            (
                                coverage,
                                coverage_indel_cat,
                                unmasked_coverage,
                                unmasked_coverage_indel_cat,
                            ),
                            coverage_leftovers,
                            coverage_leftover_gstart,
                            reference_mat_start,
                        )
                    (
                        flush_rows,
                        coverage_leftovers,
                        coverage_leftover_gstart,
                    ) = _cut_coverage_leftover(
                        (
                            coverage,
                            coverage_indel_cat,
                            unmasked_coverage,
                            unmasked_coverage_indel_cat,
                        ),
                        reference_mat_start,
                        batch_min_start,
                        chromNow == reference_mat_chrom,
                    )
                    non_zero_positions = np.nonzero(
                        coverage[0:flush_rows].sum(axis=1)
                        + coverage_indel_cat[0:flush_rows].sum(axis=1)
                    )
                    for pos in non_zero_positions[0].tolist():
                        current_pos = pos + reference_mat_start
                        if current_pos <= max_flushed_pos:
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
                        total_coverage_indel_cat += coverage_indel_cat[pos]
                        total_unmasked_coverage_indel_cat += (
                            unmasked_coverage_indel_cat[pos]
                        )
                if chromNow != reference_mat_chrom:
                    max_flushed_pos = -1
                reference_mat_chrom = chromNow
                reference_mat_start = batch_min_start
                try:
                    region_end = region[2]
                except IndexError:
                    region_end = 10e10
                contig_end = contig_len_now
                reference_mat_end = min(
                    reference_mat_start + 1000000,
                    max(region_end, max([seq.reference_end for seq in readSet])),
                    contig_end,
                )
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
                hp_np = load_repeat_context(
                    hp_h5[reference_mat_chrom],
                    str_h5[reference_mat_chrom],
                    reference_mat_start,
                    reference_mat_end,
                )
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
                    nm_mask,
                ) = prepare_reference_mats(
                    reference_mat_chrom,
                    reference_mat_start,
                    reference_mat_end,
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
                cov_rows = _coverage_window_rows(
                    coverage_leftovers, coverage_leftover_gstart, reference_mat_start
                )
                coverage = np.zeros((cov_rows, 4))
                coverage_indel_cat = np.zeros((cov_rows, 16))
                unmasked_coverage = np.zeros((cov_rows, 4))
                unmasked_coverage_indel_cat = np.zeros((cov_rows, 16))
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
                neg_lr_sbs=neg_lr_sbs,
                neg_lr_indel=neg_lr_indel,
                power_tally=power_tally,
                coverage_totals=coverage_totals,
            ):
                duplex_count += 1
    """
    Calling block ends
    """
    if isLearn:
        # Whole-number counts: cast to int so output formatting is the same
        # for -p 1 and multiprocess runs.
        return (
            mismatch_mat.astype(int),
            hp_alt_mat.astype(int),
            str_alt_mat.astype(int),
            mismatch_dmg_mat.astype(int),
            hp_dmg_mat.astype(int),
            str_dmg_mat.astype(int),
            sbs_alt_bq_hist_mat,
        )

    ### SNV depth extraction, batched (a few sequential pileup scans per
    ### BAM). Round 1 extracts PASS candidates only; round 2
    ### (coverage_only) only Caller.py's deferred set (round-1 masked
    ### candidates whose raw LR clears the final channel threshold).
    deferred_depth_keys = params.get("deferred_depth_keys", frozenset())
    mut_dict = dict()
    snv_candidate_keys = []
    # Ids of every founding family per candidate; their reads count toward
    # tumor depth regardless of base quality (see _collect_call_barcode).
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
        # Depth-check failures are labeled, not dropped: Caller.py prunes
        # them after merging round 2 unless --rescue (a round-2 record only
        # feeds that merge, so dropping it here would have no effect).
        if eligible:
            if ta == 0:
                mut["filter"] = "no_good_alt_read"
            elif ta / tdp > params["maxAF"]:
                mut["filter"] = "duplex_vaf"
            elif normalBams:
                if ndp < params["minNdepth"]:
                    mut["filter"] = "n_cov_mask"
                elif na / ndp > params["normalVAF"]:
                    mut["filter"] = "normal_vaf"
        mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
        mut_pass_filter.append(mut)
    if power_tally is not None and not params.get("no_masked_output"):
        _early_deferred_depth(
            muts, tumorBam, normalBams, params, extractDepthBatchSnv, indel=False
        )

    ### Indel depth extraction: same batching and eligibility as SNVs.
    muts_indels_dict = dict()
    indel_candidate_keys = []
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
            # ni: summed normal-BAM other-indel count, for the ni > 0 filter.
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
        # Labeled as for SNVs; ni > 0 (other-indel evidence in the normal)
        # is an unconditional drop, even with --rescue.
        if eligible:
            if ta == 0:
                mut["filter"] = "no_good_alt_read"
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
    if power_tally is not None and not params.get("no_masked_output"):
        _early_deferred_depth(
            muts_indels,
            tumorBam,
            normalBams,
            params,
            extractDepthBatchIndel,
            indel=True,
        )

    ### DBS depth extraction (round 1 only). _detect_dbs_pairs emits only
    ### PASS candidates, so there is no deferred set; failures of the
    ### maxAF/normal checks are dropped, ta == 0 is labeled as for SNVs.
    muts_dbs_pass_filter = []
    if not coverage_only:
        muts_dbs_dict = dict()
        dbs_candidate_keys = []
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
            if ta == 0:
                mut["filter"] = "no_good_alt_read"
            elif ta / tdp > params["maxAF"]:
                continue
            elif normalBams:
                if ndp < params["minNdepth"]:
                    continue
                if na / ndp > params["normalVAF"]:
                    continue
            mut["samples"] = [[ta, tr, tdp], [na, nr, ndp]]
            muts_dbs_pass_filter.append(mut)

    # Final flush of the remaining buffered coverage.
    _final_flush_t0 = time.time()
    if "coverage" in locals():
        if coverage_leftovers is not None:
            _carry_in_coverage(
                (
                    coverage,
                    coverage_indel_cat,
                    unmasked_coverage,
                    unmasked_coverage_indel_cat,
                ),
                coverage_leftovers,
                coverage_leftover_gstart,
                reference_mat_start,
            )
        non_zero_positions = np.nonzero(
            coverage.sum(axis=1) + coverage_indel_cat.sum(axis=1)
        )
        for pos in non_zero_positions[0].tolist():
            current_pos = pos + reference_mat_start
            if current_pos <= max_flushed_pos:
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
            total_coverage_indel_cat += coverage_indel_cat[pos]
            total_unmasked_coverage_indel_cat += unmasked_coverage_indel_cat[pos]

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

    # duplex_read_num_dict_trinuc values are raw (64, 4) matrices, not
    # folded by reverse complement (that happens at estimation time).
    return (
        mut_pass_filter,
        coverage_totals[0],
        recCount,
        duplex_count,
        duplex_read_num_dict,
        duplex_read_num_dict_trinuc,
        muts_indels_pass_filter,
        unique_read_num,
        pass_read_num,
        FPs,
        RPs,
        coverage_totals[1],
        total_coverage_indel_cat,
        total_unmasked_coverage_indel_cat,
        duplex_read_num_dict_indel,
        muts_dbs_pass_filter,
        duplex_read_num_dict_dbs,
        # Depth-composition totals and power tables, which Caller.py
        # combines into effective coverage (Eeff) at any channel threshold
        # without rescanning the BAM.
        depth_by_trinuc,
        depth_by_hpstr,
        L,
        L_indel_1bp,
        L_indel_len,
        [np.frombuffer(lrs, dtype=np.float32) for lrs in neg_lr_sbs],
        {k: np.frombuffer(v, dtype=np.float32) for k, v in neg_lr_indel.items()},
        # --skipCoveragePass only: power tallies + quality sample.
        power_tally.summary(all_quals) if power_tally is not None else None,
    )
