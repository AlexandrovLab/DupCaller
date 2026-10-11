import numpy as np
from .indels import (
    base_codes,
    STR_NCOL,
    STR_REF_COLS,
    str_event_col,
    str_length_bin,
    str_slip_tract,
    INDEL_ALT,
    INDEL_REF,
    findIndels,
    getIndelArr,
    indel_context_fits_window,
    indel_context_index,
    indel_has_context,
    indel_passes_mask,
    left_align_indel,
)

# BQ axis for the amp-error BQ histograms below: covers SAM/BAM's full
# Phred range (0-93); bin index == literal BQ value.
MAX_BQ = 93
NUM_BQ = MAX_BQ + 1


def _sbs_strand_counts(seq_mat, qual_mat, antimask, trinuc_int):
    """One strand's SBS amp-error tallies: the (96, 4) trinuc x base counts
    and the (96, 4, NUM_BQ) BQ histogram. Counts bases inside antimask with
    qual > 0 (the caller zeroes qual < minBq) and seq < 4."""
    seq_masked = seq_mat[:, antimask]
    qual_masked = qual_mat[:, antimask]
    alt_1d = trinuc_int[antimask][None, :] + seq_masked * 96
    valid = (qual_masked > 0) & (seq_masked < 4)
    count_mat = (
        np.bincount(
            alt_1d.ravel(), weights=valid.ravel().astype(float), minlength=96 * 4
        )[0 : 4 * 96]
        .reshape([4, 96])
        .T.astype(float)
    )
    bq_valid = np.clip(qual_masked[valid].astype(int), 0, MAX_BQ)
    flat_idx = alt_1d[valid] * NUM_BQ + bq_valid
    bq_hist = (
        np.bincount(flat_idx, minlength=4 * 96 * NUM_BQ)
        .reshape([4, 96, NUM_BQ])
        .transpose(1, 0, 2)
        .astype(float)
    )
    return count_mat, bq_hist


def _indel_informative_reads(covered, hq, span_ends, active=True):
    """Per span_ends table, a reads x positions mask: True where the read
    could show REF for an indel whose context base is p -- base at p passes
    minBq (hq) and the read aligns contiguously from the anchor p - 1
    through span_end[p], as getIndelArr requires. False for p < 2 and where
    span_end[p] is outside the window; all False when active is False."""
    m, n = covered.shape
    out = {table: np.zeros((m, n), dtype=bool) for table in span_ends}
    if not active or m == 0 or n == 0:
        return out
    # gaps[:, j]: uncovered positions in [0, j).
    gaps = np.zeros((m, n + 1), dtype=np.int64)
    gaps[:, 1:] = np.cumsum(~covered, axis=1)
    p = np.arange(n)
    for table, span_end in span_ends.items():
        ok = (p >= 2) & (span_end < n)
        lo = p[ok] - 1
        hi = span_end[ok]
        out[table][:, ok] = (gaps[:, hi + 1] - gaps[:, lo] == 0) & hq[:, ok]
    return out


def profileTriNucMismatches(
    seqs, reference_start, reference_int, trinuc_int, hp_raw, str_raw, antimask, params
):
    reverse_comp = [1, 0, 3, 2]
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    base_changes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    chrom = seqs[0].reference_name
    trinuc2num = params["trinuc2num_dict"]
    # hp_*_count (10, 12): rows = run length 1-10+, cols = ref_base*3 +
    # (idLen + 1), idLen in {-1, 0, 1}; idLen=0 is the opportunity column.
    # str_*_count (5, STR_NCOL): rows = STR length bin, cols = a (del, ref,
    # ins) triple per unit group (indels.str_col). prob.indelErrorProbs
    # reads the same cells.
    hp_alt_count = np.zeros([10, 12])
    hp_dmg_count = np.zeros([10, 12])
    str_alt_count = np.zeros([5, STR_NCOL])
    str_dmg_count = np.zeros([5, STR_NCOL])
    # SBS amp-error counts per (trinuc, base, BQ); feeds estimate_sbs_srd_rates.
    sbs_alt_bq_hist = np.zeros([64, 4, NUM_BQ])

    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)

    m_F1R2 = len(F1R2)
    m_F2R1 = len(F2R1)

    # Min reads per strand for the SBS SRD (amp) and SSM (damage) tracks; a
    # family too small for a track gets an all-False antimask for it.
    srd_min_read = params.get("srdMinRead", 3)
    ssm_min_read = params.get("ssmMinRead", 3)
    # Per-position depth and summed-BQ floors (SBS damage and indel antimasks).
    min_depth = max(params.get("minRef", 3), params.get("minAlt", 3))
    min_alt_qual = params.get("minAltQual", 90)
    # Indel (HP/STR) min reads per strand, independent of the SBS thresholds.
    min_group_indel_amp = 3
    min_group_indel_dmg = 3

    # Too small for every track (amp needs one strand, damage both).
    if max(m_F1R2, m_F2R1) < min(srd_min_read, min_group_indel_amp) and min(
        m_F1R2, m_F2R1
    ) < min(ssm_min_read, min_group_indel_dmg):
        return (
            np.zeros([64, 4]),
            np.zeros([10, 12]),
            np.zeros([5, STR_NCOL]),
            np.zeros([64, 4]),
            np.zeros([10, 12]),
            np.zeros([5, STR_NCOL]),
            np.zeros([64, 4, NUM_BQ]),
        )

    ### Prepare sequence matrix and quality matrix for each strand
    n = len(reference_int)
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}
    # 4 (no base) wherever a read has no aligned base: before its start,
    # deletions, past its end. A zero would read as an "A".
    F1R2_seq_mat = np.full([m_F1R2, n], 4, dtype=int)  # reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2, n])
    F2R1_seq_mat = np.full([m_F2R1, n], 4, dtype=int)  # reads x pos
    F2R1_qual_mat = np.zeros([m_F2R1, n])
    for mm, seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = base_codes(seq.query_alignment_sequence)
        cigartuples = seq.cigartuples
        current_seq_ind = 0
        current_mat_ind = seq.reference_start - reference_start
        reference_ind = seq.reference_start - reference_start
        ref_length_plus_del = seq.reference_length
        for ct in cigartuples:
            if ct[0] == 0:
                F1R2_seq_mat[mm, current_mat_ind : current_mat_ind + ct[1]] = sequence[
                    current_seq_ind : current_seq_ind + ct[1]
                ]
                F1R2_qual_mat[
                    mm, current_mat_ind : current_mat_ind + ct[1]
                ] = qualities[current_seq_ind : current_seq_ind + ct[1]]
                current_seq_ind += ct[1]
                reference_ind += ct[1]
                current_mat_ind += ct[1]
            elif ct[0] == 1:
                current_seq_ind += ct[1]
            elif ct[0] == 2:
                F1R2_seq_mat[mm, current_mat_ind : current_mat_ind + ct[1]] = 4
                F1R2_qual_mat[mm, current_mat_ind : current_mat_ind + ct[1]] = 0
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F1R2_seq_mat[mm, current_mat_ind:n] = 4
        F1R2_qual_mat[mm, current_mat_ind:n] = 0
    for mm, seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = base_codes(seq.query_alignment_sequence)
        cigartuples = seq.cigartuples
        current_seq_ind = 0
        current_mat_ind = seq.reference_start - reference_start
        reference_ind = seq.reference_start - reference_start
        ref_length_plus_del = seq.reference_length
        for ct in cigartuples:
            if ct[0] == 0:
                F2R1_seq_mat[mm, current_mat_ind : current_mat_ind + ct[1]] = sequence[
                    current_seq_ind : current_seq_ind + ct[1]
                ]
                F2R1_qual_mat[
                    mm, current_mat_ind : current_mat_ind + ct[1]
                ] = qualities[current_seq_ind : current_seq_ind + ct[1]]
                current_seq_ind += ct[1]
                reference_ind += ct[1]
                current_mat_ind += ct[1]
            elif ct[0] == 1:
                current_seq_ind += ct[1]
            elif ct[0] == 2:
                F2R1_seq_mat[mm, current_mat_ind : current_mat_ind + ct[1]] = 4
                F2R1_qual_mat[mm, current_mat_ind : current_mat_ind + ct[1]] = 0
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F2R1_seq_mat[mm, current_mat_ind:n] = 4
        F2R1_qual_mat[mm, current_mat_ind:n] = 0

    F1R2_qual_mat[F1R2_qual_mat < params["minBq"]] = 0
    F2R1_qual_mat[F2R1_qual_mat < params["minBq"]] = 0

    F1R2_antimask = antimask.copy()
    F2R1_antimask = antimask.copy()
    # SRD: each strand's inclusion is independent.
    if m_F1R2 < srd_min_read:
        F1R2_antimask[:] = False
    if m_F2R1 < srd_min_read:
        F2R1_antimask[:] = False

    dmg_antimask = antimask.copy()
    # SSM/damage: requires BOTH strands to meet ssm_min_read.
    if m_F1R2 < ssm_min_read or m_F2R1 < ssm_min_read:
        dmg_antimask[:] = False

    F1R2_qual_mat_merged = np.zeros([4, n])
    F1R2_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F1R2_qual_mat_merged[nn, :] = F1R2_qual_mat.sum(
            axis=0, where=(F1R2_seq_mat == nn)
        )
        F1R2_count_mat[nn, :] = (F1R2_seq_mat == nn).sum(axis=0)

    F2R1_qual_mat_merged = np.zeros([4, n])
    F2R1_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F2R1_qual_mat_merged[nn, :] = F2R1_qual_mat.sum(
            axis=0, where=(F2R1_seq_mat == nn)
        )
        F2R1_count_mat[nn, :] = (F2R1_seq_mat == nn).sum(axis=0)

    dmg_antimask[
        np.logical_or(
            (F1R2_count_mat != 0).sum(axis=0) > 1, (F2R1_count_mat != 0).sum(axis=0) > 1
        )
    ] = False  # mask strand discordant
    dmg_antimask[
        np.logical_or(
            (F1R2_count_mat).sum(axis=0) < min_depth,
            (F2R1_count_mat).sum(axis=0) < min_depth,
        )
    ] = False  # mask locations with depth < min_depth
    dmg_antimask[
        np.logical_or(
            (F1R2_qual_mat_merged).sum(axis=0) < min_alt_qual,
            (F2R1_qual_mat_merged).sum(axis=0) < min_alt_qual,
        )
    ] = False  # mask locations where summed consensus quality is less than min_alt_qual

    F1R2_alleles = np.argmax(F1R2_count_mat, axis=0)
    F2R1_alleles = np.argmax(F2R1_count_mat, axis=0)
    ds_alt = np.logical_and(
        F1R2_alleles != reference_int, F2R1_alleles != reference_int
    )  # location where both strands are alt allele
    dmg_antimask[ds_alt] = False

    F1R2_dmg_alt = reference_int.copy()
    F2R1_dmg_alt = reference_int.copy()
    F1R2_dmg_alt[F1R2_alleles != reference_int] = F1R2_alleles[
        F1R2_alleles != reference_int
    ]
    F2R1_dmg_alt[F2R1_alleles != reference_int] = F2R1_alleles[
        F2R1_alleles != reference_int
    ]
    F1R2_dmg_trinuc_masked = trinuc_int[dmg_antimask]
    F2R1_dmg_trinuc_masked = trinuc_int[dmg_antimask]
    F1R2_dmg_alt_masked = F1R2_dmg_alt[dmg_antimask]
    F2R1_dmg_alt_masked = F2R1_dmg_alt[dmg_antimask]
    F1R2_dmg_trinuc_alt_1Dmap = F1R2_dmg_trinuc_masked + F1R2_dmg_alt_masked * 96
    F2R1_dmg_trinuc_alt_1Dmap = F2R1_dmg_trinuc_masked + F2R1_dmg_alt_masked * 96
    F1R2_dmg_trinuc_alt_count_mat = (
        np.bincount(F1R2_dmg_trinuc_alt_1Dmap, minlength=96 * 4).reshape([4, 96]).T
    )
    F2R1_dmg_trinuc_alt_count_mat = (
        np.bincount(F2R1_dmg_trinuc_alt_1Dmap, minlength=96 * 4).reshape([4, 96]).T
    )
    dmg_trinuc_alt_count_mat_norm = F1R2_dmg_trinuc_alt_count_mat[0:64, :] + np.vstack(
        [
            F2R1_dmg_trinuc_alt_count_mat[32:64, [1, 0, 3, 2]],
            F2R1_dmg_trinuc_alt_count_mat[:32, [1, 0, 3, 2]],
        ]
    )

    # BQ-passing counts (qual > 0 after the minBq zeroing); these, not the
    # raw count matrices, gate SRD sites.
    F1R2_hq_count_mat = np.zeros([4, n], dtype=int)
    F2R1_hq_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F1R2_hq_count_mat[nn, :] = ((F1R2_seq_mat == nn) & (F1R2_qual_mat > 0)).sum(
            axis=0
        )
        F2R1_hq_count_mat[nn, :] = ((F2R1_seq_mat == nn) & (F2R1_qual_mat > 0)).sum(
            axis=0
        )
    F1R2_hq_count_sum = F1R2_hq_count_mat.sum(axis=0)
    F2R1_hq_count_sum = F2R1_hq_count_mat.sum(axis=0)
    F1R2_hq_ref_count = F1R2_hq_count_mat[reference_int, np.ogrid[: reference_int.size]]
    F2R1_hq_ref_count = F2R1_hq_count_mat[reference_int, np.ogrid[: reference_int.size]]

    F1R2_antimask[ds_alt] = False
    # SRD site: >= srd_min_read BQ-passing bases, < 2 of them non-reference,
    # and not a single non-reference allele.
    F1R2_antimask[F1R2_hq_count_sum < srd_min_read] = False
    F1R2_antimask[F1R2_hq_count_sum - F1R2_hq_ref_count >= 2] = False
    F1R2_antimask[
        np.logical_and((F1R2_hq_count_mat >= 1).sum(axis=0) < 2, F1R2_hq_ref_count == 0)
    ] = False

    F2R1_antimask[ds_alt] = False
    F2R1_antimask[F2R1_hq_count_sum < srd_min_read] = False
    F2R1_antimask[F2R1_hq_count_sum - F2R1_hq_ref_count >= 2] = False
    F2R1_antimask[
        np.logical_and((F2R1_hq_count_mat >= 1).sum(axis=0) < 2, F2R1_hq_ref_count == 0)
    ] = False

    # Every covered base counts regardless of the read's mismatch total;
    # call.py's NM blacklist drops whole families, not individual reads.
    F1R2_trinuc_alt_count_mat, F1R2_trinuc_alt_bq_hist = _sbs_strand_counts(
        F1R2_seq_mat, F1R2_qual_mat, F1R2_antimask, trinuc_int
    )
    F1R2_trinuc_alt_count_mat_norm = F1R2_trinuc_alt_count_mat[0:64, :] + np.vstack(
        [
            F1R2_trinuc_alt_count_mat[32:64, [1, 0, 3, 2]],
            F1R2_trinuc_alt_count_mat[:32, [1, 0, 3, 2]],
        ]
    )
    F1R2_trinuc_alt_bq_hist_norm = F1R2_trinuc_alt_bq_hist[0:64, :, :] + np.vstack(
        [
            F1R2_trinuc_alt_bq_hist[32:64, [1, 0, 3, 2], :],
            F1R2_trinuc_alt_bq_hist[:32, [1, 0, 3, 2], :],
        ]
    )

    F2R1_trinuc_alt_count_mat, F2R1_trinuc_alt_bq_hist = _sbs_strand_counts(
        F2R1_seq_mat, F2R1_qual_mat, F2R1_antimask, trinuc_int
    )
    F2R1_trinuc_alt_count_mat_norm = F2R1_trinuc_alt_count_mat[0:64, :] + np.vstack(
        [
            F2R1_trinuc_alt_count_mat[32:64, [1, 0, 3, 2]],
            F2R1_trinuc_alt_count_mat[:32, [1, 0, 3, 2]],
        ]
    )
    F2R1_trinuc_alt_bq_hist_norm = F2R1_trinuc_alt_bq_hist[0:64, :, :] + np.vstack(
        [
            F2R1_trinuc_alt_bq_hist[32:64, [1, 0, 3, 2], :],
            F2R1_trinuc_alt_bq_hist[:32, [1, 0, 3, 2], :],
        ]
    )
    sbs_alt_bq_hist = F1R2_trinuc_alt_bq_hist_norm + F2R1_trinuc_alt_bq_hist_norm

    ###INDEL LEARN
    indels = set()
    for seq in F1R2:
        current_indels = findIndels(seq)
        if len(current_indels) > 1:
            continue
        indels.update(
            left_align_indel(i, reference_int, reference_start) for i in current_indels
        )
    for seq in F2R1:
        current_indels = findIndels(seq)
        if len(current_indels) > 1:
            continue
        indels.update(
            left_align_indel(i, reference_int, reference_start) for i in current_indels
        )
    # antimask/hp_raw/str_raw/reference_int all start at reference_start
    # (the family's minimum read start), which seqs[0] need not be.
    start = reference_start
    indels = list(indels)
    indels_masked = list()
    for indel in indels:
        refPos = int(indel.split(":")[0])
        indelLen = int(indel.split(":")[1])
        # Same locus mask as genotypeDSIndel at call time.
        if (
            indel_passes_mask(antimask, refPos - start, indelLen)
            and indel_has_context(refPos - start, len(reference_int))
            and indel_context_fits_window(
                indel, refPos - start, reference_int, hp_raw, str_raw
            )
        ):
            indels_masked.append(indel)

    m = len(indels_masked)
    if m >= 2:
        return (
            F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,
            np.zeros([10, 12]),
            np.zeros([5, STR_NCOL]),
            dmg_trinuc_alt_count_mat_norm,
            np.zeros([10, 12]),
            np.zeros([5, STR_NCOL]),
            sbs_alt_bq_hist,
        )
    # Amp opportunity: each strand gated on its own read count, like the
    # per-strand amp events below.
    F1R2_opp_antimask = antimask.copy()
    F2R1_opp_antimask = antimask.copy()
    if m_F1R2 < min_group_indel_amp:
        F1R2_opp_antimask[:] = False
    if m_F2R1 < min_group_indel_amp:
        F2R1_opp_antimask[:] = False
    dmg_opp_antimask = antimask.copy()
    if m_F1R2 < min_group_indel_dmg or m_F2R1 < min_group_indel_dmg:
        dmg_opp_antimask[:] = False
    # No candidate's context base sits at the window's first two positions
    # (indel_has_context); the amp weights are 0 there too.
    dmg_opp_antimask[:2] = False

    # hp_raw (hp.h5): row0 = homopolymer run length, row1 = start-of-run.
    # str_raw (str.h5): row0 = unit length, row1 = repeat count, row2 =
    # start-of-repeat. HP and STR are credited independently: a run start
    # inside an STR (the "AA" in (AAT)n) counts toward both tables.
    hp_len_arr = hp_raw[0].astype(int)
    hp_mask = hp_raw[1] == 1

    unit_len = str_raw[0].astype(int)
    repeat_count = str_raw[1].astype(int)
    total_len = unit_len * repeat_count
    is_str = unit_len >= 2
    # STR bins: 0 = not a repeat, 1 = 2-9bp, 2 = 10-24bp, 3 = 25-39bp,
    # 4 = 40+bp.
    str_bin_arr = np.zeros_like(total_len)
    str_bin_arr[is_str] = 1
    str_bin_arr[is_str & (total_len >= 10)] = 2
    str_bin_arr[is_str & (total_len >= 25)] = 3
    str_bin_arr[is_str & (total_len >= 40)] = 4
    str_mask = (str_raw[2] == 1) & is_str
    # Unit group (index 0-4 = u1-u5+) of each tract start's reference column.
    str_group_arr = np.minimum(np.maximum(unit_len, 1), 5) - 1

    hp_rc4 = [1, 0, 3, 2]  # base-complement permutation, 4-wide axis

    # Amp opportunity weight per position = the strand's reads that could
    # show REF there, spanning the HP run, the STR tract (rows 1-4), or just
    # the base (STR row 0).
    pos_idx = np.arange(n)
    F1R2_covered = F1R2_seq_mat < 4
    F2R1_covered = F2R1_seq_mat < 4
    F1R2_hq = (F1R2_qual_mat > 0) & (F1R2_seq_mat < 4)
    F2R1_hq = (F2R1_qual_mat > 0) & (F2R1_seq_mat < 4)
    span_ends = {
        "hp": pos_idx + hp_len_arr,
        "str": pos_idx + total_len,
        "str0": pos_idx,
    }
    # A strand that fails its opportunity gate credits nothing anywhere.
    F1R2_inf = _indel_informative_reads(
        F1R2_covered, F1R2_hq, span_ends, F1R2_opp_antimask.any()
    )
    F2R1_inf = _indel_informative_reads(
        F2R1_covered, F2R1_hq, span_ends, F2R1_opp_antimask.any()
    )

    # HP amp opportunity: each run start adds to its (run length, base)
    # idLen=0 cell, weighted by informative reads and RC-folded. N skipped.
    def _hp_amp_opportunity(strand_antimask, weight):
        sel = strand_antimask & hp_mask
        hp_here = hp_len_arr[sel]
        ref_here = reference_int[sel]
        w_here = weight[sel]
        valid = ref_here <= 3
        hp_here = np.minimum(hp_here[valid], 10)
        ref_here = ref_here[valid]
        local = np.zeros([10, 4])
        np.add.at(local, (hp_here - 1, ref_here), w_here[valid])
        return local + local[:, hp_rc4]

    hp_alt_count[:, [1, 4, 7, 10]] += _hp_amp_opportunity(
        F1R2_opp_antimask, F1R2_inf["hp"].sum(axis=0)
    )
    hp_alt_count[:, [1, 4, 7, 10]] += _hp_amp_opportunity(
        F2R1_opp_antimask, F2R1_inf["hp"].sum(axis=0)
    )

    # STR amp opportunity: repeat start positions, the reference column of
    # the tract's unit group.
    def _str_amp_opportunity(strand_antimask, weight):
        sel = strand_antimask & str_mask
        local = np.zeros([5, 5])
        np.add.at(local, (str_bin_arr[sel], str_group_arr[sel]), weight[sel])
        return local

    str_alt_count[:, STR_REF_COLS] += _str_amp_opportunity(
        F1R2_opp_antimask, F1R2_inf["str"].sum(axis=0)
    )
    str_alt_count[:, STR_REF_COLS] += _str_amp_opportunity(
        F2R1_opp_antimask, F2R1_inf["str"].sum(axis=0)
    )

    # STR row 0 (not a repeat) opportunity: every unmasked position,
    # including STR ones, in every length group.
    def _str_amp_opportunity_row0(strand_antimask, weight):
        return weight[strand_antimask].sum()

    str_alt_count[0, STR_REF_COLS] += _str_amp_opportunity_row0(
        F1R2_opp_antimask, F1R2_inf["str0"].sum(axis=0)
    )
    str_alt_count[0, STR_REF_COLS] += _str_amp_opportunity_row0(
        F2R1_opp_antimask, F2R1_inf["str0"].sum(axis=0)
    )

    # HP dmg opportunity, credited to both orientations.
    hp_dmg_here = hp_len_arr[dmg_opp_antimask][hp_mask[dmg_opp_antimask]]
    ref_dmg_here = reference_int[dmg_opp_antimask][hp_mask[dmg_opp_antimask]]
    valid_dmg = ref_dmg_here <= 3
    hp_dmg_here = np.minimum(hp_dmg_here[valid_dmg], 10)
    ref_dmg_here = ref_dmg_here[valid_dmg]
    hp_dmg_local = np.zeros([10, 4])
    np.add.at(hp_dmg_local, (hp_dmg_here - 1, ref_dmg_here), 1)
    hp_dmg_count[:, [1, 4, 7, 10]] += hp_dmg_local + hp_dmg_local[:, hp_rc4]

    # STR dmg opportunity: *2 for the two orientations.
    dmg_str_sel = dmg_opp_antimask & str_mask
    str_dmg_local = np.zeros([5, 5])
    np.add.at(str_dmg_local, (str_bin_arr[dmg_str_sel], str_group_arr[dmg_str_sel]), 2)
    str_dmg_count[:, STR_REF_COLS] += str_dmg_local

    # STR row 0 dmg opportunity: every unmasked position, *2.
    str_dmg_count[0, STR_REF_COLS] += np.count_nonzero(dmg_opp_antimask) * 2

    if m == 0:
        return (
            F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,
            hp_alt_count,
            str_alt_count,
            dmg_trinuc_alt_count_mat_norm,
            hp_dmg_count,
            str_dmg_count,
            sbs_alt_bq_hist,
        )

    F1R2_alt_count = np.zeros(m)
    F1R2_ref_count = np.zeros(m)
    F2R1_alt_count = np.zeros(m)
    F2R1_ref_count = np.zeros(m)
    ctx = indel_context_index(int(indels_masked[0].split(":")[0]) - start)
    # Per-read ALT flags: only ALT reads credited to the idLen=0 column are
    # moved out of it (_alt_in_opp).
    F1R2_is_alt = np.zeros(len(F1R2), dtype=bool)
    F2R1_is_alt = np.zeros(len(F2R1), dtype=bool)
    for mm, seq in enumerate(F1R2):
        seqArr, _ = getIndelArr(
            seq, indels_masked, params["minBq"], reference_int, reference_start
        )
        F1R2_alt_count += np.count_nonzero(seqArr == INDEL_ALT)
        F1R2_ref_count += np.count_nonzero(seqArr == INDEL_REF)
        F1R2_is_alt[mm] = seqArr[0] == INDEL_ALT
    for mm, seq in enumerate(F2R1):
        seqArr, _ = getIndelArr(
            seq, indels_masked, params["minBq"], reference_int, reference_start
        )
        F2R1_alt_count += np.count_nonzero(seqArr == INDEL_ALT)
        F2R1_ref_count += np.count_nonzero(seqArr == INDEL_REF)
        F2R1_is_alt[mm] = seqArr[0] == INDEL_ALT

    def _alt_in_opp(is_alt, inf, table):
        """ALT reads counted in the table's opportunity weight at ctx."""
        return int(np.count_nonzero(is_alt & inf[table][:, ctx]))

    def _ref_in_opp(is_alt, inf, table):
        """Non-ALT reads counted in the table's opportunity weight at ctx; an
        amp event is booked only where the strand credited at least one."""
        return int(np.count_nonzero(~is_alt & inf[table][:, ctx]))

    dmg_antimask = np.ones(m, dtype=bool)
    dmg_antimask[
        np.logical_or(
            (F1R2_ref_count + F1R2_alt_count) < min_depth,
            (F2R1_ref_count + F2R1_alt_count) < min_depth,
        )
    ] = False
    dmg_antimask[np.logical_and(F1R2_ref_count != 0, F1R2_alt_count != 0)] = False
    dmg_antimask[np.logical_and(F2R1_ref_count != 0, F2R1_alt_count != 0)] = False
    dmg_antimask[np.logical_and(F1R2_alt_count > 0, F2R1_alt_count > 0)] = False
    F1R2_dmg_antimask = dmg_antimask.copy()
    F1R2_dmg_antimask[F1R2_alt_count == 0] = False

    F2R1_dmg_antimask = dmg_antimask.copy()
    F2R1_dmg_antimask[F2R1_alt_count == 0] = False

    F1R2_antimask = np.ones(m, dtype=bool)
    F2R1_antimask = np.ones(m, dtype=bool)

    F1R2_antimask[F1R2_ref_count == 0] = False
    F1R2_antimask[F1R2_ref_count + F1R2_alt_count < min_depth] = False
    F1R2_antimask[F1R2_alt_count > 1] = False

    F2R1_antimask[F2R1_ref_count == 0] = False
    F2R1_antimask[F2R1_ref_count + F2R1_alt_count < min_depth] = False
    F2R1_antimask[F2R1_alt_count > 1] = False

    def _book_str_event(row, col, mm):
        """Book an STR-table event at ctx, only where the opportunity pass
        credited ctx (row 0: any unmasked position; rows 1-4: a tract
        start). Amp moves the credited ALT reads out of the group's
        reference column; damage needs ALT on exactly one strand."""
        table = "str" if row >= 1 else "str0"
        ctx_ok = row == 0 or bool(str_mask[ctx])
        if not ctx_ok:
            return
        ref_col = col - col % 3 + 1
        if (
            F1R2_antimask[mm]
            and F1R2_opp_antimask[ctx]
            and _ref_in_opp(F1R2_is_alt, F1R2_inf, table)
        ):
            str_alt_count[row, col] += F1R2_alt_count[mm]
            str_alt_count[row, ref_col] -= _alt_in_opp(F1R2_is_alt, F1R2_inf, table)
        if (
            F2R1_antimask[mm]
            and F2R1_opp_antimask[ctx]
            and _ref_in_opp(F2R1_is_alt, F2R1_inf, table)
        ):
            str_alt_count[row, col] += F2R1_alt_count[mm]
            str_alt_count[row, ref_col] -= _alt_in_opp(F2R1_is_alt, F2R1_inf, table)
        if dmg_opp_antimask[ctx] and (F1R2_dmg_antimask[mm] or F2R1_dmg_antimask[mm]):
            str_dmg_count[row, col] += 1
            str_dmg_count[row, ref_col] -= 1

    # m == 1 here (m >= 2 and m == 0 returned above).
    for mm, indel in enumerate(indels_masked):
        parts = indel.split(":")
        pos = int(parts[0]) - start
        indelLen = int(parts[1])
        # Same HP/STR context position as genotypeDSIndel.
        anchor = indel_context_index(pos)
        hp = int(hp_raw[0, anchor])
        hp_capped = min(hp, 10)
        # Same STR0 rule as genotypeDSIndel: a non-slip >=2bp indel is row 0.
        slip = str_slip_tract(
            indelLen,
            parts[2] if indelLen > 0 and len(parts) > 2 else "",
            anchor,
            reference_int,
            str_raw,
        )
        str_bin_here = 0 if slip is None else str_length_bin(*slip)
        unit_here = 0 if slip is None else slip[0]

        if indelLen > 5:
            indelLen = 5
        if indelLen < -5:
            indelLen = -5

        if indelLen == 1 or indelLen == -1:
            ref_allele = int(reference_int[anchor])
            ref_allele_rc = int(reverse_comp[ref_allele])
            if indelLen == -1:
                # The deleted base is the run's own base.
                hp_match = True
            else:
                inserted_seq = parts[2] if len(parts) > 2 else ""
                inserted_base = (
                    base2num.get(inserted_seq[0], -1) if inserted_seq else -1
                )
                hp_match = inserted_base == ref_allele

            if hp_match:
                row = hp_capped - 1
                col = ref_allele * 3 + (indelLen + 1)
                col_rc = ref_allele_rc * 3 + (indelLen + 1)
                opp_col = ref_allele * 3 + 1
                opp_col_rc = ref_allele_rc * 3 + 1
                # Book only where the opportunity pass credited ctx.
                hp_ctx = bool(hp_mask[ctx]) and ref_allele <= 3

                # Amp: move the ALT reads credited at ctx from the idLen=0
                # cell to the event's cell, then RC-fold per strand.
                F1R2_local = np.zeros(12)
                F2R1_local = np.zeros(12)
                if (
                    F1R2_antimask[mm]
                    and hp_ctx
                    and F1R2_opp_antimask[ctx]
                    and _ref_in_opp(F1R2_is_alt, F1R2_inf, "hp")
                ):
                    F1R2_local[col] += F1R2_alt_count[mm]
                    F1R2_local[opp_col] -= _alt_in_opp(F1R2_is_alt, F1R2_inf, "hp")
                if (
                    F2R1_antimask[mm]
                    and hp_ctx
                    and F2R1_opp_antimask[ctx]
                    and _ref_in_opp(F2R1_is_alt, F2R1_inf, "hp")
                ):
                    F2R1_local[col_rc] += F2R1_alt_count[mm]
                    F2R1_local[opp_col_rc] -= _alt_in_opp(F2R1_is_alt, F2R1_inf, "hp")
                F1R2_local = F1R2_local.reshape(4, 3)
                F1R2_local = F1R2_local + F1R2_local[hp_rc4, :]
                F2R1_local = F2R1_local.reshape(4, 3)
                F2R1_local = F2R1_local + F2R1_local[hp_rc4, :]
                hp_alt_count[row, :] += (F1R2_local + F2R1_local).reshape(12)

                # Dmg: F2R1 folds in via the complemented base.
                if hp_ctx and dmg_opp_antimask[ctx]:
                    if F1R2_dmg_antimask[mm]:
                        hp_dmg_count[row, col] += 1
                        hp_dmg_count[row, opp_col] -= 1
                    if F2R1_dmg_antimask[mm]:
                        hp_dmg_count[row, col_rc] += 1
                        hp_dmg_count[row, opp_col_rc] -= 1
            else:
                # Mismatched 1bp insertion: str.txt row 0, u1 insertion.
                _book_str_event(0, str_event_col(0, 0, indelLen), mm)
        else:
            # >=2bp: row = slip tract's length bin (0 if not a slip), group =
            # the tract's unit (row 0: the indel length's group).
            _book_str_event(
                str_bin_here, str_event_col(str_bin_here, unit_here, indelLen), mm
            )

    return (
        F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,
        hp_alt_count,
        str_alt_count,
        dmg_trinuc_alt_count_mat_norm,
        hp_dmg_count,
        str_dmg_count,
        sbs_alt_bq_hist,
    )


def estimate_sbs_srd_rates(
    sbs_alt_bq_hist, pseudocount, max_iter=100, tol=1e-12, fallback=False
):
    """EM-estimate the per-trinuc SBS single-read-damage (SRD) rate matrix
    from sbs_alt_bq_hist.

    Per row, a read at base quality BQ (e = 10**(-BQ/10)) is either a true
    conversion to alt b (rate p_b) read correctly (1-e), or the reference
    (rate p = 1 - sum p_b) miscalled to b (e/3).
        E: w_b = p_b(1-e) / (p_b(1-e) + p e/3),  w_r = p_b e/3
           N_b = sum hist_b*w_b + sum hist_ref*w_r
        M: p_b = (N_b + a) / (N + 4a)   (a = pseudocount, N = row total)

    Returns (64, 4) in build_trinuc64_order rows and ATCG columns (ref
    column = p); rows with no observations get 1/4 everywhere.
    fallback=True fills sparse rows' unobserved alt cells from the bundled
    fallback profile (misc.apply_sbs_low_coverage_fallback); off by default
    so AggregateProfile, which builds those profiles, fits from data only.
    """
    from .misc import build_trinuc64_order

    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    _, num2trinuc = build_trinuc64_order()
    n_bq = sbs_alt_bq_hist.shape[2]
    bq_values = np.arange(n_bq)
    e = 10 ** (-bq_values / 10)
    one_minus_e = 1 - e
    e_over_3 = e / 3

    srd = np.zeros([64, 4])
    for row in range(64):
        ref_col = base2num[num2trinuc[row][1]]
        alt_cols = [c for c in range(4) if c != ref_col]
        hist_r = sbs_alt_bq_hist[row, ref_col, :]
        hist_alt = {c: sbs_alt_bq_hist[row, c, :] for c in alt_cols}
        N = float(sbs_alt_bq_hist[row, :, :].sum())
        if N == 0:
            # Avoids 0/0 in the init; the M step would give 1/4 anyway.
            srd[row, :] = 0.25
            continue

        # Init from the naive (unweighted) empirical fraction.
        p_b = {c: float(hist_alt[c].sum()) / N for c in alt_cols}
        for _ in range(max_iter):
            p_alt_total = sum(p_b.values())
            p_ref = 1 - p_alt_total
            new_p_b = {}
            for c in alt_cols:
                denom = p_b[c] * one_minus_e + p_ref * e_over_3
                w_b = np.divide(
                    p_b[c] * one_minus_e,
                    denom,
                    out=np.zeros_like(denom),
                    where=denom > 0,
                )
                w_r = p_b[c] * e_over_3
                N_c = float(np.dot(hist_alt[c], w_b)) + float(np.dot(hist_r, w_r))
                new_p_b[c] = (N_c + pseudocount) / (N + 4 * pseudocount)
            delta = max(abs(new_p_b[c] - p_b[c]) for c in alt_cols)
            p_b = new_p_b
            if delta < tol:
                break

        srd[row, ref_col] = 1 - sum(p_b.values())
        for c in alt_cols:
            srd[row, c] = p_b[c]
    if fallback:
        from .misc import _read_fallback_counts, apply_sbs_low_coverage_fallback

        srd = apply_sbs_low_coverage_fallback(
            srd,
            sbs_alt_bq_hist.sum(axis=2),
            _read_fallback_counts(".amp.tn.srd.txt"),
        )
    return srd
