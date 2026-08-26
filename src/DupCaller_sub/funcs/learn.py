import numpy as np
from .indels import findIndels, getIndelArr, left_align_indel


def profileTriNucMismatches(
    seqs, reference_start, reference_int, trinuc_int, hp_raw, str_raw, antimask, params
):
    # fasta = params["reference"]
    reverse_comp = [1, 0, 3, 2]
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    num2base = "ATCG"
    base_changes = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    chrom = seqs[0].reference_name
    start = seqs[0].reference_start
    trinuc2num = params["trinuc2num_dict"]
    # hp_alt_count/hp_dmg_count: (10, 12) -- rows hp run length 1-10+
    # (capped), columns ref_allele*3+(idLen+1) for idLen in {-1,0,1}
    # (idLen=0 is the reference/opportunity column). str_alt_count/
    # str_dmg_count: (5, 11) -- rows STR-length bin 0="0-1" (not a real
    # repeat) / 1="2-9" / 2="10-24" / 3="25-39" / 4="40+", columns
    # idLen+5 for idLen in -5..5 (idLen=0 the opportunity column). See
    # funcs/prob.py's indelErrorProbs for the matching call-time
    # selection logic these are built to feed.
    hp_alt_count = np.zeros([10, 12])
    hp_dmg_count = np.zeros([10, 12])
    str_alt_count = np.zeros([5, 11])
    str_dmg_count = np.zeros([5, 11])

    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)

    ### Determine match length

    m_F1R2 = len(F1R2)
    m_F2R1 = len(F2R1)

    # Minimum reads-per-strand for a family to contribute to amp/damage
    # learning, independently adjustable since damage (single-strand-
    # specific errors) and amp (errors shared by both consensus reads)
    # need different amounts of depth to reliably separate signal from
    # noise. The family-level floor below (below which NOTHING can be
    # learned from this family, amp or damage) is the smaller of the two --
    # each learning track's own, possibly-stricter minimum is applied
    # separately at its own antimask further down, so a family passing
    # only the lower of the two still contributes to whichever track it
    # actually qualifies for.
    min_group_amp = params.get("minGroupAmp", 3)
    min_group_dmg = params.get("minGroupDmg", 3)
    if m_F1R2 < min(min_group_amp, min_group_dmg) or m_F2R1 < min(
        min_group_amp, min_group_dmg
    ):
        return (
            np.zeros([64, 4]),
            np.zeros([10, 12]),
            np.zeros([5, 11]),
            np.zeros([64, 4]),
            np.zeros([10, 12]),
            np.zeros([5, 11]),
        )

    ### Prepare sequence matrix and quality matrix for each strand
    n = len(reference_int)
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}
    base2num_npfunc = np.vectorize(lambda b: base2num[b])
    F1R2_seq_mat = np.zeros([m_F1R2, n], dtype=int)  # Base(ATCG) x reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2, n])
    F2R1_seq_mat = np.zeros([m_F2R1, n], dtype=int)  # Base(ATCG) x reads x pos
    F2R1_qual_mat = np.zeros([m_F2R1, n])
    for mm, seq in enumerate(F1R2):
        qualities = seq.query_alignment_qualities
        sequence = np.array(list(seq.query_alignment_sequence))
        cigartuples = seq.cigartuples
        current_seq_ind = 0
        current_mat_ind = seq.reference_start - reference_start
        reference_ind = seq.reference_start - reference_start
        ref_length_plus_del = seq.reference_length
        for ct in cigartuples:
            if ct[0] == 0:
                F1R2_seq_mat[
                    mm, current_mat_ind : current_mat_ind + ct[1]
                ] = base2num_npfunc(sequence[current_seq_ind : current_seq_ind + ct[1]])
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
                # antimask[reference_ind : reference_ind + ct[1]] = False
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F1R2_seq_mat[mm, current_mat_ind:n] = 4
        F1R2_qual_mat[mm, current_mat_ind:n] = 0
    for mm, seq in enumerate(F2R1):
        qualities = seq.query_alignment_qualities
        sequence = np.array(list(seq.query_alignment_sequence))
        cigartuples = seq.cigartuples
        current_seq_ind = 0
        current_mat_ind = seq.reference_start - reference_start
        reference_ind = seq.reference_start - reference_start
        ref_length_plus_del = seq.reference_length
        for ct in cigartuples:
            if ct[0] == 0:
                F2R1_seq_mat[
                    mm, current_mat_ind : current_mat_ind + ct[1]
                ] = base2num_npfunc(sequence[current_seq_ind : current_seq_ind + ct[1]])
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
                # antimask[reference_ind : reference_ind + ct[1]] = False
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F2R1_seq_mat[mm, current_mat_ind:n] = 4
        F2R1_qual_mat[mm, current_mat_ind:n] = 0

    F1R2_qual_mat[F1R2_qual_mat <= params["minBq"]] = 0
    F2R1_qual_mat[F2R1_qual_mat <= params["minBq"]] = 0

    F1R2_antimask = antimask.copy()
    F2R1_antimask = antimask.copy()
    if m_F1R2 < min_group_amp or m_F2R1 < min_group_amp:
        F1R2_antimask[:] = False
        F2R1_antimask[:] = False

    dmg_antimask = antimask.copy()
    if m_F1R2 < min_group_dmg or m_F2R1 < min_group_dmg:
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
            (F1R2_count_mat).sum(axis=0) < 3, (F2R1_count_mat).sum(axis=0) < 3
        )
    ] = False  # mask locations with depth < 3
    dmg_antimask[
        np.logical_or(
            (F1R2_qual_mat_merged).sum(axis=0) < 90,
            (F2R1_qual_mat_merged).sum(axis=0) < 90,
        )
    ] = False  # mask locations where consensus bases is less than 90

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

    F1R2_antimask[ds_alt] = False
    F1R2_count_sum = F1R2_count_mat.sum(axis=0)
    F1R2_antimask[F1R2_count_sum < 3] = False
    F1R2_ref_count = F1R2_count_mat[reference_int, np.ogrid[: reference_int.size]]
    F1R2_antimask[F1R2_count_sum - F1R2_ref_count >= 2] = False
    F1R2_antimask[
        np.logical_and((F1R2_count_mat >= 1).sum(axis=0) < 2, F1R2_ref_count == 0)
    ] = False
    # F1R2_antimask[(F1R2_seq_mat == 4).any(axis=0)] = False

    F2R1_antimask[ds_alt] = False
    F2R1_count_sum = F2R1_count_mat.sum(axis=0)
    F2R1_antimask[F2R1_count_sum < 3] = False
    F2R1_ref_count = F2R1_count_mat[reference_int, np.ogrid[: reference_int.size]]
    F2R1_antimask[F2R1_count_sum - F2R1_ref_count >= 2] = False
    F2R1_antimask[
        np.logical_and((F2R1_count_mat >= 1).sum(axis=0) < 2, F2R1_ref_count == 0)
    ] = False
    # F2R1_antimask[(F2R1_seq_mat == 4).any(axis=0)] = False

    F1R2_trinuc_masked = trinuc_int[F1R2_antimask]
    F1R2_antimask_positions = np.nonzero(F1R2_antimask)[0]
    F1R2_trinuc_alt_count_mat = np.zeros([96, 4])
    F1R2_trinuc_seq_err_count_mat = np.zeros([96, 4])
    for mm in range(F1R2_seq_mat.shape[0]):
        seq_masked = F1R2_seq_mat[mm, F1R2_antimask]
        ref_masked = reference_int[F1R2_antimask]
        mismatch_bool = seq_masked != ref_masked
        n_mismatch = int(mismatch_bool.sum())
        if n_mismatch == 0:
            continue
        if n_mismatch == 1:
            F1R2_trinuc_alt_1Dmap = F1R2_trinuc_masked + seq_masked * 96
            # F1R2_trinuc_alt_1Dmap = F1R2_trinuc_alt_1Dmap[F1R2_trinuc_alt_1Dmap < 4*96]
            F1R2_trinuc_alt_count_mat += (
                np.bincount(
                    F1R2_trinuc_alt_1Dmap,
                    weights=1 - 10 ** (-F1R2_qual_mat[mm, F1R2_antimask] / 10),
                    minlength=96 * 4,
                )[0 : 4 * 96]
                .reshape([4, 96])
                .T
            )
        # n_mismatch > 1: too unreliable, excluded from SBS amp
    # F1R2_trinuc_alt_count_mat_norm = F1R2_trinuc_alt_count_mat[:32,:] + F1R2_trinuc_alt_count_mat[32:64,np.array([1,0,3,2])]
    F1R2_trinuc_alt_count_mat_norm = F1R2_trinuc_alt_count_mat[0:64, :] + np.vstack(
        [
            F1R2_trinuc_alt_count_mat[32:64, [1, 0, 3, 2]],
            F1R2_trinuc_alt_count_mat[:32, [1, 0, 3, 2]],
        ]
    )

    F2R1_trinuc_alt_count_mat = np.zeros([96, 4])
    F2R1_trinuc_masked = trinuc_int[F2R1_antimask]
    F2R1_antimask_positions = np.nonzero(F2R1_antimask)[0]
    # F1R2_alt_masked = F1R2_alt_int[F1R2_antimask]
    for mm in range(F2R1_seq_mat.shape[0]):
        seq_masked = F2R1_seq_mat[mm, F2R1_antimask]
        ref_masked = reference_int[F2R1_antimask]
        mismatch_bool = seq_masked != ref_masked
        n_mismatch = int(mismatch_bool.sum())
        if n_mismatch == 0:
            continue
        if n_mismatch == 1:
            F2R1_trinuc_alt_1Dmap = F2R1_trinuc_masked + seq_masked * 96
            # F2R1_trinuc_alt_1Dmap = F2R1_trinuc_alt_1Dmap[F2R1_trinuc_alt_1Dmap < 4*96]
            F2R1_trinuc_alt_count_mat += (
                np.bincount(
                    F2R1_trinuc_alt_1Dmap,
                    weights=1 - 10 ** (-F2R1_qual_mat[mm, F2R1_antimask] / 10),
                    minlength=96 * 4,
                )[0 : 4 * 96]
                .reshape([4, 96])
                .T
            )
        # n_mismatch > 1: too unreliable, excluded from SBS amp
    # F2R1_trinuc_alt_count_mat_norm = F2R1_trinuc_alt_count_mat[:32,:] + F2R1_trinuc_alt_count_mat[32:64,np.array([1,0,3,2])]
    F2R1_trinuc_alt_count_mat_norm = F2R1_trinuc_alt_count_mat[0:64, :] + np.vstack(
        [
            F2R1_trinuc_alt_count_mat[32:64, [1, 0, 3, 2]],
            F2R1_trinuc_alt_count_mat[:32, [1, 0, 3, 2]],
        ]
    )

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
    start = seqs[0].reference_start
    indels = list(indels)
    indels_masked = list()
    for indel in indels:
        refPos = int(indel.split(":")[0]) - 1
        indelLen = int(indel.split(":")[1])
        if antimask[refPos - start]:
            # if indelLen < 0:
            # if antimask[refPos - start : refPos - start- indelLen].all():
            indels_masked.append(indel)
        # else: indels_masked.append(indel)

    m = len(indels_masked)
    if m >= 2:
        return (
            F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,
            np.zeros([10, 12]),
            np.zeros([5, 11]),
            dmg_trinuc_alt_count_mat_norm,
            np.zeros([10, 12]),
            np.zeros([5, 11]),
        )
    F1R2_antimask = antimask.copy()
    F2R1_antimask = antimask.copy()
    if m_F1R2 < min_group_amp or m_F2R1 < min_group_amp:
        F1R2_antimask[:] = False
        F2R1_antimask[:] = False
    dmg_antimask = antimask.copy()
    if m_F1R2 < min_group_dmg or m_F2R1 < min_group_dmg:
        dmg_antimask[:] = False

    # hp_raw (hp.h5): row0 = self-derived homopolymer run length, row1 =
    # start-of-run bool. str_raw (str.h5): row0 = repeat unit length,
    # row1 = number of times the unit repeats, row2 = start-of-repeat
    # bool. These are two independent sources — a position inside an
    # annotated STR interval that also starts an embedded run of
    # identical bases (e.g. the "AA" in a (AAT)n repeat) is credited to
    # BOTH the homopolymer and the STR opportunity/error tables, never
    # mutually exclusive (unchanged from before this rewrite).
    hp_len_arr = hp_raw[0].astype(int)
    hp_mask = hp_raw[1] == 1

    unit_len = str_raw[0].astype(int)
    repeat_count = str_raw[1].astype(int)
    total_len = unit_len * repeat_count
    is_str = unit_len >= 2
    # 5 bins: 0 = not a real repeat (is_str False -- never populated by
    # the opportunity pass below, only ever reached by an actual
    # mismatched-insertion event; see the per-event loop), 1="2-9",
    # 2="10-24", 3="25-39", 4="40+".
    str_bin_arr = np.zeros_like(total_len)
    str_bin_arr[is_str] = 1
    str_bin_arr[is_str & (total_len >= 10)] = 2
    str_bin_arr[is_str & (total_len >= 25)] = 3
    str_bin_arr[is_str & (total_len >= 40)] = 4
    str_mask = (str_raw[2] == 1) & is_str

    hp_rc4 = [1, 0, 3, 2]  # base-complement permutation, 4-wide axis

    # HP amp opportunity: every start-of-run position (regardless of
    # is_str) contributes its own (hp length, own reference base) as an
    # idLen=0 (column base*3+1) observation, weighted by family size and
    # self-RC-folded per strand -- same convention as the SBS/DBS amp
    # matrices above (each strand's own tally symmetrized under base-
    # complement, then the two strands summed), since amp errors have no
    # strand-of-origin directionality. reference_int can be 4 (N/
    # ambiguous) at a run-start position; those are excluded rather than
    # risking an out-of-range column.
    def _hp_amp_opportunity(strand_antimask, weight):
        hp_here = hp_len_arr[strand_antimask][hp_mask[strand_antimask]]
        ref_here = reference_int[strand_antimask][hp_mask[strand_antimask]]
        valid = ref_here <= 3
        hp_here = np.minimum(hp_here[valid], 10)
        ref_here = ref_here[valid]
        local = np.zeros([10, 4])
        np.add.at(local, (hp_here - 1, ref_here), 1)
        local *= weight
        return local + local[:, hp_rc4]

    hp_alt_count[:, [1, 4, 7, 10]] += _hp_amp_opportunity(F1R2_antimask, m_F1R2)
    hp_alt_count[:, [1, 4, 7, 10]] += _hp_amp_opportunity(F2R1_antimask, m_F2R1)

    # STR amp opportunity: only at annotated-repeat start positions,
    # column 5 (idLen=0). No base identity/RC-fold involved (a repeat's
    # unit length has no complement).
    def _str_amp_opportunity(strand_antimask, weight):
        bin_here = str_bin_arr[strand_antimask][str_mask[strand_antimask]]
        return np.bincount(bin_here, minlength=5) * weight

    str_alt_count[:, 5] += _str_amp_opportunity(F1R2_antimask, m_F1R2)
    str_alt_count[:, 5] += _str_amp_opportunity(F2R1_antimask, m_F2R1)

    # STR amp opportunity for row 0 ("not a real repeat"): unlike rows
    # 1-4 (one credit per real, cut-gated repeat tract via str_mask),
    # row 0 has no tract to dedupe by -- str_mask structurally excludes
    # every row-0 position (str_mask requires is_str), so without this,
    # row 0's own column 5 stays permanently 0 and
    # _normalize_indel_str_mat's row-sum normalization has no real
    # opportunity denominator for row 0, silently normalizing its rare
    # event counts against each other instead and producing wildly
    # inflated "probabilities" (observed: up to 0.58, when a real
    # damage/amp rate should be ~1e-5 or smaller). Every antimask-passing
    # position independently counts here, INCLUDING real-STR-annotated
    # ones -- not gated by ~is_str. A row-0-type event (a mismatched
    # insertion/arbitrary indel unrelated to any specific repeat) is a
    # real, independent opportunity even at a position that also has a
    # real, different-unit STR annotation (e.g. an arbitrary "TAA"
    # insertion inside an ATGATGATG tract isn't a repeat of the ATG unit,
    # so it's still a genuine row-0 event there). Matches
    # misc.py's indel100_reference_bucket_indices flat rep1/rep0 credit
    # and call.py's L_indel_len[...,0,X] lookup, both of which likewise
    # credit row 0 unconditionally now, not real_str-excluded.
    def _str_amp_opportunity_row0(strand_antimask, weight):
        return np.count_nonzero(strand_antimask) * weight

    str_alt_count[0, 5] += _str_amp_opportunity_row0(F1R2_antimask, m_F1R2)
    str_alt_count[0, 5] += _str_amp_opportunity_row0(F2R1_antimask, m_F2R1)

    # HP dmg opportunity: a single combined pass (dmg_antimask doesn't
    # distinguish which strand might show damage), credited to BOTH
    # orientations from the same position data -- direct (F1R2 frame)
    # and base-complemented (F2R1/bottom-strand frame) via the same
    # self-fold operation as the amp case, mirroring how the per-event
    # dmg branch below folds F2R1's contribution onto F1R2's frame
    # instead of self-folding each independently.
    hp_dmg_here = hp_len_arr[dmg_antimask][hp_mask[dmg_antimask]]
    ref_dmg_here = reference_int[dmg_antimask][hp_mask[dmg_antimask]]
    valid_dmg = ref_dmg_here <= 3
    hp_dmg_here = np.minimum(hp_dmg_here[valid_dmg], 10)
    ref_dmg_here = ref_dmg_here[valid_dmg]
    hp_dmg_local = np.zeros([10, 4])
    np.add.at(hp_dmg_local, (hp_dmg_here - 1, ref_dmg_here), 1)
    hp_dmg_count[:, [1, 4, 7, 10]] += hp_dmg_local + hp_dmg_local[:, hp_rc4]

    # STR dmg opportunity: same *2 (both orientations, no base identity
    # to complement) as the two separate-weighted amp strand passes.
    str_dmg_here = str_bin_arr[dmg_antimask][str_mask[dmg_antimask]]
    str_dmg_count[:, 5] += np.bincount(str_dmg_here, minlength=5) * 2

    # STR dmg opportunity for row 0 -- same reasoning as the amp case
    # above (unconditional, not ~is_str-excluded), same *2 (both
    # orientations, no base identity to complement) as the two
    # separate-weighted amp strand passes.
    str_dmg_count[0, 5] += np.count_nonzero(dmg_antimask) * 2

    if m == 0:
        return (
            F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,
            hp_alt_count,
            str_alt_count,
            dmg_trinuc_alt_count_mat_norm,
            hp_dmg_count,
            str_dmg_count,
        )

    F1R2_alt_qual = np.zeros(m)
    F1R2_ref_qual = np.zeros(m)
    F2R1_alt_qual = np.zeros(m)
    F2R1_ref_qual = np.zeros(m)
    F1R2_alt_count = np.zeros(m)
    F1R2_ref_count = np.zeros(m)
    F2R1_alt_count = np.zeros(m)
    F2R1_ref_count = np.zeros(m)
    for seq in F1R2:
        seqArr, qualArr = getIndelArr(seq, indels_masked)
        ac = np.count_nonzero(seqArr == 1)
        rc = np.count_nonzero(seqArr == 0)
        aq = np.sum(qualArr[seqArr == 1])
        rq = np.sum(qualArr[seqArr == 0])
        F1R2_alt_qual += aq
        F1R2_ref_qual += rq
        F1R2_alt_count += ac
        F1R2_ref_count += rc
    for seq in F2R1:
        seqArr, qualArr = getIndelArr(seq, indels_masked)
        ac = np.count_nonzero(seqArr == 1)
        rc = np.count_nonzero(seqArr == 0)
        aq = np.sum(qualArr[seqArr == 1])
        rq = np.sum(qualArr[seqArr == 0])
        F2R1_alt_qual += aq
        F2R1_ref_qual += rq
        F2R1_alt_count += ac
        F2R1_ref_count += rc

    dmg_antimask = np.ones(m, dtype=bool)
    dmg_antimask[
        np.logical_or(
            (F1R2_ref_count + F1R2_alt_count) < 3, (F2R1_ref_count + F2R1_alt_count) < 3
        )
    ] = False
    dmg_antimask[np.logical_and(F1R2_ref_count != 0, F1R2_alt_count != 0)] = False
    dmg_antimask[np.logical_and(F2R1_ref_count != 0, F2R1_alt_count != 0)] = False
    dmg_antimask[F1R2_ref_qual + F1R2_alt_qual < 90] = False
    dmg_antimask[F2R1_ref_qual + F2R1_alt_qual < 90] = False
    dmg_antimask[np.logical_and(F1R2_alt_count > 0, F2R1_alt_count > 0)] = False
    F1R2_dmg_antimask = dmg_antimask.copy()
    F1R2_dmg_antimask[F1R2_alt_count == 0] = False

    F2R1_dmg_antimask = dmg_antimask.copy()
    F2R1_dmg_antimask[F2R1_alt_count == 0] = False

    F1R2_antimask = np.ones(m, dtype=bool)
    F2R1_antimask = np.ones(m, dtype=bool)

    F1R2_antimask[F1R2_ref_count == 0] = False
    F1R2_antimask[F1R2_ref_count + F1R2_alt_count < 3] = False
    F1R2_antimask[F1R2_alt_count > 1] = False

    F2R1_antimask[F2R1_ref_count == 0] = False
    F2R1_antimask[F2R1_ref_count + F2R1_alt_count < 3] = False
    F2R1_antimask[F2R1_alt_count > 1] = False

    # m == 1 here (m>=2 and m==0 both returned above), so this loop runs
    # exactly once -- kept as a loop (rather than indexing indels_masked[0]
    # directly) only to mirror the surrounding code's style.
    for mm, indel in enumerate(indels_masked):
        parts = indel.split(":")
        pos = int(parts[0]) - start
        indelLen = int(parts[1])
        # anchor is the reference base immediately after the insertion
        # point / the deleted base itself -- same position for both
        # directions, matching funcs/prob.py's genotypeDSIndel exactly, so
        # learn-time and call-time always agree.
        anchor = pos + 1
        hp = int(hp_raw[0, anchor])
        hp_capped = min(hp, 10)
        unit_len_here = int(str_raw[0, anchor])
        repeat_count_here = int(str_raw[1, anchor])
        if unit_len_here >= 2:
            total_len_here = unit_len_here * repeat_count_here
            if total_len_here >= 40:
                str_bin_here = 4
            elif total_len_here >= 25:
                str_bin_here = 3
            elif total_len_here >= 10:
                str_bin_here = 2
            else:
                str_bin_here = 1
        else:
            str_bin_here = 0

        if indelLen > 5:
            indelLen = 5
        if indelLen < -5:
            indelLen = -5

        if indelLen == 1 or indelLen == -1:
            ref_allele = int(reference_int[anchor])
            ref_allele_rc = int(reverse_comp[ref_allele])
            if indelLen == -1:
                # Deletion: the deleted base is trivially the
                # homopolymer's own base -- always a match.
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

                # Amp: reconcile against the opportunity credit already
                # given to this exact (hp_len, base) cell above (this
                # position was covered and its hp/base context counted
                # there without knowing an indel would land here), then
                # self-RC-fold each strand's own delta before summing --
                # same convention as the opportunity pass.
                F1R2_local = np.zeros(12)
                F2R1_local = np.zeros(12)
                if F1R2_antimask[mm]:
                    F1R2_local[col] += F1R2_alt_count[mm]
                    F1R2_local[opp_col] -= F1R2_alt_count[mm]
                if F2R1_antimask[mm]:
                    F2R1_local[col_rc] += F2R1_alt_count[mm]
                    F2R1_local[opp_col_rc] -= F2R1_alt_count[mm]
                F1R2_local = F1R2_local.reshape(4, 3)
                F1R2_local = F1R2_local + F1R2_local[hp_rc4, :]
                F2R1_local = F2R1_local.reshape(4, 3)
                F2R1_local = F2R1_local + F2R1_local[hp_rc4, :]
                hp_alt_count[row, :] += (F1R2_local + F2R1_local).reshape(12)

                # Dmg: both strands fold directly into the same cell
                # (F2R1 via the complemented base), same as the dmg
                # opportunity pass and the SBS/DBS dmg convention.
                if F1R2_dmg_antimask[mm]:
                    hp_dmg_count[row, col] += 1
                    hp_dmg_count[row, opp_col] -= 1
                if F2R1_dmg_antimask[mm]:
                    hp_dmg_count[row, col_rc] += 1
                    hp_dmg_count[row, opp_col_rc] -= 1
            else:
                # Inserted base doesn't match the flanking homopolymer:
                # not a real slippage event -- str.txt row 0, +-1
                # column. Row 0 now does receive opportunity credit (see
                # the row-0 opportunity pass above), so this position's
                # own event needs the same reconciliation subtraction
                # rows 1-4 get below -- otherwise it would be double-
                # counted as both "no event" (col 5) and "this event"
                # (col). No base identity to RC-fold either way.
                col = indelLen + 5
                if F1R2_antimask[mm]:
                    str_alt_count[0, col] += F1R2_alt_count[mm]
                    str_alt_count[0, 5] -= F1R2_alt_count[mm]
                if F2R1_antimask[mm]:
                    str_alt_count[0, col] += F2R1_alt_count[mm]
                    str_alt_count[0, 5] -= F2R1_alt_count[mm]
                if dmg_antimask[mm]:
                    str_dmg_count[0, col] += 1
                    str_dmg_count[0, 5] -= 1
        else:
            # Length >=2: always STR-context, keyed by this position's
            # real STR-length bin (0 if not actually annotated -- no
            # hp-length fallback). Row 0 receives opportunity credit same
            # as rows 1-4 (see the row-0 opportunity pass above), so the
            # reconciliation subtraction applies uniformly regardless of
            # row.
            row = str_bin_here
            col = indelLen + 5
            if F1R2_antimask[mm]:
                str_alt_count[row, col] += F1R2_alt_count[mm]
                str_alt_count[row, 5] -= F1R2_alt_count[mm]
            if F2R1_antimask[mm]:
                str_alt_count[row, col] += F2R1_alt_count[mm]
                str_alt_count[row, 5] -= F2R1_alt_count[mm]
            if dmg_antimask[mm]:
                str_dmg_count[row, col] += 1
                str_dmg_count[row, 5] -= 1

    return (
        F1R2_trinuc_alt_count_mat_norm + F2R1_trinuc_alt_count_mat_norm,
        hp_alt_count,
        str_alt_count,
        dmg_trinuc_alt_count_mat_norm,
        hp_dmg_count,
        str_dmg_count,
    )
