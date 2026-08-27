import numpy as np
from .indels import findIndels, getIndelArr, left_align_indel


def log10(mat):
    return np.log10(np.where(mat > 0, mat, np.finfo(float).eps))


def power10(mat):
    return 10 ** np.where(mat >= np.log10(np.finfo(float).eps), mat, -np.inf)


def log2(mat):
    return np.log2(np.where(mat > 0, mat, np.finfo(float).eps))


def power2(mat):
    return np.where(mat >= -100, 2**mat, 0)


def log(mat):
    return np.log(np.where(mat > 0, mat, np.finfo(float).eps))


def exp(mat):
    return np.exp(np.where(mat >= np.log10(np.finfo(float).eps), mat, -np.inf))


def _logsumexp4(terms):
    """Numerically stable log(sum(exp(terms), axis=0)) for a fixed 4-row
    stack — equivalent to scipy.special.logsumexp(terms, axis=0), but
    calculateDSPosterior calls this on tiny arrays hundreds of thousands
    of times (both per-family and while precomputing the indel power grid
    in call.py), and scipy's generic array-API dispatch/dtype-promotion
    layer (xp_promote/isdtype/_preprocess_dtype) dominates cost at that
    call volume. Since the input is always exactly 4 rows, a direct
    numpy implementation of the same max-subtraction trick skips that
    overhead entirely while returning bit-for-bit the same values,
    including scipy's -inf (not nan) result when every row is -inf.
    """
    m = np.max(terms, axis=0)
    finite_m = np.where(np.isfinite(m), m, 0)
    with np.errstate(divide="ignore"):
        out = np.log(np.sum(np.exp(terms - finite_m), axis=0)) + finite_m
    return np.where(np.isfinite(m), out, m)


def calculateDSPosterior(Pt, P_rev_t, Pb, P_rev_b, PAt, PAb, PBt, PBb):
    PA_At = PAt + log(1 - Pt)
    PA_Ab = PAb + log(1 - Pb)
    PA_Bt = PBt + log(Pt)
    PA_Bb = PBb + log(Pb)
    PB_At = PAt + log(P_rev_t)
    PB_Ab = PAb + log(P_rev_b)
    PB_Bt = PBt + log(1 - P_rev_t)
    PB_Bb = PBb + log(1 - P_rev_b)

    # ll1 = log10(power10(PA_At+PA_Ab) + power10(PA_Bt+PA_Ab) + power10(PA_At+PA_Bb) + power10(PA_Bt+PA_Bb))
    # ll2 = log10(power10(PB_At+PB_Ab) + power10(PB_Bt+PB_Ab) + power10(PB_At+PB_Bb) + power10(PB_Bt+PB_Bb))
    ll1 = _logsumexp4(
        np.vstack((PA_At + PA_Ab, PA_Bt + PA_Ab, PA_At + PA_Bb, PA_Bt + PA_Bb))
    ) / np.log(10)
    ll2 = _logsumexp4(
        np.vstack((PB_At + PB_Ab, PB_Bt + PB_Ab, PB_At + PB_Bb, PB_Bt + PB_Bb))
    ) / np.log(10)
    return ll1, ll2


def calculateSSPosterior(P, P_rev, bin_seq, Pseq):  # countb1, countb2, Pb1, Pb2):
    # Pseq is ln(base-call error rate implied by BQ): exp(Pseq) is the
    # probability of ANY miscall, but the mixture below treats it as the
    # probability of specifically miscalling to the ONE alternate allele
    # under test here, so it's divided by 3 (the number of possible wrong
    # bases) -- skipped for the "no real quality info" sentinel (Pseq
    # coming in as exactly 0, e.g. a >2bp indel insertion, see
    # funcs/indels.py's getIndelArr) so that case stays exactly the
    # uninformative 50/50 log(0.5) below, not further scaled by it.
    zero_mask = Pseq == 0
    Pseq[zero_mask] = log(0.5)
    Pseq[~zero_mask] -= log(3)
    bin_seq = bin_seq.astype(bool, copy=False)

    expP = np.exp(Pseq)

    # precompute the two mixture terms
    A = (1 - P) * (1 - expP) + P * expP
    B = (1 - P) * expP + P * (1 - expP)
    A_rev = (1 - P_rev) * (1 - expP) + P_rev * expP
    B_rev = (1 - P_rev) * expP + P_rev * (1 - expP)
    # precompute logs
    logA = np.log(A)
    logB = np.log(B)
    logA_rev = np.log(A_rev)
    logB_rev = np.log(B_rev)

    # precompute exp
    delta_fwd = logA - logB
    delta_rev = logB_rev - logA_rev

    mask0 = ~bin_seq  # sparse

    rows, cols = np.nonzero(mask0)

    prob1 = logA.sum(axis=0)
    np.add.at(prob1, cols, -delta_fwd[rows, cols])

    prob2 = logB_rev.sum(axis=0)
    np.add.at(prob2, cols, -delta_rev[rows, cols])
    return prob1, prob2


def genotypeDSSnv(
    seqs,
    reference_start,
    reference_int,
    trinuc_int,
    prior_mat,
    antimask,
    mut_antimask_scope,
    params,
    L=None,
):
    """antimask and mut_antimask_scope both start as boolean masks over
    the same window, refined identically by the structural/quality checks
    below (trinuc validity, zero-qual fraction, 3+-allele ambiguity) --
    but antimask is the narrow one (still respects n_cov_mask/nm_mask/
    trim, whatever the caller passed in) used only to gate cov_mat, while
    mut_antimask_scope is the wide one (only include_mask ever blocks it)
    used to gate mut_antimask/candidate detection. They're deliberately
    decoupled so a --rescue-eligible candidate blocked only by a
    rescuable mask (n_cov/nm/trim/indel_mask) can still get a real LR
    score without also making low-reliability positions contribute to
    cov_mat's coverage/opportunity accounting -- see call.py's rescue
    branch for how the resulting LR gets used.
    """
    prob_amp_mat = params["ampmat"]
    prob_amp_mat_rev = params["ampmat_rev"]
    prob_dmg_mat_top = params["dmgmat_top"]
    prob_dmg_mat_rev_top = params["dmgmat_rev_top"]
    prob_dmg_mat_bot = params["dmgmat_bot"]
    prob_dmg_mat_rev_bot = params["dmgmat_rev_bot"]
    trinuc_convert_np = params["trinuc_convert"]
    antimask[trinuc_int >= 64] = False
    mut_antimask_scope[trinuc_int >= 64] = False
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

    F1R2_qual_mat_0_count = np.count_nonzero(F1R2_qual_mat == 0, axis=0)
    F2R1_qual_mat_0_count = np.count_nonzero(F2R1_qual_mat == 0, axis=0)

    F1R2_qual_mat_merged = np.zeros([4, n])
    F1R2_count_mat = np.zeros([4, n], dtype=int)

    for nn in range(0, 4):
        F1R2_qual_mat_merged[nn, :] = F1R2_qual_mat.sum(
            axis=0, where=(F1R2_seq_mat == nn)
        )
        F1R2_count_mat[nn, :] = (
            np.logical_and(F1R2_seq_mat == nn, F1R2_qual_mat != 0)
        ).sum(axis=0)
    F2R1_qual_mat_merged = np.zeros([4, n])
    F2R1_count_mat = np.zeros([4, n], dtype=int)
    for nn in range(0, 4):
        F2R1_qual_mat_merged[nn, :] = F2R1_qual_mat.sum(
            axis=0, where=(F2R1_seq_mat == nn)
        )
        F2R1_count_mat[nn, :] = (F2R1_seq_mat == nn).sum(axis=0)
        F2R1_count_mat[nn, :] = (
            np.logical_and(F2R1_seq_mat == nn, F2R1_qual_mat != 0)
        ).sum(axis=0)
    total_count_mat = F1R2_count_mat + F2R1_count_mat
    zero_qual_frac_fail = (
        (F1R2_qual_mat_0_count + F2R1_qual_mat_0_count) / (m_F1R2 + m_F2R1)
    ) >= params["maxZeroQualFrac"]
    ambiguous_allele_fail = (total_count_mat >= 1).sum(axis=0) > 2
    antimask[zero_qual_frac_fail] = False
    antimask[ambiguous_allele_fail] = False
    mut_antimask_scope[zero_qual_frac_fail] = False
    mut_antimask_scope[ambiguous_allele_fail] = False
    base1_int = np.argmax(total_count_mat, axis=0)
    total_count_without_base1 = total_count_mat.copy()
    total_count_without_base1[base1_int, np.ogrid[:n]] = -1
    base2_int = np.argmax(total_count_without_base1, axis=0)
    base2_int[base1_int != reference_int] = reference_int[base1_int != reference_int]
    # Empirical duplex coverage matrix [n, 4]: L lookup per unmasked position per alt base
    cov_mat = np.zeros([n, 4])
    if L is not None:
        n_top_L = np.minimum(F1R2_count_mat.sum(axis=0), 9).astype(int)
        n_bot_L = np.minimum(F2R1_count_mat.sum(axis=0), 9).astype(int)
        valid_idx_L = np.nonzero(np.logical_and(antimask, trinuc_int < 64))[0]
        if valid_idx_L.size > 0:
            cov_mat[valid_idx_L] = L[
                n_top_L[valid_idx_L], n_bot_L[valid_idx_L], trinuc_int[valid_idx_L]
            ]

    mut_antimask = np.logical_and(mut_antimask_scope, base1_int != reference_int)
    if not mut_antimask.any():
        return (
            cov_mat,
            np.zeros(0),
            np.zeros(0),
            mut_antimask,
            base1_int,
            antimask,
            F1R2_count_mat,
            F2R1_count_mat,
        )

    F1R2_masked_qual_mat = F1R2_qual_mat[:, mut_antimask]
    F2R1_masked_qual_mat = F2R1_qual_mat[:, mut_antimask]
    F1R2_masked_seq_mat = F1R2_seq_mat[:, mut_antimask]
    F2R1_masked_seq_mat = F2R1_seq_mat[:, mut_antimask]
    F1R2_prob = -F1R2_masked_qual_mat / 10
    F2R1_prob = -F2R1_masked_qual_mat / 10
    base1_int_masked = base1_int[mut_antimask]
    base2_int_masked = base2_int[mut_antimask]
    F1R2_bin_seq_mat = F1R2_masked_seq_mat == base1_int_masked
    F2R1_bin_seq_mat = F2R1_masked_seq_mat == base1_int_masked
    trinuc_converted_masked = trinuc_convert_np[
        trinuc_int[mut_antimask], base1_int_masked
    ]
    ref_int_masked = reference_int[mut_antimask]
    base2_int_masked[
        np.logical_and(
            base1_int_masked == ref_int_masked,
            total_count_mat[:, mut_antimask][
                base2_int_masked, np.ogrid[: base2_int_masked.size]
            ]
            == 0,
        )
    ] = 4
    Pamp = prob_amp_mat[trinuc_converted_masked, base2_int_masked]
    Pamp_rev = prob_amp_mat_rev[trinuc_converted_masked, base2_int_masked]
    Pdmg_t = prob_dmg_mat_top[trinuc_converted_masked, base2_int_masked]
    Pdmg_rev_t = prob_dmg_mat_rev_top[trinuc_converted_masked, base2_int_masked]
    Pdmg_b = prob_dmg_mat_bot[trinuc_converted_masked, base2_int_masked]
    Pdmg_rev_b = prob_dmg_mat_rev_bot[trinuc_converted_masked, base2_int_masked]

    Pdmg_t[Pdmg_t == 0] = 1e-9
    Pdmg_rev_t[Pdmg_rev_t == 0] = 1e-9
    Pdmg_b[Pdmg_b == 0] = 1e-9
    Pdmg_rev_b[Pdmg_rev_b == 0] = 1e-9

    F1R2_count_b1 = F1R2_count_mat[:, mut_antimask][
        base1_int_masked, np.ogrid[: base1_int_masked.size]
    ]
    F1R2_count_b2 = np.zeros(F1R2_count_b1.size)
    alt_pos = base2_int_masked != 4
    F1R2_count_b2[alt_pos] = F1R2_count_mat[:, mut_antimask][:, alt_pos][
        base2_int_masked[alt_pos], np.ogrid[: np.count_nonzero(alt_pos)]
    ]
    F2R1_count_b1 = F2R1_count_mat[:, mut_antimask][
        base1_int_masked, np.ogrid[: base1_int_masked.size]
    ]
    F2R1_count_b2 = np.zeros(F2R1_count_b1.size)
    alt_pos = base2_int_masked != 4
    F2R1_count_b2[alt_pos] = F2R1_count_mat[:, mut_antimask][:, alt_pos][
        base2_int_masked[alt_pos], np.ogrid[: np.count_nonzero(alt_pos)]
    ]
    ln10 = np.log(10)
    F1R2_b1_prob, F1R2_b2_prob = calculateSSPosterior(
        Pamp,
        Pamp_rev,
        F1R2_bin_seq_mat,
        F1R2_prob * ln10,
    )
    F2R1_b1_prob, F2R1_b2_prob = calculateSSPosterior(
        Pamp,
        Pamp_rev,
        F2R1_bin_seq_mat,
        F2R1_prob * ln10,
    )
    ln10 = np.log(10)
    LL_B1, LL_B2 = calculateDSPosterior(
        Pdmg_t,
        Pdmg_rev_t,
        Pdmg_b,
        Pdmg_rev_b,
        F1R2_b1_prob,
        F2R1_b1_prob,
        F1R2_b2_prob,
        F2R1_b2_prob,
    )
    LR_masked = LL_B1 - LL_B2
    LR_max = (
        log10(1 - Pdmg_t) + log10(1 - Pdmg_b) - log10(Pdmg_rev_t) - log10(Pdmg_rev_b)
    )
    # A negative LR_masked means the evidence actually favors base2 over
    # the majority-vote base1 -- not a real mutation candidate, just noise
    # in the vote. genotypeDSIndel has always excluded these via its own
    # `take = LR_masked >= 0` gate, independent of CS; mirror that here so
    # SNVs get the same non-CS-based sanity filter.
    keep = LR_masked >= 0
    mut_antimask[mut_antimask] = keep
    LR_masked = LR_masked[keep]
    LR_max = LR_max[keep]
    return (
        cov_mat,
        LR_masked,
        LR_max,
        mut_antimask,
        base1_int,
        antimask,
        F1R2_count_mat,
        F2R1_count_mat,
    )


def indelErrorProbs(
    hps,
    strs,
    idLen,
    ref_allele,
    inserted_base,
    prob_amp_hp,
    prob_dmg_hp,
    prob_amp_str,
    prob_dmg_str,
):
    """Select Pamp/Pdmg values for a candidate indel from the new split
    error matrices: prob_*_hp is (10, 12) -- rows hp run length 1-10+
    (capped), columns ref_allele*3+(idLen+1) for idLen in {-1,0,1} (the
    idLen=0 column is the reference/opportunity count); prob_*_str is
    (5, 11) -- rows STR-length bin 0="0-1" (i.e. not a real repeat) / 1=
    "2-9" / 2="10-24" / 3="25-39" / 4="40+", columns idLen+5 for idLen in
    -5..5 (idLen=0 again the opportunity column). See Index.py-style
    write-out in Caller.py/Learn.py and funcs/learn.py's accumulation for
    how these are built.

    Routing rule for +-1bp indels: a deletion always removes a base that
    was actually there, so it trivially belongs to that base's own
    homopolymer (hp matrix) regardless of hp length (even hp_len==1, an
    "isolated" base, counts). An insertion only belongs to the hp matrix
    if `inserted_base` (the base actually being inserted -- for real
    candidates, read off the event itself; for context/coverage
    enumeration where every base is by construction the position's own,
    pass inserted_base == ref_allele) matches `ref_allele` (the reference
    base immediately after the insertion point, i.e. the base whose run
    would be extended); otherwise it isn't a real homopolymer-slippage
    event and falls to str.txt's row 0 ("not a repeat") at the +-1
    column. Indels of length >=2 are always STR-context (no hp.txt
    column exists for them), keyed directly by `strs` (the STR-length
    bin already computed by the caller) regardless of any hp annotation
    -- a length>=2 indel not inside an annotated STR has no hp-length
    fallback here; it still routes through str.txt (row 0, "not a
    repeat").

    No separate "_rev" matrices, unlike the SBS side: Pamp_rev/Pdmg_rev
    are the same hp/str matrix as Pamp/Pdmg, just indexed at the column
    for -idLen instead of idLen (the opposite-direction column).
    """
    if idLen == 0:
        raise ValueError("idLen must be nonzero")
    if abs(idLen) >= 2:
        row = strs
        col = -idLen + 5
        col_rev = idLen + 5
        Pamp = prob_amp_str[row, col]
        Pamp_rev = prob_amp_str[row, col_rev]
        Pdmg = prob_dmg_str[row, col]
        Pdmg_rev = prob_dmg_str[row, col_rev]
        Pdmg_bot = Pdmg
        Pdmg_rev_bot = Pdmg_rev
    elif idLen == -1:
        # Deletion: always matches its own homopolymer.
        rc = [1, 0, 3, 2]
        ref_allele_rc = rc[ref_allele]
        row = min(hps, 10) - 1
        col = ref_allele * 3 + 0
        col_rev = ref_allele * 3 + 2
        col_bot = ref_allele_rc * 3 + 0
        col_rev_bot = ref_allele_rc * 3 + 2
        Pamp = prob_amp_hp[row, col]
        Pamp_rev = prob_amp_hp[row, col_rev]
        Pdmg = prob_dmg_hp[row, col]
        Pdmg_rev = prob_dmg_hp[row, col_rev]
        Pdmg_bot = prob_dmg_hp[row, col_bot]
        Pdmg_rev_bot = prob_dmg_hp[row, col_rev_bot]
    else:  # idLen == 1
        if inserted_base == ref_allele:
            rc = [1, 0, 3, 2]
            ref_allele_rc = rc[ref_allele]
            row = min(hps, 10) - 1
            col = ref_allele * 3 + 2
            col_rev = ref_allele * 3 + 0
            col_bot = ref_allele_rc * 3 + 2
            col_rev_bot = ref_allele_rc * 3 + 0
            Pamp = prob_amp_hp[row, col]
            Pamp_rev = prob_amp_hp[row, col_rev]
            Pdmg = prob_dmg_hp[row, col]
            Pdmg_rev = prob_dmg_hp[row, col_rev]
            Pdmg_bot = prob_dmg_hp[row, col_bot]
            Pdmg_rev_bot = prob_dmg_hp[row, col_rev_bot]
        else:
            row = 0
            col = 1 + 5
            col_rev = -1 + 5
            Pamp = prob_amp_str[row, col]
            Pamp_rev = prob_amp_str[row, col_rev]
            Pdmg = prob_dmg_str[row, col]
            Pdmg_rev = prob_dmg_str[row, col_rev]
            Pdmg_bot = Pdmg
            Pdmg_rev_bot = Pdmg_rev
    if Pamp == 0:
        Pamp = 1e-9
    return Pamp, Pamp_rev, Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot


def indelMaxLR(Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot):
    """Theoretical ceiling of masked LR for an indel context: log10(1-Pdmg)
    + log10(1-Pdmg_bot) - log10(Pdmg_rev) - log10(Pdmg_rev_bot). Same shape
    as genotypeDSSnv's LR_max/"LM" field; depends only on context (hps/
    strs/idLen/ref_allele via indelErrorProbs), never on read depth, so it
    doubles as the per-context calling-threshold ceiling in call.py/
    Caller.py.
    """
    return log10(1 - Pdmg) + log10(1 - Pdmg_bot) - log10(Pdmg_rev) - log10(Pdmg_rev_bot)


def genotypeDSIndel(
    seqs,
    reference_start,
    reference_end,
    reference_int,
    antimask,
    hp_raw,
    str_raw,
    params,
):
    prob_amp_hp = params["ampmat_hp"]
    prob_dmg_hp = params["dmgmat_hp"]
    prob_amp_str = params["ampmat_str"]
    prob_dmg_str = params["dmgmat_str"]
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)
    chrom = seqs[0].reference_name
    start = reference_start
    end = reference_end
    indels = set()
    ### Geonotype indel for all found indels
    for seq in F1R2:
        indels.update(
            left_align_indel(i, reference_int, reference_start) for i in findIndels(seq)
        )
    for seq in F2R1:
        indels.update(
            left_align_indel(i, reference_int, reference_start) for i in findIndels(seq)
        )
    # start = seqs[0].reference_start
    indels = list(indels)
    indels_masked = list()
    pos_masked = list()
    indelLen_masked = list()
    for indel in indels:
        refPos = int(indel.split(":")[0])
        indelLen = int(indel.split(":")[1])
        if indelLen < 0:
            indelLen_abs = -indelLen
        else:
            indelLen_abs = 0
        if antimask[refPos - start : refPos - start + indelLen_abs + 1].all():
            indels_masked.append(indel)
            pos_masked.append(refPos)
            indelLen_masked.append(indelLen)
    pos_masked = np.array(pos_masked)
    indelLen_masked = np.array(indelLen_masked, dtype=int)
    if len(indels_masked) != 0:
        pos_arg_sorted = np.argsort(pos_masked)
        pos_sorted = pos_masked[pos_arg_sorted]
        pos_take = np.ones(pos_masked.size, dtype=bool)
        pos_take[np.ediff1d(pos_sorted, to_begin=1) == 0] = False
        pos_take[np.ediff1d(pos_sorted, to_end=1) == 0] = False
        pos_arg_masked = pos_arg_sorted[pos_take]
        pos_masked = pos_masked[pos_arg_masked]
        indels_masked = [indels_masked[_] for _ in pos_arg_masked]
        indelLen_masked = indelLen_masked[pos_arg_masked]
        indelLen_masked[indelLen_masked > 5] = 5
        indelLen_masked[indelLen_masked < -5] = -5

    m = len(indels_masked)
    if m == 0:  # or m >= 2:
        return [np.zeros(0)] * 10
    n_f1r2 = len(F1R2)
    n_f2r1 = len(F2R1)
    mask_multiallele = np.ones(m, dtype=bool)
    f1r2_seq = np.zeros([n_f1r2, m])
    f2r1_seq = np.zeros([n_f2r1, m])
    f1r2_prob = np.zeros([n_f1r2, m])
    f2r1_prob = np.zeros([n_f2r1, m])

    f1r2_alt_count = np.zeros(m)
    f1r2_ref_count = np.zeros(m)
    f2r1_alt_count = np.zeros(m)
    f2r1_ref_count = np.zeros(m)

    for nn, seq in enumerate(F1R2):
        seqArr, qualArr = getIndelArr(seq, indels_masked)
        mask_multiallele[seqArr == -1] = 0
        f1r2_seq[nn, :] = seqArr > 0
        f1r2_prob[nn, :] = qualArr
        f1r2_alt_count += (seqArr == 1).astype(int)
        f1r2_ref_count += (seqArr == 0).astype(int)
    for nn, seq in enumerate(F2R1):
        seqArr, qualArr = getIndelArr(seq, indels_masked)
        mask_multiallele[seqArr == -1] = 0
        f2r1_seq[nn, :] = seqArr > 0
        f2r1_prob[nn, :] = qualArr
        f2r1_alt_count += (seqArr == 1).astype(int)
        f2r1_ref_count += (seqArr == 0).astype(int)
    f1r2_prob = -f1r2_prob / 10
    f2r1_prob = -f2r1_prob / 10

    offset = -indelLen_masked
    offset[offset < 0] = 0
    hps = np.zeros(pos_masked.size, dtype=int)
    strs = np.zeros(pos_masked.size, dtype=int)
    Pamp = np.zeros(pos_masked.size)
    Pamp_rev = np.zeros(pos_masked.size)
    Pdmg = np.zeros(pos_masked.size)
    Pdmg_rev = np.zeros(pos_masked.size)
    Pdmg_bot = np.zeros(pos_masked.size)
    Pdmg_rev_bot = np.zeros(pos_masked.size)
    # True for every deletion (always matches its own homopolymer) and
    # for a homopolymer-matching insertion; False only for an insertion
    # that indelErrorProbs routed to str.txt row 0. Saved per-candidate
    # (unlike Pamp/Pdmg, which only need the aggregate probability) so
    # Caller.py's channel/rate-table classification can route a real
    # PASS call the same way indelErrorProbs itself did, instead of
    # approximating every +-1bp call as HP-matched.
    hp_match_arr = np.ones(pos_masked.size, dtype=bool)
    for nn in range(pos_masked.size):
        # hps is always the self-derived homopolymer run length (hp_raw
        # row0), taken as the max over the indel's reference span so an
        # anchor landing just outside the repeat still picks up the run
        # it abuts — matches main's `hps[nn] = np.max(hp_int[0, ...])`.
        # strs is a separate, independent lookup against the BED-derived
        # STR annotation (str_raw), not a fallback gated by unit_len — an
        # embedded homopolymer inside an STR unit (e.g. the "AA" in a
        # (AAT)n repeat) still gets its true run length here, and is
        # never mutually exclusive with the STR classification.
        anchor = pos_masked[nn] - start + offset[nn] + 1
        hps[nn] = np.max(
            hp_raw[
                0,
                max(0, pos_masked[nn] - start) : pos_masked[nn]
                + offset[nn]
                - start
                + 2,
            ]
        )
        unit_len_here = int(str_raw[0, anchor])
        repeat_count_here = int(str_raw[1, anchor])
        # 5 STR-length bins now (0="0-1"/not a real repeat, 1="2-9",
        # 2="10-24", 3="25-39", 4="40+"), matching str.txt's 5 rows --
        # was 4 bins (0-3, with 0 meaning "<10" rather than "not a real
        # repeat") before this rewrite.
        if unit_len_here >= 2:
            total_len = unit_len_here * repeat_count_here
            if total_len >= 40:
                strs[nn] = 4
            elif total_len >= 25:
                strs[nn] = 3
            elif total_len >= 10:
                strs[nn] = 2
            else:
                strs[nn] = 1
        else:
            strs[nn] = 0
        idLen = indelLen_masked[nn]
        pos = pos_masked[nn]
        if hps[nn] > 10:
            hps[nn] = 10
        # Only evaluated lazily for +-1bp indels, matching the original
        # inline code, so this never indexes reference_int out of bounds for
        # longer indels near the edge of the window. ref_allele is the
        # reference base immediately after the insertion point / the
        # deleted base itself -- same anchor (pos-start+1) for both
        # directions, matching funcs/learn.py's accumulation exactly (the
        # old learn-time code used a different, inconsistent anchor for
        # deletions; this rewrite aligns the two).
        if idLen == 1 or idLen == -1:
            ref_allele = int(reference_int[pos - start + 1])
        else:
            ref_allele = 0
        # inserted_base: the actual base being inserted, read off the
        # event's own sequence (indels_masked entries are "pos:len:seq"
        # for insertions) -- only meaningful/used when idLen==1, where
        # indelErrorProbs compares it against ref_allele to decide
        # whether this is a real homopolymer-extending insertion (hp.txt)
        # or an unrelated single-base insertion (str.txt row 0).
        if idLen == 1:
            inserted_seq = indels_masked[nn].split(":")[2]
            inserted_base = base2num.get(inserted_seq[0], -1) if inserted_seq else -1
            hp_match_arr[nn] = inserted_base == ref_allele
        else:
            inserted_base = ref_allele
        (
            Pamp[nn],
            Pamp_rev[nn],
            Pdmg[nn],
            Pdmg_rev[nn],
            Pdmg_bot[nn],
            Pdmg_rev_bot[nn],
        ) = indelErrorProbs(
            hps[nn],
            strs[nn],
            idLen,
            ref_allele,
            inserted_base,
            prob_amp_hp,
            prob_dmg_hp,
            prob_amp_str,
            prob_dmg_str,
        )
    ln10 = np.log(10)
    F1R2_alt_prob, F1R2_ref_prob = calculateSSPosterior(
        Pamp,
        Pamp_rev,
        # f1r2_alt_count,
        # f1r2_ref_count,
        f1r2_seq,
        f1r2_prob * ln10,
    )
    F2R1_alt_prob, F2R1_ref_prob = calculateSSPosterior(
        Pamp,
        Pamp_rev,
        # f1r2_alt_count,
        # f1r2_ref_count,
        f2r1_seq,
        f2r1_prob * ln10,
    )
    ln10 = np.log(10)
    LL_B1, LL_B2 = calculateDSPosterior(
        Pdmg,
        Pdmg_rev,
        Pdmg_bot,
        Pdmg_rev_bot,
        F1R2_alt_prob,
        F2R1_alt_prob,
        F1R2_ref_prob,
        F2R1_ref_prob,
    )
    LR_masked = LL_B1 - LL_B2
    LR_max = indelMaxLR(Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot)
    take = LR_masked >= 0
    take[mask_multiallele == 0] = 0
    return (
        LR_masked[take],
        LR_max[take],
        [indels_masked[nn] for nn in range(len(take)) if take[nn]],
        hps[take],
        strs[take],
        f1r2_ref_count[take].astype("int"),
        f1r2_alt_count[take].astype("int"),
        f2r1_ref_count[take].astype("int"),
        f2r1_alt_count[take].astype("int"),
        hp_match_arr[take],
    )
