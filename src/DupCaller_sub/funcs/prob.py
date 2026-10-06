import numpy as np
from scipy.stats import poisson_binom
from .indels import (
    base_codes,
    str_length_bin,
    str_slip_tract,
    INDEL_ALT,
    INDEL_CONFLICT,
    INDEL_REF,
    findIndels,
    getIndelArr,
    indel_context_fits_window,
    indel_context_index,
    indel_has_context,
    indel_passes_mask,
    left_align_indel,
)


# |log10 LR| below this is rounding noise around a tie and is snapped to 0.
LR_ZERO_TOL = 1e-9


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
    # Pseq is ln(e), e = 10^(-BQ/10). The read emission is conditioned on
    # the read showing one of the two modelled alleles: P(the other
    # allele) = (e/3) / (1 - 2e/3), P(the true allele) = 1 minus that.
    # Pseq == 0 (no BQ for this read) stays uninformative at 0.5.
    bin_seq = bin_seq.astype(bool, copy=False)
    e = np.exp(Pseq)
    expP = (e / 3) / (1 - 2 * e / 3)
    expP = np.where(Pseq == 0, 0.5, expP)

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


def strand_evidence(qual_mat, alt_mat, p_extra):
    """Per-position strand-test inputs for one strand: (BQs of the counted
    reads, k = number of them that are non-alt, p_extra). A read counts
    where qual > 0. qual_mat and alt_mat are reads x positions; p_extra
    is per position (Pamp)."""
    counted = qual_mat > 0
    k = np.logical_and(counted, ~alt_mat).sum(axis=0)
    return [
        (qual_mat[counted[:, j], j], int(k[j]), float(p_extra[j]))
        for j in range(qual_mat.shape[1])
    ]


def _nonalt_upper_tail(probs, k):
    """P(X >= k), X Poisson-binomial over per-read non-alt probabilities
    `probs`. Computed as 1 - cdf(k - 1), so tails below ~1e-16 round to 0."""
    if k == 0:
        return 1.0
    return min(max(1.0 - poisson_binom(probs).cdf(k - 1), 0.0), 1.0)


def strand_pvalue(quals, k, p_extra):
    """SBS strand p-value: P(X >= k) under H0 "clean mutant strand", each
    counted read showing a non-alt base independently with probability
    10^(-BQ/10) + p_extra."""
    probs = np.minimum(10 ** (-np.asarray(quals, dtype=float) / 10) + p_extra, 1.0)
    return _nonalt_upper_tail(probs, k)


def indel_strand_evidence(n_alt, n_ref, rev_rate):
    """Indel strand-test inputs for one strand: (informative reads, REF
    reads, per-read reversion rate). Informative = ALT or REF in
    getIndelArr, which already applies the minBq gate."""
    return (int(n_alt) + int(n_ref), int(n_ref), float(rev_rate))


def indel_strand_pvalue(n, k, rev_rate):
    """Indel strand p-value: P(X >= k REF reads among n informative reads)
    under H0 "clean mutant strand", each read showing the reference
    independently with the learned reversion rate (indelReversionRate),
    which already includes sequencing error. NaN if the rate is
    unavailable."""
    if not np.isfinite(rev_rate):
        return float("nan")
    return _nonalt_upper_tail(np.full(n, min(rev_rate, 1.0)), k)


def strand_nonalt_pvalues(qual_mat, alt_mat, p_extra):
    """strand_pvalue for every column of a reads x positions matrix."""
    return np.array(
        [strand_pvalue(*ev) for ev in strand_evidence(qual_mat, alt_mat, p_extra)]
    )


def genotypeDSSnv(
    seqs,
    reference_start,
    reference_int,
    trinuc_int,
    antimask,
    mut_antimask_scope,
    params,
    L=None,
):
    """Genotype every position of one duplex family's window.

    antimask and mut_antimask_scope get the same checks below (trinuc
    validity, more than one non-reference base). antimask (all of the
    caller's masks) gates cov_mat; mut_antimask_scope (include_mask only)
    gates candidate detection, so a candidate blocked only by a rescuable
    mask still gets an LR for --rescue.
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
    F1R2_seq_mat = np.zeros([m_F1R2, n], dtype=int)  # Base(ATCG) x reads x pos
    F1R2_qual_mat = np.zeros([m_F1R2, n])
    F2R1_seq_mat = np.zeros([m_F2R1, n], dtype=int)  # Base(ATCG) x reads x pos
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
                # antimask[reference_ind : reference_ind + ct[1]] = False
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
                # antimask[reference_ind : reference_ind + ct[1]] = False
                reference_ind += ct[1]
                current_mat_ind += ct[1]
                ref_length_plus_del += ct[1]
        F2R1_seq_mat[mm, current_mat_ind:n] = 4
        F2R1_qual_mat[mm, current_mat_ind:n] = 0

    # BQs for the strand independence test: no minBq cut, N/uncovered = 0.
    F1R2_strand_qual_mat = np.where(F1R2_seq_mat == 4, 0, F1R2_qual_mat)
    F2R1_strand_qual_mat = np.where(F2R1_seq_mat == 4, 0, F2R1_qual_mat)
    F1R2_qual_mat[F1R2_qual_mat < params["minBq"]] = 0
    F2R1_qual_mat[F2R1_qual_mat < params["minBq"]] = 0

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
    # base1 is the position's non-reference base, base2 the reference. A
    # position with more than one distinct non-reference base (counted
    # reads only) is masked.
    ref_valid = reference_int < 4
    alt_count_mat = total_count_mat.copy()
    alt_count_mat[reference_int[ref_valid], np.nonzero(ref_valid)[0]] = 0
    n_alt_alleles = (alt_count_mat >= 1).sum(axis=0)
    ambiguous_allele_fail = n_alt_alleles > 1
    antimask[ambiguous_allele_fail] = False
    mut_antimask_scope[ambiguous_allele_fail] = False
    base1_int = np.where(
        n_alt_alleles >= 1, np.argmax(alt_count_mat, axis=0), reference_int
    )
    base2_int = reference_int.copy()
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
    # Non-candidates are reported as reference.
    base1_int[~mut_antimask] = reference_int[~mut_antimask]
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
            [],
            [],
            np.zeros(0, dtype=int),
            np.zeros(0, dtype=int),
            np.zeros(0),
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
    # LR < 0: the reads favor the reference over the mismatch base.
    LR_masked[np.abs(LR_masked) < LR_ZERO_TOL] = 0.0
    keep = LR_masked >= 0
    # Strand independence test inputs for the LR >= 0 candidates; Caller.py
    # computes the p-values (strand_pvalue) for PASS calls only.
    ev_F1R2 = strand_evidence(
        F1R2_strand_qual_mat[:, mut_antimask][:, keep],
        F1R2_bin_seq_mat[:, keep],
        Pamp[keep],
    )
    ev_F2R1 = strand_evidence(
        F2R1_strand_qual_mat[:, mut_antimask][:, keep],
        F2R1_bin_seq_mat[:, keep],
        Pamp[keep],
    )
    # LR < 0 positions count as reference coverage; their position,
    # mismatch base and LR are returned for the per-channel mu solve.
    rejected = np.nonzero(mut_antimask)[0][~keep]
    neg_alt = base1_int[rejected].copy()
    neg_LR = LR_masked[~keep]
    base1_int[rejected] = reference_int[rejected]
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
        ev_F1R2,
        ev_F2R1,
        rejected,
        neg_alt,
        neg_LR,
    )


_RC = [1, 0, 3, 2]


def _indel_error_cells(hps, strs, idLen, ref_allele, inserted_base, strs_mut=None):
    """Matrix cells of the two amplification/damage events that matter for
    a candidate indel, each as (matrix, row, col, col_bot) with matrix
    "hp" or "str" and col_bot the bottom-strand (base-complemented) column:

    reversion: alt->ref, the event that undoes the indel in the mutant
        molecule's own context.
    forward: ref->alt, the indel itself in the reference context.

    hp rows are run length 1-10+ (rows cap at 10; hps is the uncapped run
    length, so the mutant run of a deletion from an 11+ run is still
    10+), columns base*3 + (idLen + 1); str rows are STR bins, columns
    idLen + 5.
      1bp deletion from a run of L: forward = -1 in a run of L; reversion
        = +1 in a run of L-1, or str row 0 +1 when L == 1 (no run left).
      Run-extending 1bp insertion into a run of L: forward = +1 in a run
        of L; reversion = -1 in a run of L+1.
      Mismatched 1bp insertion of base x: forward = str row 0 +1;
        reversion = -1 of x as a run of 1 (str row 0 -1 if x is not ACGT).
      Longer indels: forward = idLen in the reference tract's bin (strs);
        reversion = -idLen in the mutant tract's bin (strs_mut, default
        strs when the tract length isn't known, e.g. per-context tables).
    """
    if abs(idLen) >= 2:
        if strs_mut is None:
            strs_mut = strs
        return ("str", strs_mut, -idLen + 5, -idLen + 5), (
            "str",
            strs,
            idLen + 5,
            idLen + 5,
        )
    L_raw = max(1, int(hps))
    L = min(L_raw, 10)
    if idLen == -1:
        b = ref_allele
        forward = ("hp", L - 1, b * 3 + 0, _RC[b] * 3 + 0)
        if L_raw >= 2:
            reversion = ("hp", min(L_raw - 1, 10) - 1, b * 3 + 2, _RC[b] * 3 + 2)
        else:
            reversion = ("str", 0, 6, 6)
    elif inserted_base == ref_allele:
        b = ref_allele
        forward = ("hp", L - 1, b * 3 + 2, _RC[b] * 3 + 2)
        reversion = ("hp", min(L_raw + 1, 10) - 1, b * 3 + 0, _RC[b] * 3 + 0)
    else:
        forward = ("str", 0, 6, 6)
        x = inserted_base
        if 0 <= x <= 3:
            reversion = ("hp", 0, x * 3 + 0, _RC[x] * 3 + 0)
        else:
            reversion = ("str", 0, 4, 4)
    return reversion, forward


def _cell_value(cell, mat_hp, mat_str, bottom=False):
    kind, row, col, col_bot = cell
    mat = mat_hp if kind == "hp" else mat_str
    return mat[row, col_bot if bottom else col]


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
    strs_mut=None,
):
    """Amplification and damage rates for a candidate indel, in the slots
    calculateSSPosterior/calculateDSPosterior expect (same convention as
    genotypeDSSnv): Pamp/Pdmg/Pdmg_bot are alt->ref (the reversion in the
    mutant context), Pamp_rev/Pdmg_rev/Pdmg_rev_bot are ref->alt (the
    indel in the reference context). See _indel_error_cells.

    ref_allele is the reference base at indel_context_index; for 1bp
    insertions, inserted_base == ref_allele selects the run-extending
    (hp) case, anything else the mismatched (str row 0) case.
    """
    if idLen == 0:
        raise ValueError("idLen must be nonzero")
    reversion, forward = _indel_error_cells(
        hps, strs, idLen, ref_allele, inserted_base, strs_mut
    )
    Pamp = _cell_value(reversion, prob_amp_hp, prob_amp_str)
    Pamp_rev = _cell_value(forward, prob_amp_hp, prob_amp_str)
    Pdmg = _cell_value(reversion, prob_dmg_hp, prob_dmg_str)
    Pdmg_rev = _cell_value(forward, prob_dmg_hp, prob_dmg_str)
    Pdmg_bot = _cell_value(reversion, prob_dmg_hp, prob_dmg_str, bottom=True)
    Pdmg_rev_bot = _cell_value(forward, prob_dmg_hp, prob_dmg_str, bottom=True)
    if Pamp == 0:
        Pamp = 1e-9
    return Pamp, Pamp_rev, Pdmg, Pdmg_rev, Pdmg_bot, Pdmg_rev_bot


def indelReversionRate(
    hps,
    strs,
    idLen,
    ref_allele,
    inserted_base,
    prob_amp_hp,
    prob_amp_str,
    strs_mut=None,
):
    """Per-read probability that a read from a molecule carrying this indel
    shows the reference: the amplification reversion of _indel_error_cells
    (indelErrorProbs' Pamp without its zero floor). NaN when the matrices
    aren't loaded."""
    if prob_amp_hp is None or prob_amp_str is None:
        return np.nan
    reversion, _ = _indel_error_cells(
        hps, strs, idLen, ref_allele, inserted_base, strs_mut
    )
    return float(_cell_value(reversion, prob_amp_hp, prob_amp_str))


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
        # Same filter as profileTriNucMismatches at learn time.
        if (
            indel_passes_mask(antimask, refPos - start, indelLen)
            and indel_has_context(refPos - start, len(reference_int))
            and indel_context_fits_window(
                indel, refPos - start, reference_int, hp_raw, str_raw
            )
        ):
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
        return [np.zeros(0)] * 11
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
        seqArr, qualArr = getIndelArr(
            seq, indels_masked, params["minBq"], reference_int, reference_start
        )
        mask_multiallele[seqArr == INDEL_CONFLICT] = 0
        f1r2_seq[nn, :] = seqArr == INDEL_ALT
        f1r2_prob[nn, :] = qualArr
        f1r2_alt_count += (seqArr == INDEL_ALT).astype(int)
        f1r2_ref_count += (seqArr == INDEL_REF).astype(int)
    for nn, seq in enumerate(F2R1):
        seqArr, qualArr = getIndelArr(
            seq, indels_masked, params["minBq"], reference_int, reference_start
        )
        mask_multiallele[seqArr == INDEL_CONFLICT] = 0
        f2r1_seq[nn, :] = seqArr == INDEL_ALT
        f2r1_prob[nn, :] = qualArr
        f2r1_alt_count += (seqArr == INDEL_ALT).astype(int)
        f2r1_ref_count += (seqArr == INDEL_REF).astype(int)
    f1r2_prob = -f1r2_prob / 10
    f2r1_prob = -f2r1_prob / 10

    hps = np.zeros(pos_masked.size, dtype=int)
    strs = np.zeros(pos_masked.size, dtype=int)
    Pamp = np.zeros(pos_masked.size)
    Pamp_rev = np.zeros(pos_masked.size)
    Pdmg = np.zeros(pos_masked.size)
    Pdmg_rev = np.zeros(pos_masked.size)
    Pdmg_bot = np.zeros(pos_masked.size)
    Pdmg_rev_bot = np.zeros(pos_masked.size)
    # Per-read reversion rate for the strand independence test.
    rev_rate = np.zeros(pos_masked.size)
    # False only for a 1bp insertion whose base differs from the next
    # reference base (routed to str.txt row 0); Caller.py uses it to pick
    # the call's channel.
    hp_match_arr = np.ones(pos_masked.size, dtype=bool)
    for nn in range(pos_masked.size):
        # HP and STR context, read independently at indel_context_index
        # (same position as learning).
        anchor = indel_context_index(pos_masked[nn] - start)
        # Uncapped run length for the error cells (a deletion from an 11+
        # run leaves a 10+ run); hps (the HP INFO field / channel) caps at 10.
        hp_run = int(hp_raw[0, anchor])
        hps[nn] = min(hp_run, 10)
        # STR bins (str.txt rows): 0 = not a repeat, 1 = 2-9bp,
        # 2 = 10-24bp, 3 = 25-39bp, 4 = 40+bp.
        idLen = indelLen_masked[nn]
        pos = pos_masked[nn]
        indel_parts = indels_masked[nn].split(":")
        raw_len = int(indel_parts[1])
        # A >=2bp indel that isn't a slip of the tract at its context base
        # (str_slip_tract) is STR0 ("not a repeat").
        slip = str_slip_tract(
            raw_len,
            indel_parts[2] if raw_len > 0 else "",
            anchor,
            reference_int,
            str_raw,
        )
        if slip is not None:
            unit_len_here, total_len = slip
            strs[nn] = str_length_bin(unit_len_here, total_len)
            # Mutant tract (reference tract plus the uncapped indel
            # length): its bin holds the reversion rate.
            strs_mut = str_length_bin(unit_len_here, total_len + raw_len)
        else:
            strs[nn] = 0
            strs_mut = 0
        # ref_allele (1bp indels only): the deleted base, or the reference
        # base right after the insertion point.
        if idLen == 1 or idLen == -1:
            ref_allele = int(reference_int[anchor])
        else:
            ref_allele = 0
        # inserted_base (1bp insertions): compared with ref_allele to pick
        # hp.txt (run-extending) or str.txt row 0 (mismatched).
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
            hp_run,
            strs[nn],
            idLen,
            ref_allele,
            inserted_base,
            prob_amp_hp,
            prob_dmg_hp,
            prob_amp_str,
            prob_dmg_str,
            strs_mut,
        )
        rev_rate[nn] = indelReversionRate(
            hp_run,
            strs[nn],
            idLen,
            ref_allele,
            inserted_base,
            prob_amp_hp,
            prob_amp_str,
            strs_mut,
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
    LR_masked[np.abs(LR_masked) < LR_ZERO_TOL] = 0.0
    # Negative-LR candidates are returned too, for the per-channel mu solve.
    take = mask_multiallele != 0
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
        rev_rate[take],
    )
