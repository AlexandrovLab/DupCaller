import numpy as np
from .indels import findIndels, getIndelArr
import scipy.special as sp


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
    ll1 = sp.logsumexp(
        np.vstack((PA_At + PA_Ab, PA_Bt + PA_Ab, PA_At + PA_Bb, PA_Bt + PA_Bb)), axis=0
    ) / np.log(10)
    ll2 = sp.logsumexp(
        np.vstack((PB_At + PB_Ab, PB_Bt + PB_Ab, PB_At + PB_Bb, PB_Bt + PB_Bb)), axis=0
    ) / np.log(10)
    return ll1, ll2


def calculateSSPosterior(P, P_rev, bin_seq, Pseq):  # countb1, countb2, Pb1, Pb2):
    Pseq[Pseq == 0] = log(0.5)
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
    seqs, reference_start, reference_int, trinuc_int, prior_mat, antimask, params
):
    prob_amp_mat = params["ampmat"]
    prob_amp_mat_rev = params["ampmat_rev"]
    prob_dmg_mat_top = params["dmgmat_top"]
    prob_dmg_mat_rev_top = params["dmgmat_rev_top"]
    prob_dmg_mat_bot = params["dmgmat_bot"]
    prob_dmg_mat_rev_bot = params["dmgmat_rev_bot"]
    trinuc_convert_np = params["trinuc_convert"]
    antimask[trinuc_int >= 64] = False
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
    antimask[
        ((F1R2_qual_mat_0_count + F2R1_qual_mat_0_count) / (m_F1R2 + m_F2R1))
        >= params["maxZeroQualFrac"]
    ] = False
    antimask[(total_count_mat >= 1).sum(axis=0) > 2] = False
    base1_int = np.argmax(total_count_mat, axis=0)
    total_count_without_base1 = total_count_mat.copy()
    total_count_without_base1[base1_int, np.ogrid[:n]] = -1
    base2_int = np.argmax(total_count_without_base1, axis=0)
    base2_int[base1_int != reference_int] = reference_int[base1_int != reference_int]
    F1R2_masked_qual_mat = F1R2_qual_mat[:, antimask]
    F2R1_masked_qual_mat = F2R1_qual_mat[:, antimask]
    F1R2_masked_seq_mat = F1R2_seq_mat[:, antimask]
    F2R1_masked_seq_mat = F2R1_seq_mat[:, antimask]
    F1R2_prob = -F1R2_masked_qual_mat / 10
    F2R1_prob = -F2R1_masked_qual_mat / 10
    base1_int_masked = base1_int[antimask]
    base2_int_masked = base2_int[antimask]
    F1R2_bin_seq_mat = F1R2_masked_seq_mat == base1_int_masked
    F2R1_bin_seq_mat = F2R1_masked_seq_mat == base1_int_masked
    trinuc_converted_masked = trinuc_convert_np[trinuc_int[antimask], base1_int_masked]
    ref_int_masked = reference_int[antimask]
    base2_int_masked[
        np.logical_and(
            base1_int_masked == ref_int_masked,
            total_count_mat[:, antimask][
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

    F1R2_count_b1 = F1R2_count_mat[:, antimask][
        base1_int_masked, np.ogrid[: base1_int_masked.size]
    ]
    F1R2_count_b2 = np.zeros(F1R2_count_b1.size)
    alt_pos = base2_int_masked != 4
    F1R2_count_b2[alt_pos] = F1R2_count_mat[:, antimask][:, alt_pos][
        base2_int_masked[alt_pos], np.ogrid[: np.count_nonzero(alt_pos)]
    ]
    F2R1_count_b1 = F2R1_count_mat[:, antimask][
        base1_int_masked, np.ogrid[: base1_int_masked.size]
    ]
    F2R1_count_b2 = np.zeros(F2R1_count_b1.size)
    alt_pos = base2_int_masked != 4
    F2R1_count_b2[alt_pos] = F2R1_count_mat[:, antimask][:, alt_pos][
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
    LR_diff = LR_max - LR_masked
    # LR_diff[LR_diff < 0] = 0
    LR_diff[LR_diff < 0] = np.finfo(float).eps
    LR_score = -log10(LR_diff / LR_max)
    LR_abs = np.zeros(n)
    LR_abs[antimask] = LR_masked
    LR = np.zeros(n)
    LR[antimask] = LR_score
    LR[LR_abs <= 0] = 0
    return (
        LR,
        base1_int,
        antimask,
        F1R2_count_mat,
        F2R1_count_mat,
    )


def genotypeDSIndel(
    seqs, reference_start, reference_end, bam, antimask, hp_int, params
):
    prob_amp = params["ampmat_indel"]
    prob_amp_rev = params["ampmat_indel_rev"]
    prob_dmg = params["dmgmat_indel"]
    prob_dmg_rev = params["dmgmat_indel_rev"]
    F1R2 = []
    F2R1 = []
    for seq in seqs:
        if (seq.is_read1 and seq.is_forward) or (seq.is_read2 and seq.is_reverse):
            F1R2.append(seq)
        if (seq.is_read2 and seq.is_forward) or (seq.is_read1 and seq.is_reverse):
            F2R1.append(seq)
    chrom = seqs[0].reference_name
    # start = seqs[0].reference_start
    # end = seqs[0].reference_end
    start = reference_start
    end = reference_end
    indels = set()
    ### Geonotype indel for all found indels
    for seq in F1R2:
        indels.update(findIndels(seq))
    for seq in F2R1:
        indels.update(findIndels(seq))
    # start = seqs[0].reference_start
    indels = list(indels)
    indels_masked = list()
    pos_masked = list()
    indelLen_masked = list()
    for indel in indels:
        refPos = int(indel.split(":")[0])
        indelLen = int(indel.split(":")[1])
        if antimask[refPos - start]:
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
    if m == 0:
        return [np.zeros(0)] * 7
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
    for nn in range(pos_masked.size):
        hps[nn] = np.max(
            hp_int[pos_masked[nn] - start : pos_masked[nn] + 1 + offset[nn] - start,]
        )
    hps[hps >= 40] = 40
    Pamp = prob_amp[hps - 1, indelLen_masked + 5]
    Pamp_rev = prob_amp_rev[hps - 1, indelLen_masked + 5]
    Pdmg = prob_dmg[hps - 1, indelLen_masked + 5]
    Pdmg_rev = prob_dmg_rev[hps - 1, indelLen_masked + 5]

    Pdmg[Pdmg == 0] = 1e-9
    Pdmg_rev[Pdmg_rev == 0] = 1e-9

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
        Pdmg,
        Pdmg_rev,
        F1R2_alt_prob,
        F2R1_alt_prob,
        F1R2_ref_prob,
        F2R1_ref_prob,
    )
    LR_masked = LL_B1 - LL_B2
    LR_max = log10(1 - Pdmg) + log10(1 - Pdmg) - log10(Pdmg_rev) - log10(Pdmg_rev)
    LR_diff = LR_max - LR_masked
    LR_score = -log10(LR_diff / LR_max)
    take = LR_masked >= 0
    take[mask_multiallele == 0] = 0
    return (
        # LR_diff[take],
        LR_score[take],
        [indels_masked[nn] for nn in range(len(take)) if take[nn]],
        hps[take],
        f1r2_ref_count[take].astype("int"),
        f1r2_alt_count[take].astype("int"),
        f2r1_ref_count[take].astype("int"),
        f2r1_alt_count[take].astype("int"),
    )
