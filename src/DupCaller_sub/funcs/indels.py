import numpy as np

_BASE2NUM = {"A": 0, "T": 1, "C": 2, "G": 3}
_NUM2BASE = "ATCG"


def left_align_indel(indel, reference_int, reference_start):
    """Shift a raw CIGAR-derived indel string (findIndels output) to its
    canonical leftmost representation, the same normalization bcftools
    norm/vt normalize apply as a post-alignment step -- BWA (like most
    aligners) does not left-align indels itself, so within a repeat run,
    different reads covering the identical true event can report it at
    different raw CIGAR positions. Left-aligning here, before any indel
    ever enters a set()-based consensus (genotypeDSIndel, funcs/prob.py;
    profileTriNucMismatches, funcs/learn.py), makes reads of the same
    physical event converge onto one canonical string instead of being
    treated as distinct indels, and guarantees the "before" side of any
    downstream repeat-count scan (Estimate.py's classify_indel_record/
    classify_indel_channel) is always empty.

    Standard single-base-at-a-time shift: for a deletion, shifting the
    del_len-base deleted window one position left produces an identical
    resulting sequence exactly when the base newly excluded from the
    window (the old last deleted base) equals the base newly included
    (one before the old anchor) -- ref[anchor] == ref[anchor+del_len].
    For an insertion, the analogous condition is ref[anchor] == S[-1]
    (the last inserted base), with the inserted sequence rotating by one
    base each shift (S -> S[-1]+S[:-1]) to keep representing the same net
    change.

    indel: raw indel string from findIndels -- "{pos}:{len}:{seq}" for an
        insertion or "{pos}:-{len}" for a deletion. pos is the 0-based
        genomic anchor (last matched base before the event, VCF POS
        convention).
    reference_int: base2num-encoded reference array (A=0,T=1,C=2,G=3;
        anything else treated as non-matching/ambiguous) for the window
        containing this indel, starting at reference_start.
    reference_start: genomic position reference_int[0] corresponds to.

    Shifting simply stops early if it would run past the left edge of
    reference_int, or (for a deletion) past the right edge when reading
    reference_int[anchor + del_len] -- identical in spirit to how repeat-context
    computations elsewhere in this codebase (e.g. call.py's
    last_cut_valid) already accept degraded results right at a
    processing-window boundary, since real repeat runs are far shorter
    than a window.
    """
    parts = indel.split(":")
    pos = int(parts[0])
    length = int(parts[1])
    if length < 0:
        del_len = -length
        anchor = pos - reference_start
        while (
            anchor >= 0
            and anchor + del_len < len(reference_int)
            and 0 <= reference_int[anchor] <= 3
            and 0 <= reference_int[anchor + del_len] <= 3
            and reference_int[anchor] == reference_int[anchor + del_len]
        ):
            anchor -= 1
        return f"{anchor + reference_start}:{length}"
    else:
        seq_nums = [_BASE2NUM.get(b, -1) for b in parts[2]]
        anchor = pos - reference_start
        while (
            anchor >= 0
            and seq_nums[-1] != -1
            and 0 <= reference_int[anchor] <= 3
            and reference_int[anchor] == seq_nums[-1]
        ):
            seq_nums = [seq_nums[-1]] + seq_nums[:-1]
            anchor -= 1
        new_seq = "".join(_NUM2BASE[n] if 0 <= n <= 3 else "N" for n in seq_nums)
        return f"{anchor + reference_start}:{length}:{new_seq}"


def findIndels(seq):
    refPos = seq.reference_start
    readPos = 0
    indels = list()
    for cigar in seq.cigartuples:
        if cigar[0] == 0:
            refPos += cigar[1]
            readPos += cigar[1]
        if cigar[0] in [3, 4]:
            readPos += cigar[1]
        if cigar[0] == 1:
            sequence = seq.query_sequence[readPos : readPos + cigar[1]]
            indels.append(f"{refPos-1}:{cigar[1]}:{sequence}")
            readPos += cigar[1]
        if cigar[0] == 2:
            indels.append(f"{refPos-1}:-{cigar[1]}")
            refPos += cigar[1]
    return indels


def getIndelArr(seq, indels):
    refPosList = seq.get_reference_positions(full_length=True)
    refPosListNoNone = [_ if _ else -1 for _ in refPosList]
    reference_positions = np.array(refPosListNoNone, dtype=int)
    """
    refQualArr = np.zeros(len(indels))
    altQualArr = np.zeros(len(indels))
    refCountArr = np.zeros(len(indels))
    altCountArr = np.zeros(len(indels))
    """
    seqArr = np.zeros(len(indels), dtype=int)
    qualArr = np.zeros(len(indels))
    for nn, indel in enumerate(indels):
        refPos = int(indel.split(":")[0])
        indelLen = int(indel.split(":")[1])
        if refPos >= seq.reference_end or refPos < seq.reference_start:
            seqArr[nn] = -1
            continue
        readPos = np.where(reference_positions == refPos)[0]
        if len(readPos) == 0 or readPos >= seq.query_length - 1:
            continue
        readPos = readPos[0]
        if indelLen > 0:
            if (reference_positions[readPos + 1 : readPos + indelLen + 1] == -1).all():
                seqArr[nn] = 1
                qualArr[nn] = np.average(
                    seq.query_qualities[readPos + 1 : readPos + indelLen + 1]
                )
            elif reference_positions[readPos + 1] - reference_positions[readPos] != -1:
                seqArr[nn] = 0
                qualArr[nn] = seq.query_qualities[readPos + 1]
            else:
                seqArr[nn] = -1
        if indelLen < 0 and reference_positions.size > readPos:
            if (
                reference_positions[readPos + 1] != -1
                and (reference_positions[readPos + 1] - reference_positions[readPos])
                == -indelLen + 1
            ):
                seqArr[nn] = 1
                qualArr[nn] = seq.query_qualities[readPos + 1]
            elif (
                reference_positions[readPos + 1] != -1
                and (reference_positions[readPos + 1] - reference_positions[readPos])
                == 1
            ):
                seqArr[nn] = 0
                qualArr[nn] = np.average(
                    seq.query_qualities[readPos + 1 : readPos - indelLen + 1]
                )
            else:
                seqArr[nn] = -1
    return seqArr, qualArr
