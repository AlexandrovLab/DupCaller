import numpy as np

_BASE2NUM = {"A": 0, "T": 1, "C": 2, "G": 3}
_NUM2BASE = "ATCG"


# A/T/C/G -> 0-3, anything else (N) -> 4, for whole read sequences.
_BASE_CODE_LUT = np.full(256, 4, dtype=np.int64)
for _b, _c in (("A", 0), ("T", 1), ("C", 2), ("G", 3)):
    _BASE_CODE_LUT[ord(_b)] = _c


def base_codes(seq):
    """Base codes (A,T,C,G = 0-3, anything else 4) of a sequence string."""
    return _BASE_CODE_LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def str_unit_multiple(indel_len, unit_len):
    """True iff an indel of indel_len (signed, uncapped) is a whole number
    of repeat units of an STR with this unit length. 1bp indels always
    count (they're routed by homopolymer context, not STR bin)."""
    return abs(indel_len) < 2 or abs(indel_len) % unit_len == 0


def str_length_bin(unit_len, total_len):
    """STR bin (str.txt row) of a tract: 0 = not a repeat (unit < 2 or
    fewer than two units), 1 = 2-9bp, 2 = 10-24bp, 3 = 25-39bp, 4 = 40+bp."""
    if unit_len < 2 or total_len < 2 * unit_len:
        return 0
    if total_len >= 40:
        return 4
    if total_len >= 25:
        return 3
    if total_len >= 10:
        return 2
    return 1


def str_slip_tract(indel_len, inserted_seq, ctx, reference_int, str_raw):
    """(unit_len, total_len) of the STR tract a >=2bp indel is a slip of, or
    None if it isn't a slip (an STR0 event).

    A slip is a whole number of the tract's units, at the tract start (a
    left-aligned slip always lands there), and made of the tract itself: a
    deletion that stays inside the tract, or an insertion of copies of the
    tract's unit in phase. "GG" inserted into (AC)n, or a deletion running
    out of the tract, has the right length but is not a slip.

    ctx: indel_context_index of the indel; reference_int and str_raw (unit,
    repeat count, start flag rows) share its local frame. Shared by learning
    and calling."""
    if abs(indel_len) < 2 or ctx >= str_raw.shape[1]:
        return None
    unit_len = int(str_raw[0, ctx])
    if unit_len < 2 or not str_raw[2, ctx] or abs(indel_len) % unit_len:
        return None
    total_len = unit_len * int(str_raw[1, ctx])
    if indel_len < 0:
        return (unit_len, total_len) if -indel_len <= total_len else None
    if ctx + unit_len > len(reference_int) or len(inserted_seq) != indel_len:
        return None
    unit = reference_int[ctx : ctx + unit_len]
    inserted = [_BASE2NUM.get(b, -1) for b in inserted_seq.upper()]
    for i, b in enumerate(inserted):
        if b != unit[i % unit_len]:
            return None
    return unit_len, total_len


def left_align_indel(indel, reference_int, reference_start):
    """Shift a raw findIndels indel to its leftmost equivalent position
    (the bcftools norm convention), so reads reporting the same event at
    different places inside a repeat agree on one string.

    A deletion shifts left while ref[anchor] == ref[anchor + del_len]; an
    insertion shifts left while ref[anchor] == its last inserted base,
    rotating the inserted sequence by one base per shift. Shifting stops
    at the edges of reference_int.

    indel: "{pos}:{len}:{seq}" (insertion) or "{pos}:-{len}" (deletion);
        pos is the 0-based anchor (last reference base before the event).
    reference_int: base2num-encoded reference (A=0,T=1,C=2,G=3, other=4)
        starting at reference_start.
    """
    parts = indel.split(":")
    pos = int(parts[0])
    length = int(parts[1])
    if length < 0:
        del_len = -length
        anchor = pos - reference_start
        while (
            anchor > 0
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
            anchor > 0
            and anchor < len(reference_int)
            and seq_nums[-1] != -1
            and 0 <= reference_int[anchor] <= 3
            and reference_int[anchor] == seq_nums[-1]
        ):
            seq_nums = [seq_nums[-1]] + seq_nums[:-1]
            anchor -= 1
        new_seq = "".join(_NUM2BASE[n] if 0 <= n <= 3 else "N" for n in seq_nums)
        return f"{anchor + reference_start}:{length}:{new_seq}"


def indel_mask_span(local_pos, indel_len):
    """[lo, hi) local interval that must be unmasked for an indel anchored
    at local_pos (0-based, relative to the window's antimask): the anchor,
    plus every deleted base for a deletion or the context base (the first
    base after the insertion point) for an insertion -- the same bases the
    Eeff site masks require (funcs/misc.py's indel_eeff_site_masks). Shared
    by learning (funcs/learn.py) and calling (funcs/prob.py, funcs/call.py)
    so both apply the identical locus mask."""
    return local_pos, local_pos + max(-indel_len, 1) + 1


def indel_passes_mask(antimask, local_pos, indel_len):
    """True iff the whole indel_mask_span lies inside antimask and is
    unmasked. A span reaching outside the window fails rather than
    wrapping (negative index) or silently truncating (empty slice)."""
    lo, hi = indel_mask_span(local_pos, indel_len)
    if lo < 0 or hi > len(antimask):
        return False
    return bool(antimask[lo:hi].all())


def indel_context_index(local_pos):
    """Local index whose hp_raw/str_raw/reference_int value classifies an
    indel's HP/STR context: anchor + 1, i.e. the first deleted base of a
    deletion or the first reference base after an insertion. hp.h5 and
    str.h5 store the whole run's/tract's values at every position in it,
    and for a left-aligned event anchor + 1 is the run/tract start, where
    learning credits the opportunity. Shared by learning and calling."""
    return local_pos + 1


def indel_context_fits_window(indel, local_pos, reference_int, hp_raw, str_raw):
    """True iff the homopolymer run (1bp indels routed to hp.txt) or STR
    tract (a slip, str_slip_tract) an indel is scored in ends inside the
    window, the base after it included -- the runs/tracts the Eeff site
    masks count (hp_repeat_valid, str_tract_valid). Reads cut inside the
    run/tract can't show its length. Other indels (STR0) always fit.

    indel: left-aligned "{pos}:{len}[:{seq}]"; local_pos: its anchor in
    the window frame shared by reference_int, hp_raw and str_raw."""
    parts = indel.split(":")
    indel_len = int(parts[1])
    ctx = indel_context_index(local_pos)
    window_len = len(reference_int)
    if ctx >= window_len:
        return False
    if abs(indel_len) == 1:
        if indel_len == 1:
            inserted = _BASE2NUM.get(parts[2][:1].upper(), -1) if len(parts) > 2 else -1
            if inserted != reference_int[ctx]:
                return True
        return ctx + int(hp_raw[0, ctx]) < window_len
    slip = str_slip_tract(
        indel_len, parts[2] if indel_len > 0 else "", ctx, reference_int, str_raw
    )
    return slip is None or ctx + slip[1] < window_len


def findIndels(seq):
    refPos = seq.reference_start
    readPos = 0
    indels = list()
    for cigar in seq.cigartuples:
        if cigar[0] == 0:
            refPos += cigar[1]
            readPos += cigar[1]
        if cigar[0] == 3:
            # N: reference skip, unlike S below consumes ref not query.
            refPos += cigar[1]
        if cigar[0] == 4:
            readPos += cigar[1]
        if cigar[0] == 1:
            sequence = seq.query_sequence[readPos : readPos + cigar[1]]
            indels.append(f"{refPos-1}:{cigar[1]}:{sequence}")
            readPos += cigar[1]
        if cigar[0] == 2:
            indels.append(f"{refPos-1}:-{cigar[1]}")
            refPos += cigar[1]
    return indels


# Per-read evidence states for one candidate indel (getIndelArr).
INDEL_ALT = 1
INDEL_REF = 0
# The read carries a different indel inside the candidate's locus
# (genotypeDSIndel drops the candidate).
INDEL_CONFLICT = -1
# Informative read whose post-anchor BQ doesn't clear min_bq.
INDEL_LOW_BQ = -2
# The read can't establish either allele: it doesn't span the locus, ends
# or is soft-clipped inside it, or its only matching evidence is a soft
# clip / an unanchored end-of-read insertion.
INDEL_UNINFORMATIVE = -3


def _read_indel_events(seq):
    """findIndels, paired with whether each event has an aligned (M) block
    on both sides. Events at either end of the alignment are unanchored."""
    events = findIndels(seq)
    ops = [op for op, _ in seq.cigartuples if op in (0, 1, 2)]
    anchored = []
    for k, op in enumerate(ops):
        if op in (1, 2):
            anchored.append(0 in ops[:k] and 0 in ops[k + 1 :])
    return list(zip(events, anchored))


def _indel_right_end(local_pos, indel_len, inserted_seq, reference_int):
    """Last local reference position a read must align through
    (contiguously from the anchor) to show the reference allele: the base
    after the event's rightmost equivalent placement in a repeat. None if
    that position is outside reference_int."""
    n = len(reference_int)
    anchor = local_pos
    if indel_len < 0:
        del_len = -indel_len
        while (
            anchor + 1 + del_len < n
            and 0 <= reference_int[anchor + 1] <= 3
            and reference_int[anchor + 1] == reference_int[anchor + 1 + del_len]
        ):
            anchor += 1
        end = anchor + del_len + 1
    else:
        seq_nums = [_BASE2NUM.get(b, -1) for b in inserted_seq]
        while (
            anchor + 1 < n
            and seq_nums[0] != -1
            and 0 <= reference_int[anchor + 1] <= 3
            and reference_int[anchor + 1] == seq_nums[0]
        ):
            seq_nums = seq_nums[1:] + seq_nums[:1]
            anchor += 1
        end = anchor + 1
    return end if end < n else None


def _event_overlaps(event, lo, hi):
    """True iff raw indel `event` touches the open reference interval
    (lo, hi) between a candidate's anchor lo and its REF span end hi."""
    parts = event.split(":")
    a, length = int(parts[0]), int(parts[1])
    if length < 0:
        # deleted bases a+1..a-length
        return a + 1 < hi and a - length > lo
    # junction between a and a+1
    return lo <= a < hi


def getIndelArr(seq, indels, min_bq, reference_int, reference_start):
    """Per-candidate evidence of one read: returns (seqArr, qualArr), seqArr
    holding INDEL_ALT / INDEL_REF / INDEL_CONFLICT / INDEL_LOW_BQ /
    INDEL_UNINFORMATIVE per candidate in `indels` (left-aligned
    "pos:len[:seq]" strings), qualArr the representative BQ for ALT/REF
    reads (0 otherwise).

    ALT needs the read's own anchored CIGAR indel, left-aligned against
    the same reference_int, to equal the candidate. A different anchored
    indel inside the candidate's locus is a CONFLICT. REF needs the read
    to align contiguously from the anchor through _indel_right_end.
    Anything else (non-spanning read, soft clip, read ending inside the
    locus) is UNINFORMATIVE rather than a guess at either allele."""
    ref_pos = np.array(
        [-1 if p is None else p for p in seq.get_reference_positions(full_length=True)],
        dtype=int,
    )
    quals = seq.query_qualities
    own = {}
    for raw, anchored in _read_indel_events(seq):
        if anchored:
            own[left_align_indel(raw, reference_int, reference_start)] = raw

    def _median_bq(anchor, window_len):
        idx = np.nonzero(ref_pos == anchor)[0]
        if idx.size == 0:
            return 0.0
        i = idx[0]
        return float(np.median(quals[i + 1 : i + 1 + window_len]))

    seqArr = np.full(len(indels), INDEL_UNINFORMATIVE, dtype=int)
    qualArr = np.zeros(len(indels))
    for nn, indel in enumerate(indels):
        parts = indel.split(":")
        pos = int(parts[0])
        indel_len = int(parts[1])
        inserted_seq = parts[2] if len(parts) > 2 else ""
        if indel in own:
            state = INDEL_ALT
            bq_anchor = int(own[indel].split(":")[0])
        else:
            local_end = _indel_right_end(
                pos - reference_start, indel_len, inserted_seq, reference_int
            )
            span_end = (
                reference_start + local_end
                if local_end is not None
                else reference_start + len(reference_int)
            )
            if any(
                _event_overlaps(raw, pos, span_end)
                for key, raw in own.items()
                if key != indel
            ):
                seqArr[nn] = INDEL_CONFLICT
                continue
            if local_end is None:
                continue
            idx = np.nonzero(ref_pos == pos)[0]
            if idx.size == 0:
                continue
            i = idx[0]
            k = span_end - pos
            if i + k >= ref_pos.size or not np.array_equal(
                ref_pos[i : i + k + 1], np.arange(pos, span_end + 1)
            ):
                continue
            state = INDEL_REF
            bq_anchor = pos

        # BQ: median over the |indel_len| read bases after the anchor (the
        # inserted bases, or the bases after the deletion point), for ALT
        # and REF alike.
        median_bq = _median_bq(bq_anchor, abs(indel_len))
        if median_bq <= min_bq:
            seqArr[nn] = INDEL_LOW_BQ
            continue
        seqArr[nn] = state
        qualArr[nn] = median_bq
    return seqArr, qualArr


def indel_has_context(local_pos, window_len):
    """True iff the indel's anchor is not the window's first base (where
    left_align_indel may have been stopped by the window edge) and its
    context base (indel_context_index) is inside the window."""
    return local_pos >= 1 and indel_context_index(local_pos) < window_len
