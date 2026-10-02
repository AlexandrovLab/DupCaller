"""Indel error learning (profileTriNucMismatches) gates opportunity and
events the same way:

- a read counts as amp opportunity at a run start only if it could show
  the reference there: base above minBq and aligned through the run plus
  one base (getIndelArr's REF rule), so reads ending inside the run don't;
- an event is booked only where the opportunity pass credited its context
  base, and only ALT reads that were credited are moved out of the
  idLen=0 column (no negative cells);
- an STR-table damage event needs an ALT read on one strand.
"""

import numpy as np
import pysam

from DupCaller_sub.funcs.learn import _indel_informative_reads, profileTriNucMismatches

CHROM = "chr1"
WIN_START = 100
_B2N = {"A": 0, "T": 1, "C": 2, "G": 3}
# G C [A A A] G C G G ...: the only A/T homopolymer of length 3 starts at 102.
REF = "GC" + "AAA" + "GCGGCCGCGGCCGCGGCCGCGGCC"
RUN_START = 102
HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
)
A = _B2N["A"]


def _ref_int(seq):
    return np.array([_B2N[b] for b in seq])


def _hp_raw(seq):
    ref = _ref_int(seq)
    cut = np.ones(len(ref), dtype=bool)
    cut[1:] = ref[1:] != ref[:-1]
    run_id = np.cumsum(cut) - 1
    return np.vstack((np.bincount(run_id)[run_id], cut)).astype(float)


def _read(name, is_read1, end=None, ins_after=None, ins_base="A"):
    """Read from WIN_START to `end` (exclusive, default the window end);
    ins_after inserts ins_base after that reference position."""
    end = WIN_START + len(REF) if end is None else end
    ref_seq = REF[: end - WIN_START]
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    if ins_after is None:
        rec.query_sequence = ref_seq
        rec.cigar = [(0, len(ref_seq))]
    else:
        lead = ins_after - WIN_START + 1
        rec.query_sequence = ref_seq[:lead] + ins_base + ref_seq[lead:]
        cigar = [(0, lead), (1, 1)]
        if len(ref_seq) > lead:
            cigar.append((0, len(ref_seq) - lead))
        rec.cigar = cigar
    rec.query_qualities = pysam.qualitystring_to_array(
        chr(40 + 33) * len(rec.query_sequence)
    )
    rec.reference_id = 0
    rec.reference_start = WIN_START
    rec.mapping_quality = 60
    rec.is_paired = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = False
    return rec


def _learn(family, antimask=None):
    return profileTriNucMismatches(
        seqs=family,
        reference_start=WIN_START,
        reference_int=_ref_int(REF),
        trinuc_int=np.zeros(len(REF), dtype=int),
        hp_raw=_hp_raw(REF),
        str_raw=np.zeros([3, len(REF)]),
        antimask=np.ones(len(REF), dtype=bool) if antimask is None else antimask,
        params={"trinuc2num_dict": {}, "minBq": 10, "minRef": 1, "minAlt": 1},
    )


def _hp3_a(result):
    hp_alt_count = result[1]
    return tuple(hp_alt_count[2, A * 3 : A * 3 + 3])


def test_informative_reads_need_the_whole_span():
    covered = np.array([[True] * 6, [True, True, True, True, False, False]])
    hq = np.ones_like(covered)
    span_end = np.arange(6) + np.array([0, 0, 2, 0, 0, 0])
    inf = _indel_informative_reads(covered, hq, {"t": span_end})["t"]
    # Position 2 needs 1..4 covered: read 0 yes, read 1 no (4 uncovered).
    assert inf[:, 2].tolist() == [True, False]
    # Positions 0 and 1 are never informative (anchor at the window start).
    assert not inf[:, :2].any()


def test_read_ending_inside_run_is_not_opportunity():
    # F1R2: 3 full reads + 1 ending at 104 (inside the run 102-104, so it
    # can't show the base after the run); F2R1: 3 full reads.
    family = [_read(f"t{i}", True) for i in range(3)]
    family.append(_read("t3", True, end=RUN_START + 2))
    family += [_read(f"b{i}", False) for i in range(3)]
    assert _hp3_a(_learn(family)) == (0, 6, 0)


def test_event_at_uncredited_context_is_not_booked():
    # +A into the run, left-aligned to anchor 101 / context 102; context 102
    # is masked (e.g. trimmed), the anchor isn't, so the candidate exists
    # but its context was never credited as opportunity.
    family = [_read(f"t{i}", True) for i in range(3)]
    family.append(_read("t3", True, ins_after=RUN_START + 1))
    family += [_read(f"b{i}", False) for i in range(3)]
    antimask = np.ones(len(REF), dtype=bool)
    antimask[RUN_START - WIN_START] = False
    assert _hp3_a(_learn(family, antimask)) == (0, 0, 0)


def test_credited_insertion_moves_its_read_out_of_reference_column():
    family = [_read(f"t{i}", True) for i in range(3)]
    family.append(_read("t3", True, ins_after=RUN_START + 1))
    family += [_read(f"b{i}", False) for i in range(3)]
    # 4 + 3 reads span the run; the ALT read moves to the +1 column.
    assert _hp3_a(_learn(family)) == (0, 6, 1)


def test_str_damage_needs_an_alt_read():
    # An unanchored end-of-read insertion of T (mismatched, str.txt row 0)
    # makes a candidate, but every informative read is REF: no damage
    # event.
    anchor = WIN_START + 15
    family = [_read(f"t{i}", True) for i in range(3)]
    family.append(_read("t3", True, end=anchor + 1, ins_after=anchor, ins_base="T"))
    family += [_read(f"b{i}", False) for i in range(3)]
    str_dmg_count = _learn(family)[5]
    assert str_dmg_count[0, 6] == 0
    # The opportunity column is the only thing credited.
    assert str_dmg_count[0, 5] > 0


# ------------------------------------------- STR0 rule and SBS read starts

# G C [CA x6] G G T C ...: a (CA)6 tract (unit 2, 12bp, bin 2) at 102-113.
STR_REF = "GC" + "CA" * 6 + "GGTCGGTCGGTCGGTCGGTCGGTCGG"
TRACT_START = 102


def _str_raw(seq, start, unit, count):
    raw = np.zeros([3, len(seq)])
    lo = start - WIN_START
    raw[0, lo : lo + unit * count] = unit
    raw[1, lo : lo + unit * count] = count
    raw[2, lo] = 1
    return raw


def _str_read(name, is_read1, del_at=None, del_len=0):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    if del_at is None:
        rec.query_sequence = STR_REF
        rec.cigar = [(0, len(STR_REF))]
    else:
        lead = del_at - WIN_START
        rec.query_sequence = STR_REF[:lead] + STR_REF[lead + del_len :]
        rec.cigar = [(0, lead), (2, del_len), (0, len(STR_REF) - lead - del_len)]
    rec.query_qualities = pysam.qualitystring_to_array(
        chr(40 + 33) * len(rec.query_sequence)
    )
    rec.reference_id = 0
    rec.reference_start = WIN_START
    rec.mapping_quality = 60
    rec.is_paired = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = False
    return rec


def _learn_str(family):
    return profileTriNucMismatches(
        seqs=family,
        reference_start=WIN_START,
        reference_int=_ref_int(STR_REF),
        trinuc_int=np.zeros(len(STR_REF), dtype=int),
        hp_raw=_hp_raw(STR_REF),
        str_raw=_str_raw(STR_REF, TRACT_START, 2, 6),
        antimask=np.ones(len(STR_REF), dtype=bool),
        params={"trinuc2num_dict": {}, "minBq": 10, "minRef": 1, "minAlt": 1},
    )


def test_whole_unit_str_deletion_is_learned_in_its_bin():
    family = [_str_read(f"t{i}", True) for i in range(3)]
    family.append(_str_read("t3", True, del_at=TRACT_START + 4, del_len=2))
    family += [_str_read(f"b{i}", False) for i in range(3)]
    str_alt_count = _learn_str(family)[2]
    assert str_alt_count[2, -2 + 5] == 1
    assert str_alt_count[0, -2 + 5] == 0


def test_non_unit_str_deletion_is_learned_as_str0():
    # A 3bp deletion inside (CA)6 isn't a slip of the dinucleotide repeat.
    family = [_str_read(f"t{i}", True) for i in range(3)]
    family.append(_str_read("t3", True, del_at=TRACT_START + 5, del_len=3))
    family += [_str_read(f"b{i}", False) for i in range(3)]
    str_alt_count = _learn_str(family)[2]
    assert str_alt_count[0, -3 + 5] == 1
    assert str_alt_count[2, -3 + 5] == 0


def _short_str_read(name, is_read1, end):
    rec = pysam.AlignedSegment(HEADER)
    rec.query_name = name
    rec.query_sequence = STR_REF[: end - WIN_START]
    rec.cigar = [(0, end - WIN_START)]
    rec.query_qualities = pysam.qualitystring_to_array(chr(40 + 33) * (end - WIN_START))
    rec.reference_id = 0
    rec.reference_start = WIN_START
    rec.mapping_quality = 60
    rec.is_paired = True
    rec.is_read1 = is_read1
    rec.is_read2 = not is_read1
    rec.is_reverse = False
    return rec


def test_no_event_without_a_credited_reference_read():
    # +GG right before the (CA)6 tract: a whole-unit length, so STR bin 2,
    # but it doesn't slide through the tract, so reads ending inside the
    # tract are REF for it while none of them spans the tract (no STR
    # opportunity). Without a credited non-ALT read the event isn't booked.
    family = [_short_str_read(f"t{i}", True, TRACT_START + 8) for i in range(3)]
    ins = pysam.AlignedSegment(HEADER)
    ins.query_name = "t3"
    ins.query_sequence = STR_REF[:2] + "GG" + STR_REF[2:]
    ins.cigar = [(0, 2), (1, 2), (0, len(STR_REF) - 2)]
    ins.query_qualities = pysam.qualitystring_to_array(
        chr(40 + 33) * (len(STR_REF) + 2)
    )
    ins.reference_id = 0
    ins.reference_start = WIN_START
    ins.mapping_quality = 60
    ins.is_paired = True
    ins.is_read1 = True
    ins.is_read2 = False
    ins.is_reverse = False
    family.append(ins)
    family += [_short_str_read(f"b{i}", False, TRACT_START + 8) for i in range(3)]
    str_alt_count = _learn_str(family)[2]
    assert str_alt_count[2, 2 + 5] == 0
    assert (str_alt_count >= 0).all()


def test_positions_before_a_read_start_are_not_counted_as_a():
    # C/G-only reference; on each strand two reads start at the window
    # start and two start 15bp later. The first 15bp still have 2 real
    # reads per strand and must keep their SBS amp opportunity.
    ref = "CGGC" * 10
    win = len(ref)

    def read(name, is_read1, start):
        rec = pysam.AlignedSegment(HEADER)
        rec.query_name = name
        rec.query_sequence = ref[start:]
        rec.cigar = [(0, win - start)]
        rec.query_qualities = pysam.qualitystring_to_array(chr(40 + 33) * (win - start))
        rec.reference_id = 0
        rec.reference_start = WIN_START + start
        rec.mapping_quality = 60
        rec.is_paired = True
        rec.is_read1 = is_read1
        rec.is_read2 = not is_read1
        rec.is_reverse = False
        return rec

    family = [read(f"t{i}", True, 0) for i in range(2)]
    family += [read(f"t{i + 2}", True, 15) for i in range(2)]
    family += [read(f"b{i}", False, 0) for i in range(2)]
    family += [read(f"b{i + 2}", False, 15) for i in range(2)]
    ref_int = _ref_int(ref)
    trinuc = np.zeros(win, dtype=int)
    for i in range(1, win - 1):
        trinuc[i] = 16 * ref_int[i - 1] + 4 * ref_int[i] + ref_int[i + 1]
    out = profileTriNucMismatches(
        seqs=family,
        reference_start=WIN_START,
        reference_int=ref_int,
        trinuc_int=trinuc,
        hp_raw=_hp_raw(ref),
        str_raw=np.zeros([3, win]),
        antimask=np.ones(win, dtype=bool),
        params={
            "trinuc2num_dict": {},
            "minBq": 10,
            "minRef": 1,
            "minAlt": 1,
            "srdMinRead": 2,
            "ssmMinRead": 2,
        },
    )
    # SRD counts: every covered base of every read, on both strands, and
    # each count appears twice after the reverse-complement fold. The first
    # 15 positions have 2 reads per strand, the rest 4. Counting the
    # later-starting reads as "A" there used to drop those 15 positions.
    assert out[0].sum() == 2 * 2 * (2 * 15 + 4 * 25)
