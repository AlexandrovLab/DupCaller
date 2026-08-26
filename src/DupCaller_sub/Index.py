#!/usr/bin/env python3
import h5py
import gzip
import numpy as np
import os
from Bio import SeqIO

from .funcs.misc import INDEX_SCHEMA_VERSION


def _iter_repeat_tsv_by_chrom(path):
    """Stream a PERF-format repeat tsv (chrom, start, end, motif, length,
    strand, num_units, motif_repeat), grouped by chromosome, yielding
    (chrom, rows) with rows a list of (start, end, unit_len, repeat_count)
    tuples in file order.

    Relies on same-chromosome rows being contiguous in the file — true of
    PERF's own output, since it processes and concatenates one FASTA record
    at a time (see PERF/rep_utils.py's fasta_ssrs) — so this only needs a
    single forward pass, no sorting or indexing. A chromosome with zero
    repeats simply contributes no rows and is never yielded; do_index's
    merge loop below handles that by treating a non-matching/absent
    chromosome as "no STR annotation here" rather than erroring.

    unit_len is read directly from len(motif) and repeat_count directly
    from the num_units column — no span/unit_len division involved, since
    PERF already reports the exact count.
    Entries with unit_len==1 (homopolymers) are skipped entirely: those are
    self-derived from the reference sequence instead (see
    homopolymer_run_and_cut below), never taken from this tsv.
    """
    handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    try:
        current_chrom = None
        current_rows = []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if chrom != current_chrom:
                if current_chrom is not None:
                    yield current_chrom, current_rows
                current_chrom = chrom
                current_rows = []
            unit_length = len(fields[3])
            if unit_length <= 1:
                continue
            start = int(fields[1])
            end = int(fields[2])
            count = int(fields[6])
            current_rows.append((start, end, unit_length, count))
        if current_chrom is not None:
            yield current_chrom, current_rows
    finally:
        handle.close()


def _dinuc_int_for_chrom(reference_int, sentinel=16):
    """Per-position raw dinucleotide class (0-15 = 4*base_at_i +
    base_at_i+1, over the same A/T/C/G->0-3 encoding ref.h5 uses),
    pairing every position with its successor — mirrors
    Estimate.py's calculate_ref_dbs exactly (same `4*base1+base2`,
    same "either base >3 (N/ambiguous) is invalid" rule), just
    materialized per-position instead of summed into an aggregate
    count. The last position in a chromosome has no successor and gets
    `sentinel` (16, one past the valid 0-15 range), matching how
    Index.py's trinuc index uses a boundary "N" at contig edges.
    """
    n = reference_int.shape[0]
    dinuc_int = np.full(n, sentinel, dtype=np.uint8)
    if n >= 2:
        base1 = reference_int[:-1]
        base2 = reference_int[1:]
        valid = (base1 <= 3) & (base2 <= 3)
        dinuc_int[:-1][valid] = (4 * base1[valid] + base2[valid]).astype(np.uint8)
    return dinuc_int


def do_index_dbs(args):
    """Add a .dbs.h5 (per-position raw-16 reference dinucleotide class)
    to an already-indexed reference, derived read-only from its existing
    .ref.h5 — does not open .ref.h5/.tn.h5/.hp.h5/.str.h5 for writing and
    never touches them. Safe to run against a reference other `call`/
    `estimate` jobs currently have open, unlike re-running the full
    `index` command (which rewrites all four files and would race any
    live reader — see Caller.py's do_call comments on this exact hazard).
    Same per-position encoding do_index's main loop below also now
    produces for any *new* full index build (`_dinuc_int_for_chrom`).
    """
    ref_h5_path = args.reference + ".ref.h5"
    if not os.path.exists(ref_h5_path):
        raise FileNotFoundError(
            f"{ref_h5_path} not found -- run `DupCaller.py index` first"
        )
    ref_h5 = h5py.File(ref_h5_path, "r")
    dbs_h5 = h5py.File(args.reference + ".dbs.h5", "w")
    dbs_h5.attrs["index_schema_version"] = INDEX_SCHEMA_VERSION
    chroms = list(ref_h5.keys())
    for i, chrom in enumerate(chroms):
        print(f"[{i + 1}/{len(chroms)}] Creating dinucleotide (DBS) index: {chrom}")
        dinuc_int = _dinuc_int_for_chrom(ref_h5[chrom][:])
        dbs_h5.create_dataset(chrom, data=dinuc_int)
    ref_h5.close()
    dbs_h5.close()
    print(f"DBS index written to {args.reference}.dbs.h5")


def _iter_fasta_records(path):
    """Stream FASTA records one at a time instead of loading every contig
    into memory up front (SeqIO.to_dict keeps the whole genome resident for
    the life of the process)."""
    handle = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    try:
        yield from SeqIO.parse(handle, "fasta")
    finally:
        handle.close()


def do_index(args):
    ### Check if reference file exists
    if not os.path.exists(args.reference):
        raise FileNotFoundError(f"Reference file not found: {args.reference}")

    ### Define trinuc order, including N bases
    trinuc2num = dict()
    num2trinuc = list()
    trinuc_order = 0
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                trinuc = minus_base + ref_base + plus_base
                trinuc2num[trinuc] = trinuc_order
                num2trinuc.append(trinuc)
                trinuc_order += 1
    for plus_base in ["T", "A", "G", "C"]:
        for ref_base in ["G", "A"]:
            for minus_base in ["T", "A", "G", "C"]:
                trinuc = minus_base + ref_base + plus_base
                trinuc2num[trinuc] = trinuc_order
                num2trinuc.append(trinuc)
                trinuc_order += 1
    for ref_base in ["C", "T"]:
        for plus_base in ["A", "T", "C", "G"]:
            trinuc = "N" + ref_base + plus_base
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            trinuc = minus_base + ref_base + "N"
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    for ref_base in ["G", "A"]:
        for plus_base in ["T", "A", "G", "C"]:
            trinuc = "N" + ref_base + plus_base
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    for minus_base in ["T", "A", "G", "C"]:
        for ref_base in ["G", "A"]:
            trinuc = minus_base + ref_base + "N"
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3, "a": 0, "t": 1, "c": 2, "g": 3}

    # Whole-genome lookup tables, built once: O(1) numpy fancy-indexing for
    # per-base/per-trinuc encoding, avoiding a per-base Python string list
    # (a new 3-character str object allocated for every base in the
    # chromosome) or a np.vectorize call (a disguised per-base Python
    # loop). Any byte combination not explicitly set below keeps the
    # table's default, matching dict.get(key, default) fallback semantics
    # — including for ambiguity codes other than "N".
    base2num_lut = np.full(256, 4, dtype=np.uint8)
    for base, code in base2num.items():
        base2num_lut[ord(base)] = code
    trinuc_lut = np.full((256, 256, 256), 96, dtype=np.uint8)
    for trinuc, code in trinuc2num.items():
        minus_base, ref_base, plus_base = trinuc
        trinuc_lut[ord(minus_base), ord(ref_base), ord(plus_base)] = code
    boundary_n = ord("N")

    ### Generate reference and index matrix into h5py file
    ref_h5 = h5py.File(args.reference + ".ref.h5", "w")
    tn_h5 = h5py.File(args.reference + ".tn.h5", "w")
    hp_h5 = h5py.File(args.reference + ".hp.h5", "w")
    str_h5 = h5py.File(args.reference + ".str.h5", "w")
    dbs_h5 = h5py.File(args.reference + ".dbs.h5", "w")
    for h5_file in (ref_h5, tn_h5, hp_h5, str_h5, dbs_h5):
        h5_file.attrs["index_schema_version"] = INDEX_SCHEMA_VERSION
    repeat_tsv_iter = _iter_repeat_tsv_by_chrom(args.repeatTsv)

    def _next_or_none(it):
        try:
            return next(it)
        except StopIteration:
            return None

    # One item of lookahead: (chrom, rows) for whichever chromosome the tsv
    # stream is currently sitting on, or None once the tsv is exhausted.
    pending_tsv_chrom = _next_or_none(repeat_tsv_iter)

    def homopolymer_run_and_cut(reference_int):
        """Self-derive per-position homopolymer run length and run-start
        markers directly from the reference sequence — no BED input
        involved. cut[i] marks the start of every maximal run of identical
        bases. run_id[i] is a running count of run starts up to and
        including i, so every position in the same run shares one id;
        bincount over run_id then gives each run's length, indexed by that
        id. Clipping run_len to the uint8 range before the final gather
        (hp = run_len[run_id]) is equivalent to computing the full-width
        run length and clipping afterward — clip is elementwise and
        commutes with a gather — but avoids ever materializing a
        chromosome-length int64 array."""
        n = reference_int.shape[0]
        cut = np.ones(n, dtype=bool)
        cut[1:] = reference_int[1:] != reference_int[:-1]
        run_id = np.cumsum(cut, dtype=np.int32) - 1
        run_len = np.clip(np.bincount(run_id), 0, 127).astype(np.uint8)
        hp = run_len[run_id]
        return hp, cut

    for record in _iter_fasta_records(args.reference):
        chrom = record.id
        print(f"Currently processing:{chrom}")
        print(f"Creating reference sequence index")
        seq_bytes = np.frombuffer(
            str(record.seq).upper().encode("ascii", errors="replace"),
            dtype=np.uint8,
        )
        ref_len = seq_bytes.shape[0]
        reference_int = base2num_lut[seq_bytes]

        print(f"Creating trinucleotide index")
        minus_bytes = np.empty(ref_len, dtype=np.uint8)
        minus_bytes[0] = boundary_n
        minus_bytes[1:] = seq_bytes[:-1]
        plus_bytes = np.empty(ref_len, dtype=np.uint8)
        plus_bytes[:-1] = seq_bytes[1:]
        plus_bytes[-1] = boundary_n
        trinuc_int = trinuc_lut[minus_bytes, seq_bytes, plus_bytes]
        del minus_bytes, plus_bytes

        print(f"Creating dinucleotide (DBS) index")
        dinuc_int = _dinuc_int_for_chrom(reference_int)

        print(f"Creating homopolymer index")
        hp_run, hp_cut = homopolymer_run_and_cut(reference_int)
        # hp_lens rows: 0 = homopolymer run length (same value at every
        # position within a given run), 1 = cut (start-of-run bool)
        hp_lens = np.vstack((hp_run, hp_cut)).astype(np.uint8, copy=False)

        print(f"Creating STR repeat index")
        str_unit_len = np.zeros(ref_len, dtype=np.uint8)
        str_hp = np.zeros(ref_len, dtype=np.uint8)
        str_cut = np.zeros(ref_len, dtype=np.uint8)
        # Paint every (start, end, unit_len, count) row PERF found for this
        # chromosome directly onto the whole-chromosome arrays — unit_len
        # and count are taken as-is from the tsv (no span/unit_len
        # division), and PERF's own num_units is trusted exactly rather
        # than re-derived. If the tsv's next chromosome doesn't match this
        # one, it means PERF found zero repeats here (see
        # _iter_repeat_tsv_by_chrom's docstring) — leave the arrays at 0
        # and don't advance the tsv stream.
        if pending_tsv_chrom is not None and pending_tsv_chrom[0] == chrom:
            for start, end, unit_length, count in pending_tsv_chrom[1]:
                # PERF's reported span can overextend past the tract's true
                # length (start + num_units*unit_length) at the trailing
                # boundary -- clamp so any such trailing bases aren't
                # painted with the full tract's unit_len/count, which would
                # contaminate HP/STR opportunity counting and error-rate
                # learning at those positions.
                true_end = min(end, start + unit_length * count)
                str_unit_len[start:true_end] = unit_length
                str_cut[start] = 1
                # Clipped to the array's uint8 range before storing — since
                # every position is written at most once, this is
                # identical to storing the unclipped count and clipping
                # the whole array afterward, just without the wider dtype
                # in between.
                str_hp[start:true_end] = min(count, 127)
            pending_tsv_chrom = _next_or_none(repeat_tsv_iter)
        # str_lens rows: 0 = repeat unit length, 1 = number of times the
        # repeat unit repeats (constant across the whole repeat interval),
        # 2 = cut (start-of-repeat bool)
        str_lens = np.vstack((str_unit_len, str_hp, str_cut)).astype(
            np.uint8, copy=False
        )

        idx = ref_h5.create_dataset(chrom, data=reference_int)
        tn = tn_h5.create_dataset(chrom, data=trinuc_int)
        hp = hp_h5.create_dataset(chrom, data=hp_lens)
        st = str_h5.create_dataset(chrom, data=str_lens)
        db = dbs_h5.create_dataset(chrom, data=dinuc_int)
    ref_h5.close()
    tn_h5.close()
    hp_h5.close()
    str_h5.close()
    dbs_h5.close()
