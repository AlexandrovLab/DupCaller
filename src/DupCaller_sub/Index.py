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
    from the num_units column — no span/unit_len division involved, unlike
    the old BED-derived path, since PERF already reports the exact count.
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

    # Whole-genome lookup tables, built once: turn per-base/per-trinuc
    # encoding into O(1) numpy fancy-indexing instead of a Python-level
    # per-base string list (previously the single largest memory consumer —
    # a new 3-character str object allocated for every base in the
    # chromosome) and a np.vectorize call (a disguised per-base Python
    # loop). Any byte combination not explicitly set below keeps the
    # table's default, reproducing the old dict.get(key, default) fallback
    # exactly — including for ambiguity codes other than "N".
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
    for h5_file in (ref_h5, tn_h5, hp_h5, str_h5):
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
                str_unit_len[start:end] = unit_length
                str_cut[start] = 1
                # Clipped to the array's uint8 range before storing — since
                # every position is written at most once, this is
                # identical to storing the unclipped count and clipping
                # the whole array afterward, just without the wider dtype
                # in between.
                str_hp[start:end] = min(count, 127)
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
    ref_h5.close()
    tn_h5.close()
    hp_h5.close()
    str_h5.close()
