#!/usr/bin/env python3
from gzip import open as gzopen
from itertools import zip_longest
import os
import re

_EOF = object()
_MATE_SUFFIX = re.compile(r"/[12]$")


def _open_fastq(path):
    """Open a FASTQ as text, detecting gzip from the file's own magic bytes
    (so each mate is detected independently, whatever its extension)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzopen(path, "rt")
    return open(path)


def _validate_pattern(pattern, barcodes=None):
    if not pattern or set(pattern) - {"N", "X", "B"}:
        raise ValueError(
            f"Invalid barcode pattern '{pattern}': it must be non-empty and use "
            "only N (barcode base), X (skipped base) and B (one barcode from "
            "--barcode-list, of any length)."
        )
    n_b = pattern.count("B")
    if n_b > 1:
        raise ValueError(f"Invalid barcode pattern '{pattern}': at most one B.")
    if n_b == 1 and barcodes is None:
        raise ValueError(
            f"Barcode pattern '{pattern}' has B, which needs --barcode-list."
        )
    if n_b == 0 and barcodes is not None:
        raise ValueError(
            f"--barcode-list was given but barcode pattern '{pattern}' has no B."
        )


def load_barcode_list(path):
    """Barcodes from a text file, one per line (blank lines and lines
    starting with # are ignored). Returns them in file order."""
    barcodes = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            bc = line.strip().upper()
            if not bc or bc.startswith("#"):
                continue
            if set(bc) - set("ACGT"):
                raise ValueError(
                    f"{path}: line {line_no}: barcode {bc!r} may only contain A, C, G, T."
                )
            barcodes.append(bc)
    if not barcodes:
        raise ValueError(f"{path}: no barcodes found.")
    dup = sorted({b for b in barcodes if barcodes.count(b) > 1})
    if dup:
        raise ValueError(f"{path}: duplicate barcodes: {', '.join(dup)}")
    return barcodes


NO_MATCH = None
AMBIGUOUS = ""


def match_barcode(seq, start, barcodes, max_mismatch):
    """Match the read bases from `start` against every listed barcode at
    once, walking one position at a time.

    Each barcode still in the running is compared with the read base at the
    current position and its mismatch count updated; it is dropped once that
    count exceeds max_mismatch, and finishes when the walk reaches its
    length (provided the read is long enough). The walk stops when no
    barcode is left running. Of the finished barcodes the one with the
    fewest mismatches wins, then the shortest; two different barcodes still
    tied after that make the read AMBIGUOUS. Returns (barcode, length),
    (NO_MATCH, 0) or (AMBIGUOUS, 0).

    Cost: at most len(barcodes) * max barcode length base comparisons per
    read, whatever max_mismatch is (dropping only makes it cheaper). Used
    for max_mismatch > 1; below that, match_barcode_indexed gives the same
    result from lookup tables.
    """
    lens = [len(b) for b in barcodes]
    alive = range(len(barcodes))
    mism = [0] * len(barcodes)
    finished = []
    n = len(seq)
    i = 0
    while alive and start + i < n:
        base = seq[start + i]
        nxt = []
        for k in alive:
            if barcodes[k][i] != base:
                mism[k] += 1
                if mism[k] > max_mismatch:
                    continue
            if i + 1 == lens[k]:
                finished.append((mism[k], lens[k], k))
            else:
                nxt.append(k)
        alive = nxt
        i += 1
    if not finished:
        return NO_MATCH, 0
    best_mm, best_len, best_k = min(finished)
    if sum(1 for mm, ln, _ in finished if mm == best_mm and ln == best_len) > 1:
        return AMBIGUOUS, 0
    return barcodes[best_k], best_len


def _read_base_name(header, path, record_no):
    """Read name without the leading '@', any comment, or a trailing /1 or /2."""
    header = header.rstrip("\r\n")
    if not header.startswith("@") or len(header) < 2:
        raise ValueError(
            f"{path}: record {record_no + 1} header is not a FASTQ header: {header!r}"
        )
    return _MATE_SUFFIX.sub("", header[1:].split(None, 1)[0])


def _check_read(name, seq, qual, pattern_len, path, record_no):
    if len(seq) != len(qual):
        raise ValueError(
            f"{path}: read {name} (pair {record_no + 1}) has sequence length "
            f"{len(seq)} but quality length {len(qual)}."
        )
    if len(seq) < pattern_len:
        raise ValueError(
            f"{path}: read {name} (pair {record_no + 1}) is {len(seq)} bp, shorter "
            f"than the {pattern_len} bp barcode pattern."
        )


def trim(readPair, pattern):
    """Move each mate's barcode bases (pattern N positions) into the read
    name and DB tag and clip the pattern off the read. Records are
    [header, seq, qual] with or without trailing newlines; mate headers must
    name the same fragment (a trailing /1 or /2 is dropped) and both mates
    get the same rewritten name."""
    adapterLen = len(pattern)
    (name1, seq1, qual1), (name2, seq2, qual2) = readPair
    seq1, seq2 = seq1.rstrip("\r\n"), seq2.rstrip("\r\n")
    qual1, qual2 = qual1.rstrip("\r\n"), qual2.rstrip("\r\n")
    base_name = _read_base_name(name1, "read 1", 0)
    if _read_base_name(name2, "read 2", 0) != base_name:
        raise ValueError(f"Mate names differ: {name1.strip()!r} vs {name2.strip()!r}")
    bc1 = "".join(
        [base for nn, base in enumerate(seq1[0:adapterLen]) if pattern[nn] == "N"]
    )
    bc2 = "".join(
        [base for nn, base in enumerate(seq2[0:adapterLen]) if pattern[nn] == "N"]
    )
    if len(bc1) > 0:
        namenew = f"@{base_name}_{bc1}+{bc2} DB:Z:{bc1}-{bc2}\n"
    else:
        namenew = f"@{base_name}\n"
    seqnew1 = seq1[adapterLen:] + "\n"
    seqnew2 = seq2[adapterLen:] + "\n"
    qualnew1 = qual1[adapterLen:] + "\n"
    qualnew2 = qual2[adapterLen:] + "\n"
    return [namenew, seqnew1, qualnew1], [namenew, seqnew2, qualnew2]


HASH_MAX_MISMATCH = 1


def build_barcode_index(barcodes, max_mismatch):
    """Lookup tables for max_mismatch <= HASH_MAX_MISMATCH: for each
    barcode length, every listed barcode and (if max_mismatch == 1) every
    1-mismatch variant of it, mapped to the (mismatches, barcode) entries it
    is within max_mismatch of. Size is about len(barcodes) * length * 3,
    linear in the list (variants within m mismatches grow exponentially in
    m, hence the walk above HASH_MAX_MISMATCH)."""
    assert max_mismatch <= HASH_MAX_MISMATCH
    by_len = {}
    for bc in barcodes:
        table = by_len.setdefault(len(bc), {})
        table.setdefault(bc, []).append((0, bc))
        if max_mismatch == 1:
            for p in range(len(bc)):
                for c in "ACGTN":
                    if c != bc[p]:
                        table.setdefault(bc[:p] + c + bc[p + 1 :], []).append((1, bc))
    return sorted(by_len.items())


def match_barcode_indexed(seq, start, index):
    """Same result as match_barcode (fewest mismatches, then shortest;
    same-length tie -> AMBIGUOUS) using build_barcode_index's tables: one
    lookup per distinct barcode length, shortest first. An exact hit ends
    the search -- nothing can have fewer mismatches, and any later hit is
    longer -- while a 1-mismatch hit keeps looking for an exact longer one."""
    n = len(seq)
    best = None  # (mismatches, length, [barcodes])
    for length, table in index:
        if start + length > n:
            break
        hits = table.get(seq[start : start + length])
        if not hits:
            continue
        mm = min(h[0] for h in hits)
        if best is None or mm < best[0]:
            best = (mm, length, [h[1] for h in hits if h[0] == mm])
        if mm == 0:
            break
    if best is None:
        return NO_MATCH, 0
    if len(best[2]) > 1:
        return AMBIGUOUS, 0
    return best[2][0], best[1]


def make_barcode_matcher(barcodes, max_mismatch):
    """match(seq, start) -> (barcode, length) | (NO_MATCH/AMBIGUOUS, 0):
    lookup tables for max_mismatch <= 1, the position walk above that."""
    if max_mismatch <= HASH_MAX_MISMATCH:
        index = build_barcode_index(barcodes, max_mismatch)
        return lambda seq, start: match_barcode_indexed(seq, start, index)
    return lambda seq, start: match_barcode(seq, start, barcodes, max_mismatch)


def _extract_listed(seq, prefix, suffix, matcher):
    """Barcode and clip length for one read under a pattern with B:
    prefix and suffix are the pattern's fixed parts before and after B.
    The barcode is the prefix's N bases + the matched listed barcode + the
    suffix's N bases. Returns (barcode, clip) or (NO_MATCH/AMBIGUOUS, 0)."""
    p = len(prefix)
    listed, b_len = matcher(seq, p)
    if not listed:
        return listed, 0
    clip = p + b_len + len(suffix)
    if len(seq) < clip:
        return NO_MATCH, 0
    head = "".join(base for nn, base in enumerate(seq[:p]) if prefix[nn] == "N")
    tail = "".join(
        base for nn, base in enumerate(seq[p + b_len : clip]) if suffix[nn] == "N"
    )
    return head + listed + tail, clip


def trim_listed(readPair, pattern, matcher):
    """trim() for a pattern with B: each mate's B part is matched against
    the barcode list by matcher (make_barcode_matcher). Returns the two output records,
    or (None, reason) with reason one of "no_match_r1", "no_match_r2",
    "ambiguous" if the pair is dropped."""
    prefix, suffix = pattern.split("B")
    (name1, seq1, qual1), (name2, seq2, qual2) = readPair
    seq1, seq2 = seq1.rstrip("\r\n"), seq2.rstrip("\r\n")
    qual1, qual2 = qual1.rstrip("\r\n"), qual2.rstrip("\r\n")
    base_name = _read_base_name(name1, "read 1", 0)
    if _read_base_name(name2, "read 2", 0) != base_name:
        raise ValueError(f"Mate names differ: {name1.strip()!r} vs {name2.strip()!r}")
    bc1, clip1 = _extract_listed(seq1, prefix, suffix, matcher)
    bc2, clip2 = _extract_listed(seq2, prefix, suffix, matcher)
    if bc1 == AMBIGUOUS or bc2 == AMBIGUOUS:
        return None, "ambiguous"
    if bc1 is NO_MATCH:
        return None, "no_match_r1"
    if bc2 is NO_MATCH:
        return None, "no_match_r2"
    namenew = f"@{base_name}_{bc1}+{bc2} DB:Z:{bc1}-{bc2}\n"
    return (
        [namenew, seq1[clip1:] + "\n", qual1[clip1:] + "\n"],
        [namenew, seq2[clip2:] + "\n", qual2[clip2:] + "\n"],
    )


def do_trim(args):
    # Check if input files exist
    if not os.path.exists(args.fq):
        raise FileNotFoundError(f"Input fastq file not found: {args.fq}")
    if not os.path.exists(args.fq2):
        raise FileNotFoundError(f"Input fastq file not found: {args.fq2}")
    barcode_list = getattr(args, "barcode_list", None)
    barcodes = load_barcode_list(barcode_list) if barcode_list else None
    max_mismatch = getattr(args, "max_mismatch", 1)
    if max_mismatch is None:
        max_mismatch = 1
    if max_mismatch < 0:
        raise ValueError(f"--max-mismatch must be >= 0, got {max_mismatch}")
    _validate_pattern(args.pattern, barcodes)
    if barcodes is None:
        pattern_len = len(args.pattern)
    else:
        # shortest possible pattern: fixed parts + the shortest barcode
        pattern_len = len(args.pattern) - 1 + min(len(b) for b in barcodes)
    dropped = {"no_match_r1": 0, "no_match_r2": 0, "ambiguous": 0}
    matcher = make_barcode_matcher(barcodes, max_mismatch) if barcodes else None

    fq1 = _open_fastq(args.fq)
    fq2 = _open_fastq(args.fq2)
    with open(args.output + "_1.fastq", "w") as out1:
        out1.write("")
    with open(args.output + "_2.fastq", "w") as out2:
        out2.write("")
    fq1Out = open(args.output + "_1.fastq", "a")
    fq2Out = open(args.output + "_2.fastq", "a")
    lineIndex = 0
    record_no = 0
    for line1, line2 in zip_longest(fq1, fq2, fillvalue=_EOF):
        if line1 is _EOF or line2 is _EOF:
            if lineIndex != 0:
                raise ValueError(
                    f"{args.fq} and {args.fq2} both ended mid-FASTQ-record "
                    f"(read pair {record_no}, incomplete record); at least "
                    "one file is truncated or corrupt."
                )
            if line1 is not line2:
                shorter, longer = (
                    (args.fq, args.fq2) if line1 is _EOF else (args.fq2, args.fq)
                )
                raise ValueError(
                    f"{shorter} has fewer reads than {longer} "
                    f"({record_no} complete read pairs processed before "
                    "mismatch)."
                )
            break
        if lineIndex == 0:
            name1 = line1
            name2 = line2
            lineIndex = 1
        elif lineIndex == 1:
            seq1 = line1.rstrip("\r\n")
            seq2 = line2.rstrip("\r\n")
            lineIndex = 2
        elif lineIndex == 2:
            lineIndex = 3
        else:
            qual1 = line1.rstrip("\r\n")
            qual2 = line2.rstrip("\r\n")
            base1 = _read_base_name(name1, args.fq, record_no)
            base2 = _read_base_name(name2, args.fq2, record_no)
            if base1 != base2:
                raise ValueError(
                    f"Mate names differ at read pair {record_no + 1}: "
                    f"{base1!r} ({args.fq}) vs {base2!r} ({args.fq2}); the "
                    "FASTQs are not in the same order."
                )
            _check_read(base1, seq1, qual1, pattern_len, args.fq, record_no)
            _check_read(base2, seq2, qual2, pattern_len, args.fq2, record_no)
            readPair = [[name1, seq1, qual1], [name2, seq2, qual2]]
            lineIndex = 0
            record_no += 1
            if barcodes is None:
                read1, read2 = trim(readPair, args.pattern)
            else:
                read1, read2 = trim_listed(readPair, args.pattern, matcher)
                if read1 is None:
                    dropped[read2] += 1
                    continue
            fq1Out.write(read1[0] + read1[1] + "+\n" + read1[2])
            fq2Out.write(read2[0] + read2[1] + "+\n" + read2[2])
    fq1Out.close()
    fq2Out.close()
    fq1.close()
    fq2.close()
    if barcodes is not None:
        kept = record_no - sum(dropped.values())
        print(
            f"trim: {record_no} read pairs, {kept} written; dropped "
            f"{dropped['no_match_r1']} (read 1 barcode not in list), "
            f"{dropped['no_match_r2']} (read 2 barcode not in list), "
            f"{dropped['ambiguous']} (ambiguous barcode); "
            f"max mismatches {max_mismatch}"
        )
