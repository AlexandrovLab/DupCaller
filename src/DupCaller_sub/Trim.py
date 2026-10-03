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


def _validate_pattern(pattern):
    if not pattern or set(pattern) - {"N", "X"}:
        raise ValueError(
            f"Invalid barcode pattern '{pattern}': it must be non-empty and use "
            "only N (barcode base) and X (skipped base)."
        )


def _check_pattern_pair(pattern, pattern2):
    """Read 1 and read 2 patterns must extract barcodes of the same length:
    each end of a molecule is read by read 1 on one strand and by read 2 on
    the other, and call pairs the strands by swapping (bc1, bc2), so a
    different N count would split every duplex into two single-strand
    families."""
    n1, n2 = pattern.count("N"), pattern2.count("N")
    if n1 != n2:
        raise ValueError(
            f"Barcode patterns '{pattern}' (read 1) and '{pattern2}' (read 2) have "
            f"{n1} and {n2} N (barcode) bases. They must match: each molecule end's "
            "barcode is read by read 1 on one strand and by read 2 on the other, so "
            "different lengths would stop the two strands of a molecule from being "
            "grouped into one duplex family. Only the X (skipped) bases may differ."
        )


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


def trim(readPair, pattern, pattern2=None):
    """Move each mate's barcode bases (pattern N positions) into the read
    name and DB tag and clip the pattern off the read. Read 2 uses pattern2
    when given, else pattern. Records are [header, seq, qual] with or
    without trailing newlines; mate headers must name the same fragment (a
    trailing /1 or /2 is dropped) and both mates get the same rewritten
    name."""
    if pattern2 is None:
        pattern2 = pattern
    adapterLen = len(pattern)
    adapterLen2 = len(pattern2)
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
        [base for nn, base in enumerate(seq2[0:adapterLen2]) if pattern2[nn] == "N"]
    )
    if len(bc1) > 0:
        namenew = f"@{base_name}_{bc1}+{bc2} DB:Z:{bc1}-{bc2}\n"
    else:
        namenew = f"@{base_name}\n"
    seqnew1 = seq1[adapterLen:] + "\n"
    seqnew2 = seq2[adapterLen2:] + "\n"
    qualnew1 = qual1[adapterLen:] + "\n"
    qualnew2 = qual2[adapterLen2:] + "\n"
    return [namenew, seqnew1, qualnew1], [namenew, seqnew2, qualnew2]


def do_trim(args):
    # Check if input files exist
    if not os.path.exists(args.fq):
        raise FileNotFoundError(f"Input fastq file not found: {args.fq}")
    if not os.path.exists(args.fq2):
        raise FileNotFoundError(f"Input fastq file not found: {args.fq2}")
    _validate_pattern(args.pattern)
    pattern2 = getattr(args, "pattern2", None) or args.pattern
    _validate_pattern(pattern2)
    _check_pattern_pair(args.pattern, pattern2)
    pattern_len = len(args.pattern)
    pattern2_len = len(pattern2)

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
            _check_read(base2, seq2, qual2, pattern2_len, args.fq2, record_no)
            readPair = [[name1, seq1, qual1], [name2, seq2, qual2]]
            read1, read2 = trim(readPair, args.pattern, pattern2)
            lineIndex = 0
            record_no += 1
            fq1Out.write(read1[0] + read1[1] + "+\n" + read1[2])
            fq2Out.write(read2[0] + read2[1] + "+\n" + read2[2])
    fq1Out.close()
    fq2Out.close()
    fq1.close()
    fq2.close()
