"""Derive META-CS-style reads from the synthetic mock reads to test trim's B
pattern letter (barcode from a list, any length).

Each mock read starts with a 3-base barcode and 4 spacer bases (NNNXXXX).
Here that 7-base start is replaced by a variable-length code (10-13 nt) for
the 3-base barcode, followed by the 19-nt Tn5 ME, as in META-CS. Every 4th
read gets one mismatch in its code. Codes are chosen so that, with the ME
after them, every code is at least 3 mismatches from every other code's
read prefix, so one sequencing error can never make another code closer.

trim -p B + 19 X -bl OUT_LIST then recovers the codes, and the pipeline must
give the same duplex families as the plain run; OUT_MAP (code<TAB>barcode)
maps the codes in the outputs back to the original barcodes.

Usage: python make_barcode_list_reads.py IN_R1 IN_R2 OUT_R1 OUT_R2 OUT_LIST OUT_MAP
"""

import itertools
import random
import sys

ME = "AGATGTGTATAAGAGACAG"
SKIP = 7  # len("NNNXXXX")
MIN_DIST = 3


def _dist_as_read(code, other):
    """Mismatches between `other` and the first len(other) bases of a read
    carrying `code` followed by the ME."""
    read = code + ME
    return sum(a != b for a, b in zip(other, read[: len(other)]))


def make_codes():
    rng = random.Random(20261003)
    codes = {}
    for i, bc in enumerate("".join(p) for p in itertools.product("ACGT", repeat=3)):
        length = 10 + i % 4
        while True:
            cand = "".join(rng.choice("ACGT") for _ in range(length))
            if all(
                _dist_as_read(cand, c) >= MIN_DIST
                and _dist_as_read(c, cand) >= MIN_DIST
                for c in codes.values()
            ):
                codes[bc] = cand
                break
    return codes


def rewrite(src, dst, codes, counter):
    lines = open(src).read().splitlines()
    out = []
    for i in range(0, len(lines), 4):
        name, seq, plus, qual = lines[i : i + 4]
        code = codes[seq[:3]]
        if counter[0] % 4 == 3:
            p = counter[0] % len(code)
            code = code[:p] + ("A" if code[p] != "A" else "C") + code[p + 1 :]
        counter[0] += 1
        prefix = code + ME
        out += [name, prefix + seq[SKIP:], plus, qual[0] * len(prefix) + qual[SKIP:]]
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    in1, in2, out1, out2, out_list, out_map = sys.argv[1:7]
    codes = make_codes()
    counter = [0]
    rewrite(in1, out1, codes, counter)
    rewrite(in2, out2, codes, counter)
    with open(out_list, "w") as f:
        f.write("\n".join(codes.values()) + "\n")
    with open(out_map, "w") as f:
        f.write("".join(f"{c}\t{bc}\n" for bc, c in codes.items()))


if __name__ == "__main__":
    main()
