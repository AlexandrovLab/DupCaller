"""Derive a read-2 FASTQ whose reads are N bases shorter at the 3' end.

After trim -p NNNXXXX, read 2 is then N bases shorter than read 1, so each
duplex family's downstream mates (read 1 from one strand, read 2 from the
other) start N bases apart and must still be merged into one family by
call's rugged-mate handling. The genomic span covered by the molecule is
unchanged at its 5' ends, so the pipeline must reproduce expected/.

Usage: python make_short_read2.py IN_R2.fastq OUT_R2.fastq [N=2]
"""

import sys


def main():
    src, dst = sys.argv[1:3]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    lines = open(src).read().splitlines()
    out = []
    for i in range(0, len(lines), 4):
        name, seq, plus, qual = lines[i : i + 4]
        out += [name, seq[:-n], plus, qual[:-n]]
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
