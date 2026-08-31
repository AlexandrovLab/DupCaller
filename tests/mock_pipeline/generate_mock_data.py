#!/usr/bin/env python3
"""Generator for the synthetic dataset in tests/mock_pipeline/data/.

Not run by CI -- its output (reference.fa, repeats.tsv, mock_1.fastq,
mock_2.fastq) is committed as static data. Re-run only when deliberately
changing the planted mutations/layout below, and regenerate
tests/mock_pipeline/expected/ to match.
"""
import random
import textwrap

random.seed(20260830)

CONTIG = "mockchr1"
RAW_READ_LEN = 100  # length of each raw (pre-trim) read
BARCODE_PATTERN = "NNNXXXX"  # 3bp UMI + 4bp constant spacer (NanoSeq-style)
BC_REGION = len(BARCODE_PATTERN)  # 7
INSERT_READ_LEN = RAW_READ_LEN - BC_REGION  # 93, post-trim read length
FRAGMENT_LEN = 200  # genomic span of one duplex molecule
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def random_barcode():
    return "".join(random.choice("ACGT") for _ in range(3))


def random_skip():
    return "".join(random.choice("ACGT") for _ in range(4))


# Reference sequence, with a homopolymer run and an STR tract planted at
# fixed offsets.
CONTIG_LEN = 3200
seq = list(random_seq(CONTIG_LEN))

HP_START = 2050  # 8-base A homopolymer, for the planted deletion
HP_LEN = 8
seq[HP_START : HP_START + HP_LEN] = list("A" * HP_LEN)

STR_START = 2900  # (AC) x 8 dinucleotide STR, index-only (no reads)
STR_UNIT = "AC"
STR_UNITS = 8
seq[STR_START : STR_START + len(STR_UNIT) * STR_UNITS] = list(STR_UNIT * STR_UNITS)

reference_seq = "".join(seq)

# Somatic SBS at MUT_POS, confirmed by two independent duplex molecules,
# and a somatic 1bp deletion inside the HP run, confirmed by one.
MUT_POS = 1000  # 0-based; reference base must not already be ALT
REF_BASE = reference_seq[MUT_POS]
ALT_BASE = {"A": "C", "C": "A", "G": "T", "T": "G"}[REF_BASE]

DEL_POS = HP_START + 3  # 0-based first deleted reference base


def make_molecule(start, bc_a, bc_b, mutate=None, n_copies=3):
    """A duplex molecule spanning [start, start+FRAGMENT_LEN): two strand
    families of PCR-duplicate read pairs, each replicated n_copies times,
    tagged with the same barcode pair in swapped order (DupCaller's trim/
    DB-tag convention -- see README Steps 2-4). `mutate` is an optional
    (rel_pos, ref, alt) edit applied to the genomic sequence before both
    families are drawn from it, so every copy of both strands carries it.
    `alt` may be "" for a deletion or longer than `ref` for an insertion.
    """
    frag = list(reference_seq[start : start + FRAGMENT_LEN])
    if mutate is not None:
        rel_pos, ref, alt = mutate
        assert "".join(frag[rel_pos : rel_pos + len(ref)]) == ref
        frag[rel_pos : rel_pos + len(ref)] = list(alt)
    frag = "".join(frag)

    r1_insert = frag[0:INSERT_READ_LEN]
    r2_insert = revcomp(frag[len(frag) - INSERT_READ_LEN :])

    records = []
    for group, (bc1, bc2, ins1, ins2) in {
        "top": (bc_a, bc_b, r1_insert, r2_insert),
        "bottom": (bc_b, bc_a, r2_insert, r1_insert),
    }.items():
        for copy in range(n_copies):
            name = f"{CONTIG}_{start}_{group}_{copy}"
            seq1 = bc1 + random_skip() + ins1
            seq2 = bc2 + random_skip() + ins2
            qual1 = "F" * len(seq1)
            qual2 = "F" * len(seq2)
            records.append((name, seq1, qual1, seq2, qual2))
    return records


molecules = []
# Background (reference-only) coverage, tiled across the contig.
for start in (100, 500, 1600, 2400, 2700):
    molecules += make_molecule(start, random_barcode(), random_barcode())

# SBS somatic call: two independent duplex molecules confirming the same
# alt allele at MUT_POS.
molecules += make_molecule(
    MUT_POS - 30,
    random_barcode(),
    random_barcode(),
    mutate=(30, REF_BASE, ALT_BASE),
)
molecules += make_molecule(
    MUT_POS - 60,
    random_barcode(),
    random_barcode(),
    mutate=(60, REF_BASE, ALT_BASE),
)

# Indel somatic call: one duplex molecule with a 1bp deletion in the HP run.
molecules += make_molecule(
    DEL_POS - 25,
    random_barcode(),
    random_barcode(),
    mutate=(25, reference_seq[DEL_POS], ""),
)

random.shuffle(molecules)

with open("reference.fa", "w") as fh:
    fh.write(f">{CONTIG}\n")
    fh.write("\n".join(textwrap.wrap(reference_seq, 70)))
    fh.write("\n")

with open("repeats.tsv", "w") as fh:
    fh.write(
        "\t".join(
            [
                CONTIG,
                str(STR_START),
                str(STR_START + len(STR_UNIT) * STR_UNITS),
                STR_UNIT,
                str(len(STR_UNIT)),
                "+",
                str(STR_UNITS),
                f"({STR_UNIT}){STR_UNITS}",
            ]
        )
        + "\n"
    )

with open("mock_1.fastq", "w") as f1, open("mock_2.fastq", "w") as f2:
    for name, seq1, qual1, seq2, qual2 in molecules:
        f1.write(f"@{name}\n{seq1}\n+\n{qual1}\n")
        f2.write(f"@{name}\n{seq2}\n+\n{qual2}\n")

with open("truth.txt", "w") as fh:
    fh.write(f"contig={CONTIG} length={CONTIG_LEN}\n")
    fh.write(f"SBS: pos(0-based)={MUT_POS} ref={REF_BASE} alt={ALT_BASE} molecules=2\n")
    fh.write(
        f"INDEL: pos(0-based, first deleted base)={DEL_POS} "
        f"ref={reference_seq[DEL_POS - 1:DEL_POS + 1]} "
        f"alt={reference_seq[DEL_POS - 1]} molecules=1 (in {HP_LEN}xA homopolymer "
        f"starting at {HP_START})\n"
    )
    fh.write(
        f"STR: {STR_UNIT} x{STR_UNITS} at [{STR_START},"
        f"{STR_START + len(STR_UNIT) * STR_UNITS})\n"
    )
    fh.write(f"read pairs: {len(molecules)}\n")

print("done")
