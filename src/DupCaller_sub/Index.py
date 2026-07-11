#!/usr/bin/env python3
import h5py
import gzip
import numpy as np
import os
from Bio import SeqIO
from pysam import TabixFile as TABIX


def do_index(args):
    ### Check if reference file exists
    if not os.path.exists(args.reference):
        raise FileNotFoundError(f"Reference file not found: {args.reference}")

    ### Load Fasta
    if args.reference.endswith(".gz"):
        with gzip.open(args.reference, "rt") as handle:
            fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        fasta = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

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
    base2num_npfunc = np.vectorize(lambda b: base2num.get(b, 4))

    ### Generate reference and index matrix into h5py file
    ref_h5 = h5py.File(args.reference + ".ref.h5", "w")
    tn_h5 = h5py.File(args.reference + ".tn.h5", "w")
    hp_h5 = h5py.File(args.reference + ".hp.h5", "w")
    repeat_bed = TABIX(args.repeatBed, "r")
    str_bed = TABIX(args.strbed, "r")

    def apply_repeat_bed(bed, chrom, unit_len, is_start, count_down):
        if chrom not in bed.contigs:
            return
        for line in bed.fetch(chrom):
            fields = line.split("\t")
            start = int(fields[1])
            end = int(fields[2])
            unit_length = int(fields[3])
            unit_len[start:end] = unit_length
            is_start[start] = 1
            bases_to_end = end - np.arange(start, end)
            count_down[start:end] = -(-bases_to_end // unit_length)

    for chrom in fasta.keys():
        print(f"Currently processing:{chrom}")
        print(f"Creating reference sequence index")
        reference_seq = list(fasta[chrom].seq.upper())
        reference_int = base2num_npfunc(np.array(reference_seq)).astype(np.uint8)
        print(f"Creating trinucleotide index")
        reference_minus = [
            "N",
        ] + reference_seq[:-1]
        reference_plus = reference_seq[1:] + [
            "N",
        ]
        trinucs = [
            a + b + c for a, b, c in zip(reference_minus, reference_seq, reference_plus)
        ]
        trinuc_int = np.array([trinuc2num.get(_, 96) for _ in trinucs], dtype="uint8")
        print(f"Creating repeat index")
        ref_len = len(reference_seq)
        unit_len = np.zeros(ref_len, dtype="uint8")
        is_start = np.zeros(ref_len, dtype="uint8")
        count_down = np.zeros(ref_len, dtype=int)
        # Repeats <=12bp (including homopolymers) come from --repeatBed,
        # repeats >12bp come from --strbed
        apply_repeat_bed(repeat_bed, chrom, unit_len, is_start, count_down)
        apply_repeat_bed(str_bed, chrom, unit_len, is_start, count_down)

        hp_lens = np.vstack((unit_len, is_start, count_down))
        hp_lens[hp_lens > 127] = 127
        hp_lens = hp_lens.astype(np.uint8)
        idx = ref_h5.create_dataset(chrom, data=reference_int)
        tn = tn_h5.create_dataset(chrom, data=trinuc_int)
        hp = hp_h5.create_dataset(chrom, data=hp_lens)
    ref_h5.close()
    tn_h5.close()
    hp_h5.close()
