#!/usr/bin/env python3
import h5py

import numpy as np
from Bio import SeqIO


def do_index(args):
    ### Load Fasta
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
            trinuc = 'N' + ref_base + plus_base
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            trinuc = minus_base + ref_base + 'N'
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    for ref_base in ["G", "A"]:
        for plus_base in ["T", "A", "G", "C"]:
            trinuc = 'N' + ref_base + plus_base
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    for minus_base in ["T", "A", "G", "C"]:
        for ref_base in ["G", "A"]:
            trinuc = minus_base + ref_base + 'N'
            trinuc2num[trinuc] = trinuc_order
            num2trinuc.append(trinuc)
            trinuc_order += 1
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3,'a':0,'t':1,'c':2,'g':3}
    base2num_npfunc = np.vectorize(lambda b:base2num.get(b,4))

    ### Generate reference and index matrix into h5py file
    ref_h5 = h5py.File(args.reference+".ref.h5","w")
    tn_h5 = h5py.File(args.reference+".tn.h5","w")
    hp_h5 = h5py.File(args.reference+".hp.h5","w")
    for chrom in fasta.keys():
        print(f"Currently processing:{chrom}")
        print(f"Creating reference sequence index")
        reference_seq =  list(fasta[chrom].seq.upper())
        reference_int = base2num_npfunc(np.array(reference_seq)).astype(np.uint8)
        print(f"Creating trinucleotide index")
        reference_minus = ["N",] + reference_seq[:-1]
        reference_plus = reference_seq[1:] + ["N",]
        trinucs = [a+b+c for a,b,c in zip(reference_minus,reference_seq,reference_plus)]
        trinuc_int = np.array([trinuc2num.get(_,96) for _ in trinucs],dtype="uint8")
        print(f"Creating homopolymer index")
        hp_lens_forw = np.ones(len(reference_seq),dtype=int)
        hp_lens_rev = np.ones(len(reference_seq),dtype=int)
        for nn,b in enumerate(reference_int[1:]):
            if reference_int[nn] == reference_int[nn+1]:
                hp_lens_forw[nn+1] = hp_lens_forw[nn] + 1
        ref_len = reference_int.size
        for nn,b in enumerate(reference_int[:-1]):
            if reference_int[ref_len - nn - 1] == reference_int[ref_len - nn - 2]:
                hp_lens_rev[ref_len - nn - 2] = hp_lens_rev[ref_len - nn - 1] + 1
        hp_lens = np.vstack((hp_lens_forw,hp_lens_rev))
        hp_lens[hp_lens>127] = 127
        hp_lens = hp_lens.astype(np.uint8)
        idx = ref_h5.create_dataset(chrom,data=reference_int)
        tn = tn_h5.create_dataset(chrom,data=trinuc_int)
        hp = hp_h5.create_dataset(chrom,data=hp_lens)
    ref_h5.close()
    tn_h5.close()
    hp_h5.close()



