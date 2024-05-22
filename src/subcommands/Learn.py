#!/usr/bin/env python3
import os
import time
from multiprocessing import Pool
import errno

import pysam
import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt


from .funcs.call import callBam
from .funcs.misc import createVcfStrings
from .funcs.misc import splitBamRegions
from .funcs.misc import getAlignmentObject as BAM


# if __name__ == "__main__":
def do_learn(args):
    params = {
        "tumorBam": args.bam,
        "normalBam": None,
        "germline": args.germline,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regions,
        "threads": args.threads,
        "mutRate": 10e-7,
        "pcutoff": 2,
        "amperr": 1e-5,
        "amperr_file": None,
        "amperri": 1e-6,
        "amperri_file": None,
        "dmgerr": 1e-5,
        "dmgerri": 1e-6,
        "dmgerr_file": None,
        "dmgerri_file": None,
        "mapq": args.mapq,
        "noise": args.noise,
        "indel_bed": args.indelbed,
        "trim5": args.trimF,
        "trim3": args.trimR,
        "minNdepth": args.minNdepth,
        "germline_cutoff": args.germlineAfCutoff,
        "minBq": args.minBq,
        "minAltQual": args.minAltQual,
        "maxNM": args.nmflt,
        "step": args.windowSize,
        "minRef": args.minRef,
        "minAlt": args.minAlt,
    }
    """
    Initialze run
    """
    print("..............Loading reference genome.....................")
    startTime = time.time()
    if not os.path.exists("tmp"):
        try:
            os.mkdir("tmp")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    if not os.path.exists(params["output"]):
        try:
            os.mkdir(params["output"])
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    bamObject = BAM(args.bam, "rb")

    """
    Execulte variant calling
    """
    if args.threads <= 2:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        # contigs = [(r.strip('\n'),) for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        paramsNow = params
        # paramsNow["reference"] = fasta
        paramsNow["isLearn"] = True
        regions = params["regions"]
        paramsNow["regions"] = [
            (chrom, 0, bamObject.get_reference_length(chrom) - 1) for chrom in regions
        ]
        (
            mismatch_profile,
            indelerr_profile,
            mismatch_dmg_profile,
            indelerr_dmg_profile,
        ) = callBam(paramsNow, 0)
    else:
        """
        Multi-thread execution
        """
        args.threads = args.threads - 1
        # contigs = [r.strip('\n') for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        print(
            "...........Spliting genomic regions for parallel execution................"
        )
        # print(args.threads)
        # if args.normalBam:
        cutSites, chunkSize, contigs = splitBamRegions(
            [args.bam], args.threads, contigs, args.windowSize
        )
        # else:
        # cutSites, chunkSize, contigs = splitBamRegions(
        # [args.bam], args.threads, contigs, args.windowSize
        # )
        # print(cutSites,chunkSize,contigs)# Split the whole genome for parallel execution
        regionSequence = []
        currentContigIndex = 0
        usedTime = (time.time() - startTime) / 60
        print(f"....Genomic regions splitted in {usedTime} minutes...")
        """
        Determine regions for each process
        """

        for nn, site in enumerate(cutSites[1:]):
            pSite = cutSites[nn]
            if site[0] == pSite[0]:
                regionSequence.append((contigs[site[0]], pSite[1], site[1]))
            else:
                if pSite[1] != 0:
                    regionSequence.append((contigs[pSite[0]], pSite[1]))
                else:
                    regionSequence.append((contigs[pSite[0]],))
                for ii in range(pSite[0] + 1, site[0]):
                    regionSequence.append((contigs[ii],))
                regionSequence.append((contigs[site[0]], 0, site[1]))
        regionSequence.append((contigs[site[0]], site[1]))
        for ii in range(site[0] + 1, len(contigs)):
            regionSequence.append((contigs[ii],))
        print(
            "............Completed region splitting in "
            + str((time.time() - startTime) / 60)
            + " minutes,loading reference genome.................."
        )

        """
        Start variant calling 
        """

        callArguments = []
        startTime2 = time.time()
        print(".........Starting variant calling..............")
        pool = Pool()
        for nn in range(args.threads):
            regions = []
            while len(regionSequence) != 0:
                if len(regionSequence[0]) != 3:
                    regions.append(regionSequence.pop(0))
                else:
                    regions.append(regionSequence.pop(0))
                    break
            chroms = [region[0] for region in regions]
            paramsNow = params.copy()
            paramsNow["regions"] = regions
            paramsNow["isLearn"] = True
            callArgument = (paramsNow, nn)
            callArguments.append(callArgument)
            regions = []
        results = pool.starmap(
            callBam, callArguments
        )  # each result return three list: number of duplex reads, effective lengths, list of mutations
        print(
            "..............Completed bam calling in "
            + str((time.time() - startTime2) / 60)
            + " minutes,merging results................."
        )
        pool.close()
        pool.terminate()
        pool.join()

        mismatch_profile = sum([_[0] for _ in results]).astype(int)
        indelerr_profile = sum([_[1] for _ in results]).astype(int)
        mismatch_dmg_profile = sum([_[2] for _ in results]).astype(int)
        # mismatch_F2R1_dmg_profile = sum([_[3] for _ in results]).astype(int)
        indelerr_dmg_profile = sum([_[3] for _ in results]).astype(int)

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
    amp_tn_pd = pd.DataFrame(
        mismatch_profile, columns=["A", "T", "C", "G"], index=num2trinuc
    )
    dmg_tn_pd = pd.DataFrame(
        mismatch_dmg_profile, columns=["A", "T", "C", "G"], index=num2trinuc
    )
    # np.savetxt(params["output"] + "/" + args.output + ".amp.tn.txt",np.hstack([trinuc_cols[0:32],mismatch_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")
    amp_tn_pd.to_csv(params["output"] + "/" + args.output + ".amp.tn.txt", sep="\t")
    np.savetxt(
        params["output"] + "/" + args.output + ".amp.id.txt",
        indelerr_profile,
        delimiter="\t",
        fmt="%d",
    )
    dmg_tn_pd.to_csv(params["output"] + "/" + args.output + ".dmg.tn.txt", sep="\t")
    # np.savetxt(params["output"] + "/" + args.output + ".dmg.tn.txt",np.hstack([trinuc_cols,mismatch_dmg_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")
    np.savetxt(
        params["output"] + "/" + args.output + ".dmg.id.txt",
        indelerr_dmg_profile,
        delimiter="\t",
        fmt="%d",
    )
    print(
        "..............Completed error learning "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )
