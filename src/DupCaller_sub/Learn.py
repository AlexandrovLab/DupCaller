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
from .funcs.indels import STR_COL_LABELS
from .funcs.learn import NUM_BQ, estimate_sbs_srd_rates
from .funcs.misc import createVcfStrings
from .funcs.misc import splitBamRegions
from .funcs.misc import drop_empty_regions
from .funcs.misc import getAlignmentObject as BAM
from .funcs.misc import check_mate_cigar_tags


def do_learn(args):
    params = {
        "tumorBam": args.bam,
        "normalBams": args.normalBams,
        "germline": args.germline,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regions,
        "region_file": None,
        "threads": args.threads,
        # callBam reads seed and the barcode keys even when isLearn is set.
        "seed": int(np.random.SeedSequence().generate_state(1)[0]),
        "barcodeTag": "DB",
        "barcodeNormalize": True,
        "barcodeSep": "-",
        "rescue": False,
        "pcutoff": 1,
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
        "minMeanASXS": args.minMeanASXS,
        "step": args.windowSize,
        "minRef": args.minRef,
        "minAlt": args.minAlt,
        "pseudocount": args.pseudocount,
    }
    """
    Initialze run
    """
    print("..............Loading reference genome.....................")
    startTime = time.time()
    if not os.path.exists(params["output"]):
        try:
            os.mkdir(params["output"])
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    tmp_dir = os.path.join(params["output"], "tmp")
    params["tmp_dir"] = tmp_dir
    os.makedirs(tmp_dir, exist_ok=True)
    # Same ERROR/ layout as Caller.py's auto-learn.
    error_dir = os.path.join(params["output"], "ERROR")
    os.makedirs(error_dir, exist_ok=True)
    error_prefix = os.path.join(error_dir, args.output)
    bamObject = BAM(args.bam, "rb", args.reference)

    """
    Execulte variant calling
    """
    # Only the multi-thread path needs MC tags to hand reads across chunk
    # boundaries.
    if args.threads > 2:
        check_mate_cigar_tags(args.bam, args.reference)
    if args.threads <= 2:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        paramsNow = params
        paramsNow["isLearn"] = True
        regions = params["regions"]
        paramsNow["regions"] = [
            (chrom, 0, bamObject.get_reference_length(chrom)) for chrom in regions
        ]
        (
            mismatch_profile,
            hp_alt_profile,
            str_alt_profile,
            mismatch_dmg_profile,
            hp_dmg_profile,
            str_dmg_profile,
            sbs_alt_bq_hist,
        ) = callBam(paramsNow, 0)
        mismatch_profile = mismatch_profile.astype(int)
        hp_alt_profile = hp_alt_profile.astype(int)
        str_alt_profile = str_alt_profile.astype(int)
        mismatch_dmg_profile = mismatch_dmg_profile.astype(int)
        hp_dmg_profile = hp_dmg_profile.astype(int)
        str_dmg_profile = str_dmg_profile.astype(int)
        sbs_alt_bq_hist = sbs_alt_bq_hist.astype(int)
    else:
        """
        Multi-thread execution
        """
        args.threads = args.threads - 1
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        print(
            "...........Spliting genomic regions for parallel execution................"
        )
        cutSites, chunkSize, contigs = splitBamRegions(
            [args.bam],
            args.threads,
            contigs,
            args.windowSize,
            min_chunk_length=getattr(args, "minChunkLength", 10000),
        )
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
        site = cutSites[-1]
        regionSequence.append((contigs[site[0]], site[1]))
        for ii in range(site[0] + 1, len(contigs)):
            regionSequence.append((contigs[ii],))
        regionSequence = drop_empty_regions(regionSequence, bamObject)
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
        pool = Pool(args.threads)
        for nn in range(args.threads):
            regions = []
            while len(regionSequence) != 0:
                if len(regionSequence[0]) != 3:
                    regions.append(regionSequence.pop(0))
                else:
                    regions.append(regionSequence.pop(0))
                    break
            if len(regions) == 0:
                # Fewer chunks than workers; callBam can't take an empty list.
                continue
            chroms = [region[0] for region in regions]
            paramsNow = params.copy()
            paramsNow["regions"] = regions
            paramsNow["isLearn"] = True
            callArgument = (paramsNow, nn)
            callArguments.append(callArgument)
            regions = []
        results = pool.starmap(callBam, callArguments)
        print(
            "..............Completed bam calling in "
            + str((time.time() - startTime2) / 60)
            + " minutes,merging results................."
        )
        pool.close()
        pool.terminate()
        pool.join()

        mismatch_profile = sum([_[0] for _ in results]).astype(int)
        hp_alt_profile = sum([_[1] for _ in results]).astype(int)
        str_alt_profile = sum([_[2] for _ in results]).astype(int)
        mismatch_dmg_profile = sum([_[3] for _ in results]).astype(int)
        hp_dmg_profile = sum([_[4] for _ in results]).astype(int)
        str_dmg_profile = sum([_[5] for _ in results]).astype(int)
        sbs_alt_bq_hist = sum([_[6] for _ in results]).astype(int)

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
    amp_tn_pd.to_csv(error_prefix + ".amp.tn.txt", sep="\t")
    dmg_tn_pd.to_csv(error_prefix + ".dmg.tn.txt", sep="\t")

    # Indel hp/str row/column labels match Caller.py's auto-learn output.
    hp_row_labels = [f"HP{i}" for i in range(1, 10)] + ["HP10+"]
    hp_col_labels = [f"{b}_{d}" for b in ["A", "T", "C", "G"] for d in (-1, 0, 1)]
    str_row_labels = ["0-1", "2-9", "10-24", "25-39", "40+"]
    str_col_labels = STR_COL_LABELS
    amp_hp_pd = pd.DataFrame(hp_alt_profile, columns=hp_col_labels, index=hp_row_labels)
    dmg_hp_pd = pd.DataFrame(hp_dmg_profile, columns=hp_col_labels, index=hp_row_labels)
    amp_str_pd = pd.DataFrame(
        str_alt_profile, columns=str_col_labels, index=str_row_labels
    )
    dmg_str_pd = pd.DataFrame(
        str_dmg_profile, columns=str_col_labels, index=str_row_labels
    )
    amp_hp_pd.to_csv(error_prefix + ".amp.hp.txt", sep="\t")
    dmg_hp_pd.to_csv(error_prefix + ".dmg.hp.txt", sep="\t")
    amp_str_pd.to_csv(error_prefix + ".amp.str.txt", sep="\t")
    dmg_str_pd.to_csv(error_prefix + ".dmg.str.txt", sep="\t")

    # SBS amp-error counts per error-matrix cell per base quality
    # (64 x 4 x NUM_BQ; bin index == BQ). QC diagnostic, so it goes to
    # tmp_dir rather than ERROR/.
    bq_values = np.arange(NUM_BQ)
    np.savez(
        os.path.join(tmp_dir, os.path.basename(args.output) + ".amp.tn.bqhist.npz"),
        hist=sbs_alt_bq_hist,
        row_labels=np.array(num2trinuc),
        col_labels=np.array(["A", "T", "C", "G"]),
        bq_values=bq_values,
    )

    # SBS SRD rate matrix, EM-estimated from the BQ histogram; this (not
    # amp.tn.txt's raw counts) is what calling loads as the amp error.
    srd_mat = estimate_sbs_srd_rates(
        sbs_alt_bq_hist, params["pseudocount"], fallback=True
    )
    srd_pd = pd.DataFrame(srd_mat, columns=["A", "T", "C", "G"], index=num2trinuc)
    srd_pd.to_csv(error_prefix + ".amp.tn.srd.txt", sep="\t")

    print(
        "..............Completed error learning "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )
