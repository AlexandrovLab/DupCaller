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
from .funcs.learn import NUM_BQ, estimate_sbs_srd_rates
from .funcs.misc import createVcfStrings
from .funcs.misc import splitBamRegions
from .funcs.misc import drop_empty_regions
from .funcs.misc import getAlignmentObject as BAM


# if __name__ == "__main__":
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
        # callBam's isLearn path shares the round-0/round-1/round-2 code
        # with call.py (funcs/call.py's threshold_rng Monte Carlo seeding,
        # _compute_read_label's NanoSeq-bam check) -- both keys are read
        # unconditionally regardless of isLearn, so they must be present
        # even though a standalone `learn` run has no call-specific CLI
        # flags for them.
        "seed": int(np.random.SeedSequence().generate_state(1)[0]),
        "nanoSeqBam": False,
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
    # Learned error matrices get their own ERROR/ subfolder, matching
    # Caller.py's internal auto-learn convention.
    error_dir = os.path.join(params["output"], "ERROR")
    os.makedirs(error_dir, exist_ok=True)
    error_prefix = os.path.join(error_dir, args.output)
    bamObject = BAM(args.bam, "rb", args.reference)

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
            hp_alt_profile,
            str_alt_profile,
            mismatch_dmg_profile,
            hp_dmg_profile,
            str_dmg_profile,
            sbs_alt_bq_hist,
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
                # Fewer real chunks than -p requested (regionSequence
                # already ran out) -- just use fewer effective parallel
                # chunks for this pass instead of handing callBam an
                # empty regions list, which it can't index into
                # (regions[0]).
                continue
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
    # np.savetxt(params["output"] + "/" + args.output + ".amp.tn.txt",np.hstack([trinuc_cols[0:32],mismatch_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")
    amp_tn_pd.to_csv(error_prefix + ".amp.tn.txt", sep="\t")
    dmg_tn_pd.to_csv(error_prefix + ".dmg.tn.txt", sep="\t")
    # np.savetxt(params["output"] + "/" + args.output + ".dmg.tn.txt",np.hstack([trinuc_cols,mismatch_dmg_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")

    # Indel hp/str: see Caller.py's internal auto-learn write-out for the
    # column/row convention (same one here).
    hp_row_labels = [f"HP{i}" for i in range(1, 10)] + ["HP10+"]
    hp_col_labels = [f"{b}_{d}" for b in ["A", "T", "C", "G"] for d in (-1, 0, 1)]
    str_row_labels = ["0-1", "2-9", "10-24", "25-39", "40+"]
    str_col_labels = [str(d) for d in range(-5, 6)]
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

    # Amp-error base-quality histogram (SBS only): not weighted by BQ
    # like mismatch_profile above -- instead records, per error-matrix
    # cell, how many raw amp-error observations were seen at each base
    # quality. Too large to flatten sensibly into a 2D CSV (64x4x
    # {NUM_BQ}), so saved as .npz with the row/column/BQ axis labels
    # alongside the counts -- bq_values (0..{NUM_BQ}-1) doubles as both
    # the BQ axis labels and, since a bin's index equals its BQ value,
    # the histogram's per-slice BQ values. A debugging/QC diagnostic, not
    # a deliverable, so it goes to tmp_dir rather than the final ERROR/
    # output dir -- same split as Caller.py's _indel_rate_by_hp_str.txt/
    # _sbs96_rate_by_trinuc.txt.
    bq_values = np.arange(NUM_BQ)
    np.savez(
        os.path.join(tmp_dir, os.path.basename(args.output) + ".amp.tn.bqhist.npz"),
        hist=sbs_alt_bq_hist,
        row_labels=np.array(num2trinuc),
        col_labels=np.array(["A", "T", "C", "G"]),
        bq_values=bq_values,
    )

    # SBS SRD (single-read-damage) rate matrix, EM-estimated from the BQ
    # histogram above -- this, not amp.tn.txt's raw counts, is what
    # params["amperr_file"] now points at and load_error_matrices
    # (funcs/misc.py) reads directly for calling, replacing the old
    # in-situ row-normalization of raw counts.
    srd_mat = estimate_sbs_srd_rates(sbs_alt_bq_hist, params["pseudocount"])
    srd_pd = pd.DataFrame(srd_mat, columns=["A", "T", "C", "G"], index=num2trinuc)
    srd_pd.to_csv(error_prefix + ".amp.tn.srd.txt", sep="\t")

    print(
        "..............Completed error learning "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )
