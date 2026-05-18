#!/usr/bin/env python3
import os
import sys
import time
import shutil
from collections import OrderedDict
from multiprocessing import Pool
import errno
import h5py
import subprocess
from Bio import bgzf


import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# import pysam
from matplotlib import pyplot as plt

# from pysam import AlignmentFile as BAM

from .funcs.call import callBam  # , output_masked_mutations
from .funcs.misc import createVcfStrings
from .funcs.misc import splitBamRegions
from .funcs.misc import getAlignmentObject as BAM
from .funcs.misc import check_h5_usable
from pysam import TabixFile as BED
import pysam

# from heapq import nlargest


def check_input_files_exist(args):
    """
    Check if all required input files exist and exit gracefully if not.
    """
    missing_files = []

    # Required files
    if not os.path.exists(args.bam):
        missing_files.append(f"Tumor BAM file: {args.bam}")

    if not os.path.exists(args.reference):
        missing_files.append(f"Reference genome: {args.reference}")

    # Check for associated reference files (h5 files)
    ref_base = os.path.splitext(args.reference)[0]
    ref_h5_files = [f"{ref_base}.h5", f"{args.reference}.ref.h5"]
    tn_h5_files = [f"{ref_base}.tn.h5", f"{args.reference}.tn.h5"]
    hp_h5_files = [f"{ref_base}.hp.h5", f"{args.reference}.hp.h5"]

    bad_h5_files = []
    for label, candidates, expected_ndim in [
        ("Reference", ref_h5_files, 1),
        ("Trinucleotide", tn_h5_files, 1),
        ("Homopolymer", hp_h5_files, 2),
    ]:
        found = next((f for f in candidates if os.path.exists(f)), None)
        if found is None:
            missing_files.append(f"{label} h5 file: {candidates[0]} or {candidates[1]}")
        else:
            ok, msg = check_h5_usable(found, expected_ndim=expected_ndim)
            if not ok:
                bad_h5_files.append(msg)

    if bad_h5_files:
        print(
            "ERROR: One or more reference index files are unreadable or incompatible:"
        )
        for msg in bad_h5_files:
            print(f"  - {msg}")
        print(
            "\nPlease re-index the reference genome with the current version of DupCaller:\n"
            f"  DupCaller.py index -f {args.reference} -s <str_regions.bed.gz>"
        )
        sys.exit(1)

    # Optional files that should be checked if provided
    if args.normalBams:
        for normal_bam in args.normalBams:
            if not os.path.exists(normal_bam):
                missing_files.append(f"Normal BAM file: {normal_bam}")

    if args.germline and not os.path.exists(args.germline):
        missing_files.append(f"Germline VCF file: {args.germline}")

    if args.regionfile and not os.path.exists(args.regionfile):
        missing_files.append(f"Region file: {args.regionfile}")

    if args.noise:
        for noise_bedfile in args.noise:
            if not os.path.exists(noise_bedfile):
                missing_files.append(f"Noise mask BED file: {noise_bedfile}")

    if args.indelbed and not os.path.exists(args.indelbed):
        missing_files.append(f"Indel BED file: {args.indelbed}")

    # Check optional error profile files
    if args.errprefix:
        for suffix, label in [
            (".amp.tn.txt", "Amplification SBS error file"),
            (".amp.id.txt", "Amplification indel error file"),
            (".dmg.tn.txt", "Damage SBS error file"),
            (".dmg.id.txt", "Damage indel error file"),
        ]:
            path = args.errprefix + suffix
            if not os.path.exists(path):
                missing_files.append(f"{label} (from -e): {path}")

    if args.amperrfile and not os.path.exists(args.amperrfile):
        missing_files.append(f"Amplification SBS error file: {args.amperrfile}")

    if args.amperrfileindel and not os.path.exists(args.amperrfileindel):
        missing_files.append(f"Amplification indel error file: {args.amperrfileindel}")

    if args.dmgerrfile and not os.path.exists(args.dmgerrfile):
        missing_files.append(f"Damage SBS error file: {args.dmgerrfile}")

    if args.dmgerrfileindel and not os.path.exists(args.dmgerrfileindel):
        missing_files.append(f"Damage indel error file: {args.dmgerrfileindel}")

    # If any files are missing, print error and exit
    if missing_files:
        print("ERROR: The following required input files are missing:")
        for missing_file in missing_files:
            print(f"  - {missing_file}")
        print("\nPlease ensure all input files exist before running DupCaller.")
        sys.exit(1)


# if __name__ == "__main__":
def do_call(args):
    # Check if all input files exist before proceeding
    check_input_files_exist(args)
    output_dir = args.output
    sample_name = os.path.basename(args.output)
    os.makedirs(output_dir, exist_ok=True)
    args.output = os.path.join(output_dir, sample_name)
    params = {
        "tumorBam": args.bam,
        "normalBams": args.normalBams,
        "germline": args.germline,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regions,
        "region_file": args.regionfile,
        "threads": args.threads,
        # "amperr": args.amperrs,
        # "amperri": args.amperri,
        "amperr_file": args.output + ".amp.tn.txt",
        "amperri_file": args.output + ".amp.id.txt",
        "dmgerr_file": args.output + ".dmg.tn.txt",
        "dmgerri_file": args.output + ".dmg.id.txt",
        "mutRate": args.mutRate,
        "pcutoff": args.thresholdSnv,
        "pcutoffi": args.thresholdIndel,
        "cscutoff": args.scoreSnv,
        "cscutoffi": args.scoreIndel,
        "mapq": args.mapq,
        "noise": args.noise,
        "indel_bed": args.indelbed,
        "trim5": args.trimF,
        "trim3": args.trimR,
        "minNdepth": args.minNdepth,
        "minBq": args.minBq,
        "maxAF": args.maxAF,
        "germline_cutoff": args.germlineAfCutoff,
        "maxNM": args.nmflt,
        "step": args.windowSize,
        "minMeanASXS": args.minMeanASXS,
        "isLearn": None,
        "normalVAF": args.naf,
        "rescue": args.rescue,
        "maxZeroQualFrac": args.maxZeroQualFrac,
        "maxDepth": args.maxPileupDepth,
    }
    if args.errprefix:
        params["amperr_file"] = args.errprefix + ".amp.tn.txt"
        params["amperri_file"] = args.errprefix + ".amp.id.txt"
        params["dmgerr_file"] = args.errprefix + ".dmg.tn.txt"
        params["dmgerri_file"] = args.errprefix + ".dmg.id.txt"
    if args.amperrfile:
        params["amperr_file"] = args.amperrfile
    if args.amperrfileindel:
        params["amperri_file"] = args.amperrfileindel
    if args.dmgerrfile:
        params["dmgerr_file"] = args.dmgerrfile
    if args.dmgerrfileindel:
        params["dmgerri_file"] = args.dmgerrfileindel

    # Write parameter log to the results folder
    log_path = args.output + "_call_params.log"
    with open(log_path, "w") as log:
        log.write(f"DupCaller call — parameter log\n")
        log.write(f"Run time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"Command:  {' '.join(sys.argv)}\n")
        log.write("\n--- All parameters (resolved values) ---\n")
        arg_dict = vars(args)
        for key, value in sorted(arg_dict.items()):
            log.write(f"  {key}: {value}\n")
        log.write("\n--- Resolved error file paths ---\n")
        log.write(f"  amperr_file:  {params['amperr_file']}\n")
        log.write(f"  amperri_file: {params['amperri_file']}\n")
        log.write(f"  dmgerr_file:  {params['dmgerr_file']}\n")
        log.write(f"  dmgerri_file: {params['dmgerri_file']}\n")

    params_learn = {
        "tumorBam": args.bam,
        "normalBams": None,
        "germline": None,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regionst,
        "region_file": None,
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
        "noise": None,
        "indel_bed": None,
        "trim5": args.trimF,
        "trim3": args.trimR,
        "minMeanASXS": args.minMeanASXS,
        "germline_cutoff": args.germlineAfCutoff,
        "minBq": args.minBq,
        "minAltQual": args.minAltQual,
        "maxNM": args.nmflt,
        "step": args.windowSize,
        "minRef": args.minRef,
        "minAlt": args.minAlt,
        "isLearn": True,
        "rescue": False,
    }
    same_regions_flag = False
    if not params_learn["regions"]:
        params_learn["regions"] = params["regions"]
        same_regions_flag = True
    if not params["normalBams"]:
        print(
            f"A matched normal is not used. \
            The maximum allele fraction to call a somatic mutation is set to be {args.maxAF}"
        )
    else:
        print(
            f"Matched normal: {args.normalBams}. \
            The maximum allele fraction to call a somatic mutation is set to be {args.maxAF}"
        )
    """
    Initialze run
    """
    # print("..............Loading reference genome.....................")
    # fasta = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    startTime = time.time()
    results_dir = os.path.dirname(args.output)
    tmp_dir = os.path.join(results_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    params["tmp_dir"] = tmp_dir
    params_learn["tmp_dir"] = tmp_dir
    bamObject = BAM(args.bam, "rb", args.reference)
    if params["region_file"]:
        regionfile = params["region_file"]
    else:
        regionfile = None

    """
    Execulte variant calling
    """

    """
    Learn
    """
    learn_flag = False
    if not (
        os.path.exists(params["amperr_file"])
        and os.path.exists(params["amperri_file"])
        and os.path.exists(params["dmgerr_file"])
        and os.path.exists(params["dmgerri_file"])
    ):
        learn_flag = True
        if args.threads == 1:
            """
            Single-thread execution
            """
            print(".........Starting estimating error rates..............")
            # contigs = [(r.strip('\n'),) for r in open(args.regions,'r').readlines()] # Only process contigs in region file
            paramsNow = params_learn
            # paramsNow["reference"] = fasta
            # paramsNow["isLearn"] = True
            regions = params_learn["regions"]
            paramsNow["regions"] = [
                (chrom, 0, bamObject.get_reference_length(chrom) - 1)
                for chrom in regions
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
            contigs = params_learn["regions"]
            contigLengths = [
                bamObject.get_reference_length(contig) for contig in contigs
            ]
            print(
                "...........Spliting genomic regions for parallel execution (error estimation)................"
            )
            # print(args.threads)
            # if args.normalBam:
            cutSites, chunkSize, contigs = splitBamRegions(
                [args.bam],
                args.threads,
                contigs,
                args.windowSize,
                args.reference,  # , regionfile
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
                + " minutes............"
            )

            """
            Start estimating error rate
            """

            callArguments = []
            startTime2 = time.time()
            print(".........Starting estimating error rate.............")
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
                paramsNow = params_learn.copy()
                paramsNow["regions"] = regions
                # paramsNow["isLearn"] = True
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

            mismatch_profile = sum([_[0] for _ in results])
            indelerr_profile = sum([_[1] for _ in results]).astype(int)
            mismatch_dmg_profile = sum([_[2] for _ in results]).astype(int)
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
        amp_tn_pd.to_csv(args.output + ".amp.tn.txt", sep="\t")
        np.savetxt(
            args.output + ".amp.id.txt",
            indelerr_profile,
            delimiter="\t",
            fmt="%d",
        )
        dmg_tn_pd.to_csv(args.output + ".dmg.tn.txt", sep="\t")
        # np.savetxt(params["output"] + "/" + args.output + ".dmg.tn.txt",np.hstack([trinuc_cols,mismatch_dmg_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")
        np.savetxt(
            args.output + ".dmg.id.txt",
            indelerr_dmg_profile,
            delimiter="\t",
            fmt="%d",
        )
        print(
            "..............Completed error estimation in "
            + str((time.time() - startTime) / 60)
            + " minutes..............."
        )

    if args.threads == 1:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        paramsNow = params
        # paramsNow["reference"] = fasta
        regions = params["regions"]
        paramsNow["regions"] = [
            (chrom, 0, bamObject.get_reference_length(chrom) - 1) for chrom in regions
        ]
        (
            mutsAll,
            coverage,
            rec_num,
            duplex_num,
            duplex_read_num_single,
            duplex_read_num_trinuc_single,
            indelsAll,
            coverage_indel,
            unique_read_num,
            pass_read_num,
            FPAll,
            RPAll,
            unmasked_coverage,
            # unmasked_duplex_read_num_dict_trinuc,
            unmasked_coverage_indel,
        ) = callBam(paramsNow, 0)
        muts_positions = [
            mut["chrom"] + str(mut["pos"]) + mut["ref"] + mut["alt"] for mut in mutsAll
        ]
        muts_dict = dict()
        take_ind = list()
        muts_num = len(mutsAll)
        indels_positions = [
            indel["chrom"] + str(indel["pos"]) + indel["ref"] + ":" + indel["alt"]
            for indel in indelsAll
        ]
        indels_num = len(indelsAll)
        duplex_combinations = list(duplex_read_num_single.keys())
        duplex_combinations.sort()
        duplex_read_num = OrderedDict(
            {num: duplex_read_num_single[num][0] for num in duplex_combinations}
        )
        duplex_coverage_by_group = OrderedDict(
            {num: duplex_read_num_single[num][1] for num in duplex_combinations}
        )
        total_family_num = OrderedDict(
            {num: duplex_read_num_single[num][2] for num in duplex_combinations}
        )
        duplex_read_num_trinuc = OrderedDict(
            {num: duplex_read_num_trinuc_single[num] for num in duplex_combinations}
        )

    else:
        """
        Multi-thread execution
        """
        # args.threads = args.threads - 1
        regions_list = list()
        regionSequence = []
        if (not same_regions_flag) or (not learn_flag):
            contigs = params["regions"]
            contigLengths = [
                bamObject.get_reference_length(contig) for contig in contigs
            ]
            print("....Splitting genomic regions for parallel execution.....")
            if args.normalBams:
                cutSites, chunkSize, contigs = splitBamRegions(
                    [args.bam],
                    args.threads,
                    contigs,
                    args.windowSize,
                    args.reference,  # , regionfile
                )
            else:
                cutSites, chunkSize, contigs = splitBamRegions(
                    [args.bam],
                    args.threads,
                    contigs,
                    args.windowSize,
                    args.reference,  # , regionfile
                )
            currentContigIndex = 0
            usedTime = (time.time() - startTime) / 60
            print(f"....Genomic regions splitted in {usedTime} minutes...")
        """
        Determine regions for each process
        """
        # print(cutSites)
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
            # print(regions)
            if len(regions) == 0:
                raise Exception(
                    f"Window size {args.windowSize} is proabably too large. Change -w to a smaller value. An ideal value will be targetSize/(threads * 100)"
                )
            paramsNow = params.copy()
            regions_list.append(regions)
            # paramsNow["reference"] = fastaNow
            paramsNow["regions"] = regions
            callArgument = (paramsNow, nn)
            callArguments.append(callArgument)
        results = pool.starmap(callBam, callArguments)

        muts = [_[0] for _ in results]
        coverages = [_[1] for _ in results]
        rec_nums = [_[2] for _ in results]
        duplex_nums = [_[3] for _ in results]
        duplex_read_nums = [_[4] for _ in results]
        duplex_read_nums_trinuc = [_[5] for _ in results]
        indels = [_[6] for _ in results]
        coverages_indels = [_[7] for _ in results]
        unique_read_nums = [_[8] for _ in results]
        pass_read_nums = [_[9] for _ in results]
        FPs = [_[10] for _ in results]
        RPs = [_[11] for _ in results]
        unmasked_coverages = [_[12] for _ in results]
        # unmasked_duplex_read_nums_trinuc = [_[13] for _ in results]
        unmasked_coverages_indels = [_[13] for _ in results]
        print(
            "..............Completed bam calling in "
            + str((time.time() - startTime2) / 60)
            + " minutes,estimating mutation rates................."
        )
        pool.close()
        pool.terminate()
        pool.join()
        mutsAll = sum(muts, [])
        muts_positions = [
            mut["chrom"] + str(mut["pos"]) + mut["ref"] + mut["alt"] for mut in mutsAll
        ]
        muts_dict = dict()
        take_ind = list()
        muts_num = len(mutsAll)
        coverage = sum(coverages)
        coverage_indel = sum(coverages_indels)
        unmasked_coverage = sum(unmasked_coverages)
        unmasked_coverage_indel = sum(unmasked_coverages_indels)
        rec_num = sum(rec_nums)
        duplex_num = sum(duplex_nums)
        unique_read_num = sum(unique_read_nums)
        pass_read_num = sum(pass_read_nums)
        indelsAll = sum(indels, [])
        indels_num = len(indelsAll)
        indels_positions = [
            indel["chrom"] + str(indel["pos"]) + indel["ref"] + ":" + indel["alt"]
            for indel in indelsAll
        ]

        duplex_combinations = list(
            set.union(*[set(d.keys()) for d in duplex_read_nums])
        )
        duplex_combinations.sort()
        duplex_read_num = OrderedDict(
            {
                num: sum([d.get(num, [0, 0, 0])[0] for d in duplex_read_nums])
                for num in duplex_combinations
            }
        )
        duplex_coverage_by_group = OrderedDict(
            {
                num: sum([d.get(num, [0, 0, 0])[1] for d in duplex_read_nums])
                for num in duplex_combinations
            }
        )
        total_family_num = OrderedDict(
            {
                num: sum([d.get(num, [0, 0, 0])[2] for d in duplex_read_nums])
                for num in duplex_combinations
            }
        )

        duplex_read_num_trinuc = OrderedDict(
            {
                num: sum([d.get(num, np.zeros(32)) for d in duplex_read_nums_trinuc])
                for num in duplex_combinations
            }
        )
        """
        # Combine unmasked trinucleotide data
        unmasked_duplex_read_num_dict_trinuc = OrderedDict(
            {
                num: sum([d.get(num, np.zeros(32)) for d in unmasked_duplex_read_nums_trinuc])
                for num in duplex_combinations
            }
        )
        """

        FPAll = sum(FPs, [])
        RPAll = sum(RPs, [])

    tBam = BAM(args.bam, "rb", args.reference)
    contigs = tBam.references
    # print(contigs)
    chromDict = {contig: tBam.get_reference_length(contig) for contig in contigs}
    infoDict = {
        "F1R2": [1, "Integer", "Number of F1R2 read(s) in the read bundle"],
        "F2R1": [1, "Integer", "Number of F2R1 read(s) in the read bundle"],
        # "TLR": [1, "Float", "Alt/Ref log likelihood ratio of top strand"],
        # "BLR": [1, "Float", "Alt/Ref log likelihood ratio of bottom strand"],
        "CS": [1, "Float", "confidence score of major base over minor base"],
        "LR": [1, "Float", "Log-Likelihood ratio of major base over minor base"],
        "LM": [
            1,
            "Float",
            "maximum log-Likelihood ratio of major base over minor base",
        ],
        "TC": [4, "Integer", "Top strand base count"],
        "BC": [4, "Float", "Bottom strand base count"],
        "DF": [1, "Integer", "Distance from fragment end"],
        "DR": [1, "Integer", "Distance from read end"],
        "TAG1": [1, "String", "Barcode of top strand 5 prime"],
        "TAG2": [1, "String", "Barcode of bottom strand 5 prime"],
        "SP": [1, "Integer", "Read family reference start position"],
        "TN": [1, "String", "trinucleotide context"],
        "HP": [1, "Integer", "Homopolymer length. Always 0 for SBS"],
        "STR": [
            1,
            "Integer",
            "reference allele length bin of short tandem repeats. 0: no STR or less than 10bp. 1: 10-24bp. 2: 25-39bp. 3: 40bp+. Always 0 for SBS",
        ],
    }
    formatDict = {
        "AC": [1, "Integer", "Count of alt allele"],
        "RC": [1, "Integer", "Count of ref allele"],
        "DP": [1, "Integer", "Depth at the location"],
    }
    filterDict = {
        "PASS": "All filter Passed",
        "masked": "Mutation filtered by noise mask",
        "underpowered": "Mutation that are underpowered under set threshold",
        "high_nm": "Mutation filtered by high NM tag filter",
    }

    # Separate mutations by filter type
    pass_muts = [mut for mut in mutsAll if mut.get("filter", "PASS") == "PASS"]
    masked_muts = [mut for mut in mutsAll if mut.get("filter") == "masked"]
    pass_indels = [
        indel for indel in indelsAll if indel.get("filter", "PASS") == "PASS"
    ]
    masked_indels = [indel for indel in indelsAll if indel.get("filter") == "masked"]

    # Create VCF with all mutations (PASS and masked)
    vcfLines = createVcfStrings(chromDict, infoDict, formatDict, filterDict, mutsAll)
    with open(args.output + "_snv.vcf", "w") as vcf:
        vcf.write(vcfLines)

    vcfLines = createVcfStrings(chromDict, infoDict, formatDict, filterDict, indelsAll)
    with open(args.output + "_indel.vcf", "w") as vcf:
        vcf.write(vcfLines)

    burden_naive = muts_num / coverage if coverage > 0 else 0.0
    indel_burden = indels_num / coverage if coverage > 0 else 0.0
    efficiency = duplex_num / rec_num if rec_num > 0 else 0.0
    pass_duprate = unique_read_num / pass_read_num if pass_read_num > 0 else 0.0

    # Build muts_by_duplex_group over all keys (including 0+n / n+0)
    muts_by_duplex_group = {k: 0 for k in duplex_read_num.keys()}
    for mut in mutsAll:
        TC_total = int(mut["infos"]["F1R2"])
        BC_total = int(mut["infos"]["F2R1"])
        key_fwd = str(TC_total) + "+" + str(BC_total)
        key_rev = str(BC_total) + "+" + str(TC_total)
        if key_fwd in muts_by_duplex_group:
            muts_by_duplex_group[key_fwd] += 1
        elif key_rev in muts_by_duplex_group:
            muts_by_duplex_group[key_rev] += 1

    # All keys, including 0+n and n+0 (read_set_count may be 0 for duplex-only counts
    # but the read set itself exists); include every key present in duplex_read_num
    all_keys = sorted(
        duplex_read_num.keys(),
        key=lambda s: (int(s.split("+")[0]), int(s.split("+")[1])),
    )
    with open(args.output + "_duplex_family_strand_composition.txt", "w") as f:
        f.write(
            "duplex_group_strand_composition\tread_family_number\tduplex_group_number\teffective_coverage\n"
        )
        for read_num in all_keys:
            f.write(
                f"{read_num}\t{total_family_num[read_num]}\t{duplex_read_num[read_num]}\t{duplex_coverage_by_group[read_num]}\n"
            )

    # Heatmap: x = F1R2 count, y = F2R1 count, value = proportion of read sets
    max_val = max(
        (max(int(k.split("+")[0]), int(k.split("+")[1])) for k in all_keys),
        default=0,
    )
    total_read_sets = sum(total_family_num.values())
    # Determine n such that families with either strand > n are < 1% of the library
    n = 0
    if total_read_sets > 0:
        for candidate in range(max_val + 1):
            outlier_count = sum(
                v
                for k, v in total_family_num.items()
                if max(int(k.split("+")[0]), int(k.split("+")[1])) > candidate
            )
            if outlier_count / total_read_sets < 0.01:
                n = candidate
                break
        else:
            n = max_val
    heatmap_data = np.zeros((n + 1, n + 1), dtype=float)
    for k in all_keys:
        x, y = int(k.split("+")[0]), int(k.split("+")[1])
        if x <= n and y <= n:
            heatmap_data[y, x] = total_family_num[k] / total_read_sets
    fig, ax = plt.subplots(figsize=(max(6, n + 2), max(5, n + 1)))
    im = ax.imshow(heatmap_data, aspect="auto", origin="lower", cmap="YlOrRd")
    ax.set_xlabel("F1R2 read count (top strand)")
    ax.set_ylabel("F2R1 read count (bottom strand)")
    ax.set_title("Proportion of read sets by duplex group composition")
    ax.set_xticks(range(n + 1))
    ax.set_yticks(range(n + 1))
    plt.colorbar(im, ax=ax, label="Proportion of read sets")
    plt.tight_layout()
    fig.savefig(args.output + "_duplex_family_strand_composition_heatmap.pdf")
    plt.close(fig)

    non_zero_keys = [k for k in all_keys if duplex_read_num[k] != 0]
    trinuc_list = list()
    trinuc2num = dict()
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                trinuc2num[minus_base + ref_base + plus_base] = len(trinuc_list)
                trinuc_list.append(minus_base + ref_base + plus_base)
    duplex_read_num_trinuc = {_: duplex_read_num_trinuc[_] for _ in non_zero_keys}
    trinuc_by_duplex_group = pd.DataFrame(duplex_read_num_trinuc)
    trinuc_by_duplex_group.insert(0, "", trinuc_list)
    trinuc_by_duplex_group.to_csv(
        args.output + "_trinuc_by_duplex_group.txt",
        sep="\t",
        index=False,
    )

    # Also output unmasked trinucleotide counts
    # unmasked_trinuc_by_duplex_group = pd.DataFrame(unmasked_duplex_read_num_dict_trinuc)
    # unmasked_trinuc_by_duplex_group.insert(0, "", trinuc_list)
    # unmasked_trinuc_by_duplex_group.to_csv(
    # args.output + "_unmasked_trinuc_by_duplex_group.txt",
    # sep="\t",
    # index=False,
    # )

    with open(args.output + "_stats.txt", "w") as f:
        f.write(f"Number of Read Families\t{unique_read_num}\n")
        f.write(f"Number of Pass-filter Reads\t{pass_read_num}\n")
        f.write(f"Number of Effective Read Families\t{duplex_num}\n")
        f.write(f"Effective Coverage\t{coverage}\n")
        f.write(f"Unmasked Coverage\t{unmasked_coverage}\n")
        f.write(f"Effective Indel Coverage\t{coverage_indel}\n")
        f.write(f"Unmasked Indel Coverage\t{unmasked_coverage_indel}\n")
        f.write(
            f"Per Read Family Coverage \t{coverage/duplex_num if duplex_num > 0 else 0.0}\n"
        )
        f.write(
            f"Pass-filter Duplication Rate\t\
        {1-unique_read_num/pass_read_num if pass_read_num > 0 else 0.0}\n"
        )
        f.write(f"Efficiency\t{efficiency}\n")

    print(
        "..............Completed variant calling "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )

    # Merge and combine coverage files after multi-threaded calling
    mergeStartTime = time.time()
    sample_name = os.path.basename(args.output)
    sample_dir = params["tmp_dir"]
    if args.threads > 1:
        merge_and_combine_coverage_files(sample_name, sample_dir, args.threads)
    else:
        # Single-threaded: rename the single coverage file and index it
        single_cov = os.path.join(sample_dir, f"{sample_name}_0_coverage.bed.gz")
        merged_cov = os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz")
        subprocess.run(f"mv {single_cov} {merged_cov}", shell=True, check=True)
        subprocess.run(f"tabix -f -p bed {merged_cov}", shell=True, check=True)
    subprocess.run(
        f"mv {os.path.join(sample_dir, f'{sample_name}_coverage.bed.gz')} {os.path.dirname(args.output)}/",
        shell=True,
        check=True,
    )
    subprocess.run(
        f"mv {os.path.join(sample_dir, f'{sample_name}_coverage.bed.gz.tbi')} {os.path.dirname(args.output)}/",
        shell=True,
        check=True,
    )
    shutil.rmtree(params["tmp_dir"])
    print(
        "..............Completed coverage merging in "
        + str((time.time() - mergeStartTime) / 60)
        + " minutes..............."
    )


def merge_and_combine_coverage_files(sample_name, sample_dir, nprocess):
    """
    1. Merge adjacent next_region and prev_region files by summing columns 4 and 5
    2. Combine all files in the correct order using cat
    """
    print("Merging adjacent region files and combining coverage files...")

    # Step 1: Create overlap files by merging adjacent regions
    for n in range(nprocess - 1):
        next_file = os.path.join(
            sample_dir, f"{sample_name}_{n}_coverage_next_region.tmp.bed.gz"
        )
        prev_file = os.path.join(
            sample_dir, f"{sample_name}_{n+1}_coverage_prev_region.tmp.bed.gz"
        )
        overlap_file = os.path.join(
            sample_dir, f"{sample_name}_{n}_{n+1}_overlap_coverage.tmp.bed.gz"
        )

        merge_adjacent_bed_files(next_file, prev_file, overlap_file)

    # Step 2: Create list of files in the correct order
    files_to_combine = []

    for n in range(nprocess):
        # Add main coverage file
        main_file = os.path.join(sample_dir, f"{sample_name}_{n}_coverage.bed.gz")
        if not os.path.exists(main_file):
            raise FileNotFoundError(f"Expected coverage file not found: {main_file}")
        files_to_combine.append(main_file)

        # Add overlap file (except for the last process)
        if n < nprocess - 1:
            overlap_file = os.path.join(
                sample_dir, f"{sample_name}_{n}_{n+1}_overlap_coverage.tmp.bed.gz"
            )
            if not os.path.exists(overlap_file):
                raise FileNotFoundError(
                    f"Expected overlap coverage file not found: {overlap_file}"
                )
            files_to_combine.append(overlap_file)

    # Step 3: Combine files using cat command
    if files_to_combine:
        final_output = os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz")

        # Use cat to combine files
        cmd = f"cat {' '.join(files_to_combine)} > {final_output}"

        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"Combined coverage file created: {final_output}")

            # Index the combined bed file with tabix
            index_cmd = f"tabix -f -p bed {final_output}"
            try:
                subprocess.run(index_cmd, shell=True, check=True)
                print(f"Tabix index created for: {final_output}")
            except subprocess.CalledProcessError as e:
                print(f"Warning: Could not create tabix index: {e}")

        except subprocess.CalledProcessError as e:
            print(f"Error combining files: {e}")

        # Clean up temporary files
        cleanup_temp_files(sample_name, sample_dir, nprocess)
    else:
        print("No coverage files found to combine")


def merge_adjacent_bed_files(next_file, prev_file, output_file):
    """
    Merge two adjacent bed files by summing columns 4 and 5 for matching positions
    """
    coverage_dict = {}

    # Read next_region file
    if not os.path.exists(next_file):
        raise FileNotFoundError(f"Expected file for merging not found: {next_file}")
    with bgzf.open(next_file, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split("\t")
                if len(parts) >= 5:
                    chrom, start, end, cov1, cov2 = parts[:5]
                    key = (chrom, start, end)
                    coverage_dict[key] = [int(cov1), int(cov2)]

    # Read prev_region file and merge
    if not os.path.exists(prev_file):
        raise FileNotFoundError(f"Expected file for merging not found: {prev_file}")
    with bgzf.open(prev_file, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split("\t")
                if len(parts) >= 5:
                    chrom, start, end, cov1, cov2 = parts[:5]
                    key = (chrom, start, end)
                    if key in coverage_dict:
                        coverage_dict[key][0] += int(cov1)
                        coverage_dict[key][1] += int(cov2)
                    else:
                        coverage_dict[key] = [int(cov1), int(cov2)]

    # Write merged result only if there's data
    with bgzf.open(output_file, "wt") as f:
        for (chrom, start, end), (cov1, cov2) in sorted(coverage_dict.items()):
            f.write(f"{chrom}\t{start}\t{end}\t{cov1}\t{cov2}\n")


def cleanup_temp_files(sample_name, sample_dir, nprocess):
    """
    Clean up temporary files
    """
    for n in range(nprocess):
        temp_files = [
            os.path.join(
                sample_dir, f"{sample_name}_{n}_coverage_prev_region.tmp.bed.gz"
            ),
            os.path.join(
                sample_dir, f"{sample_name}_{n}_coverage_next_region.tmp.bed.gz"
            ),
            os.path.join(sample_dir, f"{sample_name}_{n}_coverage.bed.gz"),
        ]

        if n < nprocess - 1:
            temp_files.append(
                os.path.join(
                    sample_dir, f"{sample_name}_{n}_{n+1}_overlap_coverage.tmp.bed.gz"
                )
            )

        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
