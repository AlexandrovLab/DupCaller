#!/usr/bin/env python3
import os
import sys
import time
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

from . import __version__
from .funcs.call import callBam  # , output_masked_mutations
from .funcs.learn import NUM_BQ, estimate_sbs_srd_rates
from .funcs.misc import simulate_power_grid, load_error_matrices
from .funcs.misc import init_refine_worker, refine_channel_task
from .funcs.prob import indelErrorProbs, indelMaxLR
from .funcs.misc import createVcfStrings
from .funcs.misc import splitBamRegions
from .funcs.misc import drop_empty_regions
from .funcs.misc import getAlignmentObject as BAM
from .funcs.misc import check_h5_usable
from .funcs.misc import build_trinuc64_order, build_trinuc192_labels
from .funcs.misc import _build_sbs96_labels
from .funcs.misc import _REVCOMP
from .funcs.misc import build_indel100_labels
from .funcs.misc import build_dbs_raw144_labels
from .funcs.misc import _ensure_type_subdirs
from .funcs.misc import load_repeat_context
from pysam import TabixFile as BED
import pysam

# from heapq import nlargest
# refine_channel_task/init_refine_worker (the per-channel FDR-refinement
# Pool it and _refine_channel run inside of) live in funcs/call.py -- see
# below in do_call.


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
    str_h5_files = [f"{ref_base}.str.h5", f"{args.reference}.str.h5"]

    bad_h5_files = []
    for label, candidates, expected_ndim in [
        ("Reference", ref_h5_files, 1),
        ("Trinucleotide", tn_h5_files, 1),
        ("Homopolymer", hp_h5_files, 2),
        ("STR repeat", str_h5_files, 2),
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
            (".amp.tn.srd.txt", "Amplification SBS error file (SRD rate matrix)"),
            (".amp.hp.txt", "Amplification indel (homopolymer) error file"),
            (".amp.str.txt", "Amplification indel (STR) error file"),
            (".dmg.tn.txt", "Damage SBS error file"),
            (".dmg.hp.txt", "Damage indel (homopolymer) error file"),
            (".dmg.str.txt", "Damage indel (STR) error file"),
        ]:
            path = args.errprefix + suffix
            if not os.path.exists(path):
                missing_files.append(f"{label} (from -E): {path}")

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
    # Resolve the Monte Carlo base seed once, up front: if the user didn't
    # pass --seed, generate one now and write it back onto args so it (a)
    # gets threaded through params below and (b) shows up in the
    # _call_params.log's "resolved values" dump (sourced from vars(args))
    # for later reuse -- reproducing a run with no explicit --seed still
    # requires knowing what seed it actually used.
    if args.seed is None:
        args.seed = int(np.random.SeedSequence().generate_state(1)[0])
    output_dir = args.output
    sample_name = os.path.basename(args.output)
    os.makedirs(output_dir, exist_ok=True)
    args.output = os.path.join(output_dir, sample_name)
    # Learned error matrices (.amp/.dmg .tn/.hp/.str.txt) get their own
    # ERROR/ subfolder, same idea as SBS/INDEL/DBS -- keeps the sample's
    # top level to just the shared, whole-sample outputs.
    error_dir = os.path.join(output_dir, "ERROR")
    os.makedirs(error_dir, exist_ok=True)
    error_prefix = os.path.join(error_dir, sample_name)
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
        "amperr_file": error_prefix + ".amp.tn.srd.txt",
        "amperri_file": error_prefix + ".amp.id.txt",
        "dmgerr_file": error_prefix + ".dmg.tn.txt",
        "dmgerri_file": error_prefix + ".dmg.id.txt",
        "mutRate": args.mutRate,
        # ts/ti removed: every candidate with LR>0 is treated as a
        # candidate mutation (see the >0 gates in funcs/call.py); FDR
        # refinement (-fdr) is the only thing that can raise this baseline
        # per-channel afterward.
        "pcutoff": 0,
        "pcutoffi": 0,
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
        "minGroupAmp": args.minGroupAmp,
        "minGroupDmg": args.minGroupDmg,
        "nanoSeqBam": args.NanoSeqBam,
        "seed": args.seed,
    }
    if args.errprefix:
        params["amperr_file"] = args.errprefix + ".amp.tn.srd.txt"
        params["amperri_file"] = args.errprefix + ".amp.id.txt"
        params["dmgerr_file"] = args.errprefix + ".dmg.tn.txt"
        params["dmgerri_file"] = args.errprefix + ".dmg.id.txt"

    # Write parameter log to the results folder
    log_path = args.output + "_call_params.log"
    with open(log_path, "w") as log:
        log.write(f"DupCaller call — parameter log\n")
        log.write(f"Version:  {__version__}\n")
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
        # Round 0 (error-rate learning) shares callBam with rounds 1/2,
        # which requires params["seed"] (funcs/call.py's threshold_rng-based
        # Monte Carlo seeding). Same base seed as params so a single -seed
        # value seeds the whole run reproducibly, not just rounds 1/2.
        "seed": args.seed,
        # _compute_read_label (funcs/call.py) is shared unconditionally
        # across all rounds and reads params.get("nanoSeqBam"); without
        # this key here, a --NanoSeqBam run falls back to the DB tag
        # during round 0 even though NanoSeq-format bams only carry MB/RB.
        "nanoSeqBam": args.NanoSeqBam,
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
    amperri_hp_file = params["amperri_file"].replace(".amp.id.txt", ".amp.hp.txt")
    amperri_str_file = params["amperri_file"].replace(".amp.id.txt", ".amp.str.txt")
    dmgerri_hp_file = params["dmgerri_file"].replace(".dmg.id.txt", ".dmg.hp.txt")
    dmgerri_str_file = params["dmgerri_file"].replace(".dmg.id.txt", ".dmg.str.txt")
    if not (
        os.path.exists(params["amperr_file"])
        and os.path.exists(amperri_hp_file)
        and os.path.exists(amperri_str_file)
        and os.path.exists(params["dmgerr_file"])
        and os.path.exists(dmgerri_hp_file)
        and os.path.exists(dmgerri_str_file)
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
                args.reference,
                params_learn["region_file"],
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
                + " minutes............"
            )

            """
            Start estimating error rate
            """

            callArguments = []
            startTime2 = time.time()
            print(".........Starting estimating error rate.............")
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
                    # already ran out) -- just use fewer effective
                    # parallel chunks for this pass instead of handing
                    # callBam an empty regions list, which it can't index
                    # into (regions[0]).
                    continue
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

        # Indel hp/str: see funcs/learn.py's profileTriNucMismatches for
        # the amp/dmg accumulation and funcs/prob.py's indelErrorProbs for
        # the matching call-time selection. hp.txt: 10 rows (HP1..HP9,
        # HP10+), 12 columns (base x idLen in {-1,0,1}, idLen=0 the
        # reference/opportunity count). str.txt: 5 rows (STR-length bin
        # "0-1"=not a real repeat through "40+"), 11 columns (idLen -5..5,
        # idLen=0 the opportunity count).
        hp_row_labels = [f"HP{i}" for i in range(1, 10)] + ["HP10+"]
        hp_col_labels = [f"{b}_{d}" for b in ["A", "T", "C", "G"] for d in (-1, 0, 1)]
        str_row_labels = ["0-1", "2-9", "10-24", "25-39", "40+"]
        str_col_labels = [str(d) for d in range(-5, 6)]
        amp_hp_pd = pd.DataFrame(
            hp_alt_profile, columns=hp_col_labels, index=hp_row_labels
        )
        dmg_hp_pd = pd.DataFrame(
            hp_dmg_profile, columns=hp_col_labels, index=hp_row_labels
        )
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

        # Amp-error BQ histogram (SBS only) -- see Learn.py's matching
        # write-out for the shape/convention (same one here). Too large
        # for a 2D CSV, so saved as .npz with the row/column/BQ axis
        # labels alongside the counts. A debugging/QC diagnostic, not a
        # deliverable, so it goes to tmp_dir rather than the final ERROR/
        # output dir -- same split as _indel_rate_by_hp_str.txt/
        # _sbs96_rate_by_trinuc.txt below.
        bq_values = np.arange(NUM_BQ)
        np.savez(
            os.path.join(params["tmp_dir"], sample_name + ".amp.tn.bqhist.npz"),
            hist=sbs_alt_bq_hist,
            row_labels=np.array(num2trinuc),
            col_labels=np.array(["A", "T", "C", "G"]),
            bq_values=bq_values,
        )

        # SBS SRD (single-read-damage) rate matrix, EM-estimated from the
        # BQ histogram above -- this, not amp.tn.txt's raw counts, is what
        # params["amperr_file"] now points at and load_error_matrices
        # (funcs/misc.py) reads directly for calling, replacing the old
        # in-situ row-normalization of raw counts.
        srd_mat = estimate_sbs_srd_rates(sbs_alt_bq_hist, args.pseudocount)
        srd_pd = pd.DataFrame(srd_mat, columns=["A", "T", "C", "G"], index=num2trinuc)
        srd_pd.to_csv(error_prefix + ".amp.tn.srd.txt", sep="\t")
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
        # Single "worker" region list, so round 2 below (always dispatched
        # via Pool.starmap over regions_list) works the same way regardless
        # of thread count -- Pool with one task is just a slightly less
        # direct callBam(..., 0) call.
        regions_list = [paramsNow["regions"]]
        (
            mutsAll,
            coverage,
            rec_num,
            duplex_num,
            duplex_read_num_single,
            duplex_read_num_trinuc_single,
            indelsAll,
            unique_read_num,
            pass_read_num,
            FPAll,
            RPAll,
            unmasked_coverage,
            # unmasked_duplex_read_num_dict_trinuc,
            coverage_indel_cat,
            unmasked_coverage_indel_cat,
            duplex_read_num_indel_single,
            dbsAll,
            duplex_read_num_dbs_single,
            depth_by_trinuc,
            depth_by_hpstr,
            L_power,
            L_indel_1bp,
            L_indel_len,
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
        duplex_read_num_indel = OrderedDict(
            {num: duplex_read_num_indel_single[num] for num in duplex_combinations}
        )
        duplex_read_num_dbs = OrderedDict(
            {num: duplex_read_num_dbs_single[num] for num in duplex_combinations}
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
                    args.reference,
                    regionfile,
                )
            else:
                cutSites, chunkSize, contigs = splitBamRegions(
                    [args.bam],
                    args.threads,
                    contigs,
                    args.windowSize,
                    args.reference,
                    regionfile,
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
        regionSequence = drop_empty_regions(regionSequence, bamObject)
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
            # print(regions)
            if len(regions) == 0:
                if nn == 0:
                    raise Exception(
                        f"Window size {args.windowSize} is proabably too large. Change -w to a smaller value. An ideal value will be targetSize/(threads * 100)"
                    )
                # Fewer real chunks than -p requested -- e.g. splitBamRegions
                # skipped a duplicate chunk boundary (see its own comment),
                # or the target region just doesn't have enough distinct
                # windows for this many threads. Neither is an error: just
                # continue with fewer effective parallel chunks.
                print(
                    f"....Only {nn} parallel chunk(s) available (requested -p "
                    f"{args.threads}); continuing with fewer effective threads "
                    "for this pass...."
                )
                break
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
        unique_read_nums = [_[7] for _ in results]
        pass_read_nums = [_[8] for _ in results]
        FPs = [_[9] for _ in results]
        RPs = [_[10] for _ in results]
        unmasked_coverages = [_[11] for _ in results]
        coverages_indels_cat = [_[12] for _ in results]
        unmasked_coverages_indels_cat = [_[13] for _ in results]
        duplex_read_nums_indel = [_[14] for _ in results]
        dbs_muts = [_[15] for _ in results]
        duplex_read_nums_dbs = [_[16] for _ in results]
        depth_by_trinuc = sum(_[17] for _ in results)
        depth_by_hpstr = sum(_[18] for _ in results)
        # L/L_indel_1bp/L_indel_len are deterministic given params up to
        # Monte Carlo noise (each process re-simulates its own copy from
        # the same ampmat/dmgmat_indel); averaging across processes just
        # reduces that noise, it isn't combining different information.
        L_power = np.mean(np.stack([_[19] for _ in results]), axis=0)
        L_indel_1bp = np.mean(np.stack([_[20] for _ in results]), axis=0)
        L_indel_len = np.mean(np.stack([_[21] for _ in results]), axis=0)
        dbsAll = sum(dbs_muts, [])
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
        coverage_indel_cat = sum(coverages_indels_cat)
        unmasked_coverage = sum(unmasked_coverages)
        unmasked_coverage_indel_cat = sum(unmasked_coverages_indels_cat)
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

        # duplex_combinations/duplex_read_num/duplex_coverage_by_group/
        # total_family_num/duplex_read_num_trinuc/_indel/_dbs are real,
        # position-weighted coverage data (cov_mat-derived) that round 1 no
        # longer computes (coverage_only gating in call.py) -- rebuilt
        # further down from round 2's results, once final thresholds are
        # known.
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
        "SNPM": [
            1,
            "Integer",
            "1 if this position falls in snp_mask. A candidate blocked only by snp_mask/noise_mask (see NOISEM) is still fully evaluated (LR + depth) and feeds the unmasked-burden numerator, but stays filter=masked rather than PASS -- never set on a PASS record. SBS only, always 0 for indel.",
        ],
        "NOISEM": [
            1,
            "Integer",
            "1 if this position (or, for indels, any base of its span) falls in noise_mask. A candidate blocked only by snp_mask/noise_mask is still fully evaluated (LR + depth) and feeds the unmasked-burden numerator, but stays filter=masked rather than PASS -- never set on a PASS record",
        ],
        "FDR": [
            1,
            "Float",
            "Local false discovery rate for this call: 1/(1+rawLR*mu), mu being this call's channel's final (FDR-refined) mutation rate estimate. Computed and reported regardless of filter, so a fail.vcf record's own FDR explains why it didn't clear its channel's threshold. 1.0 if this call's channel has no determinable mu (e.g. an indel with no HP/STR channel at all).",
        ],
    }
    dbsInfoDict = {
        "F1R2": infoDict["F1R2"],
        "F2R1": infoDict["F2R1"],
        "TAG1": infoDict["TAG1"],
        "TAG2": infoDict["TAG2"],
        "SP": infoDict["SP"],
        "TL": [1, "Integer", "Template length of the read pair"],
        "FDR": infoDict["FDR"],
    }
    formatDict = {
        "AC": [1, "Integer", "Count of alt allele"],
        "RC": [1, "Integer", "Count of ref allele"],
        "DP": [1, "Integer", "Depth at the location"],
    }
    filterDict = {
        "PASS": "All filters passed, including snp_mask/noise_mask (see SNPM/NOISEM INFO fields, both 0 here) -- a fully unmasked call",
        "masked": "Blocked by snp_mask/noise_mask only (SNPM/NOISEM INFO fields mark which); other mask types never reach this label (see the *_rescued reasons below, only emitted under --rescue). A masked candidate whose raw LR clears its channel's final refined threshold gets real depth extracted and feeds the unmasked-burden numerator, unless that real depth then fails a post-hoc sanity check (see no_good_alt_read/duplex_vaf/normal_vaf/n_cov_mask below), in which case it's relabeled to that specific reason and, like any other reject reason, only kept under --rescue. A masked candidate that never got real depth extracted (LR below threshold, or --skipCoveragePass) is dropped entirely and never appears in the fail vcf unless --rescue is on",
        "underpowered": "Passed the default calling threshold but failed the FDR-refined per-channel threshold",
        "high_nm": "Read family kept under --rescue despite failing the NM/blacklist filter",
        "low_mapq": "Read family kept under --rescue despite failing the mapq filter",
        "low_ASXS": "Read family kept under --rescue despite failing the AS-XS filter",
        "ncov_rescued": "Mutation kept under --rescue despite failing the coverage-depth mask; no coverage/depth extracted",
        "nm_rescued": "Mutation kept under --rescue despite failing the per-family NM mask at this position; no coverage/depth extracted",
        "trim_rescued": "Mutation kept under --rescue despite falling in the read-end trim zone; no coverage/depth extracted",
        "indelregion_rescued": "Indel kept under --rescue despite failing indel_mask; no coverage/depth extracted",
        "masked_rescued": "Mutation kept under --rescue despite failing a mask not covered by a more specific *_rescued reason above; no coverage/depth extracted",
        "no_good_alt_read": "Kept under --rescue despite real depth-extraction finding zero alt-supporting reads in the tumor BAM pileup (LR cleared its channel's threshold on family-consensus evidence alone) -- real AC/RC/DP attached",
        "duplex_vaf": "Kept under --rescue despite the extracted tumor allele fraction (AC/DP) exceeding --maxAF -- real AC/RC/DP attached",
        "normal_vaf": "Kept under --rescue despite the extracted matched-normal allele fraction exceeding --naf (likely germline or a systematic artifact) -- real AC/RC/DP attached",
        "n_cov_mask": "Kept under --rescue despite matched-normal depth at this position falling below --minNdepth once real depth was extracted -- real AC/RC/DP attached",
    }

    # Each mutation type gets its own SBS/INDEL/DBS subfolder (matching
    # where Estimate.py writes every downstream per-type output) rather
    # than the sample's top-level output directory. VCF writing itself
    # (PASS-only + new _fail.vcf) happens once, at the very end, after the
    # re-filter step and round 2's deferred depth-extraction have both
    # landed -- no need for the early draft-then-conditionally-rewrite
    # dance the single-pass architecture used to need.
    sbs_dir, indel_dir, dbs_dir = _ensure_type_subdirs(output_dir)

    # muts_by_duplex_group/all_keys/duplex_family_strand_composition.txt/
    # the composition heatmap/non_zero_keys all depend on duplex_read_num*
    # (real coverage data) -- moved to after round 2, below.
    # duplex_read_num_trinuc[key] is a raw (64, 4) matrix: 64 trinuc contexts
    # (un-folded by reverse complement) x 4 alt bases, with the context's own
    # reference-base column always 0 (see call.py's L-table self-skip).
    # Select the 192 valid (context, alt) cells, in the same order the
    # labels are built in, rather than flattening all 256.
    trinuc2num_64, num2trinuc_64 = build_trinuc64_order()
    label2num_192, index_192 = build_trinuc192_labels(num2trinuc_64)
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}

    # load_error_matrices(params) runs here so the initial, pre-refinement
    # per-channel Eeff estimate can use the exact same fresh-simulation
    # machinery as every later refinement iteration, instead of round-1
    # coverage (which round 1 no longer computes -- see call.py's
    # coverage_only gating). It sets ampmat/dmgmat_top/etc/trinuc_convert
    # on params; all_quals_fdr feeds simulate_power_grid (each refinement
    # worker builds its own rng, see init_refine_worker).
    load_error_matrices(params)
    all_quals_fdr = []
    reads_sampled_fdr = 0
    for read in tBam.fetch():
        if not read.is_unmapped and read.query_alignment_qualities is not None:
            quals = np.array(read.query_alignment_qualities, dtype=float)
            all_quals_fdr.extend(quals[quals > params["minBq"]].tolist())
            reads_sampled_fdr += 1
            if reads_sampled_fdr >= 1000:
                break
    if not all_quals_fdr:
        all_quals_fdr = [30]
    all_quals_fdr = np.array(all_quals_fdr, dtype=float)

    def _indel_ctx_threshold(hps_c, strs_c, id_len, pool=None):
        # Uniform default cutoff, capped per-context by the theoretical
        # ceiling of masked LR (same formula real calling now uses -- see
        # indelContextThreshold in call.py / indelMaxLR). For id_len +-1
        # with no STR annotation, Pdmg is also ref_allele-dependent (see
        # indelErrorProbs); pool ("C"=C/G, "T"=A/T, matching
        # classify_indel_channel's own base pooling) restricts this to
        # the 2 bases real calls in that specific HP channel span, so the
        # min across just those 2 (not all 4) is this channel's
        # representative (most restrictive) threshold -- matching this
        # same min-across-bases approximation _channel_eeff_at_threshold's
        # HP-channel Eeff also uses (funcs/call.py). hps_c==0 (never a
        # real hp length, which always starts at 1) is a sentinel for the
        # mismatched-insertion background channel (str.txt row 0) -- no
        # base axis there, so just force the mismatch branch once instead
        # of looping bases.
        pcutoffi = params["pcutoffi"]

        def _maxlr(ref_allele, inserted_base):
            # inserted_base=ref_allele forces the genuine same-base
            # homopolymer-matching scenario (mirrors
            # _channel_eeff_at_threshold's "hp" kind and the L_indel_1bp
            # coverage grid); a distinct inserted_base forces the row-0
            # mismatch scenario instead.
            _, _, Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c = indelErrorProbs(
                hps_c,
                strs_c,
                id_len,
                ref_allele,
                inserted_base,
                params["ampmat_hp"],
                params["dmgmat_hp"],
                params["ampmat_str"],
                params["dmgmat_str"],
            )
            return float(indelMaxLR(Pdmg_c, Pdmg_rev_c, Pdmg_bot_c, Pdmg_rev_bot_c))

        if hps_c == 0 and strs_c == 0 and id_len == 1:
            maxLR_ctx = _maxlr(0, 1)
        elif id_len in (1, -1) and strs_c == 0:
            bases = _HP_POOL_BASES[pool] if pool is not None else range(4)
            maxLR_ctx = min(_maxlr(base, base) for base in bases)
        else:
            maxLR_ctx = _maxlr(0, 0)
        return min(pcutoffi, maxLR_ctx)

    # _channel_eeff_at_threshold (funcs/call.py) is a single, picklable,
    # top-level dispatch on channel kind, callable from inside the
    # per-channel refinement Pool below -- closures over depth_by_trinuc/
    # depth_by_hpstr/all_quals_fdr/rng_fdr/params can't be pickled across
    # a Pool, so each per-channel Eeff/refine step can run in parallel here.

    # trinuc_by_duplex_group.txt (row_idx_192/col_idx_192/duplex_read_num_
    # trinuc-derived) also moved to after round 2, alongside the other
    # duplex-group report files.

    # SBS-96 mutation rate per channel (mutnum/Eeff/mutation_rate_mle),
    # written out below alongside indel_rate_by_hp_str once sbs_results is
    # available -- same columns/source (the per-channel brentq/pseudocount
    # solve in refine_channel_task) as the indel table, for consistency.
    startTime_refine = time.time()
    raw_lr_192_n1 = [[] for _ in range(192)]
    for mut in mutsAll:
        if mut.get("filter", "PASS") != "PASS":
            continue
        if min(int(mut["infos"]["F1R2"]), int(mut["infos"]["F2R1"])) < 1:
            continue
        idx_192 = label2num_192[f"{mut['infos']['TN']}>{mut['alt']}"]
        raw_lr_192_n1[idx_192].append(10 ** mut["infos"]["LR"])
    # Same raw-192 -> canonical-96 reverse-complement pairing as
    # combine_raw192_to_sbs96 (misc.py), applied here to lists of raw LR
    # values instead of scalar counts, so each of the 96 canonical classes
    # gets the pooled raw LR list of both RC-paired raw-192 classes. Eeff/
    # mutation_rate_mle per channel are filled in below, after the
    # refinement pool dispatch -- building the (name, kind, ctx_key,
    # raw_lr, threshold0, fdr) job spec here is all that's cheap enough to
    # do inline.
    _, sbs96_labels_for_rc = _build_sbs96_labels()
    raw_lr_96_n1 = []
    sbs_jobs = []
    for k, label in enumerate(sbs96_labels_for_rc):
        trinuc = label[0] + label[2] + label[6]
        alt = label[4]
        rc_trinuc = "".join(_REVCOMP[b] for b in reversed(trinuc))
        rc_alt = _REVCOMP[alt]
        raw_lr_96_n1.append(
            raw_lr_192_n1[label2num_192[f"{trinuc}>{alt}"]]
            + raw_lr_192_n1[label2num_192[f"{rc_trinuc}>{rc_alt}"]]
        )
        t_fwd, b_fwd = trinuc2num_64[trinuc], base2num[alt]
        t_rc = trinuc2num_64[rc_trinuc]
        sbs_jobs.append(
            (
                f"SBS96:{label}",
                "sbs96",
                (t_fwd, b_fwd, t_rc),
                raw_lr_96_n1[k],
                params["pcutoff"],
                args.fdrThreshold,
                args.pseudocount,
            )
        )

    # HP/STR indel mutation rate, analogous to the SBS-96 block above but
    # for indels and unrestricted by strand-depth stratum (Eeff below is
    # already an exact detection-power-weighted opportunity count, not an
    # n1-only approximation). rawLR list per (hp length 1-10 capped, indel
    # length) and (STR bin 1-3, indel length) context comes from PASS
    # indel calls; Eeff (effective opportunity) for the same context is
    # computed directly from depth_by_hpstr (read-family depth-composition
    # counts, base-aware for hp -- 40 buckets: base A/T/C/G x length 1-10,
    # then 3 STR buckets -- see funcs/call.py's _accumulate_depth_matrix)
    # weighted by the matching L_indel_1bp/L_indel_len detection-power
    # table cell, exactly the same power tables genotypeDSIndel itself
    # uses to classify these calls (real hp length, not the fixed hps=1
    # context indel100/cov_mat_indel fall back to for length>=2 opportunity).
    # Summing a context's Eeff/rawLR over all 4 hp bases is what "pools
    # reverse complementary bases" (A<->T, C<->G) amounts to here: the
    # output grid has no base axis, so forward and RC-partner base
    # contributions are always combined, not just for indel length +1.
    # HP channels only ever exist for id_len in {-1, 1} now -- there's no
    # hp.txt column for multi-bp lengths, so a length>=2 indel is always
    # STR-context (4 real bins, 1="2-9" through 4="40+"). Bin 0 ("not a
    # real repeat") gets its own channel too, but only for id_len=1 --
    # only a mismatched insertion can land there (a deletion always
    # matches its own homopolymer, so it's never a row-0 event); its
    # per-candidate routing uses the "HM" VCF field (0/1, set by
    # genotypeDSIndel/funcs/call.py at call time) rather than
    # re-deriving the match check here, since Caller.py only has the
    # final mutation record, not the reference array indelErrorProbs
    # compared against.
    hp_indel_lengths = [-1, 1]
    str_indel_lengths = [-5, -4, -3, -2, 2, 3, 4, 5]
    # Pool an HP channel's 1bp events by ref/alt base -- "C" for C/G,
    # "T" for A/T -- the same two-way split classify_indel_channel/ID83
    # labeling already uses for HP indels (Cinshp/Tinshp vs Cdelhp/
    # Tdelhp), matching the fact that the underlying learned error
    # matrices (ampmat_hp/dmgmat_hp) are indexed per-base (ref_allele*3+
    # ...) rather than pooled.
    _HP_POOL_BASES = {"C": (2, 3), "T": (0, 1)}  # base2num order A,T,C,G

    def _hp_pool(base_char):
        return "C" if base_char in ("C", "G") else "T"

    def _cap_id_len(id_len):
        # STR-length channels only exist for |id_len|<=5 (matching
        # ampmat_str/dmgmat_str's 11 columns, idLen+5 for idLen -5..5 --
        # funcs/learn.py clamps to this same range when building those
        # matrices, so an indel longer than 5bp was already folded into
        # the id_len=+-5 column at learning time).
        if id_len > 5:
            return 5
        if id_len < -5:
            return -5
        return id_len

    raw_lr_hp = {
        (hp, idl, pool): []
        for hp in range(1, 11)
        for idl in hp_indel_lengths
        for pool in ("C", "T")
    }
    raw_lr_str = {(sb, idl): [] for sb in range(0, 5) for idl in str_indel_lengths}
    raw_lr_str0 = []
    for indel in indelsAll:
        if indel.get("filter", "PASS") != "PASS":
            continue
        id_len = len(indel["alt"]) - len(indel["ref"])
        if id_len == 0:
            continue
        hp_len = int(indel["infos"]["HP"])
        str_bin = int(indel["infos"]["STR"])
        raw_lr = 10 ** indel["infos"]["LR"]
        # Matches genotypeDSIndel/indelErrorProbs's own branch selection
        # exactly: the STR branch applies at |id_len|>=2 regardless of
        # str_bin (0 is "not a real repeat," same as row 0 elsewhere in
        # str.txt -- see learn.py's str_alt_count/str_dmg_count routing,
        # which accumulates these events the same way); a +-1bp event
        # routes to the HP channel unless it's a mismatched insertion
        # (HM==0), which goes to the row-0 channel instead.
        if abs(id_len) >= 2:
            raw_lr_str[(str_bin, _cap_id_len(id_len))].append(raw_lr)
        elif abs(id_len) == 1 and indel["infos"].get("HM", 1) == 0:
            raw_lr_str0.append(raw_lr)
        elif hp_len >= 1:
            # id_len==-1 (deletion): ref is anchor+deleted-base, so the
            # deleted base is ref[1]. id_len==1 reaching this branch
            # (HM!=0) means the inserted base matches the flanking run,
            # so alt[1] (the inserted base) is that same run's base.
            base_char = indel["ref"][1] if id_len == -1 else indel["alt"][1]
            pool = _hp_pool(base_char)
            raw_lr_hp[(hp_len, id_len, pool)].append(raw_lr)

    # threshold0 is cheap (indelErrorProbs + indelMaxLR, no Monte Carlo) --
    # computed here in the main process, same as always. Only the Eeff
    # simulation (simulate_power_grid) and the refinement loop itself move
    # into the pool below.
    indel_jobs = [
        (
            "STR0_len1",
            "str",
            (0, 1),
            raw_lr_str0,
            _indel_ctx_threshold(0, 0, 1),
            args.fdrThreshold,
            args.pseudocount,
        )
    ]
    for hp_len in range(1, 11):
        for id_len in hp_indel_lengths:
            for pool in ("C", "T"):
                threshold0 = _indel_ctx_threshold(hp_len, 0, id_len, pool)
                indel_jobs.append(
                    (
                        f"HP{hp_len}_len{id_len}_{pool}",
                        "hp",
                        (hp_len, id_len, pool),
                        raw_lr_hp[(hp_len, id_len, pool)],
                        threshold0,
                        args.fdrThreshold,
                        args.pseudocount,
                    )
                )
    for str_bin in range(0, 5):
        for id_len in str_indel_lengths:
            threshold0 = _indel_ctx_threshold(1, str_bin, id_len)
            indel_jobs.append(
                (
                    f"STR{str_bin}_len{id_len}",
                    "str",
                    (str_bin, id_len),
                    raw_lr_str[(str_bin, id_len)],
                    threshold0,
                    args.fdrThreshold,
                    args.pseudocount,
                )
            )

    # Per-channel FDR control: for every channel above (96 SBS96 classes +
    # 90 HP combos + 24 STR combos, ~210 total), mu0 (this channel's MLE
    # mixture weight at its default threshold0, from a direct brentq solve
    # regularized by args.pseudocount) is taken as the channel's mutation
    # rate outright, with no iterative re-simulation/re-MLE step-up. The LR
    # threshold whose local fdr 1/(1+LR*mu0) equals args.fdrThreshold is
    # solved for directly and used unconditionally as this channel's round
    # 2 calling threshold -- see funcs/call.py's _refine_channel. The
    # pseudocount regularization guarantees mu0 has a real root in (0, 1)
    # for every channel, so there is no exclusion case to handle here.
    #
    # Each channel's Eeff0/mu0/threshold solve is independent of every
    # other channel's, and each Eeff evaluation is a Monte Carlo
    # simulate_power_grid call -- the dominant cost of this phase. One task
    # per channel is dispatched across a worker Pool
    # (refine_channel_task/init_refine_worker, funcs/call.py), capped at
    # args.threads to match what -p promises the user (same as round 1/
    # round 2's pools -- all four Pool() calls in this file are capped at
    # args.threads rather than left to default to os.cpu_count()).
    refine_pool = Pool(
        args.threads,
        initializer=init_refine_worker,
        initargs=(
            depth_by_trinuc,
            depth_by_hpstr,
            all_quals_fdr,
            params["ampmat"],
            params["ampmat_rev"],
            params["dmgmat_top"],
            params["dmgmat_rev_top"],
            params["dmgmat_bot"],
            params["dmgmat_rev_bot"],
            params["trinuc_convert"],
            params["ampmat_hp"],
            params["dmgmat_hp"],
            params["ampmat_str"],
            params["dmgmat_str"],
            params["seed"],
        ),
    )
    refine_results = refine_pool.map(refine_channel_task, sbs_jobs + indel_jobs)
    refine_pool.close()
    refine_pool.terminate()
    refine_pool.join()
    sbs_results = refine_results[: len(sbs_jobs)]
    indel_results = refine_results[len(sbs_jobs) :]

    channel_thresholds = {}
    # Every channel's mu0 (its round-1 MLE mixture weight). Used below to
    # compute each PASS call's own local FDR (1/(1+rawLR*mu0_channel)) for
    # the per-type "Total FDR" stats.txt lines.
    channel_final_mu = {}
    # Every channel gets a real new_threshold from _refine_channel's direct
    # FDR-at-mu0 solve; round 2 re-simulates calling from real per-position
    # data at that threshold.

    sbs96_rate_rows = []
    for k, (name, kind, ctx_key, Eeff0, mu0, new_threshold) in enumerate(sbs_results):
        channel_final_mu[name] = mu0
        channel_thresholds[name] = new_threshold
        sbs96_rate_rows.append(
            {
                "context": name,
                "mutnum": len(raw_lr_96_n1[k]),
                "Eeff": Eeff0,
                "mutation_rate_mle": mu0,
            }
        )
    sbs96_rate_n1 = pd.DataFrame(sbs96_rate_rows)
    sbs96_rate_n1.to_csv(
        os.path.join(
            params["tmp_dir"],
            os.path.basename(args.output) + "_sbs96_rate_n1.txt",
        ),
        sep="\t",
        index=False,
    )

    indel_rate_rows = []
    indel_results_iter = iter(indel_results)
    # STR0_len1 (mismatched-insertion background) was prepended first in
    # indel_jobs above -- consume it here, in the same order.
    name, kind, ctx_key, Eeff0, mu0, new_threshold = next(indel_results_iter)
    indel_rate_rows.append(
        {
            "context": "STR0",
            "indel_length": 1,
            "mutnum": len(raw_lr_str0),
            "Eeff": Eeff0,
            "mutation_rate_mle": mu0,
        }
    )
    channel_final_mu[name] = mu0
    channel_thresholds[name] = new_threshold
    for hp_len in range(1, 11):
        for id_len in hp_indel_lengths:
            for pool in ("C", "T"):
                raw_lr_list = raw_lr_hp[(hp_len, id_len, pool)]
                name, kind, ctx_key, Eeff0, mu0, new_threshold = next(
                    indel_results_iter
                )
                channel_final_mu[name] = mu0
                indel_rate_rows.append(
                    {
                        "context": f"HP{hp_len}_{pool}",
                        "indel_length": id_len,
                        "mutnum": len(raw_lr_list),
                        "Eeff": Eeff0,
                        "mutation_rate_mle": mu0,
                    }
                )
                channel_thresholds[name] = new_threshold
    for str_bin in range(0, 5):
        for id_len in str_indel_lengths:
            raw_lr_list = raw_lr_str[(str_bin, id_len)]
            name, kind, ctx_key, Eeff0, mu0, new_threshold = next(indel_results_iter)
            channel_final_mu[name] = mu0
            indel_rate_rows.append(
                {
                    "context": f"STR{str_bin}",
                    "indel_length": id_len,
                    "mutnum": len(raw_lr_list),
                    "Eeff": Eeff0,
                    "mutation_rate_mle": mu0,
                }
            )
            channel_thresholds[name] = new_threshold
    indel_rate_by_hp_str = pd.DataFrame(indel_rate_rows)
    indel_rate_by_hp_str.to_csv(
        os.path.join(
            params["tmp_dir"],
            os.path.basename(args.output) + "_indel_rate_by_hp_str.txt",
        ),
        sep="\t",
        index=False,
    )
    print(
        "..............Completed per-channel FDR refinement ("
        f"{len(sbs_jobs) + len(indel_jobs)} channels) in "
        + str((time.time() - startTime_refine) / 60)
        + " minutes,filtering calls................."
    )

    # indel_by_duplex_group.txt/dbs_by_duplex_group.txt/_stats.txt's
    # coverage-derived lines/coverage.bed.gz itself all need real
    # (round-2) coverage too -- moved down, after round 2 runs.
    sample_name = os.path.basename(args.output)
    sample_dir = params["tmp_dir"]

    # Build the final per-channel threshold override dicts, covering every
    # channel _refine_channel processed (not just ones that needed
    # raising) -- every channel now gets its own FDR-at-mu0-derived
    # threshold applied unconditionally for round 2, always defined
    # (possibly empty) so round 2 below works uniformly.
    pcutoff_sbs_override = {}
    indel_1bp_threshold_override = {}
    indel_len_threshold_override = {}
    for name, new_threshold in channel_thresholds.items():
        if name.startswith("SBS96:"):
            label = name.split(":", 1)[1]
            trinuc = label[0] + label[2] + label[6]
            alt = label[4]
            rc_trinuc = "".join(_REVCOMP[b] for b in reversed(trinuc))
            rc_alt = _REVCOMP[alt]
            t_fwd, b_fwd = trinuc2num_64[trinuc], base2num[alt]
            t_rc, b_rc = trinuc2num_64[rc_trinuc], base2num[rc_alt]
            pcutoff_sbs_override[(t_fwd, b_fwd)] = new_threshold
            pcutoff_sbs_override[(t_rc, b_rc)] = new_threshold
        elif name.startswith("HP"):
            # HP jobs only ever have id_len in {-1, 1} (see
            # hp_indel_lengths above), so this is always a 1bp override.
            # Name is "HP{hp_len}_len{id_len}_{pool}" -- pool ("C"/"T")
            # is the base-pool axis.
            hp_len_str, rest = name[2:].split("_len")
            id_len_str, pool = rest.split("_")
            hp_len, id_len = int(hp_len_str), int(id_len_str)
            sign_idx = 1 if id_len == 1 else 0
            indel_1bp_threshold_override[(hp_len, sign_idx, pool)] = new_threshold
        elif name.startswith("STR"):
            # Row = str_bin directly (0-4, matching L_indel_len's row
            # axis) -- no more "9+" offset now that STR has its own
            # dedicated 5-row table instead of sharing rows with hp.txt.
            str_bin_str, id_len_str = name[3:].split("_len")
            str_bin, id_len = int(str_bin_str), int(id_len_str)
            indel_len_threshold_override[(str_bin, id_len)] = new_threshold

    # SBS/ID/DBS re-filter: a candidate's LR is a property of the observed
    # reads, not the calling threshold, so this is pure reclassification
    # against each channel's final threshold, no rescan needed.
    # _refine_channel clamps its FDR-at-mu0 solve to max(threshold0, ...),
    # so the final threshold is still guaranteed >= threshold0 -- never
    # lower -- same invariant the old iterative step-up loop preserved,
    # just for a different reason now (round 1's bam scan only ever
    # records a candidate whose LR already clears threshold0, so nothing
    # below it exists in mutsAll/indelsAll to newly promote here). A
    # record already failing the default threshold can therefore never
    # newly clear a higher one -- only currently-PASS records are eligible
    # to flip, and they flip down to "underpowered" (distinct from
    # structural "masked") rather than being relabeled "masked" outright.
    # Structurally-masked candidates (antimask-failed, never LR-gated to
    # begin with) get a separate check: if their raw LR clears this
    # channel's final threshold, they're flagged here for depth-extraction
    # in round 2 below (still reported filter="masked", just with real
    # depth instead of none, never promoted to PASS). filter=="masked" is
    # only ever set by call.py for candidates blocked solely by snp_mask/
    # noise_mask (see call.py's unmasked_pass_bool branch), so SNPM or
    # NOISEM is always set on every record reaching this check -- the
    # eligibility rule is purely LR>=threshold and (SNPM or NOISEM),
    # independent of --rescue. --rescue is a separate, earlier concern in
    # call.py: it only controls whether *_rescued-labeled records (a
    # different filter value entirely, for candidates blocked by n_cov/nm/
    # trim/indel_mask) are emitted into the output VCF at all -- those
    # never get real depth either way, regardless of rescue, per rescue's
    # own "no depth" semantics.
    def _mut_key(m):
        return (m["chrom"], m["pos"], m["ref"], m["alt"])

    deferred_depth_keys = set()
    for mut in mutsAll:
        old_filter = mut.get("filter", "PASS")
        if old_filter not in ("PASS", "masked"):
            continue
        ctx_key = (trinuc2num_64.get(mut["infos"]["TN"]), base2num.get(mut["alt"]))
        threshold = pcutoff_sbs_override.get(ctx_key, params["pcutoff"])
        if old_filter == "PASS":
            if mut["infos"]["LR"] < threshold:
                mut["filter"] = "underpowered"
        elif mut["infos"]["LR"] >= threshold and (
            mut["infos"].get("SNPM") or mut["infos"].get("NOISEM")
        ):
            deferred_depth_keys.add(_mut_key(mut))

    for indel in indelsAll:
        old_filter = indel.get("filter", "PASS")
        if old_filter not in ("PASS", "masked"):
            continue
        id_len = len(indel["alt"]) - len(indel["ref"])
        hp_len = int(indel["infos"]["HP"])
        str_bin = int(indel["infos"]["STR"])
        # Multi-bp indels are always STR-context now (row = str_bin, 0-4
        # -- 0 is "not a real repeat"), independent of hp_len. +-1bp
        # indels route by hp_len, UNLESS this specific call is a
        # mismatched insertion (HM==0, set by genotypeDSIndel at call
        # time), in which case it uses the STR0_len1 override instead --
        # same channel indelErrorProbs itself scored it against.
        if abs(id_len) >= 2:
            threshold = indel_len_threshold_override.get(
                (str_bin, _cap_id_len(id_len)), params["pcutoffi"]
            )
        elif abs(id_len) == 1 and indel["infos"].get("HM", 1) == 0:
            threshold = indel_len_threshold_override.get((0, 1), params["pcutoffi"])
        elif hp_len >= 1:
            sign_idx = 1 if id_len == 1 else 0
            base_char = indel["ref"][1] if id_len == -1 else indel["alt"][1]
            pool = _hp_pool(base_char)
            threshold = indel_1bp_threshold_override.get(
                (hp_len, sign_idx, pool), params["pcutoffi"]
            )
        else:
            threshold = params["pcutoffi"]
        if old_filter == "PASS":
            if indel["infos"]["LR"] < threshold:
                indel["filter"] = "underpowered"
        elif indel["infos"]["LR"] >= threshold and indel["infos"].get("NOISEM"):
            deferred_depth_keys.add(_mut_key(indel))

    # A DBS call requires both constituent SNVs to individually PASS (see
    # _detect_dbs_pairs in call.py); if either flipped above, the DBS call
    # needs re-validating too. DBS is never "masked" by construction (only
    # ever "PASS" or dropped entirely), so only the downgrade path applies
    # -- no deferred set for DBS.
    pass_positions = {
        (mut["chrom"], mut["pos"])
        for mut in mutsAll
        if mut.get("filter", "PASS") == "PASS"
    }
    for dbs in dbsAll:
        if dbs.get("filter", "PASS") != "PASS":
            continue
        new_pass = (dbs["chrom"], dbs["pos"]) in pass_positions and (
            dbs["chrom"],
            dbs["pos"] + 1,
        ) in pass_positions
        if not new_pass:
            dbs["filter"] = "underpowered"

    efficiency = duplex_num / rec_num if rec_num > 0 else 0.0

    # --skipCoveragePass: only round 1 calling + FDR-threshold
    # determination run. Round 2 (below) is what computes real,
    # position-weighted duplex depth/coverage -- coverage.bed.gz, the
    # duplex-family composition file/heatmap, and the three
    # *_by_duplex_group.txt files all come from it, so skipping it means
    # none of those get produced (and masked-but-now-qualifying
    # candidates from the re-filter above get no real depth either,
    # same as any other structurally-masked record). VCFs and _stats.txt
    # (minus the coverage-derived lines) are unaffected -- both are
    # sourced from round 1 + the re-filter, not round 2.
    coverage_pass_enabled = not args.skipCoveragePass
    if coverage_pass_enabled:
        # Round 2: coverage-only pass, reusing round 1's region split
        # (regions_list) directly instead of recomputing it, with every
        # channel's final threshold baked in as the actual calling threshold
        # this time -- not layered as an override on top of round 1's uniform
        # default. Builds L/L_indel_1bp/L_indel_len and real coverage.bed
        # once more, and performs the deferred depth-extraction for the
        # masked-but-now-qualifying set identified above.
        print("..............Starting coverage pass (round 2)..............")
        startTime3 = time.time()
        round2_args = []
        for nn, r2_regions in enumerate(regions_list):
            paramsNow2 = params.copy()
            paramsNow2["regions"] = r2_regions
            paramsNow2["isLearn"] = False
            paramsNow2["coverage_only"] = True
            paramsNow2["pcutoff_sbs_override"] = pcutoff_sbs_override
            paramsNow2["indel_1bp_threshold_override"] = indel_1bp_threshold_override
            paramsNow2["indel_len_threshold_override"] = indel_len_threshold_override
            paramsNow2["deferred_depth_keys"] = deferred_depth_keys
            round2_args.append((paramsNow2, nn))
        round2_pool = Pool(args.threads)
        round2_results = round2_pool.starmap(callBam, round2_args)
        round2_pool.close()
        round2_pool.terminate()
        round2_pool.join()
        print(
            "..............Completed coverage pass in "
            + str((time.time() - startTime3) / 60)
            + " minutes,merging results................."
        )

        # Merge round 2's deferred depth-extraction (a small subset -- only
        # masked candidates whose LR now clears their channel's final
        # threshold) back into round 1's mutsAll/indelsAll, by position key.
        # Carries filter along with samples: call.py relabels a deferred
        # candidate to a reject-reason filter (no_good_alt_read/duplex_vaf/
        # normal_vaf/n_cov_mask) when real depth fails a post-hoc sanity
        # check, or leaves it "masked" (with real depth) when it clears
        # everything -- either way the region-local round-2 record is the
        # source of truth, since round 1 never saw real depth for these.
        deferred_mut_updates = {
            _mut_key(m): (m["samples"], m["filter"])
            for r in round2_results
            for m in r[0]
            if _mut_key(m) in deferred_depth_keys
        }
        deferred_indel_updates = {
            _mut_key(m): (m["samples"], m["filter"])
            for r in round2_results
            for m in r[6]
            if _mut_key(m) in deferred_depth_keys
        }
        for mut in mutsAll:
            key = _mut_key(mut)
            if key in deferred_mut_updates:
                mut["samples"], mut["filter"] = deferred_mut_updates[key]
        for indel in indelsAll:
            key = _mut_key(indel)
            if key in deferred_indel_updates:
                indel["samples"], indel["filter"] = deferred_indel_updates[key]

        # Everything below that needs real position-weighted coverage now
        # sources it from round2_results (final thresholds) instead of round
        # 1's results -- mirrors exactly how round 1's own results were merged
        # near the top of this function.
        coverages = [_[1] for _ in round2_results]
        duplex_read_nums = [_[4] for _ in round2_results]
        duplex_read_nums_trinuc = [_[5] for _ in round2_results]
        unmasked_coverages = [_[11] for _ in round2_results]
        coverages_indels_cat = [_[12] for _ in round2_results]
        unmasked_coverages_indels_cat = [_[13] for _ in round2_results]
        duplex_read_nums_indel = [_[14] for _ in round2_results]
        duplex_read_nums_dbs = [_[16] for _ in round2_results]

        coverage = sum(coverages)
        coverage_indel_cat = sum(coverages_indels_cat)
        unmasked_coverage = sum(unmasked_coverages)
        unmasked_coverage_indel_cat = sum(unmasked_coverages_indels_cat)

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
                num: sum(
                    [d.get(num, np.zeros((64, 4))) for d in duplex_read_nums_trinuc]
                )
                for num in duplex_combinations
            }
        )
        duplex_read_num_indel = OrderedDict(
            {
                num: sum([d.get(num, np.zeros(100)) for d in duplex_read_nums_indel])
                for num in duplex_combinations
            }
        )
        duplex_read_num_dbs = OrderedDict(
            {
                num: sum([d.get(num, np.zeros(144)) for d in duplex_read_nums_dbs])
                for num in duplex_combinations
            }
        )

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

        max_val = max(
            (max(int(k.split("+")[0]), int(k.split("+")[1])) for k in all_keys),
            default=0,
        )
        total_read_sets = sum(total_family_num.values())
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
        row_idx_192 = np.array(
            [trinuc2num_64[label.split(">")[0]] for label in index_192]
        )
        col_idx_192 = np.array([base2num[label.split(">")[1]] for label in index_192])
        duplex_read_num_trinuc = {_: duplex_read_num_trinuc[_] for _ in non_zero_keys}
        data_192 = {
            k: v[row_idx_192, col_idx_192] for k, v in duplex_read_num_trinuc.items()
        }
        trinuc_by_duplex_group = pd.DataFrame(data_192, index=index_192)
        trinuc_by_duplex_group.to_csv(
            os.path.join(sbs_dir, sample_name + "_trinuc_by_duplex_group.txt"),
            sep="\t",
            index=True,
        )

        _, index_indel100 = build_indel100_labels()
        duplex_read_num_indel = {_: duplex_read_num_indel[_] for _ in non_zero_keys}
        data_indel100 = {k: v for k, v in duplex_read_num_indel.items()}
        indel_by_duplex_group = pd.DataFrame(data_indel100, index=index_indel100)
        indel_by_duplex_group.to_csv(
            os.path.join(indel_dir, sample_name + "_indel_by_duplex_group.txt"),
            sep="\t",
            index=True,
        )

        _, index_dbs144 = build_dbs_raw144_labels()
        duplex_read_num_dbs = {_: duplex_read_num_dbs[_] for _ in non_zero_keys}
        data_dbs144 = {k: v for k, v in duplex_read_num_dbs.items()}
        dbs_by_duplex_group = pd.DataFrame(data_dbs144, index=index_dbs144)
        dbs_by_duplex_group.to_csv(
            os.path.join(dbs_dir, sample_name + "_dbs_by_duplex_group.txt"),
            sep="\t",
            index=True,
        )
    else:
        print(
            "..............Skipping coverage pass (round 2) -- duplex depth, "
            "per-locus coverage.bed.gz, duplex-family composition, and "
            "by-duplex-group output files will not be generated "
            "(--skipCoveragePass).............."
        )

    # Records tagged with one of these reject-reason filters cleared LR
    # (round 1 PASS, or a round-2 masked-deferred candidate) but real
    # depth-extraction then failed a post-hoc sanity check (see call.py's
    # SNV/indel filtering loops). They're only worth reporting under
    # --rescue's "junk mutations that may carry real biological signal"
    # opt-in; otherwise drop them here so they never reach either VCF, same
    # as any other non-rescued reject. This is the single, centralized
    # keep/drop decision for these four labels -- call.py itself never
    # drops them, since a round-2 record is only a region-local copy that
    # feeds the merge above, not the final mutsAll/indelsAll directly.
    if not params["rescue"]:
        _reject_reasons = {"no_good_alt_read", "duplex_vaf", "normal_vaf", "n_cov_mask"}
        mutsAll = [m for m in mutsAll if m.get("filter", "PASS") not in _reject_reasons]
        indelsAll = [
            m for m in indelsAll if m.get("filter", "PASS") not in _reject_reasons
        ]

    # "masked" (SNPM/NOISEM-blocked) records that never went through
    # round 2's deferred depth-extraction -- raw LR below the channel's
    # final refined threshold, or --skipCoveragePass -- still carry
    # round 1's zero-depth placeholder in "samples". A masked record that
    # *did* get real depth can never land back here with an all-zero
    # tumor triple: call.py's eligible branch relabels ta==0 to
    # "no_good_alt_read" instead of leaving it "masked" (see its
    # eligible-branch comment). So filter=="masked" + all-zero samples
    # unambiguously means "never depth-extracted", and those shouldn't
    # reach the fail vcf without --rescue -- same opt-in as every other
    # no-real-depth reason above.
    if not params["rescue"]:

        def _masked_no_depth(m):
            return m.get("filter") == "masked" and m.get("samples") == [
                [0, 0, 0],
                [0, 0, 0],
            ]

        mutsAll = [m for m in mutsAll if not _masked_no_depth(m)]
        indelsAll = [m for m in indelsAll if not _masked_no_depth(m)]

    # Per-call FDR: the same local-fdr formula every channel already uses
    # internally (1/(1+rawLR*mu0), mu0 being that call's own channel's
    # round-1 MLE mutation rate -- see _refine_channel). Computed and
    # stamped onto every SBS/indel record's own "FDR" INFO field (not
    # just PASS -- a fail.vcf record's own FDR is exactly what explains
    # why it didn't clear its channel), and separately averaged over just
    # the PASS population per type for the "Total X FDR" stats.txt lines
    # below. channel_final_mu has an entry (mu0) for every channel
    # _refine_channel processed, so this is a straight per-call lookup, no
    # extra Monte Carlo. A DBS call's own FDR isn't
    # independently modeled -- it's derived from its two constituent SNVs
    # (call.py's _detect_dbs_pairs) -- so its local fdr is the probability
    # at least one of those two is wrong, 1-(1-fdr_pos1)*(1-fdr_pos2),
    # reusing the per-position fdrs already computed for SBS just below.
    # Records whose channel has no determinable mu (indels with no HP/STR
    # channel at all, hp_len==0 and not STR) get FDR=1.0 -- "not modeled",
    # not a claim of zero confidence.
    def _sbs_channel_name(trinuc, alt):
        ref = trinuc[1]
        if ref in ("C", "T"):
            return f"SBS96:{trinuc[0]}[{ref}>{alt}]{trinuc[2]}"
        rc_trinuc = "".join(_REVCOMP[b] for b in reversed(trinuc))
        rc_alt = _REVCOMP[alt]
        return f"SBS96:{rc_trinuc[0]}[{rc_trinuc[1]}>{rc_alt}]{rc_trinuc[2]}"

    sbs_local_fdrs = []
    pos_to_local_fdr = {}
    for mut in mutsAll:
        name = _sbs_channel_name(mut["infos"]["TN"], mut["alt"])
        mu = channel_final_mu.get(name)
        local_fdr = (
            1.0 / (1.0 + (10 ** mut["infos"]["LR"]) * mu) if mu is not None else 1.0
        )
        mut["infos"]["FDR"] = local_fdr
        pos_to_local_fdr[(mut["chrom"], mut["pos"])] = local_fdr
        if mu is not None and mut.get("filter", "PASS") == "PASS":
            sbs_local_fdrs.append(local_fdr)
    total_sbs_fdr = float(np.mean(sbs_local_fdrs)) if sbs_local_fdrs else float("nan")

    indel_local_fdrs = []
    for indel in indelsAll:
        id_len = len(indel["alt"]) - len(indel["ref"])
        hp_len = int(indel["infos"]["HP"])
        str_bin = int(indel["infos"]["STR"])
        if abs(id_len) >= 2:
            name = f"STR{str_bin}_len{_cap_id_len(id_len)}"
        elif abs(id_len) == 1 and indel["infos"].get("HM", 1) == 0:
            name = "STR0_len1"
        elif hp_len >= 1:
            base_char = indel["ref"][1] if id_len == -1 else indel["alt"][1]
            name = f"HP{hp_len}_len{id_len}_{_hp_pool(base_char)}"
        else:
            name = None
        mu = channel_final_mu.get(name)
        local_fdr = (
            1.0 / (1.0 + (10 ** indel["infos"]["LR"]) * mu) if mu is not None else 1.0
        )
        indel["infos"]["FDR"] = local_fdr
        if mu is not None and indel.get("filter", "PASS") == "PASS":
            indel_local_fdrs.append(local_fdr)
    total_indel_fdr = (
        float(np.mean(indel_local_fdrs)) if indel_local_fdrs else float("nan")
    )

    dbs_local_fdrs = []
    for dbs in dbsAll:
        fdr1 = pos_to_local_fdr.get((dbs["chrom"], dbs["pos"]))
        fdr2 = pos_to_local_fdr.get((dbs["chrom"], dbs["pos"] + 1))
        local_fdr = (
            1.0 - (1.0 - fdr1) * (1.0 - fdr2)
            if fdr1 is not None and fdr2 is not None
            else 1.0
        )
        dbs["infos"]["FDR"] = local_fdr
        if (
            fdr1 is not None
            and fdr2 is not None
            and dbs.get("filter", "PASS") == "PASS"
        ):
            dbs_local_fdrs.append(local_fdr)
    total_dbs_fdr = float(np.mean(dbs_local_fdrs)) if dbs_local_fdrs else float("nan")

    with open(args.output + "_stats.txt", "w") as f:
        f.write(f"Number of Read Families\t{unique_read_num}\n")
        f.write(f"Number of Pass-filter Reads\t{pass_read_num}\n")
        f.write(f"Number of Effective Read Families\t{duplex_num}\n")
        if coverage_pass_enabled:
            f.write(f"Effective Coverage\t{coverage}\n")
            f.write(f"Unmasked Coverage\t{unmasked_coverage}\n")
            # Single aggregate total, not broken out per
            # INDEL_COVERAGE_CATEGORY_LABELS category -- every consumer
            # (Estimate.py's unmasked_indel_cov) summed the per-category
            # breakdown back into one number anyway, and the per-category
            # "Effective"/"Unmasked Indel Coverage (label)" lines this used
            # to write alongside it were never read by anything.
            f.write(f"Unmasked Indel Coverage\t{unmasked_coverage_indel_cat.sum()}\n")
            f.write(
                f"Per Read Family Coverage \t{coverage/duplex_num if duplex_num > 0 else 0.0}\n"
            )
        f.write(
            f"Pass-filter Duplication Rate\t\
        {1-unique_read_num/pass_read_num if pass_read_num > 0 else 0.0}\n"
        )
        f.write(f"Efficiency\t{efficiency}\n")
        f.write(f"Total SBS FDR\t{total_sbs_fdr}\n")
        f.write(f"Total Indel FDR\t{total_indel_fdr}\n")
        f.write(f"Total DBS FDR\t{total_dbs_fdr}\n")

    if coverage_pass_enabled:
        # Merge and combine round 2's per-worker coverage.bed.gz files (same
        # naming/numbering as round 1's, since round 2 reused regions_list --
        # round 2's real data simply overwrote round 1's empty placeholders at
        # the same paths) into the final args.output + "_coverage.bed.gz".
        # Round 2 is always a fresh whole-region write, not a partial-region
        # patch.
        # Use the actual number of dispatched chunks, not args.threads --
        # regionSequence can legitimately run out before -p chunks are
        # filled (empty/duplicate-boundary chunks dropped upstream), in
        # which case fewer real *_<n>_coverage*.tmp.bed.gz files exist
        # than -p would imply.
        merge_and_combine_coverage_files(sample_name, sample_dir, len(regions_list))
        os.replace(
            os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz"),
            args.output + "_coverage.bed.gz",
        )
        os.replace(
            os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz.tbi"),
            args.output + "_coverage.bed.gz.tbi",
        )

    # Final VCF write, once: PASS-only in the main per-type files,
    # everything else (masked/underpowered/any rescue-reason label) in a
    # parallel _fail.vcf per type.
    pass_mutsAll = [m for m in mutsAll if m.get("filter", "PASS") == "PASS"]
    fail_mutsAll = [m for m in mutsAll if m.get("filter", "PASS") != "PASS"]
    pass_indelsAll = [m for m in indelsAll if m.get("filter", "PASS") == "PASS"]
    fail_indelsAll = [m for m in indelsAll if m.get("filter", "PASS") != "PASS"]
    pass_dbsAll = [m for m in dbsAll if m.get("filter", "PASS") == "PASS"]
    fail_dbsAll = [m for m in dbsAll if m.get("filter", "PASS") != "PASS"]

    vcfLines = createVcfStrings(
        chromDict, infoDict, formatDict, filterDict, pass_mutsAll
    )
    with open(os.path.join(sbs_dir, sample_name + "_sbs.vcf"), "w") as vcf:
        vcf.write(vcfLines)
    vcfLines = createVcfStrings(
        chromDict, infoDict, formatDict, filterDict, fail_mutsAll
    )
    with open(os.path.join(sbs_dir, sample_name + "_sbs_fail.vcf"), "w") as vcf:
        vcf.write(vcfLines)

    vcfLines = createVcfStrings(
        chromDict, infoDict, formatDict, filterDict, pass_indelsAll
    )
    with open(os.path.join(indel_dir, sample_name + "_indel.vcf"), "w") as vcf:
        vcf.write(vcfLines)
    vcfLines = createVcfStrings(
        chromDict, infoDict, formatDict, filterDict, fail_indelsAll
    )
    with open(os.path.join(indel_dir, sample_name + "_indel_fail.vcf"), "w") as vcf:
        vcf.write(vcfLines)

    vcfLines = createVcfStrings(
        chromDict, dbsInfoDict, formatDict, filterDict, pass_dbsAll
    )
    with open(os.path.join(dbs_dir, sample_name + "_dbs.vcf"), "w") as vcf:
        vcf.write(vcfLines)
    vcfLines = createVcfStrings(
        chromDict, dbsInfoDict, formatDict, filterDict, fail_dbsAll
    )
    with open(os.path.join(dbs_dir, sample_name + "_dbs_fail.vcf"), "w") as vcf:
        vcf.write(vcfLines)

    print(
        "..............Completed variant calling "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )


def merge_and_combine_coverage_files(sample_name, sample_dir, nprocess):
    """
    1. Merge adjacent next_region and prev_region files by summing all numeric columns
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
    Merge two adjacent bed files by summing every numeric column (everything
    after chrom/start/end) for matching positions. Column count is read from
    the data rather than hardcoded, so this stays correct as the coverage
    bed schema grows (currently: cov_A/T/C/G and 6 indel category coverage
    columns).
    """
    coverage_dict = {}

    def accumulate(path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Expected file for merging not found: {path}")
        with bgzf.open(path, "rt") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                chrom, start, end = parts[:3]
                values = [float(v) for v in parts[3:]]
                key = (chrom, start, end)
                if key in coverage_dict:
                    coverage_dict[key] = [
                        existing + new
                        for existing, new in zip(coverage_dict[key], values)
                    ]
                else:
                    coverage_dict[key] = values

    accumulate(next_file)
    accumulate(prev_file)

    # Write merged result only if there's data
    with bgzf.open(output_file, "wt") as f:
        for (chrom, start, end), values in sorted(coverage_dict.items()):
            f.write("\t".join([chrom, start, end] + [str(v) for v in values]) + "\n")


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
