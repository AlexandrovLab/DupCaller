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

from matplotlib import pyplot as plt

from . import __version__
from .funcs.call import callBam, contract_worker_tally
from .funcs.learn import NUM_BQ, estimate_sbs_srd_rates
from .funcs.misc import simulate_power_grid, load_error_matrices
from .funcs.misc import init_refine_worker, refine_channel_task, mu_at_min_lr
from .funcs.misc import indel_mu_channel, str_channel_name, parse_str_channel_name
from .funcs.misc import str_profile_is_current
from .funcs.indels import STR_COL_LABELS
from .funcs.prob import indelErrorProbs, indelMaxLR, indel_strand_pvalue, strand_pvalue
from .funcs.misc import createVcfStrings
from .funcs.misc import splitBamRegions
from .funcs.misc import drop_empty_regions
from .funcs.misc import previous_region_coverage_start
from .funcs.misc import check_mate_cigar_tags
from .funcs.misc import getAlignmentObject as BAM
from .funcs.misc import check_h5_usable
from .funcs.misc import build_trinuc64_order, build_trinuc192_labels
from .funcs.misc import _build_sbs96_labels
from .funcs.misc import _REVCOMP
from .funcs.misc import build_indel100_labels
from .funcs.misc import build_dbs_raw144_labels
from .funcs.misc import _ensure_type_subdirs
from .funcs.misc import load_repeat_context
from .funcs.misc import parse_barcode_option
from pysam import TabixFile as BED
import pysam


# Written to every row of _indel_rate_by_hp_str.txt ("eeff_sites" column) to
# identify the HP/STR Eeff site set, so -mr refuses a table whose mu is on a
# different scale: one site per homopolymer run / STR tract, deletion and
# insertion sites counted separately, STR channels per repeat unit
# (STR{bin}_u{U}_{del|ins}).
INDEL_EEFF_SITE_SET = "per_unit_signed"


def _load_mutation_rate_override(prefix):
    """Load {prefix}_sbs96_rate_n1.txt and {prefix}_indel_rate_by_hp_str.txt
    (as written by do_call's refine step) into a {channel_name: mu0} dict
    keyed by do_call's channel job names, for -mr/--muterateprefix: a
    small-region call can reuse a genome-wide run's mu0 estimates.

    SBS "context" is already the job name; indel names are rebuilt from
    "context"/"indel_length" ("HP{hp_len}_len{id_len}_{pool}", or
    str_channel_name for STR).
    """
    sbs_table = pd.read_csv(prefix + "_sbs96_rate_n1.txt", sep="\t")
    indel_table = pd.read_csv(prefix + "_indel_rate_by_hp_str.txt", sep="\t")
    if (
        "eeff_sites" not in indel_table.columns
        or (indel_table["eeff_sites"] != INDEL_EEFF_SITE_SET).any()
    ):
        raise ValueError(
            f"-mr/--muterateprefix table {prefix}_indel_rate_by_hp_str.txt was "
            "written by a DupCaller version whose HP/STR Eeff counted a different "
            f"site set (expected eeff_sites={INDEL_EEFF_SITE_SET!r}); its indel "
            "mutation rates are on a different scale. Regenerate it from a full "
            "run with the current DupCaller version."
        )
    override = dict(zip(sbs_table["context"], sbs_table["mutation_rate_mle"]))
    for context, indel_length, mu0 in zip(
        indel_table["context"],
        indel_table["indel_length"],
        indel_table["mutation_rate_mle"],
    ):
        if context.startswith("HP"):
            hp_len_str, pool = context[2:].rsplit("_", 1)
            name = f"HP{hp_len_str}_len{indel_length}_{pool}"
        elif context.startswith("STR"):
            name = str_channel_name(int(context[3:].split("_u")[0]), int(indel_length))
        else:
            raise ValueError(
                f"unrecognized indel mutation-rate-table context {context!r} "
                f"in {prefix}_indel_rate_by_hp_str.txt"
            )
        override[name] = mu0
    return override


def _sbs_family_key(mut):
    """(chrom, pos, TAG1, TAG2, SP, TL): one read family's call at one
    position. TL is part of the family label because a fragment no longer
    than the read makes two families with the same barcodes and start."""
    infos = mut["infos"]
    return (
        mut["chrom"],
        mut["pos"],
        infos["TAG1"],
        infos["TAG2"],
        infos["SP"],
        infos["TL"],
    )


def _dbs_family_keys(dbs):
    """_sbs_family_key of a DBS record's two constituent SBS (same family)."""
    return [
        _sbs_family_key(
            {"chrom": dbs["chrom"], "pos": dbs["pos"] + k, "infos": dbs["infos"]}
        )
        for k in (0, 1)
    ]


def _relabel_dbs_members(sbs_records, dbsAll):
    """Set FILTER dbs_member on each PASS SBS that is one of the two bases
    of a PASS DBS from the same read family, so SBS and DBS are mutually
    exclusive event classes. Returns the number relabelled."""
    member_keys = {
        site
        for dbs in dbsAll
        if dbs.get("filter", "PASS") == "PASS"
        for site in _dbs_family_keys(dbs)
    }
    n = 0
    for mut in sbs_records:
        if mut.get("filter", "PASS") == "PASS" and _sbs_family_key(mut) in member_keys:
            mut["filter"] = "dbs_member"
            n += 1
    return n


def _apply_strand_filter(records, pvalue_fn, p_threshold):
    """Strand independence filter for PASS records: stamp INFO MSP from each
    record's _strand inputs (pvalue_fn(*inputs) per strand) and, when
    p_threshold > 0, set FILTER strand_independence on records with a
    strand p <= p_threshold, comparing the unrounded p-values. Records
    without _strand inputs, or with a non-finite p, keep MSP "." and are
    not tested. Returns the failed records."""
    failed = []
    for rec in records:
        strand_inputs = rec["infos"].get("_strand")
        if strand_inputs is None:
            continue
        ps = [pvalue_fn(*ev) for ev in strand_inputs]
        if not all(np.isfinite(ps)):
            continue
        rec["infos"]["MSP"] = ",".join(f"{p:.4g}" for p in ps)
        if p_threshold > 0 and min(ps) <= p_threshold:
            rec["filter"] = "strand_independence"
            failed.append(rec)
    return failed


def _fail_dbs_of_failed_sbs(dbsAll, sbs_records, failed_sbs):
    """Propagate strand_independence from SBS to DBS within a read family:
    a PASS DBS with a failed SBS of its own family fails, and so does its
    other SBS. Repeats until stable (three adjacent SBS make two DBS that
    share the middle one). Returns (failed DBS, partner SBS failed)."""
    failed_keys = {_sbs_family_key(m) for m in failed_sbs}
    pass_by_key = {
        _sbs_family_key(m): m for m in sbs_records if m.get("filter", "PASS") == "PASS"
    }
    n_dbs = 0
    n_partner = 0
    changed = True
    while changed:
        changed = False
        for dbs in dbsAll:
            if dbs.get("filter", "PASS") != "PASS":
                continue
            sites = _dbs_family_keys(dbs)
            if not any(site in failed_keys for site in sites):
                continue
            dbs["filter"] = "strand_independence"
            n_dbs += 1
            changed = True
            for site in sites:
                partner = pass_by_key.pop(site, None)
                if partner is not None:
                    partner["filter"] = "strand_independence"
                    failed_keys.add(site)
                    n_partner += 1
    return n_dbs, n_partner


def check_input_files_exist(args):
    """Exit with an error message if any required input file is missing."""
    missing_files = []

    if not os.path.exists(args.bam):
        missing_files.append(f"Tumor BAM file: {args.bam}")

    if not os.path.exists(args.reference):
        missing_files.append(f"Reference genome: {args.reference}")

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
            f"  DupCaller.py index -f {args.reference} -rt <repeats.tsv>"
        )
        sys.exit(1)

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
            elif suffix.endswith(".str.txt") and not str_profile_is_current(path):
                missing_files.append(
                    f"{label} (from -E) has the pre-1.2.12 indel-length layout; "
                    f"relearn it with this version: {path}"
                )

    if args.muterateprefix:
        for suffix, label in [
            ("_sbs96_rate_n1.txt", "SBS96 mutation-rate table"),
            ("_indel_rate_by_hp_str.txt", "Indel mutation-rate table"),
        ]:
            path = args.muterateprefix + suffix
            if not os.path.exists(path):
                missing_files.append(f"{label} (from -mr): {path}")

    if missing_files:
        print("ERROR: The following required input files are missing:")
        for missing_file in missing_files:
            print(f"  - {missing_file}")
        print("\nPlease ensure all input files exist before running DupCaller.")
        sys.exit(1)


def do_call(args):
    # Defaults for args namespaces built without these options.
    if not hasattr(args, "minChunkLength"):
        args.minChunkLength = 10000
    if not hasattr(args, "srdMinRead"):
        args.srdMinRead = 3
    if not hasattr(args, "ssmMinRead"):
        args.ssmMinRead = 3
    check_input_files_exist(args)
    # -p 1 has no chunk boundaries to reconcile via MC tags.
    if args.threads != 1:
        check_mate_cigar_tags(args.bam, args.reference)
    coverage_starts = [None]
    # Resolve the Monte Carlo seed up front and write it back onto args so
    # it reaches params and is logged in _call_params.log for reproduction.
    if args.seed is None:
        args.seed = int(np.random.SeedSequence().generate_state(1)[0])
    try:
        barcode_spec = parse_barcode_option(args.barcode)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    output_dir = args.output
    sample_name = os.path.basename(args.output)
    os.makedirs(output_dir, exist_ok=True)
    args.output = os.path.join(output_dir, sample_name)
    # Learned error matrices go in their own ERROR/ subfolder.
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
        "amperr_file": error_prefix + ".amp.tn.srd.txt",
        "amperri_file": error_prefix + ".amp.id.txt",
        "dmgerr_file": error_prefix + ".dmg.tn.txt",
        "dmgerri_file": error_prefix + ".dmg.id.txt",
        # Every candidate with LR > 0 is a candidate mutation; per-channel
        # FDR refinement is what raises the threshold afterward.
        "pcutoff": 0,
        "pcutoffi": 0,
        "minLR": args.minLR,
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
        "pseudocount": args.pseudocount,
        "normalVAF": args.naf,
        "rescue": args.rescue,
        "maxDepth": args.maxPileupDepth,
        "barcodeTag": barcode_spec["tag"],
        "barcodeNormalize": barcode_spec["normalize"],
        "barcodeSep": barcode_spec["sep"],
        "seed": args.seed,
    }
    if args.errprefix:
        params["amperr_file"] = args.errprefix + ".amp.tn.srd.txt"
        params["amperri_file"] = args.errprefix + ".amp.id.txt"
        params["dmgerr_file"] = args.errprefix + ".dmg.tn.txt"
        params["dmgerri_file"] = args.errprefix + ".dmg.id.txt"

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
        "srdMinRead": args.srdMinRead,
        "ssmMinRead": args.ssmMinRead,
        "isLearn": True,
        "rescue": False,
        # callBam requires a seed in every round; one base seed makes the
        # whole run reproducible.
        "seed": args.seed,
        # Read labelling is shared by all rounds; without these, round 0
        # would fall back to the default DB-tag barcode parsing.
        "barcodeTag": barcode_spec["tag"],
        "barcodeNormalize": barcode_spec["normalize"],
        "barcodeSep": barcode_spec["sep"],
    }
    same_regions_flag = False
    if not params_learn["regions"]:
        params_learn["regions"] = params["regions"]
        same_regions_flag = True
    if not params["normalBams"]:
        print(
            "A matched normal is not used. The maximum allele fraction to call "
            f"a somatic mutation is set to be {args.maxAF}"
        )
    else:
        print(
            f"Matched normal: {args.normalBams}. The maximum allele fraction to "
            f"call a somatic mutation is set to be {args.maxAF}"
        )
    """
    Initialze run
    """
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
    profiles_exist = (
        os.path.exists(params["amperr_file"])
        and os.path.exists(amperri_hp_file)
        and os.path.exists(amperri_str_file)
        and os.path.exists(params["dmgerr_file"])
        and os.path.exists(dmgerri_hp_file)
        and os.path.exists(dmgerri_str_file)
    )
    if profiles_exist and not (
        str_profile_is_current(amperri_str_file)
        and str_profile_is_current(dmgerri_str_file)
    ):
        # Outdated profile at the default prefix (-E ones are refused in
        # check_input_files_exist): relearn rather than fail.
        print(
            "..........STR error profiles at "
            + amperri_str_file.rsplit(".amp.", 1)[0]
            + " predate the per-repeat-unit layout; relearning............"
        )
        profiles_exist = False
    if not profiles_exist:
        learn_flag = True
        if args.threads == 1:
            """
            Single-thread execution
            """
            print(".........Starting estimating error rates..............")
            paramsNow = params_learn
            regions = params_learn["regions"]
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
            contigs = params_learn["regions"]
            contigLengths = [
                bamObject.get_reference_length(contig) for contig in contigs
            ]
            print(
                "...........Spliting genomic regions for parallel execution (error estimation)................"
            )
            cutSites, chunkSize, contigs = splitBamRegions(
                [args.bam],
                args.threads,
                contigs,
                args.windowSize,
                args.reference,
                params_learn["region_file"],
                min_chunk_length=args.minChunkLength,
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
                    # Fewer chunks than -p: callBam can't take an empty
                    # region list.
                    continue
                chroms = [region[0] for region in regions]
                paramsNow = params_learn.copy()
                paramsNow["regions"] = regions
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

        # hp.txt: rows HP1..HP9, HP10+; columns base x idLen in {-1,0,1},
        # idLen=0 being the opportunity count. str.txt: rows are STR-length
        # bins ("0-1" = not a real repeat); columns are a del/ref/ins triple
        # per repeat-unit group u1-u5+, ref being the opportunity count.
        hp_row_labels = [f"HP{i}" for i in range(1, 10)] + ["HP10+"]
        hp_col_labels = [f"{b}_{d}" for b in ["A", "T", "C", "G"] for d in (-1, 0, 1)]
        str_row_labels = ["0-1", "2-9", "10-24", "25-39", "40+"]
        str_col_labels = STR_COL_LABELS
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

        # QC diagnostic (3D histogram with axis labels), kept in tmp_dir.
        bq_values = np.arange(NUM_BQ)
        np.savez(
            os.path.join(params["tmp_dir"], sample_name + ".amp.tn.bqhist.npz"),
            hist=sbs_alt_bq_hist,
            row_labels=np.array(num2trinuc),
            col_labels=np.array(["A", "T", "C", "G"]),
            bq_values=bq_values,
        )

        # SBS SRD (single-read-damage) rate matrix, EM-estimated from the BQ
        # histogram; this is the params["amperr_file"] used for calling.
        srd_mat = estimate_sbs_srd_rates(
            sbs_alt_bq_hist, args.pseudocount, fallback=True
        )
        srd_pd = pd.DataFrame(srd_mat, columns=["A", "T", "C", "G"], index=num2trinuc)
        srd_pd.to_csv(error_prefix + ".amp.tn.srd.txt", sep="\t")
        print(
            "..............Completed error estimation in "
            + str((time.time() - startTime) / 60)
            + " minutes..............."
        )

    if getattr(args, "learnOnly", False):
        print(
            "..............--learnOnly set: stopping after error-rate estimation, "
            f"skipping variant calling (error files in {error_dir}/)................."
        )
        return

    # --skipCoveragePass: round 1 also tallies coverage opportunity by
    # power-table index (funcs/power_tally.py), contracted below once the
    # final thresholds are known, in place of round 2.
    params["power_tally"] = bool(getattr(args, "skipCoveragePass", False))
    # --doNotOutputMaskedSite: masked candidates get no depth extraction in
    # either round and are dropped before output.
    params["no_masked_output"] = bool(getattr(args, "doNotOutputMaskedSite", False))
    if args.threads == 1:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        paramsNow = params
        regions = params["regions"]
        paramsNow["regions"] = [
            (chrom, 0, bamObject.get_reference_length(chrom)) for chrom in regions
        ]
        # Single worker region list so round 2 (Pool.starmap over
        # regions_list) is the same for any thread count.
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
            neg_lr_sbs,
            neg_lr_indel,
            power_tally_single,
        ) = callBam(paramsNow, 0)
        duplex_read_nums = [duplex_read_num_single]
        power_tallies = [power_tally_single]
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
                    min_chunk_length=args.minChunkLength,
                )
            else:
                cutSites, chunkSize, contigs = splitBamRegions(
                    [args.bam],
                    args.threads,
                    contigs,
                    args.windowSize,
                    args.reference,
                    regionfile,
                    min_chunk_length=args.minChunkLength,
                )
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
                if nn == 0:
                    raise Exception(
                        f"Window size {args.windowSize} is proabably too large. Change -w to a smaller value. An ideal value will be targetSize/(threads * 100)"
                    )
                # Fewer chunks than -p (duplicate chunk boundary skipped, or
                # too few windows in the target): not an error.
                print(
                    f"....Only {nn} parallel chunk(s) available (requested -p "
                    f"{args.threads}); continuing with fewer effective threads "
                    "for this pass...."
                )
                break
            paramsNow = params.copy()
            regions_list.append(regions)
            paramsNow["regions"] = regions
            callArgument = (paramsNow, nn)
            callArguments.append(callArgument)
        print("Calculating coverage boundary extents from preceding regions...")
        coverage_starts = pool.starmap(
            previous_region_coverage_start,
            [
                (
                    args.bam,
                    regions_list[nn - 1] if nn else None,
                    args.reference,
                    params.get("region_file"),
                )
                for nn in range(len(regions_list))
            ],
        )
        for (worker_params, nn), coverage_start in zip(callArguments, coverage_starts):
            worker_params["coverage_start"] = coverage_start
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
        power_tallies = [_[24] for _ in results]
        depth_by_trinuc = sum(_[17] for _ in results)
        depth_by_hpstr = sum(_[18] for _ in results)
        # Each worker simulates the same power tables independently;
        # averaging only reduces Monte Carlo noise.
        L_power = np.mean(np.stack([_[19] for _ in results]), axis=0)
        L_indel_1bp = np.mean(np.stack([_[20] for _ in results]), axis=0)
        L_indel_len = np.mean(np.stack([_[21] for _ in results]), axis=0)
        neg_lr_sbs = [np.concatenate([_[22][k] for _ in results]) for k in range(256)]
        neg_lr_indel = {
            name: np.concatenate([_[23][name] for _ in results if name in _[23]])
            for name in set().union(*(_[23].keys() for _ in results))
        }
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

        # Per-duplex-group coverage tables (duplex_read_num*, etc.) are not
        # computed in round 1; they are built from round 2's results below.
        FPAll = sum(FPs, [])
        RPAll = sum(RPs, [])

    tBam = BAM(args.bam, "rb", args.reference)
    contigs = tBam.references
    chromDict = {contig: tBam.get_reference_length(contig) for contig in contigs}
    infoDict = {
        "F1R2": [1, "Integer", "Number of F1R2 read(s) in the read bundle"],
        "F2R1": [1, "Integer", "Number of F2R1 read(s) in the read bundle"],
        "LR": [
            1,
            "Float",
            "log10 likelihood ratio of the alt allele (SBS: the position's non-reference base; indel: the indel) over the reference allele",
        ],
        "LM": [
            1,
            "Float",
            "maximum log10 likelihood ratio the context's damage rates allow (the LR ceiling)",
        ],
        "MSP": [
            2,
            "String",
            "Strand independence hypothesis p-values (F1R2,F2R1): P(>= observed non-alt reads on that strand | clean mutant strand). SBS: Poisson-binomial over every read covering the position (actual BQ, no minBq cut) with per-read non-alt probability 10^(-BQ/10) + amplification error. Indel: binomial over the informative (ALT+REF) reads with per-read probability of the learned amplification error that reverts the indel. Small = more discordant reads than read error explains (mixed read family). Computed only for calls that are PASS going into the strand independence filter; '.' on every other record. See --p_threshold",
        ],
        "TC": [4, "Integer", "Top strand base count"],
        "BC": [4, "Float", "Bottom strand base count"],
        "DF": [1, "Integer", "Distance from fragment end"],
        "DR": [1, "Integer", "Distance from read end"],
        "TAG1": [1, "String", "Barcode of top strand 5 prime"],
        "TAG2": [1, "String", "Barcode of bottom strand 5 prime"],
        "SP": [1, "Integer", "Read family reference start position"],
        "TL": [
            1,
            "Integer",
            "Template length of the read pair (SAM TLEN of the family's first read; negative when that read is the downstream mate)",
        ],
        "TN": [1, "String", "trinucleotide context"],
        "HP": [1, "Integer", "Homopolymer length. Always 0 for SBS"],
        "STR": [
            1,
            "Integer",
            "STR bin of the indel's reference tract: 0 = no STR, or a 2bp+ indel that is not a whole number of the tract's repeat units; 1 = 2-9bp; 2 = 10-24bp; 3 = 25-39bp; 4 = 40bp+. Always 0 for SBS",
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
            "Local false discovery rate for this call: (1-mu)/(rawLR*mu+1-mu), mu being this call's channel's mutation rate estimate capped at the rate whose LR threshold is --minLR. Computed and reported regardless of filter, so a fail.vcf record's own FDR explains why it didn't clear its channel's threshold. 1.0 if this call's channel has no determinable mu (e.g. an indel with no HP/STR channel at all).",
        ],
    }
    dbsInfoDict = {
        "F1R2": infoDict["F1R2"],
        "F2R1": infoDict["F2R1"],
        "TAG1": infoDict["TAG1"],
        "TAG2": infoDict["TAG2"],
        "SP": infoDict["SP"],
        "TL": infoDict["TL"],
        "FDR": infoDict["FDR"],
    }
    formatDict = {
        "AC": [1, "Integer", "Count of alt allele"],
        "RC": [1, "Integer", "Count of ref allele"],
        "DP": [1, "Integer", "Depth at the location"],
    }
    filterDict = {
        "PASS": "All filters passed, including snp_mask/noise_mask (see SNPM/NOISEM INFO fields, both 0 here) -- a fully unmasked call",
        "masked": "Blocked by snp_mask/noise_mask only (SNPM/NOISEM INFO fields mark which); other mask types never reach this label (see the reasons below, only emitted under --rescue). A masked candidate whose raw LR clears its channel's final refined threshold gets real depth extracted and feeds the unmasked-burden numerator, unless that real depth then fails a post-hoc sanity check (see no_good_alt_read/duplex_vaf/normal_vaf/n_cov_mask below), in which case it's relabeled to that specific reason and, like any other reject reason, only kept under --rescue. A masked candidate that never got real depth extracted (LR below its channel's final threshold) is dropped entirely and never appears in the fail vcf unless --rescue is on",
        "underpowered": "Candidate with LR > 0 that failed the FDR-refined per-channel threshold",
        "zero_LR": "Candidate with no evidence for the alt over the ref (log10 LR == 0, e.g. alt and non-alt reads cancelling); reported for completeness, never PASS. A DBS inherits it from either of its own read family's SNVs",
        "high_nm": "Read family failed the NM/blacklist filter -- only reported under --rescue",
        "low_mapq": "Read family failed the mapq filter -- only reported under --rescue",
        "low_ASXS": "Read family failed the AS-XS filter -- only reported under --rescue",
        "cov_mask": "Blocked by the coverage-depth mask (n_cov_mask); no coverage/depth extracted -- only reported under --rescue",
        "nm_mask": "Blocked by the per-family NM mask at this position; no coverage/depth extracted -- only reported under --rescue",
        "trim_mask": "Falls in the read-end trim zone; no coverage/depth extracted -- only reported under --rescue",
        "indel_mask": "Indel blocked by indel_mask; no coverage/depth extracted -- only reported under --rescue",
        "other_mask": "Blocked by a mask not covered by a more specific reason above; no coverage/depth extracted -- only reported under --rescue",
        "no_good_alt_read": "Real depth-extraction found zero alt-supporting reads in the tumor BAM pileup (LR cleared its channel's threshold on family-consensus evidence alone) -- real AC/RC/DP attached; only reported under --rescue",
        "duplex_vaf": "The extracted tumor allele fraction (AC/DP) exceeds --maxAF -- real AC/RC/DP attached; only reported under --rescue",
        "normal_vaf": "The extracted matched-normal allele fraction exceeds --naf (likely germline or a systematic artifact) -- real AC/RC/DP attached; only reported under --rescue",
        "n_cov_mask": "Matched-normal depth at this position fell below --minNdepth once real depth was extracted -- real AC/RC/DP attached; only reported under --rescue",
        "dbs_member": "SBS that is one of the two bases of a PASS DBS from the same read family; reported once, as the DBS, so SBS and DBS calls (and their burdens/spectra) are mutually exclusive",
        "strand_independence": "SBS or indel call with a strand whose strand independence hypothesis p-value (INFO MSP) is <= --p_threshold: more non-alt reads than read error explains, i.e. a mixed read family. A DBS containing such an SBS (same read family) fails too, and so does that DBS's other SBS",
    }

    # Per-type SBS/INDEL/DBS subfolders (same layout Estimate.py uses).
    # VCFs are written once, at the end, after refinement and round 2's
    # depth extraction.
    sbs_dir, indel_dir, dbs_dir = _ensure_type_subdirs(output_dir)

    # duplex_read_num_trinuc[key] is a raw (64, 4) matrix: 64 trinuc contexts
    # (not RC-folded) x 4 alt bases, the reference-base column always 0.
    # index_192 selects the 192 valid (context, alt) cells in label order.
    trinuc2num_64, num2trinuc_64 = build_trinuc64_order()
    label2num_192, index_192 = build_trinuc192_labels(num2trinuc_64)
    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}

    # Sets ampmat/dmgmat_*/trinuc_convert on params, used by
    # _indel_ctx_threshold and by round 1/2 calling.
    load_error_matrices(params)

    def _indel_ctx_threshold(hps_c, strs_c, id_len, pool=None):
        # Default cutoff capped by the context's LR ceiling (indelMaxLR).
        # For +-1bp non-STR indels Pdmg depends on the base, so take the min
        # over the channel's pooled bases (pool "C"=C/G, "T"=A/T), the same
        # approximation the HP-channel Eeff uses. hps_c == 0 (real HP
        # lengths start at 1) is the sentinel for the mismatched-insertion
        # background channel (str.txt row 0), which has no base axis.
        pcutoffi = params["pcutoffi"]

        def _maxlr(ref_allele, inserted_base):
            # inserted_base == ref_allele: same-base homopolymer scenario;
            # a different inserted_base: row-0 mismatch scenario.
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
                # STR channels: id_len is the signed repeat unit (bins 1-4).
                unit_len=abs(id_len),
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

    # -mr/--muterateprefix: reuse a genome-wide run's mu0 table instead of
    # a noisy small-region solve. A channel missing from the table is an
    # error, not a silent fallback to a fresh solve.
    mu0_override_table = (
        _load_mutation_rate_override(args.muterateprefix)
        if args.muterateprefix
        else None
    )

    def _mu0_for(name):
        return mu0_override_table.get(name) if mu0_override_table is not None else None

    startTime_refine = time.time()
    # The mu solve uses every unmasked site that passes all other filters,
    # whatever the sign of its LR: round-1 PASS records (raw LRs) plus the
    # LR < 0 sites in neg_lr_sbs (float32 log10 LRs, kept apart and only
    # exponentiated inside refine_channel_task). Eeff excludes the same
    # masked positions.
    raw_lr_192_n1 = [[] for _ in range(192)]
    for mut in mutsAll:
        if mut.get("filter", "PASS") != "PASS":
            continue
        if min(int(mut["infos"]["F1R2"]), int(mut["infos"]["F2R1"])) < 1:
            continue
        idx_192 = label2num_192[f"{mut['infos']['TN']}>{mut['alt']}"]
        raw_lr_192_n1[idx_192].append(10 ** mut["infos"]["LR"])
    raw_lr_192_n1 = [np.asarray(lrs, dtype=float) for lrs in raw_lr_192_n1]
    neg_lr_192_n1 = [np.empty(0, dtype=np.float32)] * 192
    for k, lrs in enumerate(neg_lr_sbs):
        if lrs.size:
            label = f"{num2trinuc_64[k // 4]}>{'ATCG'[k % 4]}"
            neg_lr_192_n1[label2num_192[label]] = lrs
    # Pool each canonical SBS-96 class's raw LRs with its RC-paired raw-192
    # class (as combine_raw192_to_sbs96 does for counts) and build the
    # per-channel refinement job specs.
    _, sbs96_labels_for_rc = _build_sbs96_labels()
    raw_lr_96_n1 = []
    neg_lr_96_n1 = []
    sbs_jobs = []
    for k, label in enumerate(sbs96_labels_for_rc):
        trinuc = label[0] + label[2] + label[6]
        alt = label[4]
        rc_trinuc = "".join(_REVCOMP[b] for b in reversed(trinuc))
        rc_alt = _REVCOMP[alt]
        raw_lr_96_n1.append(
            np.concatenate(
                [
                    raw_lr_192_n1[label2num_192[f"{trinuc}>{alt}"]],
                    raw_lr_192_n1[label2num_192[f"{rc_trinuc}>{rc_alt}"]],
                ]
            )
        )
        neg_lr_96_n1.append(
            np.concatenate(
                [
                    neg_lr_192_n1[label2num_192[f"{trinuc}>{alt}"]],
                    neg_lr_192_n1[label2num_192[f"{rc_trinuc}>{rc_alt}"]],
                ]
            )
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
                args.lfdrThreshold,
                _mu0_for(f"SBS96:{label}"),
                neg_lr_96_n1[k],
                args.minLR,
            )
        )

    # HP/STR indel channels, analogous to SBS-96. rawLRs come from PASS
    # candidates plus neg_lr_indel; Eeff is a depth sum from depth_by_hpstr
    # (one site per homopolymer run / STR tract; see indel_eeff_site_masks).
    # HP channels exist only for id_len +-1 (hp.txt has no multi-bp
    # columns); longer indels are STR context. The mismatched-insertion
    # background channel (STR0_u1_ins) is routed by the call-time "HM"
    # INFO field, since only the record, not the reference, is available
    # here.
    hp_indel_lengths = [-1, 1]
    # STR channels: signed repeat unit 2-5+ in bins 1-4 (any whole number
    # of units), signed indel length 2-5+ in STR0 (str_channel_name).
    str_indel_lengths = [-5, -4, -3, -2, 2, 3, 4, 5]
    # HP channels pool 1bp events by base: "C" for C/G, "T" for A/T (the
    # same split as ID83's C/T homopolymer classes).
    _HP_POOL_BASES = {"C": (2, 3), "T": (0, 1)}  # base2num order A,T,C,G

    def _indel_channel(indel):
        # Channel of an indel record (None if none), routed exactly as
        # call.py routes LR < 0 candidates. Deletion: ref is anchor+deleted
        # bases, so ref[1] is the deleted base; insertion: alt[1].
        id_len = len(indel["alt"]) - len(indel["ref"])
        if id_len == 0:
            return None
        return indel_mu_channel(
            id_len,
            int(indel["infos"]["HP"]),
            int(indel["infos"]["STR"]),
            indel["infos"].get("HM", 1),
            indel["ref"][1] if id_len < 0 else indel["alt"][1],
            indel["infos"].get("SU", 0),
        )

    # Raw LRs of the PASS candidates per channel; LR < 0 candidates arrive
    # pre-routed (by the same indel_mu_channel) as float32 log10 LRs in
    # neg_lr_indel and never as records.
    indel_channel_names = (
        [str_channel_name(0, 1)]
        + [
            f"HP{hp_len}_len{id_len}_{pool}"
            for hp_len in range(1, 11)
            for id_len in hp_indel_lengths
            for pool in ("C", "T")
        ]
        + [
            str_channel_name(str_bin, id_len)
            for str_bin in range(0, 5)
            for id_len in str_indel_lengths
        ]
    )
    raw_lr_indel = {name: [] for name in indel_channel_names}
    _no_neg_lr = np.empty(0, dtype=np.float32)
    for indel in indelsAll:
        if indel.get("filter", "PASS") != "PASS":
            continue
        name = _indel_channel(indel)
        if name is not None:
            raw_lr_indel[name].append(10 ** indel["infos"]["LR"])

    # Nothing with LR < 0 is ever written to a VCF (PASS or _fail).
    indelsAll = [m for m in indelsAll if m["infos"]["LR"] >= 0]
    mutsAll = [m for m in mutsAll if m["infos"]["LR"] >= 0]

    # threshold0 is cheap (no Monte Carlo) so it is computed here; the Eeff
    # simulation and refinement run in the pool below.
    str0_ins1 = str_channel_name(0, 1)
    indel_jobs = [
        (
            str0_ins1,
            "str",
            (0, 1),
            raw_lr_indel[str0_ins1],
            _indel_ctx_threshold(0, 0, 1),
            args.lfdrThreshold,
            _mu0_for(str0_ins1),
            neg_lr_indel.get(str0_ins1, _no_neg_lr),
            args.minLR,
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
                        raw_lr_indel[f"HP{hp_len}_len{id_len}_{pool}"],
                        threshold0,
                        args.lfdrThreshold,
                        _mu0_for(f"HP{hp_len}_len{id_len}_{pool}"),
                        neg_lr_indel.get(f"HP{hp_len}_len{id_len}_{pool}", _no_neg_lr),
                        args.minLR,
                    )
                )
    for str_bin in range(0, 5):
        for id_len in str_indel_lengths:
            threshold0 = _indel_ctx_threshold(1, str_bin, id_len)
            str_name = str_channel_name(str_bin, id_len)
            indel_jobs.append(
                (
                    str_name,
                    "str",
                    (str_bin, id_len),
                    raw_lr_indel[str_name],
                    threshold0,
                    args.lfdrThreshold,
                    _mu0_for(str_name),
                    neg_lr_indel.get(str_name, _no_neg_lr),
                    args.minLR,
                )
            )

    # Per-channel FDR control (96 SBS96 + 40 HP + 40 STR + STR0_u1_ins):
    # mu0 is the channel's MLE mixture weight at threshold0 (direct brentq
    # solve), and the round-2 threshold is the LR whose local fdr
    # (1-mu0)/(LR*mu0+1-mu0) equals args.lfdrThreshold (_refine_channel).
    # A channel with no root in (0, 1) gets mu0 = 1/(Eeff+1), whose
    # threshold none of its own sites can reach.
    # Channels are independent and each Eeff is a Monte Carlo simulation,
    # so one task per channel runs in a Pool capped at args.threads.
    if mu0_override_table is not None:
        all_job_names = {job[0] for job in sbs_jobs + indel_jobs}
        missing = sorted(all_job_names - mu0_override_table.keys())
        if missing:
            raise ValueError(
                f"-mr/--muterateprefix table {args.muterateprefix!r} is missing "
                f"{len(missing)} channel(s) this run expects (e.g. "
                f"{missing[:5]}) -- likely a stale or mismatched-version table; "
                "regenerate it from a full run with the current DupCaller version."
            )
    refine_pool = Pool(
        args.threads,
        initializer=init_refine_worker,
        initargs=(depth_by_trinuc, depth_by_hpstr),
    )
    refine_results = refine_pool.map(refine_channel_task, sbs_jobs + indel_jobs)
    refine_pool.close()
    refine_pool.terminate()
    refine_pool.join()
    sbs_results = refine_results[: len(sbs_jobs)]
    indel_results = refine_results[len(sbs_jobs) :]

    channel_thresholds = {}
    # mu for the per-call local FDR and stats.txt "Total FDR": mu0 capped
    # at the rate whose threshold is --minLR, matching the channel's actual
    # threshold. The rate tables keep the uncapped mu0.
    channel_final_mu = {}
    mu_fdr_cap = mu_at_min_lr(args.lfdrThreshold, args.minLR)

    sbs96_rate_rows = []
    for k, (name, kind, ctx_key, Eeff0, mu0, new_threshold) in enumerate(sbs_results):
        channel_final_mu[name] = min(mu0, mu_fdr_cap)
        channel_thresholds[name] = new_threshold
        sbs96_rate_rows.append(
            {
                "context": name,
                "n_sites": len(raw_lr_96_n1[k]) + len(neg_lr_96_n1[k]),
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
    # Consume results in indel_jobs order (STR0_u1_ins first).
    name, kind, ctx_key, Eeff0, mu0, new_threshold = next(indel_results_iter)
    indel_rate_rows.append(
        {
            "context": "STR0_u1",
            "indel_length": 1,
            "n_sites": len(raw_lr_indel[name])
            + len(neg_lr_indel.get(name, _no_neg_lr)),
            "Eeff": Eeff0,
            "mutation_rate_mle": mu0,
        }
    )
    channel_final_mu[name] = min(mu0, mu_fdr_cap)
    channel_thresholds[name] = new_threshold
    for hp_len in range(1, 11):
        for id_len in hp_indel_lengths:
            for pool in ("C", "T"):
                name, kind, ctx_key, Eeff0, mu0, new_threshold = next(
                    indel_results_iter
                )
                channel_final_mu[name] = min(mu0, mu_fdr_cap)
                indel_rate_rows.append(
                    {
                        "context": f"HP{hp_len}_{pool}",
                        "indel_length": id_len,
                        "n_sites": len(raw_lr_indel[name])
                        + len(neg_lr_indel.get(name, _no_neg_lr)),
                        "Eeff": Eeff0,
                        "mutation_rate_mle": mu0,
                    }
                )
                channel_thresholds[name] = new_threshold
    for str_bin in range(0, 5):
        for id_len in str_indel_lengths:
            name, kind, ctx_key, Eeff0, mu0, new_threshold = next(indel_results_iter)
            channel_final_mu[name] = min(mu0, mu_fdr_cap)
            indel_rate_rows.append(
                {
                    "context": f"STR{str_bin}_u{abs(id_len)}",
                    "indel_length": id_len,
                    "n_sites": len(raw_lr_indel[name])
                    + len(neg_lr_indel.get(name, _no_neg_lr)),
                    "Eeff": Eeff0,
                    "mutation_rate_mle": mu0,
                }
            )
            channel_thresholds[name] = new_threshold
    indel_rate_by_hp_str = pd.DataFrame(indel_rate_rows)
    indel_rate_by_hp_str["eeff_sites"] = INDEL_EEFF_SITE_SET
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

    sample_name = os.path.basename(args.output)
    sample_dir = params["tmp_dir"]

    # Per-channel round-2 threshold overrides, covering every channel.
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
            # "HP{hp_len}_len{id_len}_{pool}", id_len always +-1.
            hp_len_str, rest = name[2:].split("_len")
            id_len_str, pool = rest.split("_")
            hp_len, id_len = int(hp_len_str), int(id_len_str)
            sign_idx = 1 if id_len == 1 else 0
            indel_1bp_threshold_override[(hp_len, sign_idx, pool)] = new_threshold
        elif name.startswith("STR"):
            # Row = str_bin (0-4), matching L_indel_len's row axis.
            str_bin, id_len = parse_str_channel_name(name)
            indel_len_threshold_override[(str_bin, id_len)] = new_threshold

    # Re-filter against each channel's final threshold. LR does not depend
    # on the threshold, so no rescan is needed. The final threshold is
    # >= threshold0, so only PASS records can flip (to "underpowered").
    # "masked" records (blocked only by snp_mask/noise_mask) whose LR clears
    # the threshold are deferred for real depth extraction in round 2; they
    # stay "masked" and are never promoted to PASS. Other mask labels
    # (--rescue only) never get depth.
    def _mut_key(m):
        return (m["chrom"], m["pos"], m["ref"], m["alt"])

    deferred_depth_keys = set()
    # Final threshold of each deferred key: picks the --skipCoveragePass
    # round-1 depth entry extracted with the families clearing it.
    deferred_thresholds = {}
    for mut in mutsAll:
        old_filter = mut.get("filter", "PASS")
        if old_filter not in ("PASS", "masked"):
            continue
        ctx_key = (trinuc2num_64.get(mut["infos"]["TN"]), base2num.get(mut["alt"]))
        threshold = pcutoff_sbs_override.get(
            ctx_key, max(params["pcutoff"], args.minLR)
        )
        if old_filter == "PASS":
            # LR == 0: no evidence either way.
            if mut["infos"]["LR"] <= 0:
                mut["filter"] = "zero_LR"
            elif mut["infos"]["LR"] < threshold:
                mut["filter"] = "underpowered"
        elif (
            not params["no_masked_output"]
            and mut["infos"]["LR"] >= threshold
            and (mut["infos"].get("SNPM") or mut["infos"].get("NOISEM"))
        ):
            deferred_depth_keys.add(_mut_key(mut))
            deferred_thresholds[_mut_key(mut)] = threshold

    for indel in indelsAll:
        old_filter = indel.get("filter", "PASS")
        if old_filter not in ("PASS", "masked"):
            continue
        # The channel indelErrorProbs scored this call against (multi-bp:
        # STR bin and repeat unit; 1bp: HP, or STR0_u1_ins for a mismatched
        # insertion).
        threshold = channel_thresholds.get(
            _indel_channel(indel), max(params["pcutoffi"], args.minLR)
        )
        if old_filter == "PASS":
            if indel["infos"]["LR"] <= 0:
                indel["filter"] = "zero_LR"
            elif indel["infos"]["LR"] < threshold:
                indel["filter"] = "underpowered"
        elif (
            not params["no_masked_output"]
            and indel["infos"]["LR"] >= threshold
            and indel["infos"].get("NOISEM")
        ):
            deferred_depth_keys.add(_mut_key(indel))

    # A DBS requires both of its SNVs (same read family) to PASS, so
    # revalidate after the SBS re-filter. DBS is never "masked", so there
    # is no deferred set for DBS.
    pass_keys = {
        _sbs_family_key(mut) for mut in mutsAll if mut.get("filter", "PASS") == "PASS"
    }
    no_alt_keys = {
        _sbs_family_key(mut)
        for mut in mutsAll
        if mut.get("filter") == "no_good_alt_read"
    }
    zero_lr_keys = {
        _sbs_family_key(mut) for mut in mutsAll if mut.get("filter") == "zero_LR"
    }
    for dbs in dbsAll:
        if dbs.get("filter", "PASS") != "PASS":
            continue
        sites = _dbs_family_keys(dbs)
        if any(site in no_alt_keys for site in sites):
            # Inherit the constituent SNV's rescue-only reason instead of
            # "underpowered", so it's pruned below unless --rescue.
            dbs["filter"] = "no_good_alt_read"
        elif any(site in zero_lr_keys for site in sites):
            dbs["filter"] = "zero_LR"
        elif not all(site in pass_keys for site in sites):
            dbs["filter"] = "underpowered"

    efficiency = duplex_num / rec_num if rec_num > 0 else 0.0

    # Round 2 computes position-weighted duplex coverage and coverage.bed.gz.
    # Under --skipCoveragePass it is skipped: round 1's power tallies are
    # contracted at the final thresholds (same by-duplex-group/stats
    # coverage), deferred masked candidates take the depth round 1
    # extracted for them, and coverage.bed.gz is not written.
    coverage_pass_enabled = not args.skipCoveragePass
    if coverage_pass_enabled:
        # Round 2: coverage-only pass over round 1's region split with the
        # final per-channel thresholds; also does the deferred depth
        # extraction for qualifying masked candidates.
        print("..............Starting coverage pass (round 2)..............")
        startTime3 = time.time()
        round2_args = []
        for nn, r2_regions in enumerate(regions_list):
            paramsNow2 = params.copy()
            paramsNow2["regions"] = r2_regions
            paramsNow2["coverage_start"] = coverage_starts[nn]
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

        # Merge deferred depth back by position key. The filter comes along:
        # call.py may relabel a deferred candidate to a reject reason when
        # its real depth fails a sanity check, else it stays "masked".
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
        # Coverage below comes from round 2 (final thresholds).
        coverages = [_[1] for _ in round2_results]
        duplex_read_nums = [_[4] for _ in round2_results]
        duplex_read_nums_trinuc = [_[5] for _ in round2_results]
        unmasked_coverages = [_[11] for _ in round2_results]
        coverages_indels_cat = [_[12] for _ in round2_results]
        unmasked_coverages_indels_cat = [_[13] for _ in round2_results]
        duplex_read_nums_indel = [_[14] for _ in round2_results]
        duplex_read_nums_dbs = [_[16] for _ in round2_results]
    else:
        print(
            "..............Skipping coverage pass (round 2): coverage from "
            "round 1's power tallies; per-locus coverage.bed.gz will not be "
            "generated (--skipCoveragePass).............."
        )
        startTime3 = time.time()
        params_final = params.copy()
        params_final["isLearn"] = False
        params_final["pcutoff_sbs_override"] = pcutoff_sbs_override
        params_final["indel_1bp_threshold_override"] = indel_1bp_threshold_override
        params_final["indel_len_threshold_override"] = indel_len_threshold_override
        tally_pool = Pool(args.threads)
        contracted = tally_pool.starmap(
            contract_worker_tally, [(params_final, t) for t in power_tallies]
        )
        tally_pool.close()
        tally_pool.terminate()
        tally_pool.join()
        print(
            "..............Contracted round 1's power tallies in "
            + str((time.time() - startTime3) / 60)
            + " minutes................."
        )
        duplex_read_nums_trinuc = [c[0] for c in contracted]
        duplex_read_nums_indel = [c[1] for c in contracted]
        duplex_read_nums_dbs = [c[2] for c in contracted]
        coverages = [c[3] for c in contracted]
        unmasked_coverages = [c[4] for c in contracted]

        def _deferred_entry(m):
            # The entry extracted with exactly the families clearing the
            # key's final threshold (smallest cut at or above it), matching
            # round 2. Indel entries have cut -inf.
            final = deferred_thresholds.get(_mut_key(m), -np.inf)
            entries = [e for e in m["deferred_depth"] if e[0] >= final]
            return min(entries, key=lambda e: e[0])[1] if entries else None

        deferred_mut_updates = {}
        deferred_indel_updates = {}
        for records, updates in (
            (mutsAll, deferred_mut_updates),
            (indelsAll, deferred_indel_updates),
        ):
            for m in records:
                if "deferred_depth" in m and _mut_key(m) in deferred_depth_keys:
                    entry = _deferred_entry(m)
                    if entry is not None:
                        updates[_mut_key(m)] = entry
        for m in mutsAll + indelsAll:
            m.pop("deferred_depth", None)

    for mut in mutsAll:
        key = _mut_key(mut)
        if key in deferred_mut_updates:
            mut["samples"], mut["filter"] = deferred_mut_updates[key]
    for indel in indelsAll:
        key = _mut_key(indel)
        if key in deferred_indel_updates:
            indel["samples"], indel["filter"] = deferred_indel_updates[key]

    coverage = sum(coverages)
    unmasked_coverage = sum(unmasked_coverages)

    duplex_combinations = list(set.union(*[set(d.keys()) for d in duplex_read_nums]))
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
            num: sum([d.get(num, np.zeros((64, 4))) for d in duplex_read_nums_trinuc])
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

    # Includes the 0+n / n+0 keys.
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
    row_idx_192 = np.array([trinuc2num_64[label.split(">")[0]] for label in index_192])
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
    # Records that cleared LR but failed a post-hoc depth sanity check are
    # reported only under --rescue. This is the single keep/drop decision
    # for these labels (call.py never drops them). no_good_alt_read (ta==0)
    # means the pileup re-check found no alt-supporting read at all.
    if not params["rescue"]:
        _reject_reasons = {"no_good_alt_read", "duplex_vaf", "normal_vaf", "n_cov_mask"}
        mutsAll = [m for m in mutsAll if m.get("filter", "PASS") not in _reject_reasons]
        indelsAll = [
            m for m in indelsAll if m.get("filter", "PASS") not in _reject_reasons
        ]
        dbsAll = [m for m in dbsAll if m.get("filter", "PASS") not in _reject_reasons]

    # "masked" + all-zero samples unambiguously means never depth-extracted
    # (LR below the final threshold): a masked record that did get depth
    # with ta==0 is relabelled no_good_alt_read by call.py. Drop these
    # unless --rescue.
    if not params["rescue"]:

        def _masked_no_depth(m):
            return m.get("filter") == "masked" and m.get("samples") == [
                [0, 0, 0],
                [0, 0, 0],
            ]

        mutsAll = [m for m in mutsAll if not _masked_no_depth(m)]
        indelsAll = [m for m in indelsAll if not _masked_no_depth(m)]
    # --doNotOutputMaskedSite: no masked record is output, --rescue or not.
    if params["no_masked_output"]:
        mutsAll = [m for m in mutsAll if m.get("filter") != "masked"]
        indelsAll = [m for m in indelsAll if m.get("filter") != "masked"]

    # Strand independence filter on final PASS SBS/indel calls (MSP stamped
    # only on these). A DBS fails with either of its family's SBS, and takes
    # the other SBS down with it.
    pass_sbs = [m for m in mutsAll if m.get("filter", "PASS") == "PASS"]
    failed_sbs = _apply_strand_filter(pass_sbs, strand_pvalue, args.p_threshold)
    pass_indels = [m for m in indelsAll if m.get("filter", "PASS") == "PASS"]
    failed_indels = _apply_strand_filter(
        pass_indels, indel_strand_pvalue, args.p_threshold
    )
    if args.p_threshold > 0:
        n_si_dbs, n_si_partner = _fail_dbs_of_failed_sbs(dbsAll, pass_sbs, failed_sbs)
        print(
            f"..............Strand independence hypothesis filter (--p_threshold "
            f"{args.p_threshold}): {len(failed_sbs)} of {len(pass_sbs)} PASS SBS "
            f"calls failed ({n_si_partner} more as the other SBS of {n_si_dbs} "
            f"failed DBS), {len(failed_indels)} of {len(pass_indels)} PASS indel "
            f"calls failed..............."
        )

    # Per-call local FDR, (1-mu)/(rawLR*mu+1-mu) with the channel's capped
    # mu, stamped on every record (not just PASS) and averaged over PASS
    # calls for the stats.txt "Total FDR" lines. A DBS's FDR is the chance
    # at least one of its two SNVs is wrong, 1-(1-fdr1)*(1-fdr2). Records
    # with no channel mu get FDR=1.0 ("not modeled").
    def _sbs_channel_name(trinuc, alt):
        ref = trinuc[1]
        if ref in ("C", "T"):
            return f"SBS96:{trinuc[0]}[{ref}>{alt}]{trinuc[2]}"
        rc_trinuc = "".join(_REVCOMP[b] for b in reversed(trinuc))
        rc_alt = _REVCOMP[alt]
        return f"SBS96:{rc_trinuc[0]}[{rc_trinuc[1]}>{rc_alt}]{rc_trinuc[2]}"

    def _local_fdr(log10_lr, mu):
        """(1-mu)/(rawLR*mu+1-mu); must match the formula _refine_channel
        (funcs/misc.py) inverts for its LR threshold."""
        return (1.0 - mu) / ((10**log10_lr) * mu + (1.0 - mu))

    sbs_local_fdrs = []
    key_to_local_fdr = {}
    for mut in mutsAll:
        name = _sbs_channel_name(mut["infos"]["TN"], mut["alt"])
        mu = channel_final_mu.get(name)
        local_fdr = _local_fdr(mut["infos"]["LR"], mu) if mu is not None else 1.0
        mut["infos"]["FDR"] = local_fdr
        key_to_local_fdr[_sbs_family_key(mut)] = local_fdr
        if mu is not None and mut.get("filter", "PASS") == "PASS":
            sbs_local_fdrs.append(local_fdr)
    total_sbs_fdr = float(np.mean(sbs_local_fdrs)) if sbs_local_fdrs else float("nan")

    indel_local_fdrs = []
    for indel in indelsAll:
        mu = channel_final_mu.get(_indel_channel(indel))
        local_fdr = _local_fdr(indel["infos"]["LR"], mu) if mu is not None else 1.0
        indel["infos"]["FDR"] = local_fdr
        if mu is not None and indel.get("filter", "PASS") == "PASS":
            indel_local_fdrs.append(local_fdr)
    total_indel_fdr = (
        float(np.mean(indel_local_fdrs)) if indel_local_fdrs else float("nan")
    )

    dbs_local_fdrs = []
    for dbs in dbsAll:
        fdr1, fdr2 = (key_to_local_fdr.get(site) for site in _dbs_family_keys(dbs))
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

    # SBS and DBS are mutually exclusive: a PASS DBS's two bases become
    # dbs_member. Must run after the strand filter (needs PASS members to
    # propagate failures) and the DBS FDR (built from the members' FDRs).
    n_dbs_member = _relabel_dbs_members(mutsAll, dbsAll)
    if n_dbs_member:
        print(f"{n_dbs_member} SBS reported as part of a PASS DBS (dbs_member)")
    sbs_local_fdrs = [
        mut["infos"]["FDR"]
        for mut in mutsAll
        if mut.get("filter", "PASS") == "PASS"
        and channel_final_mu.get(_sbs_channel_name(mut["infos"]["TN"], mut["alt"]))
        is not None
    ]
    total_sbs_fdr = float(np.mean(sbs_local_fdrs)) if sbs_local_fdrs else float("nan")

    with open(args.output + "_stats.txt", "w") as f:
        f.write(f"Number of Read Families\t{unique_read_num}\n")
        f.write(f"Number of Pass-filter Reads\t{pass_read_num}\n")
        f.write(f"Number of Effective Read Families\t{duplex_num}\n")
        # SNV coverage is divided by 3 (possible alt bases per locus), as
        # Estimate.py does for trinuc_cov. DBS and indel coverage need
        # reference-wide normalization, so the estimate step writes them.
        # Unmasked Coverage is N/A under --doNotOutputMaskedSite.
        if params["no_masked_output"]:
            f.write("Unmasked Coverage\tN/A\n")
        else:
            f.write(f"Unmasked Coverage\t{unmasked_coverage / 3}\n")
        f.write(
            f"Per Read Family Coverage \t{(coverage / 3) / duplex_num if duplex_num > 0 else 0.0}\n"
        )
        f.write(
            "Pass-filter Duplication Rate\t"
            f"{1-unique_read_num/pass_read_num if pass_read_num > 0 else 0.0}\n"
        )
        f.write(f"Efficiency\t{efficiency}\n")
        f.write(f"Total SBS FDR\t{total_sbs_fdr}\n")
        f.write(f"Total Indel FDR\t{total_indel_fdr}\n")
        f.write(f"Total DBS FDR\t{total_dbs_fdr}\n")

    if coverage_pass_enabled:
        # Combine round 2's per-worker coverage files. Use the number of
        # dispatched chunks, which can be fewer than -p.
        if merge_and_combine_coverage_files(sample_name, sample_dir, len(regions_list)):
            os.replace(
                os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz"),
                args.output + "_coverage.bed.gz",
            )
            os.replace(
                os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz.tbi"),
                args.output + "_coverage.bed.gz.tbi",
            )
        else:
            print(
                f"Warning: {args.output}_coverage.bed.gz was not written; "
                f"per-worker coverage files are in {sample_dir}.",
                file=sys.stderr,
            )

    # PASS calls in the main per-type VCFs, everything else in _fail.vcf.
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


def _boundary_bed_extent(path):
    """Return (chromosome, largest start) for a next_region boundary
    coverage file, or None if it doesn't exist or has no data rows.
    """
    if not os.path.exists(path):
        return None
    max_pos = None
    boundary_chrom = None
    with bgzf.open(path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start = parts[0], int(parts[1])
            if boundary_chrom is not None and chrom != boundary_chrom:
                raise ValueError(
                    f"Boundary coverage spans multiple chromosomes: {path}"
                )
            boundary_chrom = chrom
            if max_pos is None or start > max_pos:
                max_pos = start
    return (boundary_chrom, max_pos) if max_pos is not None else None


def _split_bed_by_start(path, cutoff, below_path, above_path):
    """Peel rows on cutoff's chromosome at or below its start coordinate.

    All other chromosomes and later positions remain in above_path.
    """
    boundary_chrom, max_pos = cutoff
    with bgzf.open(path, "rt") as fin, bgzf.open(below_path, "wt") as fbelow, bgzf.open(
        above_path, "wt"
    ) as fabove:
        for line in fin:
            if not line.strip():
                continue
            chrom, start, _ = line.split("\t", 2)
            in_overlap = chrom == boundary_chrom and int(start) <= max_pos
            (fbelow if in_overlap else fabove).write(line)


def _main_bed_overlaps(path, extent):
    """True if the main file's first row falls within the boundary extent.

    Workers normally route overlap rows to prev_region, so only the first
    row needs checking; if it overlaps, the caller splits the file.
    """
    with bgzf.open(path, "rt") as source:
        for line in source:
            if line.strip():
                chrom, start, _ = line.split("\t", 2)
                return chrom == extent[0] and int(start) <= extent[1]
    return False


def merge_and_combine_coverage_files(sample_name, sample_dir, nprocess):
    """Combine per-worker coverage files into {sample_name}_coverage.bed.gz.

    1. Sum adjacent next_region/prev_region files at matching positions
       (plus any overlapping rows split off the following main file).
    2. Concatenate BGZF blocks in order. Order holds as long as overflow
       stays within the following worker's chunk (minimum chunk lengths
       make that likely).

    Returns True if the combined file and its tabix index were written,
    False (per-worker files left in place) otherwise.
    """
    print("Merging adjacent region files and combining coverage files...")

    # Step 1: merge boundary files. A main file whose leading rows were
    # folded into an overlap is replaced by its remainder in step 2.
    remainder_files = {}
    try:
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
            main_next_file = os.path.join(
                sample_dir, f"{sample_name}_{n+1}_coverage.bed.gz"
            )

            overlap_extent = _boundary_bed_extent(next_file)
            extra_files = []
            if (
                overlap_extent is not None
                and os.path.exists(main_next_file)
                and _main_bed_overlaps(main_next_file, overlap_extent)
            ):
                peeled_file = os.path.join(
                    sample_dir, f"{sample_name}_{n+1}_coverage_peeled.tmp.bed.gz"
                )
                remainder_file = os.path.join(
                    sample_dir, f"{sample_name}_{n+1}_coverage_remainder.tmp.bed.gz"
                )
                _split_bed_by_start(
                    main_next_file, overlap_extent, peeled_file, remainder_file
                )
                extra_files.append(peeled_file)
                remainder_files[n + 1] = remainder_file

            merge_adjacent_bed_files(next_file, prev_file, overlap_file, extra_files)

            for f in extra_files:
                os.remove(f)
    except Exception as e:
        print(
            f"Warning: could not merge boundary coverage files (worker pair "
            f"{n}/{n+1}): {e}. Per-worker coverage files were left in place "
            "for diagnosis/retry.",
            file=sys.stderr,
        )
        return False

    # Step 2: list files in genomic order.
    files_to_combine = []

    for n in range(nprocess):
        main_file = remainder_files.get(
            n, os.path.join(sample_dir, f"{sample_name}_{n}_coverage.bed.gz")
        )
        if not os.path.exists(main_file):
            raise FileNotFoundError(f"Expected coverage file not found: {main_file}")
        files_to_combine.append(main_file)

        if n < nprocess - 1:
            overlap_file = os.path.join(
                sample_dir, f"{sample_name}_{n}_{n+1}_overlap_coverage.tmp.bed.gz"
            )
            if not os.path.exists(overlap_file):
                raise FileNotFoundError(
                    f"Expected overlap coverage file not found: {overlap_file}"
                )
            files_to_combine.append(overlap_file)

    # BGZF blocks can be concatenated directly; no decompression or recompression.
    final_output = os.path.join(sample_dir, f"{sample_name}_coverage.bed.gz")
    try:
        _concatenate_bgzf_files(files_to_combine, final_output)
        print(f"Combined coverage file created: {final_output}")
    except Exception as e:
        print(
            f"Warning: could not combine coverage files into {final_output}: {e}. "
            "Per-worker coverage files were left in place for diagnosis/retry.",
            file=sys.stderr,
        )
        return False
    # Keep worker files if indexing fails, so the output can be diagnosed/retried.
    try:
        subprocess.run(["tabix", "-f", "-p", "bed", final_output], check=True)
        print(f"Tabix index created for: {final_output}")
    except (subprocess.CalledProcessError, OSError) as e:
        print(
            f"Warning: could not create tabix index for {final_output}: {e}. "
            "The combined coverage file exists but is not indexed; per-worker "
            "files were left in place for diagnosis/retry.",
            file=sys.stderr,
        )
        return False
    cleanup_temp_files(sample_name, sample_dir, nprocess)
    for remainder_file in remainder_files.values():
        if os.path.exists(remainder_file):
            os.remove(remainder_file)
    return True


def _concatenate_bgzf_files(paths, output):
    """Copy compressed blocks, retaining just one terminal BGZF EOF marker."""
    eof = bgzf._bgzf_eof
    with open(output, "wb") as dest:
        for path in paths:
            with open(path, "rb") as source:
                size = os.fstat(source.fileno()).st_size
                if size >= len(eof):
                    source.seek(-len(eof), os.SEEK_END)
                    if source.read() == eof:
                        size -= len(eof)
                source.seek(0)
                while size:
                    block = source.read(min(size, 1024 * 1024))
                    if not block:
                        raise EOFError(f"Truncated BGZF file: {path}")
                    dest.write(block)
                    size -= len(block)
        dest.write(eof)


def merge_adjacent_bed_files(next_file, prev_file, output_file, extra_files=None):
    """Merge boundary bed files (and any extra_files) by summing every
    column after chrom/start/end at matching positions. The column count
    is taken from the data, not hardcoded.
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
    for path in extra_files or []:
        accumulate(path)

    # Sort start/end numerically (they are strings here).
    with bgzf.open(output_file, "wt") as f:
        for (chrom, start, end), values in sorted(
            coverage_dict.items(),
            key=lambda kv: (kv[0][0], int(kv[0][1]), int(kv[0][2])),
        ):
            f.write("\t".join([chrom, start, end] + [str(v) for v in values]) + "\n")


def cleanup_temp_files(sample_name, sample_dir, nprocess):
    """Remove per-worker coverage temp files."""
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
