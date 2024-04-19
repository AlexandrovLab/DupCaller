#!/usr/bin/env python3
import argparse
import os
import time
from collections import OrderedDict
from multiprocessing import Pool
import errno


import matplotlib.patches as mpatches
import numpy as np

# import pysam
from matplotlib import pyplot as plt
from pysam import AlignmentFile as BAM

from DCutils.call import callBam
from DCutils.funcs import createVcfStrings
from DCutils.funcs import splitBamRegions
from DCutils.funcs import getAlignmentObject

if __name__ == "__main__":
    """
    Parse Arguments
    """
    parser = argparse.ArgumentParser(description="Call any mutations from a bam file")
    parser.add_argument(
        "-b", "--bam", type=str, help="bam file of sample sequencing reads"
    )
    parser.add_argument(
        "-g", "--germline", type=str, help="indexed germline vcf with AF field"
    )
    parser.add_argument(
        "-gaf",
        "--germlineAfCutoff",
        type=str,
        help="minimum population af to exclude a germline mutation",
        default=0.001,
    )
    parser.add_argument("-f", "--reference", type=str, help="reference fasta file")
    parser.add_argument("-o", "--output", type=str, help="prefix of the output files")
    parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for variant calling",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX", "chrY"],
    )
    parser.add_argument(
        "-p", "--threads", type=int, help="number of threads", default=1
    )
    parser.add_argument(
        "-aes",
        "--amperrs",
        type=float,
        help="estimated polymerase substitutionerror rate",
        default=1e-5,
    )
    parser.add_argument(
        "-aei",
        "--amperri",
        type=float,
        help="estimated polymerase indel error rate",
        default=3e-7,
    )
    parser.add_argument(
        "-mr",
        "--mutRate",
        type=float,
        help="estimated somatic mutation rate per base",
        default=2.5e-7,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        help="log likelihood ratio threshold of making a mutation call",
        default=2,
    )
    parser.add_argument(
        "-mq",
        "--mapq",
        type=float,
        help="minumum mapq for an alignment to be considered",
        default=30,
    )
    # parser.add_argument('-da','--damage',type=float,default=5E-7)

    parser.add_argument(
        "-n", "--normalBam", type=str, help="bam file of matched normal"
    )
    parser.add_argument("-m", "--noise", type=str, help="noise mask")
    parser.add_argument(
        "-tt",
        "--trimF",
        type=int,
        help="ignore mutation if it is less than n bps from ends of template",
        default=8,
    )
    parser.add_argument(
        "-tr",
        "--trimR",
        type=int,
        help="ignore mutation if it is less than n bps from ends of read",
        default=8,
    )
    parser.add_argument(
        "-d",
        "--minNdepth",
        type=int,
        help="minumum coverage in normal for called variants",
        default=10,
    )
    parser.add_argument(
        "-maf",
        "--maxAF",
        type=float,
        help="maximum allele fraction to call a somatic mutation",
        default=1,
    )
    parser.add_argument(
        "-mac",
        "--maxAltCount",
        type=int,
        help="maximum allele count of alt allele",
        default=100000,
    )
    parser.add_argument(
        "-mnv",
        "--maxMNVlen",
        type=int,
        help="maximum length of MNV to be considered a real mutation",
        default=2,
    )
    parser.add_argument(
        "-x",
        "--panel",
        type=bool,
        help="enable panel region splitting. Must be set if the input data is from small targeted panel",
        default=False,
    )

    parser.add_argument(
        "-id",
        "--indelbed",
        type=str,
        help="noise bed file for indels",
        default=False,
    )
    args = parser.parse_args()

    ## If the user points to an absolute output path, it will make a giant path under tmp/
    ## Quick check if params["output"] starts is an absolute path pointing to the working directory, and changes it to a non-absolute path. 
    ## Otherwise, keeps the large path.
    cwd = os.getcwd()
    if args.output.startswith(cwd):
        abs_path = args.output
        outpath = abs_path.replace(cwd, "")
        ## Make sure we don't end up with a root path
        if outpath.startswith("/"):
            outpath = outpath[1:]
    else:
        outpath = args.output

    """
    Store Parameters
    """
    params = {
        "tumorBam": args.bam,
        "normalBam": args.normalBam,
        "germline": args.germline,
        "reference": args.reference,
        "output": outpath,
        "regions": args.regions,
        "threads": args.threads,
        "amperr": args.amperrs,
        "amperri": args.amperri,
        "mutRate": args.mutRate,
        "pcutoff": args.threshold,
        "mapq": args.mapq,
        "noise": args.noise,
        "indel_bed": args.indelbed,
        "trim5": args.trimF,
        "trim3": args.trimR,
        "minNdepth": args.minNdepth,
        "maxAltCount": args.maxAltCount,
        "maxAF": args.maxAF,
        "maxMnv": args.maxMNVlen,
        "germline_cutoff": args.germlineAfCutoff,
    }
    if not params["normalBam"]:
        print(
            f"A matched normal is not used. \
            The maximum allele fraction to call a somatic mutation is set to be {args.maxAF}"
        )
    else:
        print(
            f"Matched normal: {args.normalBam}. \
            The maximum allele fraction to call a somatic mutation is set to be {args.maxAF}"
        )
    """
    Initialze run
    """
    # print("..............Loading reference genome.....................")
    # fasta = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    startTime = time.time()
    if not os.path.exists(os.path.join("tmp", outpath)):
        try:
            os.makedirs(os.path.join("tmp", outpath))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise 
    if not os.path.exists(outpath):
        try:
            os.makedirs(outpath)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise 
    bamObject = getAlignmentObject(args.bam, args.reference)
    """
    Execulte variant calling
    """
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
            indelsAll,
            coverage_indel,
            unique_read_num,
            pass_read_num,
            FPAll,
            RPAll,
        ) = callBam(paramsNow, 0, 1000000)
        muts_positions = [
            mut["chrom"] + str(mut["pos"]) + mut["ref"] + mut["alt"] for mut in mutsAll
        ]
        muts_dict = dict()
        take_ind = list()
        """
        for nnn, mut in enumerate(muts_positions):
            if muts_dict.get(mut) is None:
                take_ind.append(nnn)
                muts_dict[mut] = 1
        mutsAll_new = [mutsAll[ind] for ind in take_ind]
        mutsAll = copy.copy(mutsAll_new)
        """
        muts_num = len(mutsAll)
        indels_positions = [
            indel["chrom"] + str(indel["pos"]) + indel["ref"] + ":" + indel["alt"]
            for indel in indelsAll
        ]
        """
        indels_dict = dict()
        take_ind = list()
        for nnn, indel in enumerate(indels_positions):
            if indels_dict.get(indel) is None:
                take_ind.append(nnn)
                indels_dict[mut] = 1
        indelsAll_new = [indelsAll[ind] for ind in take_ind]
        indelsAll = copy.copy(indelsAll_new)
        """
        indels_num = len(indelsAll)
        duplex_combinations = list(duplex_read_num_single.keys())
        duplex_combinations.sort()
        duplex_read_num = OrderedDict(
            {num: duplex_read_num_single[num][0] for num in duplex_combinations}
        )
        duplex_coverage_by_group = OrderedDict(
            {num: duplex_read_num_single[num][1] for num in duplex_combinations}
        )

    else:
        """
        Multi-thread execution
        """
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        print("....Spliting genomic regions for parallel execution.....")
        if not args.panel:
            cutSites, chunkSize = splitBamRegions(
                args.bam, args.threads, contigs
            )  # Split the whole genome for parallel execution
        else:
            cutSites, chunkSize = splitBamRegions(
                args.bam, args.threads, contigs, fast=False
            )
        regionSequence = []
        currentContigIndex = 0

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
            chroms = [region[0] for region in regions]
            # fastaNow = {
            # chrom: fasta[chrom] for chrom in chroms
            # }  # Takes partial fasta to reduce memory usage
            # paramsNow = None
            paramsNow = params.copy()
            # paramsNow["reference"] = fastaNow
            paramsNow["regions"] = regions
            callArgument = (paramsNow, nn, chunkSize)
            callArguments.append(callArgument)
            regions = []
        del bamObject
        results = pool.starmap(callBam, callArguments)
        muts = [_[0] for _ in results]
        coverages = [_[1] for _ in results]
        rec_nums = [_[2] for _ in results]
        duplex_nums = [_[3] for _ in results]
        duplex_read_nums = [_[4] for _ in results]
        indels = [_[5] for _ in results]
        coverages_indels = [_[6] for _ in results]
        unique_read_nums = [_[7] for _ in results]
        pass_read_nums = [_[8] for _ in results]
        FPs = [_[9] for _ in results]
        RPs = [_[10] for _ in results]
        print(
            "..............Completed bam calling in "
            + str((time.time() - startTime2) / 60)
            + " minutes,merging results................."
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
        """
        for nnn, mut in enumerate(muts_positions):
            if muts_dict.get(mut) is None:
                take_ind.append(nnn)
                muts_dict[mut] = 1
        mutsAll_new = [mutsAll[ind] for ind in take_ind]
        mutsAll = copy.copy(mutsAll_new)
        """
        muts_num = len(mutsAll)
        coverage = sum(coverages)
        coverage_indel = sum(coverages_indels)
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
        """
        indels_dict = dict()
        take_ind = list()
        for nnn, indel in enumerate(indels_positions):
            if indels_dict.get(indel) is None:
                take_ind.append(nnn)
                indels_dict[mut] = 1
        indelsAll_new = [indelsAll[ind] for ind in take_ind]
        indelsAll = indelsAll_new
        """

        duplex_combinations = list(
            set.union(*[set(d.keys()) for d in duplex_read_nums])
        )
        duplex_combinations.sort()
        duplex_read_num = OrderedDict(
            {
                num: sum([d.get(num, [0, 0])[0] for d in duplex_read_nums])
                for num in duplex_combinations
            }
        )
        duplex_coverage_by_group = OrderedDict(
            {
                num: sum([d.get(num, [0, 0])[1] for d in duplex_read_nums])
                for num in duplex_combinations
            }
        )
        FPAll = sum(FPs, [])
        RPAll = sum(RPs, [])

        coverage_beds = [
            f"tmp/{outpath}_{_}_coverage.bed.gz" for _ in range(args.threads)
        ]
        coverage_bed = f"{outpath}_coverage.bed.gz"

        # with gzopen(coverage_bed, "w") as outfile:
        # for bed_file in coverage_beds:
        # with gzopen(bed_file, "r") as infile:
        # shutil.copyfileobj(infile, outfile)
        # subprocess.run(["cat"]+ coverage_beds+[f"> {coverage_bed}"])
        # subprocess.run(["rm"] + coverage_beds)
        # subprocess.run(["tabix",coverage_bed])

    tBam = getAlignmentObject(args.bam, args.reference)
    contigs = tBam.references
    # print(contigs)
    chromDict = {contig: tBam.get_reference_length(contig) for contig in contigs}
    infoDict = {
        "F1R2": [1, "Integer", "Number of F1R2 read(s) in the read bundle"],
        "F2R1": [1, "Integer", "Number of F2R1 read(s) in the read bundle"],
        "TG": [1, "Float", "Alt/Ref log likelihood ratio of top strand"],
        "BG": [1, "Float", "Alt/Ref log likelihood ratio of bottom strand"],
        "TC": [4, "Integer", "Top strand base count"],
        "BC": [4, "Float", "Bottom strand base count"],
    }
    formatDict = {
        "AC": [1, "Integer", "Count of alt allele"],
        "RC": [1, "Integer", "Count of ref allele"],
        "DP": [1, "Integer", "Depth at the location"],
    }
    filterDict = {"PASS": "All filter Passed"}
    vcfLines = createVcfStrings(chromDict, infoDict, formatDict, filterDict, mutsAll)
    with open(outpath + "_snv.vcf", "w") as vcf:
        vcf.write(vcfLines)

    vcfLines = createVcfStrings(chromDict, infoDict, formatDict, filterDict, indelsAll)
    with open(outpath + "_indel.vcf", "w") as vcf:
        vcf.write(vcfLines)

    burden_naive = muts_num / (coverage)
    indel_burden = indels_num / (coverage + coverage_indel)
    efficiency = duplex_num / rec_num
    pass_duprate = unique_read_num / pass_read_num

    with open(
        outpath + "_duplex_group_stats.txt", "w"
    ) as f:
        f.write(
            "duplex_group_strand_composition\tduplex_group_number\t\
            effective_coverage\tmutation_count\n"
        )
        muts_by_duplex_group = OrderedDict()
        non_zero_keys = []
        for read_num in duplex_read_num.keys():
            if duplex_read_num[read_num] != 0:
                non_zero_keys.append(read_num)
            muts_by_duplex_group[read_num] = 0
        for mut in mutsAll:
            TC_total = int(mut["infos"]["F1R2"])
            BC_total = int(mut["infos"]["F2R1"])
            if (
                muts_by_duplex_group.get(str(TC_total) + "+" + str(BC_total))
                is not None
            ):
                muts_by_duplex_group[str(TC_total) + "+" + str(BC_total)] += 1
            else:
                muts_by_duplex_group[str(BC_total) + "+" + str(TC_total)] += 1
        for read_num in non_zero_keys:
            f.write(
                f"{read_num}\t{duplex_read_num[read_num]}\t\
                {duplex_coverage_by_group[read_num]}\t{muts_by_duplex_group[read_num]}\n"
            )

    muts_by_group = np.loadtxt(
        outpath + "_duplex_group_stats.txt",
        skiprows=1,
        dtype=float,
        delimiter="\t",
        usecols=(2, 3),
    ).transpose()
    burden_lstsq = np.linalg.lstsq(
        np.vstack([muts_by_group[0, :], np.ones(muts_by_group.shape[1])]).T,
        muts_by_group[1, :],
        rcond=None,
    )[0][0]
    bootstrap_lstsqs = []
    for _ in range(10000):
        muts_by_group_resampled = np.random.default_rng().choice(
            muts_by_group, muts_by_group.shape[1], axis=1
        )
        burden_lstsq_resampled = np.linalg.lstsq(
            np.vstack(
                [muts_by_group_resampled[0, :], np.ones(muts_by_group.shape[1])]
            ).T,
            muts_by_group_resampled[1, :],
            rcond=None,
        )[0][0]
        bootstrap_lstsqs.append(burden_lstsq_resampled)
    bootstrap_lstsqs.sort()
    burden_lstsq_lci = bootstrap_lstsqs[250]
    burden_lstsq_uci = bootstrap_lstsqs[9750]
    x = np.linspace(0, muts_by_group[0, :].max() * 1.1)
    plt.scatter(muts_by_group[0, :], muts_by_group[1, :])
    plt.plot(x, burden_lstsq * x, color="r")
    plt.plot(x, burden_lstsq_lci * x, color="r", linestyle="dashed")
    plt.plot(x, burden_lstsq_uci * x, color="r", linestyle="dashed")
    plt.plot(x, burden_naive * x, color="b")
    plt.xlabel("Coverage")
    plt.ylabel("Mutation Count")
    lgd1 = mpatches.Patch(color="red", label="Least Square")
    lgd2 = mpatches.Patch(color="blue", label="Naive")
    plt.legend(handles=[lgd1, lgd2])
    plt.savefig(
        outpath + "_burden_by_duplex_group_size.png"
    )
    if len(FPAll + RPAll) != 0:
        FPs_count = [0 for _ in range(max(FPAll + RPAll))]
        RPs_count = [0 for _ in range(max(FPAll + RPAll))]
        for nn in range(max(FPAll + RPAll)):
            FPs_count[nn] = FPAll.count(nn + 1)
            RPs_count[nn] = RPAll.count(nn + 1)
        with open(
            outpath + "_SBS_end_profile.txt", "w"
        ) as f:
            f.write("Distance\tMutations_fragment_end\tMutations_read_end\n")
            for nn in range(max(FPAll + RPAll)):
                f.write(f"{nn+1}\t{FPs_count[nn]}\t{RPs_count[nn]}\n")
        ACs = [_["samples"][0][0] for _ in mutsAll]
        """
        plt.figure()
        plt.hist(FPAll, bins=range(0, max(FPAll)))
        plt.savefig(params["output"] + "/" + args.output + "_fragment_end_distance.png")
        plt.figure()
        plt.hist(RPAll, bins=range(0, max(RPAll)))
        plt.savefig(params["output"] + "/" + args.output + "_read_end_distance.png")
        plt.figure()
        plt.hist(ACs, bins=range(0, max(ACs)))
        plt.savefig(params["output"] + "/" + args.output + "_alt_read_count.png")
        """
        ACs_clonal = [_ for _ in ACs if _ > 1]
        clonal_num = len(ACs_clonal)
    else:
        print(f"No mutations detected.")
        clonal_num = 0

    with open(outpath + "_stats.txt", "w") as f:
        f.write(f"Number of Read Families\t{unique_read_num}\n")
        f.write(f"Number of Pass-filter Reads\t{pass_read_num}\n")
        f.write(f"Number of Effective Read Families\t{duplex_num}\n")
        f.write(f"Effective Coverage\t{coverage}\n")
        f.write(f"Per Read Family Coverage \t{coverage/duplex_num}\n")
        f.write(
            f"Pass-filter Duplication Rate\t\
        {1-unique_read_num/pass_read_num}\n"
        )
        f.write(f"Efficiency\t{efficiency}\n")

    with open(outpath + "_snv_burden.txt", "w") as f:
        f.write(f"Number of Mutations\t{muts_num}\n")
        f.write(f"Number of multi-clonal Mutations\t{clonal_num}\n")
        f.write(f"Estimated Naive Burden\t{burden_naive}\n")
        f.write(f"Estimated Least-square Burden\t{burden_lstsq}\n")
        f.write(f"Least-square Burden Upper 95% CI\t{burden_lstsq_uci}\n")
        f.write(f"Least-square Burden Lower 95% CI\t{burden_lstsq_lci}\n")

    with open(outpath + "_indel_burden.txt", "w") as f:
        f.write(f"Number of Mutations\t{indels_num}\n")
        f.write(f"Effective Coverage\t{coverage+coverage_indel}\n")
        f.write(f"Estimated Naive Burden\t{indel_burden}\n")

    print(
        "..............Completed variant calling "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )
