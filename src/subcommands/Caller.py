#!/usr/bin/env python3
import argparse
import os
import time
from collections import OrderedDict
from multiprocessing import Pool
import errno
import h5py


import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# import pysam
from matplotlib import pyplot as plt
from pysam import AlignmentFile as BAM

from . funcs.call import callBam
from . funcs.misc import createVcfStrings
from . funcs.misc import splitBamRegions
from . funcs.misc import getAlignmentObject as BAM
#from heapq import nlargest

#if __name__ == "__main__":
def do_call(args):
    params = {
        "tumorBam": args.bam,
        "normalBams": args.normalBams,
        "germline": args.germline,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regions,
        "threads": args.threads,
        #"amperr": args.amperrs,
        #"amperri": args.amperri,
        "amperr_file": args.output+"/"+args.output+".amp.tn.txt",
        "amperri_file": args.output+"/"+args.output+".amp.id.txt",
        #"dmgerr": args.dmgerrs,
        #"dmgerri": args.dmgerri,
        "dmgerr_file": args.output+"/"+args.output+".dmg.tn.txt",
        "dmgerri_file": args.output+"/"+args.output+".dmg.id.txt",
        "mutRate": args.mutRate,
        "pcutoff": args.threshold,
        "mapq": args.mapq,
        "noise": args.noise,
        "indel_bed": args.indelbed,
        "trim5": args.trimF,
        "trim3": args.trimR,
        "minNdepth": args.minNdepth,
        "minBq": args.minBq,
        "maxAF": args.maxAF,
        "maxMnv": args.maxMNVlen,
        "germline_cutoff": args.germlineAfCutoff,
        "maxNM": args.nmflt,
        "step": args.windowSize,
        "isLearn": None
    }

    params_learn = {
        "tumorBam": args.bam,
        "normalBams": None,
        "germline": None,
        "reference": args.reference,
        "output": args.output,
        "regions": args.regionst,
        "threads": args.threads,
        "mutRate": 10e-7,
        "pcutoff":2,
        "amperr":1e-5,
        "amperr_file":None,
        "amperri":1e-6,
        "amperri_file":None,
        "dmgerr": 1e-5,
        "dmgerri": 1e-6,
        "dmgerr_file": None,
        "dmgerri_file": None,
        "mapq": args.mapq,
        "noise": None,
        "indel_bed": None,
        "trim5": args.trimF,
        "trim3": args.trimR,
        "germline_cutoff": args.germlineAfCutoff,
        "minBq": args.minBq,
        "minAltQual": args.minAltQual,
        "maxNM": args.nmflt,
        "step": args.windowSize,
        "minRef":args.minRef,
        "minAlt":args.minAlt,
        "isLearn":True
    }
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


    """
    Learn
    """
    if not (os.path.exists(params["amperr_file"]) and os.path.exists(params["amperri_file"]) and os.path.exists(params["dmgerr_file"]) and os.path.exists(params["dmgerri_file"])):
        if args.threads == 1:
            """
            Single-thread execution
            """
            print(".........Starting estimating error rates..............")
            # contigs = [(r.strip('\n'),) for r in open(args.regions,'r').readlines()] # Only process contigs in region file
            paramsNow = params_learn
            #paramsNow["reference"] = fasta
            #paramsNow["isLearn"] = True
            regions = params["regions"]
            paramsNow["regions"] = [
                (chrom, 0, bamObject.get_reference_length(chrom) - 1) for chrom in regions
            ]
            mismatch_profile, indelerr_profile,mismatch_dmg_profile,indelerr_dmg_profile = callBam(paramsNow, 0)
        else:
            """
            Multi-thread execution
            """
            contigs = args.regionst
            contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
            print(
                "...........Spliting genomic regions for parallel execution................"
            )
            #print(args.threads)
            #if args.normalBam:
            cutSites, chunkSize, contigs = splitBamRegions(
                [args.bam], args.threads, contigs, args.windowSize
            )
            #else:
                #cutSites, chunkSize, contigs = splitBamRegions(
                    #[args.bam], args.threads, contigs, args.windowSize
                #)            
            #print(cutSites,chunkSize,contigs)# Split the whole genome for parallel execution
            regionSequence = []
            currentContigIndex = 0
            usedTime = (time.time()-startTime)/60 
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
                #paramsNow["isLearn"] = True
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
        amp_tn_pd = pd.DataFrame(mismatch_profile,columns=["A","T","C","G"],index=num2trinuc)
        dmg_tn_pd = pd.DataFrame(mismatch_dmg_profile,columns=["A","T","C","G"],index=num2trinuc)
        #np.savetxt(params["output"] + "/" + args.output + ".amp.tn.txt",np.hstack([trinuc_cols[0:32],mismatch_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")
        amp_tn_pd.to_csv(params["output"] + "/" + args.output + ".amp.tn.txt",sep='\t')
        np.savetxt(params["output"] + "/" + args.output + ".amp.id.txt",indelerr_profile,delimiter="\t",fmt="%d")
        dmg_tn_pd.to_csv(params["output"] + "/" + args.output + ".dmg.tn.txt",sep='\t')
        #np.savetxt(params["output"] + "/" + args.output + ".dmg.tn.txt",np.hstack([trinuc_cols,mismatch_dmg_profile]),delimiter="\t",header=" \tA\tT\tC\tG\n")
        np.savetxt(params["output"] + "/" + args.output + ".dmg.id.txt",indelerr_dmg_profile,delimiter="\t",fmt="%d")
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
        duplex_read_num_trinuc = OrderedDict(
            {
            num: duplex_read_num_trinuc_single[num] for num in duplex_combinations
            }
        )

    else:
        """
        Multi-thread execution
        """
        #args.threads = args.threads - 1
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        regions_list = list()
        print("....Splitting genomic regions for parallel execution.....")
        if args.normalBams:
            cutSites, chunkSize, contigs = splitBamRegions(
                [args.bam], args.threads, contigs, args.windowSize
            )
        else:
            cutSites, chunkSize, contigs = splitBamRegions(
                [args.bam], args.threads, contigs, args.windowSize
            )    
        regionSequence = []
        currentContigIndex = 0
        usedTime = (time.time()-startTime)/60   
        print(f"....Genomic regions splitted in {usedTime} minutes...")
        """
        Determine regions for each process
        """
        #print(cutSites)
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
                raise Exception(f"Window size {args.windowSize} is proabably too large. Change -w to a smaller value. An ideal value will be targetSize/(threads * 100)")
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
        
        duplex_read_num_trinuc = OrderedDict(
            {
                num: sum([d.get(num, np.zeros(32)) for d in duplex_read_nums_trinuc])
                for num in duplex_combinations
            }
        )

        
        FPAll = sum(FPs, [])
        RPAll = sum(RPs, [])



    tBam = BAM(args.bam, "rb")
    contigs = tBam.references
    # print(contigs)
    chromDict = {contig: tBam.get_reference_length(contig) for contig in contigs}
    infoDict = {
        "F1R2": [1, "Integer", "Number of F1R2 read(s) in the read bundle"],
        "F2R1": [1, "Integer", "Number of F2R1 read(s) in the read bundle"],
        #"TLR": [1, "Float", "Alt/Ref log likelihood ratio of top strand"],
        #"BLR": [1, "Float", "Alt/Ref log likelihood ratio of bottom strand"],
        "LR": [1,"Float","Log-Likelihood ratio of major base over minor base"],
        "TC": [4, "Integer", "Top strand base count"],
        "BC": [4, "Float", "Bottom strand base count"],
        "DF": [1,"Integer","Distance from fragment end"],
        "DR": [1,"Integer","Distance from read end"],
        "TAG1":[1,"String","Barcode of top strand 5 prime"],
        "TAG2":[1,"String","Barcode of bottom strand 5 prime"],
        "TN": [1,"String","trinucleotide context"],
        "HP": [1,"Integer","Homopolymer length"]
    }
    formatDict = {
        "AC": [1, "Integer", "Count of alt allele"],
        "RC": [1, "Integer", "Count of ref allele"],
        "DP": [1, "Integer", "Depth at the location"],
    }
    filterDict = {"PASS": "All filter Passed"}
    vcfLines = createVcfStrings(chromDict, infoDict, formatDict, filterDict, mutsAll)
    with open(args.output + "/" + args.output + "_snv.vcf", "w") as vcf:
        vcf.write(vcfLines)

    vcfLines = createVcfStrings(chromDict, infoDict, formatDict, filterDict, indelsAll)
    with open(args.output + "/" + args.output + "_indel.vcf", "w") as vcf:
        vcf.write(vcfLines)

    burden_naive = muts_num / (coverage)
    indel_burden = indels_num / (coverage + coverage_indel)
    efficiency = duplex_num / rec_num
    pass_duprate = unique_read_num / pass_read_num

    with open(
        params["output"] + "/" + args.output + "_duplex_group_stats.txt", "w"
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
    trinuc_list = list()
    trinuc2num = dict()
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                trinuc2num[minus_base + ref_base + plus_base] = len(trinuc_list)
                trinuc_list.append(minus_base + ref_base + plus_base)
    duplex_read_num_trinuc = {_:duplex_read_num_trinuc[_] for _ in non_zero_keys}
    trinuc_by_duplex_group = pd.DataFrame(duplex_read_num_trinuc)
    trinuc_by_duplex_group.insert(0, "", trinuc_list)
    trinuc_by_duplex_group.to_csv(params["output"] + "/" + args.output + "_trinuc_by_duplex_group.txt",sep='\t',index=False)    
    
        

    muts_by_group = np.loadtxt(
        params["output"] + "/" + args.output + "_duplex_group_stats.txt",
        skiprows=1,
        dtype=float,
        delimiter="\t",
        usecols=(2, 3),
        ndmin=2,
    ).transpose()
    x = muts_by_group[0,:]
    y = muts_by_group[1,:]
    burden_lstsq = x.dot(y) / x.dot(x)
    bootstrap_lstsqs = []
    for _ in range(10000):
        muts_by_group_resampled = np.random.default_rng().choice(
            muts_by_group, muts_by_group.shape[1], axis=1
        )
        x = muts_by_group_resampled[0,:]
        y = muts_by_group_resampled[1,:]
        burden_lstsq_resampled = x.dot(y) / x.dot(x)
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
        params["output"] + "/" + args.output + "_burden_by_duplex_group_size.png"
    )
    """
    if len(FPAll + RPAll) != 0:
        FPs_count = [0 for _ in range(max(FPAll + RPAll))]
        RPs_count = [0 for _ in range(max(FPAll + RPAll))]
        for nn in range(max(FPAll + RPAll)):
            FPs_count[nn] = FPAll.count(nn + 1)
            RPs_count[nn] = RPAll.count(nn + 1)
        with open(
            params["output"] + "/" + args.output + "_SBS_end_profile.txt", "w"
        ) as f:
            f.write("Distance\tMutations_fragment_end\tMutations_read_end\n")
            for nn in range(max(FPAll + RPAll)):
                f.write(f"{nn+1}\t{FPs_count[nn]}\t{RPs_count[nn]}\n")
        ACs = [_["samples"][0][0] for _ in mutsAll]
        ACs_clonal = [_ for _ in ACs if _ > 1]
        clonal_num = len(ACs_clonal)
    else:
        print(f"No mutations detected.")
        clonal_num = 0
    """

    with open(params["output"] + "/" + args.output + "_stats.txt", "w") as f:
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


    print(
        "..............Completed variant calling "
        + str((time.time() - startTime) / 60)
        + " minutes..............."
    )
