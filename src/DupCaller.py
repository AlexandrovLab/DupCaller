#!/usr/bin/env python3
import argparse
import copy
from subcommands.Caller import do_call
from subcommands.Trim import do_trim
from subcommands.Summarize import do_summarize
from subcommands.Learn import do_learn
from subcommands.AggregateProfile import do_aggregate
from subcommands.Estimate import do_estimate
from subcommands.Index import do_index
#from Estimate import do_estimate
if __name__ == "__main__":
    """
    Parse Arguments
    """
    master_parser = argparse.ArgumentParser(description="DupCaller is a set of tool to call mutations and estimate mutation rate from ecNGS dataset")
    subparsers = master_parser.add_subparsers(dest="command")

    trim_parser = subparsers.add_parser("trim",help="Trim an ecNGS fastq file. The input should be either fastq or gzipped fastq.")
    trim_parser.add_argument("-i", "--fq", type=str, help="fastq file (read 1 if paired)")
    trim_parser.add_argument("-i2", "--fq2", type=str, help="read 2 fastq file")
    trim_parser.add_argument("-p", "--pattern", type=str, help="pattern of sequence barcode")
    trim_parser.add_argument(
        "-o", "--output", type=str, help="prefix of the output fastq files"
    )
    ### Call arguments
    call_parser = subparsers.add_parser("call",help="Call mutations from an aligned ecNGS bam.")
    call_parser.add_argument(
        "-b", "--bam", type=str, help="bam file of sample sequencing reads"
    )
    call_parser.add_argument(
        "-g", "--germline", type=str, help="indexed germline vcf with AF field"
    )
    call_parser.add_argument(
        "-gaf",
        "--germlineAfCutoff",
        type=str,
        help="minimum population af to exclude a germline mutation",
        default=0.001,
    )
    call_parser.add_argument("-f", "--reference", type=str, help="reference fasta file")
    call_parser.add_argument("-o", "--output", type=str, help="prefix of the output files")
    call_parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for variant calling",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX"],
    )
    call_parser.add_argument(
        "-rt",
        "--regionst",
        nargs="+",
        type=str,
        help="contigs to consider for training",
        default=["chr1"]
    )
    call_parser.add_argument(
        "-p", "--threads", type=int, help="number of threads", default=1
    )
    ### Learning file locations
    """
    call_parser.add_argument(
        "-as",
        "--amperrs",
        type=float,
        help="estimated polymerase substitutionerror rate",
        default=1e-5,
    )
    call_parser.add_argument(
        "-ai",
        "--amperri",
        type=float,
        help="estimated polymerase indel error rate",
        default=5e-7,
    )
    call_parser.add_argument(
        "-AS",
        "--amperrfile",
        type=str,
        help="amplification error matrix",
    )
    call_parser.add_argument(
        "-AI",
        "--amperrfileindel",
        type=str,
        help="amplification error matrix",
    )
    call_parser.add_argument(
        "-ds",
        "--dmgerrs",
        type=float,
        help="estimated polymerase substitutionerror rate",
        default=1e-5,
    )
    call_parser.add_argument(
        "-di",
        "--dmgerri",
        type=float,
        help="estimated polymerase indel error rate",
        default=3e-7,
    )
    call_parser.add_argument(
        "-DS",
        "--dmgerrfile",
        type=str,
        help="amplification error matrix",
    )
    call_parser.add_argument(
        "-DI",
        "--dmgerrfileindel",
        type=str,
        help="amplification error matrix",
    )
    """
    call_parser.add_argument(
        "-mr",
        "--mutRate",
        type=float,
        help="estimated somatic mutation rate per base",
        default=2.5e-7,
    )
    call_parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        help="log likelihood ratio threshold of making a mutation call",
        default=3.5,
    )
    call_parser.add_argument(
        "-mq",
        "--mapq",
        type=float,
        help="minumum mapq for an alignment to be considered",
        default=40,
    )

    call_parser.add_argument(
        "-n", "--normalBams", nargs="+",type=str, help="bam file of matched normal"
    )
    call_parser.add_argument("-m", "--noise", type=str, help="noise mask")
    call_parser.add_argument(
        "-tt",
        "--trimF",
        type=int,
        help="ignore mutation if it is less than n bps from ends of template",
        default=7,
    ) 
    call_parser.add_argument(
        "-tr",
        "--trimR",
        type=int,
        help="ignore mutation if it is less than n bps from ends of read",
        default=7,
    )
    call_parser.add_argument(
        "-d",
        "--minNdepth",
        type=int,
        help="minumum coverage in normal for called variants",
        default=10,
    )
    call_parser.add_argument(
        "-maf",
        "--maxAF",
        type=float,
        help="maximum allele fraction to call a somatic mutation",
        default=1,
    )
    call_parser.add_argument(
        "-mnv",
        "--maxMNVlen",
        type=int,
        help="maximum length of MNV to be considered a real mutation",
        default=2,
    )

    call_parser.add_argument(
        "-id",
        "--indelbed",
        type=str,
        help="noise bed file for indels",
        default=False,
    )
    call_parser.add_argument(
        "-nm",
        "--nmflt",
        type=int,
        help="if set to a number, any read group and half of reads has a higher NM will be filtered",
        default = 4,
    )
    call_parser.add_argument(
        "-w",
        "--windowSize",
        type=int,
        help="genomic window size when calculating rough coverage and split bam files into equal regions. Adjust for smaller panel",
        default = 100000
    ) 
    call_parser.add_argument(
        "-bq",
        "--minBq",
        type=int,
        help="bases with quality less than this number will be set to 6",
        default = 18
    )
    call_parser.add_argument(
        "-aq",
        "--minAltQual",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 60,
    ) 

    call_parser.add_argument(
        "--minRef",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 2,
    )
    call_parser.add_argument(
        "--minAlt",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 2,
    ) 

    ###########
    """
    Learn Arguments
    """
    ###########
    """
    learn_parser = subparsers.add_parser("learn",help="Call mutations from an aligned ecNGS bam.")
    learn_parser.add_argument(
        "-b", "--bam", type=str, help="bam file of sample sequencing reads"
    )
    learn_parser.add_argument(
        "-g", "--germline", type=str, help="indexed germline vcf with AF field"
    )
    learn_parser.add_argument(
        "-gaf",
        "--germlineAfCutoff",
        type=str,
        help="minimum population af to exclude a germline mutation",
        default=0.001,
    )
    learn_parser.add_argument("-f", "--reference", type=str, help="reference fasta file")
    learn_parser.add_argument("-o", "--output", type=str, help="prefix of the output files")
    learn_parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for variant calling",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX"],
    )
    learn_parser.add_argument(
        "-p", "--threads", type=int, help="number of threads", default=1
    )
    learn_parser.add_argument(
        "-mq",
        "--mapq",
        type=float,
        help="minumum mapq for an alignment to be considered",
        default=50,
    )
    # parser.add_argument('-da','--damage',type=float,default=5E-7)

    learn_parser.add_argument(
        "-n", "--normalBams",nargs="+",type=str, help="bam file of matched normal"
    )
    learn_parser.add_argument("-m", "--noise", type=str, help="noise mask")
    learn_parser.add_argument(
        "-tf",
        "--trimF",
        type=int,
        help="ignore mutation if it is less than n bps from ends of template",
        default=8,
    )
    learn_parser.add_argument(
        "-tr",
        "--trimR",
        type=int,
        help="ignore mutation if it is less than n bps from ends of read",
        default=8,
    )
    learn_parser.add_argument(
        "-d",
        "--minNdepth",
        type=int,
        help="minumum coverage in normal for called variants",
        default=10,
    )

    learn_parser.add_argument(
        "-id",
        "--indelbed",
        type=str,
        help="noise bed file for indels",
        default=False,
    )

    learn_parser.add_argument(
        "-rq",
        "--minBq",
        type=float,
        help="minimum base quality to be considered for training",
        default = 18,    
    )

    learn_parser.add_argument(
        "-aq",
        "--minAltQual",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 60,
    ) 

    learn_parser.add_argument(
        "-nm",
        "--nmflt",
        type=float,
        help="",
        default = 4,
    ) 

    learn_parser.add_argument(
        "-w",
        "--windowSize",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 100000,
    )
    learn_parser.add_argument(
        "--minRef",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 2,
    )
    learn_parser.add_argument(
        "--minAlt",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default = 2,
    ) 
    """
    aggregate_parser = subparsers.add_parser("aggregate",help="Aggregate learned mismatch profile from multiple samples")
    aggregate_parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=str,
        help="folder where DupCallerCall results are stored",
    )
    aggregate_parser.add_argument("-o", "--output", type=str, help="output filename")
    
    estimate_parser = subparsers.add_parser("estimate",help="Estimate mutation rate and SBS96 from results")
    estimate_parser.add_argument("-i","--prefix",type=str,help="Input prefix of results from call")
    estimate_parser.add_argument("-f","--reference",type=str,help="Fasta file of reference. Either -f or -ft should be set")
    estimate_parser.add_argument("-ft","--refTrinuc",type=str,help="Precomputed trinucleotide composition of reference genome. Either -f or -ft should be set")
    estimate_parser.add_argument("-ot","--outTrinuc",type=str,help="If ref is set, output the computed trinucleotide composition file for future use")
    estimate_parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for trinucleotide calculation",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX"],
    )

    summarize_parser = subparsers.add_parser("summarize",help="Summarize results from multiple samples")
    summarize_parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=str,
        help="folder where DupCallerCall results are stored",
    )
    summarize_parser.add_argument("-o", "--output", type=str, help="output filename")
    

    index_parser = subparsers.add_parser("index",help="Index reference genomes")
    index_parser.add_argument("-f","--reference",type=str,help="Fasta file of reference. Either -f or -ft should be set")
    args = master_parser.parse_args()

    """
    Store Parameters
    """

    if args.command == "trim":
        do_trim(args)
    elif args.command == "call":
        do_call(args)
    elif args.command == "summarize":
        do_summarize(args)
    elif args.command == "estimate":
        do_estimate(args)
    elif args.command == "aggregate":
        do_aggregate(args)
    elif args.command == "index":
        do_index(args)
    else:
        do_learn(args)