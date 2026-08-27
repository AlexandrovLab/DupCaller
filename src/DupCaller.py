#!/usr/bin/env python3
import argparse
import sys
from DupCaller_sub.Caller import do_call
from DupCaller_sub.Trim import do_trim
from DupCaller_sub.Summarize import do_summarize
from DupCaller_sub.Learn import do_learn
from DupCaller_sub.AggregateProfile import do_aggregate
from DupCaller_sub.Index import do_index, do_index_dbs

# from Estimate import do_estimate
if __name__ == "__main__":
    """
    Parse Arguments
    """
    master_parser = argparse.ArgumentParser(
        description="DupCaller is a set of tool to call mutations and estimate mutation rate from ecNGS dataset"
    )
    subparsers = master_parser.add_subparsers(dest="command")

    trim_parser = subparsers.add_parser(
        "trim",
        help="Trim an ecNGS fastq file. The input should be either fastq or gzipped fastq.",
    )
    trim_parser.add_argument(
        "-i", "--fq", type=str, required=True, help="fastq file (read 1 if paired)"
    )
    trim_parser.add_argument(
        "-i2", "--fq2", type=str, required=True, help="read 2 fastq file"
    )
    trim_parser.add_argument(
        "-p", "--pattern", type=str, required=True, help="pattern of sequence barcode"
    )
    trim_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="prefix of the output fastq files",
    )
    ### Call arguments
    call_parser = subparsers.add_parser(
        "call", help="Call mutations from an aligned ecNGS bam."
    )
    call_parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="bam file of sample sequencing reads",
    )
    call_parser.add_argument(
        "-g", "--germline", type=str, help="indexed germline vcf with AF field"
    )
    call_parser.add_argument(
        "-gaf",
        "--germlineAfCutoff",
        type=float,
        help="minimum population af to exclude a germline mutation",
        default=0.001,
    )
    call_parser.add_argument(
        "-f", "--reference", type=str, required=True, help="reference fasta file"
    )
    call_parser.add_argument(
        "-o", "--output", type=str, required=True, help="prefix of the output files"
    )
    call_parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for variant calling",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX"],
    )
    call_parser.add_argument(
        "-R",
        "--regionfile",
        type=str,
        help="an inclusive bed file",
        default=None,
    )
    call_parser.add_argument(
        "-rt",
        "--regionst",
        nargs="+",
        type=str,
        help="contigs to consider for training",
        default=None,
    )
    call_parser.add_argument(
        "-p", "--threads", type=int, help="number of threads", default=1
    )
    call_parser.add_argument(
        "-ax",
        "--minMeanASXS",
        type=int,
        help="minimum mean AS-XS for a read group to be considered for calling",
        default=50,
    )
    ### Learning file locations
    call_parser.add_argument(
        "-E",
        "--errprefix",
        type=str,
        help="prefix for all six error files ({prefix}.amp.tn.srd.txt, {prefix}.amp.hp.txt, "
        "{prefix}.amp.str.txt, {prefix}.dmg.tn.txt, {prefix}.dmg.hp.txt, "
        "{prefix}.dmg.str.txt); overrides the default (output prefix)",
    )
    call_parser.add_argument(
        "-mr",
        "--mutRate",
        type=float,
        help="estimated somatic mutation rate per base",
        default=3e-10,
    )
    call_parser.add_argument(
        "-fdr",
        "--fdrThreshold",
        type=float,
        help="target per-channel FDR (max of 1/(1+LR*mutation_rate) over that channel's PASS calls, i.e. its weakest surviving call); channels above this get their LR threshold raised and mutation rate re-estimated iteratively",
        default=0.05,
    )
    call_parser.add_argument(
        "-a",
        "--pseudocount",
        type=float,
        help="regularization pseudocount added to each channel's per-channel mixture-weight MLE solve (mu), guaranteeing a root strictly between 0 and 1 without needing a channel-exclusion fallback",
        default=0.5,
    )
    call_parser.add_argument(
        "-sc",
        "--skipCoveragePass",
        action="store_true",
        help="skip round 2 (the coverage-only pass): no duplex depth / per-locus coverage.bed.gz, "
        "no duplex-family composition or by-duplex-group output files, and masked candidates "
        "that only clear the final FDR-refined threshold get no real depth. Only round 1 calling "
        "and FDR-threshold determination run; VCFs and a trimmed _stats.txt are still written. "
        "Roughly halves total runtime.",
        default=False,
    )
    call_parser.add_argument(
        "-mq",
        "--mapq",
        type=float,
        help="minumum mapq for an alignment to be considered",
        default=40,
    )

    call_parser.add_argument(
        "-n", "--normalBams", nargs="+", type=str, help="bam file of matched normal"
    )
    call_parser.add_argument("-m", "--noise", nargs="+", type=str, help="noise mask")
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

    """
    call_parser.add_argument(
        "-mnv",
        "--maxMNVlen",
        type=int,
        help="maximum length of MNV to be considered a real mutation",
        default=2,
    )
    """

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
        default=5,
    )
    call_parser.add_argument(
        "-w",
        "--windowSize",
        type=int,
        help="genomic window size when calculating rough coverage and split bam files into equal regions. Adjust for smaller panel",
        default=100000,
    )
    call_parser.add_argument(
        "-bq",
        "--minBq",
        type=int,
        help="bases with quality less than this number will be set to 6",
        default=18,
    )
    call_parser.add_argument(
        "-aq",
        "--minAltQual",
        type=float,
        help="minimum consensus quality of alt allele, if not 0, in a read group to be considered for training",
        default=60,
    )

    call_parser.add_argument(
        "--minRef",
        type=float,
        help="minimum number of ref allele, if not 0, in a read group to be considered for training",
        default=2,
    )
    call_parser.add_argument(
        "--minAlt",
        type=float,
        help="minimum number of alt allele, if not 0, in a read group to be considered for training",
        default=2,
    )
    call_parser.add_argument(
        "--naf",
        type=float,
        help="maximum VAF in matched normal for a mutation to be called",
        default=0.01,
    )
    call_parser.add_argument(
        "--minGroupAmp",
        type=int,
        help="minimum reads on each strand (F1R2/F2R1) a duplex family needs to contribute to amplification-error learning",
        default=3,
    )
    call_parser.add_argument(
        "--minGroupDmg",
        type=int,
        help="minimum reads on each strand (F1R2/F2R1) a duplex family needs to contribute to damage-error learning",
        default=3,
    )
    call_parser.add_argument(
        "--rescue",
        "-res",
        action="store_true",
        help="output discarded variants with reason in the filter field",
        default=False,
    )
    call_parser.add_argument(
        "--maxZeroQualFrac",
        "-z",
        type=float,
        help="Maximum fraction of bases in a read family that has 0 quality",
        default=0.5,
    )
    call_parser.add_argument(
        "--maxPileupDepth",
        "-pd",
        type=float,
        help="Maximum depth for samtools mpileup",
        default=1000000,
    )
    call_parser.add_argument(
        "--NanoSeqBam",
        "-nb",
        action="store_true",
        help="bam uses NanoSeq-style per-mate RB/MB tags instead of a DB tag: the duplex "
        "barcode is read as {MB}-{RB} for read 1 and {RB}-{MB} for read 2",
        default=False,
    )
    call_parser.add_argument(
        "--seed",
        type=int,
        help="RNG seed for the Monte Carlo detection-power simulation used by mutation-burden "
        "calling, so results are reproducible across repeated runs and across different -p "
        "thread counts. If not set, a random seed is generated and recorded in the run's "
        "_call_params.log for later reuse",
        default=None,
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
        type=float,
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
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX","chrY"],
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
        "-a",
        "--pseudocount",
        type=float,
        help="regularization pseudocount added to each trinuc-context's per-base "
        "EM estimate of the SBS single-read-damage (SRD) rate matrix",
        default=0.5,
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
    aggregate_parser = subparsers.add_parser(
        "aggregate", help="Aggregate learned mismatch profile from multiple samples"
    )
    aggregate_parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=str,
        help="folder where DupCallerCall results are stored",
    )
    aggregate_parser.add_argument(
        "-f",
        "--input-file",
        type=str,
        help="file with one error file prefix per line",
    )
    aggregate_parser.add_argument("-o", "--output", type=str, help="output filename")

    estimate_parser = subparsers.add_parser(
        "estimate", help="Estimate mutation rate and SBS96 from results"
    )
    estimate_parser.add_argument(
        "-i",
        "--prefix",
        type=str,
        required=True,
        help="Input prefix of results from call",
    )
    estimate_parser.add_argument(
        "-f",
        "--reference",
        type=str,
        required=True,
        help="Fasta file of reference.",
    )
    estimate_parser.add_argument(
        "-ft",
        "--refTrinuc",
        type=str,
        help="Currently unused -- -f/--reference is always required regardless of this option.",
    )
    estimate_parser.add_argument(
        "-ot",
        "--outTrinuc",
        type=str,
        help="If ref is set, output the computed trinucleotide composition file for future use",
    )
    estimate_parser.add_argument(
        "-r",
        "--regions",
        nargs="+",
        type=str,
        help="contigs to consider for trinucleotide calculation",
        default=["chr" + str(_) for _ in range(1, 23, 1)] + ["chrX"],
    )
    estimate_parser.add_argument(
        "-c",
        "--clonal",
        action="store_true",
        help="If set, mutations detected in more than one molecule will be considered as clonal mutations",
        default=False,
    )
    estimate_parser.add_argument(
        "-d",
        "--dilute",
        action="store_true",
        help="Set when sample and matched normal are from the same starting DNA material",
        default=False,
    )
    estimate_parser.add_argument(
        "-gb",
        "--genebed",
        type=str,
        help="gene bed file, if gene coverage needs to be calculated",
        default=None,
    )
    estimate_parser.add_argument(
        "-rb",
        "--reestimatebed",
        type=str,
        help="re-estimate bed file, if burden re-estimation is needed",
        default=None,
    )
    estimate_parser.add_argument(
        "--seed",
        type=int,
        help="RNG seed for the parametric-bootstrap confidence intervals on corrected "
        "mutation burden, so results are reproducible across repeated runs. If not set, "
        "a random seed is generated and recorded in the run's _estimate_params.log for "
        "later reuse",
        default=None,
    )

    summarize_parser = subparsers.add_parser(
        "summarize", help="Summarize results from multiple samples"
    )
    summarize_parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        type=str,
        required=True,
        help="folder where DupCallerCall results are stored",
    )
    summarize_parser.add_argument(
        "-o", "--output", type=str, required=True, help="output filename"
    )

    index_parser = subparsers.add_parser("index", help="Index reference genomes")
    index_parser.add_argument(
        "-f",
        "--reference",
        type=str,
        required=True,
        help="Fasta file of reference.",
    )
    index_parser.add_argument(
        "-rt",
        "--repeatTsv",
        type=str,
        required=True,
        help="PERF-format repeat tsv (chrom, start, end, motif, length, strand, num_units, motif_repeat) for the reference, e.g. produced by `PERF.core -m 1 -M <N> -u 2 -i reference.fa`. Repeat unit length and repeat count are read directly from the motif and num_units columns. Entries with unit length 1 (homopolymers) are ignored -- those are self-derived from the reference sequence instead.",
    )

    index_dbs_parser = subparsers.add_parser(
        "index-dbs",
        help="Add a DBS (dinucleotide) index to an already-indexed reference. "
        "Derived read-only from the existing .ref.h5 -- does not touch "
        "ref/tn/hp/str.h5, so it's safe to run against a reference other "
        "call/estimate jobs currently have open.",
    )
    index_dbs_parser.add_argument(
        "-f",
        "--reference",
        type=str,
        help="Fasta file of reference (must already have a .ref.h5 from `DupCaller.py index`)",
    )
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
        # Deferred import: Estimate.py pulls in matplotlib + sigProfilerPlotting,
        # which every other subcommand (call, index, etc.) would otherwise pay
        # the import cost for on every invocation despite never using them.
        from DupCaller_sub.Estimate import do_estimate

        do_estimate(args)
    elif args.command == "aggregate":
        do_aggregate(args)
    elif args.command == "index":
        do_index(args)
    elif args.command == "index-dbs":
        do_index_dbs(args)
    elif args.command == "learn":
        do_learn(args)
    else:
        # No subcommand given (args.command is None) -- argparse's own
        # "invalid choice" error already rejects anything else before
        # reaching here, so this only fires on a bare `DupCaller.py`.
        master_parser.print_help()
        sys.exit(1)
