#!/usr/bin/env python3
from collections import OrderedDict
import h5py
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pysam import VariantFile as VCF
from scipy.stats import chi2, barnard_exact


def calculate_ref_trinuc(args):
    tn_int = h5py.File(args.reference + ".tn.h5", "r")
    trinuc_count = np.zeros(97)
    for chrom in args.regions:
        trinuc_count += np.bincount(tn_int[chrom], minlength=97)
    trinuc_count_32 = trinuc_count[0:32] + trinuc_count[32:64]
    return trinuc_count_32


def poisson_confint(k, cov, alpha=0.05):
    low = chi2.ppf(alpha / 2, 2 * k) / 2
    high = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    if k == 0:
        low = 0
    return low / cov, high / cov


def plot_96(ax, trinuc_pd):
    # fig, ax = plt.subplots(figsize=(60, 10))
    SBS96_order = sorted(trinuc_pd.index, key=lambda x: (x[2:5], x[0] + x[6]))
    palette = OrderedDict(
        {
            "C>A": "#03BDEF",
            "C>G": "#010101",
            "C>T": "#E42926",
            "T>A": "#CBCACA",
            "T>C": "#A2CF63",
            "T>G": "#ECC7C5",
        }
    )
    palette_list = ["#03BDEF", "#010101", "#E42926", "#CBCACA", "#A2CF63", "#ECC7C5"]
    palette_list = np.repeat(np.array(palette_list), 16).tolist()
    ax.bar(
        SBS96_order,
        [
            trinuc_pd.loc[trinuc_pd.index == m, trinuc_pd.columns[0]].values[0]
            for m in SBS96_order
        ],
        color=palette_list,
    )
    plt.yticks(size=60)
    ax.set(ylabel=None)
    ax.get_xaxis().get_label().set_visible(False)
    for nnn, muttype in enumerate(palette.keys()):
        rect = mpatches.Rectangle(
            (-0.5 + nnn * 16, 1),
            16,
            0.2,
            color=palette[muttype],
            zorder=0,
            transform=ax.get_xaxis_transform(),
            clip_on=False,
        )
        ax.add_patch(rect)
        ax.text(
            -0.5 + nnn * 16 + 8,
            1.065,
            muttype,
            color="#FFFFFF",
            fontsize=100,
            horizontalalignment="center",
            transform=ax.get_xaxis_transform(),
            verticalalignment="center",
        )
    labels = [label.get_text() for label in ax.get_xticklabels()]
    new_labels = []
    for label in SBS96_order:
        new_label = label[0] + label[2] + label[6]
        new_labels.append(new_label)
    ax.set_xticklabels(new_labels, fontsize=30, rotation=90, weight="bold")
    ax.set_xlabel(
        "Trinucleotide Context", fontsize=100, weight="bold", fontname="Arial"
    )
    return ax


def estimate_96(trinuc_cov_by_rf, trinuc_mut_by_rf, ref_trinuc, n):
    print("........Estimating mutation rate for each trinucleotide context.......")
    trinuc_mut_cov_by_rf = np.repeat(trinuc_cov_by_rf, 3, axis=0)
    trinuc_rate = np.zeros(96)
    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)
    burden_uncorrected = np.zeros(5)
    burden_uncorrected_ub = np.zeros(5)
    burden_uncorrected_lb = np.zeros(5)
    burden_corrected = np.zeros(5)
    burden_corrected_ub = np.zeros(5)
    burden_corrected_lb = np.zeros(5)
    hap_trinuc = np.zeros([96, 5])
    ## Calculate burden when take 5 as mininum
    trinuc_mut = trinuc_mut_by_rf[:, nmin >= 5].sum(axis=1)
    trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= 5].sum(axis=1)
    mutnum = trinuc_mut.sum()
    cov = trinuc_cov.sum() / 3
    burden_uncorrected[4] = mutnum / cov
    burden_uncorrected_lb[4], burden_uncorrected_ub[4] = poisson_confint(mutnum, cov)
    trinuc_rate = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
    mutnum_corrected = trinuc_rate.dot(np.repeat(ref_trinuc, 3))
    burden_corrected[4] = mutnum_corrected / ref_trinuc.sum()
    burden_corrected_lb[4], burden_corrected_ub[4] = poisson_confint(
        mutnum_corrected, ref_trinuc.sum()
    )
    hap_trinuc[:, 4] = np.ceil(trinuc_rate * np.repeat(ref_trinuc, 3))
    # trinuc_rate[:,9] = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
    for nn in range(4, 0, -1):
        trinuc_mut = trinuc_mut + trinuc_mut_by_rf[:, nmin == nn].sum(axis=1)
        trinuc_cov = trinuc_cov + trinuc_mut_cov_by_rf[:, nmin == nn].sum(axis=1)
        mutnum = trinuc_mut.sum()
        cov = trinuc_cov.sum() / 3
        burden_uncorrected[nn - 1] = mutnum / cov
        burden_uncorrected_lb[nn - 1], burden_uncorrected_ub[nn - 1] = poisson_confint(
            mutnum, cov
        )
        trinuc_rate = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
        trinuc_rate = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
        mutnum_corrected = trinuc_rate.dot(np.repeat(ref_trinuc, 3))
        burden_corrected[nn - 1] = mutnum_corrected / ref_trinuc.sum()
        burden_corrected_lb[nn - 1], burden_corrected_ub[nn - 1] = poisson_confint(
            mutnum_corrected, ref_trinuc.sum()
        )
        hap_trinuc[:, nn - 1] = np.ceil(trinuc_rate * np.repeat(ref_trinuc, 3))

    return (
        hap_trinuc,
        burden_corrected,
        burden_corrected_lb,
        burden_corrected_ub,
        burden_uncorrected,
        burden_uncorrected_lb,
        burden_uncorrected_ub,
        mutnum,
        mutnum_corrected,
        ref_trinuc.sum(),
    )


def estimate_id(trinuc_cov_by_rf, muts_by_rf, n):
    print("........Estimating indel rate")
    trinuc_mut_cov_by_rf = np.repeat(trinuc_cov_by_rf, 3, axis=0)
    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)
    burden_indel = np.zeros(5)
    burden_indel_ub = np.zeros(5)
    burden_indel_lb = np.zeros(5)
    ## Calculate burden when take 10 as mininum
    trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= 5].sum(axis=1)
    mutnum = muts_by_rf[:, nmin >= 5].sum(axis=1).sum()
    cov = trinuc_cov.sum() / 3
    burden_indel[4] = mutnum / cov
    burden_indel_lb[4], burden_indel_ub[4] = poisson_confint(mutnum, cov)
    # trinuc_rate[:,9] = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
    for nn in range(4, 0, -1):
        trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= nn].sum(axis=1)
        mutnum = muts_by_rf[:, nmin >= nn].sum(axis=1).sum()
        cov = trinuc_cov.sum() / 3
        burden_indel[nn - 1] = mutnum / cov
        burden_indel_lb[nn - 1], burden_indel_ub[nn - 1] = poisson_confint(mutnum, cov)

    return (burden_indel, burden_indel_lb, burden_indel_ub, mutnum)


def do_estimate(args):
    ref_trinuc = calculate_ref_trinuc(args)
    prefix = args.prefix
    if len(prefix.split("/")[-1]) == 0:
        sample = prefix.split("/")[-2]
    else:
        sample = prefix.split("/")[-1]
    trinuc_by_rf = pd.read_csv(
        prefix + "/" + sample + "_trinuc_by_duplex_group.txt", sep="\t", index_col=0
    )
    ###Estimate SNV burden
    vcf = VCF(prefix + "/" + sample + "_snv.vcf", "r")
    mut_by_rf = dict()
    trinuc_list = list()
    trinuc2num = dict()
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                trinuc2num[minus_base + ref_base + plus_base] = len(trinuc_list)
                trinuc_list.append(minus_base + ref_base + plus_base)
    trinucSbs2num = dict()
    num2trinucSbs = list()
    for minus_base in ["A", "T", "C", "G"]:
        for ref_base in ["C", "T"]:
            for plus_base in ["A", "T", "C", "G"]:
                alts = ["A", "T", "C", "G"]
                alts.remove(ref_base)
                for alt_base in alts:
                    trinucSbs2num[
                        minus_base + "[" + ref_base + ">" + alt_base + "]" + plus_base
                    ] = len(num2trinucSbs)
                    num2trinucSbs.append(
                        minus_base + "[" + ref_base + ">" + alt_base + "]" + plus_base
                    )
    # rf_sizes = np.zeros(len(trinuc_by_rf.columns),dtype=int)
    # for nn,rf in enumerate(trinuc_by_rf.columns):
    # ts = int(rf.split("+")[0])
    # bs = int(rf.split("+")[1])
    # rf_sizes[nn] = ts+bs
    # sorted_index = np.argsort(rf_sizes)
    # rf_sizes_sorted = rf_sizes[sorted_index]
    trinuc_by_rf_np = trinuc_by_rf.to_numpy()
    trinuc_mut_np = np.zeros([96, len(trinuc_by_rf.columns)], dtype=int)
    duplex_no_dict = dict()
    revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for nn, duplex_no in enumerate(trinuc_by_rf.columns):
        duplex_no_dict[duplex_no] = nn

    if args.dilute:
        vcf_out = VCF(
            args.prefix + "/" + sample + "_snv_flt.vcf", "w", header=vcf.header
        )
    rec_set = set()
    for rec in vcf.fetch():
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0]
        rec_str = chrom + str(pos) + ":" + ref + ">" + alt
        if rec_str in rec_set:
            continue
        rec_set.add(rec_str)
        F1R2 = rec.info["F1R2"]
        F2R1 = rec.info["F2R1"]
        TAC = rec.samples["TUMOR"]["AC"]
        if TAC > 1 and args.dilute:
            TDP = rec.samples["TUMOR"]["DP"]
            NAC = rec.samples["NORMAL"]["AC"]
            NDP = rec.samples["NORMAL"]["DP"]
            barnard_p = barnard_exact([[TAC, TDP - TAC], [NAC, NDP - NAC]]).pvalue
            if barnard_p <= 0.05:
                continue
        if args.dilute:
            vcf_out.write(rec)
        # duplex_no = str(min(F1R2, F2R1)) + "+" + str(max(F1R2, F2R1))
        duplex_no = str(F1R2) + "+" + str(F2R1)
        ref = rec.ref
        if ref == "C" or ref == "T":
            trinuc = rec.info["TN"]
            alt = rec.alts[0]
        else:
            trinuc_revcomp = rec.info["TN"]
            trinuc = (
                revcomp[trinuc_revcomp[2]]
                + revcomp[trinuc_revcomp[1]]
                + revcomp[trinuc_revcomp[0]]
            )
            alt = revcomp[rec.alts[0]]
        trinucSbs = trinuc[0] + "[" + trinuc[1] + ">" + alt + "]" + trinuc[2]
        trinuc_mut_np[trinucSbs2num[trinucSbs], duplex_no_dict[duplex_no]] += 1
        trinuc_by_rf_np[trinuc2num[trinuc], duplex_no_dict[duplex_no]] -= 1

    base2num = {"A": 0, "T": 1, "C": 2, "G": 3}
    print("......Estimating mutational burden and SBS96 profile........")
    (
        corrected_trinuc_num,
        burden,
        burden_lb,
        burden_ub,
        uburden,
        uburden_lb,
        uburden_ub,
        mutnum_uncorrected,
        mutnum_per_genome,
        genome_cov,
    ) = estimate_96(trinuc_by_rf_np, trinuc_mut_np, ref_trinuc, trinuc_by_rf.columns)
    corrected_trinuc_pd = pd.DataFrame(
        corrected_trinuc_num[:, [0]].astype(int),
        index=num2trinucSbs,
        columns=["number"],
    )
    corrected_trinuc_pd.to_csv(
        args.prefix + "/" + sample + "_sbs_96_corrected.txt", sep="\t"
    )
    fig, ax = plt.subplots(figsize=(60, 10))
    ax = plot_96(ax, corrected_trinuc_pd)
    fig.savefig(args.prefix + "/" + sample + "_sbs_96_corrected.png", dpi=300)
    fig, ax = plt.subplots(figsize=(20, 15))
    ax.plot(range(1, 6), uburden)
    ax.plot(range(1, 6), uburden_lb)
    ax.plot(range(1, 6), uburden_ub)
    ax.set_yscale("log")
    fig.savefig(
        args.prefix + "/" + sample + "_sbs_burden_by_min_read_group_size.png", dpi=300
    )
    with open(args.prefix + "/" + sample + "_sbs_burden.txt", "w") as f:
        f.write(f"Uncorrected burden\t{uburden[0]}\n")
        f.write(f"Uncorrected burden 95% lower\t{uburden_lb[0]}\n")
        f.write(f"Uncorrected burden 95% upper\t{uburden_ub[0]}\n")
        f.write(f"Uncorrected mutation number\t{mutnum_uncorrected}\n")
        f.write(f"Corrected burden\t{burden[0]}\n")
        f.write(f"Corrected burden 95% lower\t{burden_lb[0]}\n")
        f.write(f"Corrected burden 95% upper\t{burden_ub[0]}\n")
        f.write(f"mutation number per genome\t{mutnum_per_genome}\n")
        f.write(f"genome coverage\t{genome_cov}\n")

    table = pd.DataFrame(
        np.hstack(
            [
                _.reshape(5, 1)
                for _ in [
                    np.arange(1, 6, 1, dtype=np.int16),
                    uburden,
                    uburden_lb,
                    uburden_ub,
                    burden,
                    burden_lb,
                    burden_ub,
                ]
            ]
        ),
        columns=[
            "read number",
            "Corrected_burden",
            "Corrected_burden_lower",
            "Corrected_burden_upper",
            "Uncorrrected_burden",
            "Uncorrected_burden_lower",
            "Uncorrected_burden_upper",
        ],
    )
    table.to_csv(
        args.prefix + "/" + sample + "_sbs_burden_by_min_read_group_size.txt",
        sep="\t",
        index=False,
    )
    """
    prefix = args.prefix
    if len(prefix.split("/")[-1]) == 0:
        sample = prefix.split("/")[-2]
    else:
        sample = prefix.split("/")[-1]
    trinuc_by_rf = pd.read_csv(
        prefix + "/" + sample + "_trinuc_by_duplex_group.txt", sep="\t", index_col=0
    )
    ###Estimate INDEL burden
    vcf = VCF(prefix + "/" + sample + "_indel.vcf", "r")
    mutnum_by_rf = dict()
    for n in trinuc_by_rf.columns:
        mutnum_by_rf[n] = [0]

    if args.dilute:
        vcf_out = VCF(
            args.prefix + "/" + sample + "_indel_flt.vcf", "w", header=vcf.header
        )
    for rec in vcf.fetch():
        F1R2 = rec.info["F1R2"]
        F2R1 = rec.info["F2R1"]
        TAC = rec.samples["TUMOR"]["AC"]
        if TAC > 1 and args.dilute:
            TDP = rec.samples["TUMOR"]["DP"]
            NAC = rec.samples["NORMAL"]["AC"]
            NDP = rec.samples["NORMAL"]["DP"]
            barnard_p = barnard_exact([[TAC, TDP - TAC], [NAC, NDP - NAC]]).pvalue
            if barnard_p <= 0.05:
                continue
        if args.dilute:
            vcf_out.write(rec)
        # duplex_no = str(min(F1R2, F2R1)) + "+" + str(max(F1R2, F2R1))
        duplex_no = str(F1R2) + "+" + str(F2R1)
        mutnum_by_rf[duplex_no] = [1]
    mutnum_by_rf_np = pd.DataFrame(mutnum_by_rf).to_numpy()
    (burden, burden_lb, burden_ub, mutnum) = estimate_id(
        trinuc_by_rf_np, mutnum_by_rf_np, trinuc_by_rf.columns
    )

    with open(args.prefix + "/" + sample + "_indel_burden.txt", "w") as f:
        f.write(f"Indel burden\t{burden[0]}\n")
        f.write(f"Indel burden 95% lower\t{burden_lb[0]}\n")
        f.write(f"Indel burden 95% upper\t{burden_ub[0]}\n")
        f.write(f"Indel number\t{mutnum}\n")

    table = pd.DataFrame(
        np.hstack(
            [
                _.reshape(10, 1)
                for _ in [
                    np.arange(1, 11, 1, dtype=np.int16),
                    burden,
                    burden_lb,
                    burden_ub,
                ]
            ]
        ),
        columns=[
            "read number",
            "Indel_burden",
            "Indel_burden_lower",
            "Indel_burden_upper",
        ],
    )
    table.to_csv(
        args.prefix + "/" + sample + "_indel_burden_by_min_group_size.txt",
        sep="\t",
        index=False,
    )
    """
