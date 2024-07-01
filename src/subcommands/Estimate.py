#!/usr/bin/env python3
from collections import OrderedDict
import h5py
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pysam import VariantFile as VCF


def calculate_ref_trinuc(args):
    tn_int = h5py.File(args.reference + ".tn.h5", "r")
    trinuc_count = np.zeros(97)
    for chrom in args.regions:
        trinuc_count += np.bincount(tn_int[chrom], minlength=97)
    trinuc_count_32 = trinuc_count[0:32] + trinuc_count[32:64]
    return trinuc_count_32


def estimate_96(trinuc_cov_by_rf, trinuc_mut_by_rf, ref_trinuc, n):
    print("........Estimating mutation rate for each trinucleotide context.......")
    trinuc_mut_cov_by_rf = np.repeat(trinuc_cov_by_rf, 3, axis=0)
    trinuc_rate = np.zeros(96)
    n1 = np.array([int(_.split("+")[0]) for _ in n])
    n2 = np.array([int(_.split("+")[1]) for _ in n])
    nmin = np.vstack((n1, n2)).min(axis=0)
    trinuc_rate = np.zeros([96, nmin.max()])
    burden_uncorrected = np.zeros(nmin.max())
    mutnum_uncorrected = np.zeros(nmin.max())
    for nn in range(nmin.max()):
        trinuc_mut = trinuc_mut_by_rf[:, nmin >= nn].sum(axis=1)
        trinuc_cov = trinuc_mut_cov_by_rf[:, nmin >= nn].sum(axis=1)
        trinuc_rate[:, nn] = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
        mutnum_uncorrected[nn] = trinuc_rate[:, nn].dot(
            trinuc_mut_cov_by_rf[:, nmin >= nn].sum(axis=1)
        )
        burden_uncorrected[nn] = mutnum_uncorrected[nn] / (
            trinuc_cov_by_rf[:, nmin >= nn].sum(axis=1).sum()
        )
    corrected_mutnum = trinuc_rate.T.dot(np.repeat(ref_trinuc, 3))
    burden = corrected_mutnum / ref_trinuc.sum()
    CI95 = 1.96 * np.sqrt(corrected_mutnum)
    burden_lb = burden - CI95 / ref_trinuc.sum()
    burden_ub = burden + CI95 / ref_trinuc.sum()
    hap_trinuc = np.ceil(trinuc_rate * np.repeat(ref_trinuc, 3).reshape(96, 1))

    CI95_uncorrected = 1.96 * np.sqrt(mutnum_uncorrected)
    burden_uncorrected_lb = burden_uncorrected - CI95_uncorrected / ref_trinuc.sum()
    burden_uncorrected_ub = burden_uncorrected + CI95_uncorrected / ref_trinuc.sum()
    return (
        hap_trinuc[:, 2],
        burden[2],
        burden_lb[2],
        burden_ub[2],
        trinuc_rate[:, 2],
        burden_uncorrected[2],
        burden_uncorrected_lb[2],
        burden_uncorrected_ub[2],
    )


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
    for rec in vcf.fetch():
        F1R2 = rec.info["F1R2"]
        F2R1 = rec.info["F2R1"]
        duplex_no = str(min(F1R2, F2R1)) + "+" + str(max(F1R2, F2R1))
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
        trinuc_rate,
        uburden,
        uburden_lb,
        uburden_ub,
    ) = estimate_96(trinuc_by_rf_np, trinuc_mut_np, ref_trinuc, trinuc_by_rf.columns)
    corrected_trinuc_pd = pd.DataFrame(corrected_trinuc_num, index=num2trinucSbs)
    fig, ax = plt.subplots(figsize=(60, 10))
    SBS96_order = sorted(corrected_trinuc_pd.index, key=lambda x: (x[2:5], x[0] + x[6]))
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
            corrected_trinuc_pd.loc[
                corrected_trinuc_pd.index == m, corrected_trinuc_pd.columns[0]
            ].values[0]
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
    fig.savefig(args.prefix + "/" + sample + "_sbs_96_corrected.png", dpi=300)

    # corrected_burden = corrected_trinuc_num.sum()/ref_trinuc.sum()
    # original_burden = trinuc_rate
    # print(corrected_burden)
    """
    fig, axs = plt.subplots(16, 6)
    for nn in range(96):
        nnn = math.floor(nn / 6)
        mmm = nn - nnn * 6
        axs[nnn, mmm].scatter(
            trinuc_by_rf_np[int(np.floor(nn / 3)), :], trinuc_mut_np[nn, :]
        )
        x = np.linspace(0, trinuc_mut_np[nn, :].max() * 1.1)
        axs[nnn, mmm].plot(x, trinuc_rate[nn] * x, color="r")
        axs[nnn, mmm].set_title(num2trinucSbs[nn])
    fig.savefig(args.prefix + "/" + sample + "_sbs96_mutrates.png", dpi=300)
    """
    with open(args.prefix + "/" + sample + "_burden.txt", "w") as f:
        f.write(f"Uncorrected burden\t{uburden}\n")
        f.write(f"Uncorrected burden 95% lower\t{uburden_lb}\n")
        f.write(f"Uncorrected burden 95% upper\t{uburden_ub}\n")
        f.write(f"Corrected burden\t{burden}\n")
        f.write(f"Corrected burden 95% lower\t{burden_lb}\n")
        f.write(f"Corrected burden 95% upper\t{burden_ub}\n")
