#!/usr/bin/env python3
from collections import OrderedDict
import h5py
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pysam import VariantFile as VCF, TabixFile
from scipy.stats import chi2, barnard_exact
import pysam


def calculate_ref_trinuc(args):
    tn_int = h5py.File(args.reference + ".tn.h5", "r")
    trinuc_count = np.zeros(97)
    for chrom in args.regions:
        trinuc_count += np.bincount(tn_int[chrom], minlength=97)
    trinuc_count_32 = trinuc_count[0:32] + trinuc_count[32:64]
    return trinuc_count_32


def poisson_confint(k, cov, alpha=0.05):
    if cov == 0:
        return float("nan"), float("nan")
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
    covs = np.zeros(5)
    hap_trinuc = np.zeros([96, 5])
    ## Calculate burden when take 5 as mininum
    trinuc_mut = trinuc_mut_by_rf[:, nmin >= 5].sum(axis=1)
    trinuc_cov = trinuc_cov_by_rf[:, nmin >= 5].sum(axis=1)
    mutnum = trinuc_mut.sum()
    cov = trinuc_cov.sum()
    burden_uncorrected[4] = mutnum / cov
    burden_uncorrected_lb[4], burden_uncorrected_ub[4] = poisson_confint(mutnum, cov)
    trinuc_rate = np.where(
        np.repeat(trinuc_cov, 3) > 0, trinuc_mut / np.repeat(trinuc_cov, 3), 0
    )
    correction_ratio = (ref_trinuc / ref_trinuc.sum()) / (trinuc_cov / trinuc_cov.sum())
    correction_ratio = np.repeat(correction_ratio, 3)
    mutnum_corrected = correction_ratio * trinuc_mut
    # burden_corrected[4] = mutnum_corrected / ref_trinuc.sum()
    burden_corrected[4] = mutnum_corrected.sum() / cov
    # burden_corrected_lb[4], burden_corrected_ub[4] = poisson_confint(
    # mutnum_corrected, ref_trinuc.sum()
    # )
    burden_corrected_lb[4], burden_corrected_ub[4] = poisson_confint(
        mutnum_corrected.sum(), cov
    )
    # hap_trinuc[:, 4] = np.ceil(trinuc_rate * np.repeat(ref_trinuc, 3))
    hap_trinuc[:, 4] = trinuc_rate * np.repeat(ref_trinuc, 3)
    covs[4] = cov
    # trinuc_rate[:,9] = np.where(trinuc_cov > 0, trinuc_mut / trinuc_cov, 0)
    for nn in range(4, 0, -1):
        trinuc_mut = trinuc_mut + trinuc_mut_by_rf[:, nmin == nn].sum(axis=1)
        trinuc_cov = trinuc_cov + trinuc_cov_by_rf[:, nmin == nn].sum(axis=1)
        mutnum = trinuc_mut.sum()
        cov = trinuc_cov.sum()
        burden_uncorrected[nn - 1] = mutnum / cov
        burden_uncorrected_lb[nn - 1], burden_uncorrected_ub[nn - 1] = poisson_confint(
            mutnum, cov
        )
        trinuc_rate = np.where(
            np.repeat(trinuc_cov, 3) > 0, trinuc_mut / np.repeat(trinuc_cov, 3), 0
        )
        covs[nn - 1] = cov
        correction_ratio = (ref_trinuc / ref_trinuc.sum()) / (
            trinuc_cov / trinuc_cov.sum()
        )
        correction_ratio = np.repeat(correction_ratio, 3)
        # mutnum_corrected = trinuc_rate.dot(np.repeat(ref_trinuc, 3))
        mutnum_corrected = trinuc_mut * correction_ratio
        # burden_corrected[4] = mutnum_corrected / ref_trinuc.sum()
        burden_corrected[nn - 1] = mutnum_corrected.sum() / cov
        # burden_corrected_lb[4], burden_corrected_ub[4] = poisson_confint(
        # mutnum_corrected, ref_trinuc.sum()
        # )
        burden_corrected_lb[nn - 1], burden_corrected_ub[nn - 1] = poisson_confint(
            mutnum_corrected.sum(), cov
        )
        hap_trinuc[:, nn - 1] = trinuc_rate * np.repeat(ref_trinuc, 3)
        # hap_trinuc[:, nn - 1] = mutnum_corrected

    return (
        hap_trinuc,
        burden_corrected,
        burden_corrected_lb,
        burden_corrected_ub,
        burden_uncorrected,
        burden_uncorrected_lb,
        burden_uncorrected_ub,
        mutnum,
        mutnum_corrected.sum(),
        ref_trinuc.sum(),
        covs,
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
    if not args.reestimatebed:
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
                            minus_base
                            + "["
                            + ref_base
                            + ">"
                            + alt_base
                            + "]"
                            + plus_base
                        ] = len(num2trinucSbs)
                        num2trinucSbs.append(
                            minus_base
                            + "["
                            + ref_base
                            + ">"
                            + alt_base
                            + "]"
                            + plus_base
                        )
        trinuc_by_rf_np = trinuc_by_rf.to_numpy()
        trinuc_mut_np = np.zeros([96, len(trinuc_by_rf.columns)], dtype=int)
        duplex_no_dict = dict()
        revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        for nn, duplex_no in enumerate(trinuc_by_rf.columns):
            duplex_no_dict[duplex_no] = nn

        # Dictionary to track unique mutations and their counts
        unique_mutations = {}

        if args.dilute:
            vcf_out = VCF(
                args.prefix + "/" + sample + "_snv_flt.vcf", "w", header=vcf.header
            )
        with open(f"{sample}/{sample}_stats.txt") as st:
            lines = st.readlines()
            unmasked_cov = int(lines[4].split("\t")[1])
            unmasked_indel_cov = int(lines[6].split("\t")[1])
            indel_cov = int(lines[5].split("\t")[1])
        ###Calculate masked mutation rate
        unmasked_mut_count = 0
        for rec in vcf.fetch():
            if "PASS" not in rec.filter and "masked" not in rec.filter:
                continue
            unmasked_mut_count += 1
        unmasked_sbs_burden = (
            float(unmasked_mut_count) / float(unmasked_cov)
            if unmasked_cov > 0
            else float("nan")
        )
        unmasked_sbs_burden_lb, unmasked_sbs_burden_ub = poisson_confint(
            unmasked_mut_count, unmasked_cov
        )
        for rec in vcf.fetch():
            count_flag = 0
            if "PASS" not in rec.filter:
                continue
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref
            alt = rec.alts[0]
            TAC = rec.samples["TUMOR"]["AC"]
            TDP = rec.samples["TUMOR"]["DP"]
            mutation_key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            if not unique_mutations.get(mutation_key):
                unique_mutations[mutation_key] = [0, TAC, TDP]
                count_flag = 1
            unique_mutations[mutation_key][0] += 1
            """
            if count_flag == 0:
                continue
            """
            F1R2 = rec.info["F1R2"]
            F2R1 = rec.info["F2R1"]
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
            cov_by_minread,
        ) = estimate_96(
            trinuc_by_rf_np, trinuc_mut_np, ref_trinuc, trinuc_by_rf.columns
        )
        corrected_trinuc_pd = pd.DataFrame(
            corrected_trinuc_num[:, [0]],  # .astype(int),
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
            args.prefix + "/" + sample + "_sbs_burden_by_min_read_group_size.png",
            dpi=300,
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
            f.write(f"Unmasked burden\t{unmasked_sbs_burden}\n")
            f.write(f"Unmasked burden 95% lower\t{unmasked_sbs_burden_lb}\n")
            f.write(f"Unmasked burden 95% upper\t{unmasked_sbs_burden_ub}\n")

        table = pd.DataFrame(
            np.hstack(
                [
                    _.reshape(5, 1)
                    for _ in [
                        np.arange(1, 6, 1, dtype=np.int16),
                        uburden,
                        uburden_lb,
                        uburden_ub,
                        cov_by_minread,
                        burden,
                        burden_lb,
                        burden_ub,
                    ]
                ]
            ),
            columns=[
                "read number",
                "Uncorrected_burden",
                "Uncorrected_burden_lower",
                "Uncorrected_burden_upper",
                "Coverage base",
                "Corrrected_burden",
                "Corrected_burden_lower",
                "Corrected_burden_upper",
            ],
        )
        table.to_csv(
            args.prefix + "/" + sample + "_sbs_burden_by_min_read_group_size.txt",
            sep="\t",
            index=False,
        )
        vcf.close()
        ###Calculate indel burdens
        print("......Estimating Indel Burden.......")
        vcf_indel = VCF(f"{sample}/{sample}_indel.vcf", "r")
        unmasked_indel_count = 0
        indel_count = 0
        for rec in vcf_indel.fetch():
            if "PASS" not in rec.filter and "masked" not in rec.filter:
                continue
            unmasked_indel_count += 1
            if "PASS" not in rec.filter:
                continue
            indel_count += 1
            mutation_key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
            TAC = rec.samples["TUMOR"]["AC"]
            TDP = rec.samples["TUMOR"]["DP"]
            if not unique_mutations.get(mutation_key):
                unique_mutations[mutation_key] = [0, TAC, TDP]
            unique_mutations[mutation_key][0] += 1
        unmasked_indel_burden = (
            unmasked_indel_count / unmasked_indel_cov
            if unmasked_indel_cov > 0
            else float("nan")
        )
        unmasked_indel_burden_lb, unmasked_indel_burden_ub = poisson_confint(
            unmasked_indel_count, unmasked_indel_cov
        )
        indel_burden = indel_count / indel_cov if indel_cov > 0 else float("nan")
        indel_burden_lb, indel_burden_ub = poisson_confint(indel_count, indel_cov)
        with open(args.prefix + "/" + sample + "_indel_burden.txt", "w") as f:
            f.write(f"Indel burden\t{indel_burden}\n")
            f.write(f"Indel burden 95% lower\t{indel_burden_lb}\n")
            f.write(f"Indel burden 95% upper\t{indel_burden_ub}\n")
            f.write(f"Indel number\t{indel_count}\n")
            f.write(f"Unmasked Indel burden\t{unmasked_indel_burden}\n")
            f.write(f"Unmasked Indel burden 95% lower\t{unmasked_indel_burden_lb}\n")
            f.write(f"Unmasked Indel burden 95% upper\t{unmasked_indel_burden_ub}\n")
        # Process unique mutations to extract duplex depths and create allele counts table
        print("......Processing unique mutations for duplex allele counts.......")
        coverage_file = f"{sample}/{sample}_coverage.bed.gz"
        allele_counts_data = []

        # Open the coverage file using TabixFile
        tbx = TabixFile(coverage_file)

        for (chrom, pos, ref, alt), counts in unique_mutations.items():
            duplex_depth = 0

            try:
                # Query the coverage file for this position
                for row in tbx.fetch(chrom, pos - 1, pos):
                    parts = row.split("\t")
                    if len(parts) >= 5:
                        # For SNVs, duplex depth is 4th column (index 3)
                        # For INDELs, duplex depth is 5th column (index 4)
                        if len(ref) == 1 and len(alt) == 1:
                            # SNV case
                            duplex_depth = int(parts[3]) if parts[3].isdigit() else 0
                        else:
                            # INDEL case
                            duplex_depth = (
                                int(parts[4])
                                if len(parts) > 4 and parts[4].isdigit()
                                else 0
                            )
                        break
                gene_name = "."
                if args.genebed:
                    for g in TabixFile(args.genebed).fetch(chrom, pos - 1, pos):
                        gene_name = g.split("\t")[3].split("_")[0]
                        break

            except Exception as e:
                print(f"Warning: Could not query coverage for {chrom}:{pos}: {e}")
                duplex_depth = 0
            if duplex_depth == 0:
                print(
                    f"duplex depth for location {chrom}:{pos} is 0, skipping the mutation..."
                )
                continue
            allele_counts_data.append(
                {
                    "chromosome": chrom,
                    "position_start": pos,
                    "ref": ref,
                    "alt": alt,
                    "count": counts[0],
                    "duplex_depth": duplex_depth,
                    "bam_alt_count": counts[1],
                    "bam_depth": counts[2],
                    "duplex_vaf": float(counts[0]) / float(duplex_depth),
                    "bam_vaf": float(counts[1]) / float(counts[2]),
                    "gene": gene_name,
                }
            )
        tbx.close()

        # Create DataFrame and write to file
        allele_counts_df = pd.DataFrame(allele_counts_data)
        allele_counts_file = args.prefix + "/" + sample + "_duplex_allele_counts.txt"
        allele_counts_df.to_csv(allele_counts_file, sep="\t", index=False)
        print(f"Duplex allele counts written to: {allele_counts_file}")
        if args.genebed:
            print(f"Calculating per-gene duplex coverage")
            gene_dict = dict()
            cds_len = dict()
            for rec in TabixFile(args.genebed).fetch():
                (chrom, start, end, gene_exon) = rec.split("\t")[0:4]
                start = int(start)
                end = int(end)
                gene = gene_exon.split("_")[0]
                gene_exon_cov = 0
                for loc in TabixFile(f"{sample}/{sample}_coverage.bed.gz").fetch(
                    chrom, start, end
                ):
                    gene_exon_cov += int(loc.split("\t")[3])
                if gene_dict.get(gene):
                    gene_dict[gene] += gene_exon_cov
                    cds_len[gene] += end - start + 1
                else:
                    gene_dict[gene] = gene_exon_cov
                    cds_len[gene] = end - start + 1
            with open(args.prefix + "/" + sample + "_gene_coverage.txt", "w") as f:
                for gene in gene_dict.keys():
                    f.write(f"{gene}\t{gene_dict[gene]/cds_len[gene]}\n")

    else:
        print("Re-estimating mutational burden from provided bed file")
        ref_trinuc = calculate_ref_trinuc(args)
        prefix = args.prefix
        if len(prefix.split("/")[-1]) == 0:
            sample = prefix.split("/")[-2]
        else:
            sample = prefix.split("/")[-1]
        tn_int = h5py.File(args.reference + ".tn.h5", "r")
        trinuc_count = np.zeros(97)
        trinucSbs_count = np.zeros(96)
        cov_total = 0
        indel_cov_total = 0
        trinucSbs2num = dict()
        num2trinucSbs = list()
        for minus_base in ["A", "T", "C", "G"]:
            for ref_base in ["C", "T"]:
                for plus_base in ["A", "T", "C", "G"]:
                    alts = ["A", "T", "C", "G"]
                    alts.remove(ref_base)
                    for alt_base in alts:
                        trinucSbs2num[
                            minus_base
                            + "["
                            + ref_base
                            + ">"
                            + alt_base
                            + "]"
                            + plus_base
                        ] = len(num2trinucSbs)
                        num2trinucSbs.append(
                            minus_base
                            + "["
                            + ref_base
                            + ">"
                            + alt_base
                            + "]"
                            + plus_base
                        )
        revcomp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        vcf = VCF(prefix + "/" + sample + "_snv.vcf", "r")
        vcf_indel = VCF(prefix + "/" + sample + "_indel.vcf", "r")
        indel_count = 0
        for interval in TabixFile(args.reestimatebed, parser=pysam.asBed()).fetch():
            for loc in TabixFile(f"{sample}/{sample}_coverage.bed.gz").fetch(
                interval.contig, interval.start, interval.end
            ):
                chrom, start, end, cov, indel_cov = loc.split("\t")[0:5]
                trinuc = tn_int[chrom][int(start)]
                trinuc_count[trinuc] += int(cov)
                cov_total += int(cov)
                indel_cov_total += int(indel_cov)
            for rec in vcf_indel.fetch():
                if "PASS" not in rec.filter:
                    continue
                if rec.chrom != interval.contig:
                    continue
                if rec.pos <= interval.start or rec.pos > interval.end:
                    continue
                indel_count += 1
            for rec in vcf.fetch():
                if "PASS" not in rec.filter:
                    continue
                if rec.chrom != interval.contig:
                    continue
                if rec.pos <= interval.start or rec.pos > interval.end:
                    continue
                chrom = rec.chrom
                pos = rec.pos
                ref = rec.ref
                alt = rec.alts[0]
                F1R2 = rec.info["F1R2"]
                F2R1 = rec.info["F2R1"]
                TAC = rec.samples["TUMOR"]["AC"]
                if TAC > 1 and args.dilute:
                    TDP = rec.samples["TUMOR"]["DP"]
                    NAC = rec.samples["NORMAL"]["AC"]
                    NDP = rec.samples["NORMAL"]["DP"]
                    barnard_p = barnard_exact(
                        [[TAC, TDP - TAC], [NAC, NDP - NAC]]
                    ).pvalue
                    if barnard_p <= 0.05:
                        continue
                if args.dilute:
                    vcf_out.write(rec)
                # duplex_no = str(min(F1R2, F2R1)) + "+" + str(max(F1R2, F2R1))
                duplex_no = str(F1R2) + "+" + str(F2R1)
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
                trinucSbs_count[trinucSbs2num[trinucSbs]] += 1
        trinuc_cov_32 = trinuc_count[0:32] + trinuc_count[32:64]
        trinuc_mut_cov = np.repeat(trinuc_cov_32, 3, axis=0)
        trinucSbs_count = trinucSbs_count.astype(float)
        trinuc_mut_cov = trinuc_mut_cov.astype(float)
        trinuc_rate = np.where(trinuc_mut_cov > 0, trinucSbs_count / trinuc_mut_cov, 0)
        corrected_trinuc_num = np.repeat(ref_trinuc, 3, axis=0) * trinuc_rate
        mutnum_corrected = corrected_trinuc_num.sum()
        ref_trinuc_sum = ref_trinuc.sum()
        burden = mutnum_corrected / ref_trinuc_sum
        burden_lb, burden_ub = poisson_confint(mutnum_corrected, ref_trinuc_sum)
        mutnum = trinucSbs_count.sum()
        uburden = mutnum / float(cov_total)
        uburden_lb, uburden_ub = poisson_confint(mutnum, float(cov_total))
        genome_cov = ref_trinuc_sum

        with open(args.prefix + "/" + sample + "_sbs_burden_re_estimate.txt", "w") as f:
            f.write(f"Uncorrected burden\t{uburden}\n")
            f.write(f"Uncorrected burden 95% lower\t{uburden_lb}\n")
            f.write(f"Uncorrected burden 95% upper\t{uburden_ub}\n")
            f.write(f"Uncorrected mutation number\t{int(mutnum)}\n")
            f.write(f"Corrected burden\t{burden}\n")
            f.write(f"Corrected burden 95% lower\t{burden_lb}\n")
            f.write(f"Corrected burden 95% upper\t{burden_ub}\n")
            f.write(f"mutation number per genome\t{mutnum_corrected}\n")
            f.write(f"genome coverage\t{genome_cov}\n")
        indel_burden_re = (
            indel_count / float(indel_cov_total)
            if indel_cov_total > 0
            else float("nan")
        )
        indel_burden_re_lb, indel_burden_re_ub = poisson_confint(
            indel_count, indel_cov_total
        )
        with open(
            args.prefix + "/" + sample + "_indel_burden_re_estimate.txt", "w"
        ) as f:
            f.write(f"Indel burden\t{indel_burden_re}\n")
            f.write(f"Indel burden 95% lower\t{indel_burden_re_lb}\n")
            f.write(f"Indel burden 95% upper\t{indel_burden_re_ub}\n")
            f.write(f"Indel number\t{indel_count}\n")
            f.write(f"Indel coverage\t{indel_cov_total}\n")
        corrected_trinuc_pd = pd.DataFrame(
            corrected_trinuc_num,
            index=num2trinucSbs,
            columns=["number"],
        )
        corrected_trinuc_pd.to_csv(
            args.prefix + "/" + sample + "_sbs_96_corrected_re_estimate.txt", sep="\t"
        )
        fig, ax = plt.subplots(figsize=(60, 10))
        ax = plot_96(ax, corrected_trinuc_pd)
        fig.savefig(
            args.prefix + "/" + sample + "_sbs_96_corrected_re_estimate.png", dpi=300
        )
