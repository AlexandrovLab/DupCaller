#!/usr/bin/env python3
# if __name__ == "__main__":
import pandas as pd
import os


def do_summarize(args):
    samples = args.input
    with open(args.output, "w") as output:
        output.write(
            f"sample\tpass_filter_reads\tunique_reads\tread_families\tduplication_rate\tread_family_efficiency\t"
            f"snv_effective_coverage\tuncorrected_mutations\tuncorrected_burden\tuncorrected_burden_upper_ci\tuncorrected_burden_lower_ci\t"
            f"mutations_per_genome\tgenome_length\tcorrected_burden\tcorrected_burden_upper_ci\tcorrected_burden_lower_ci\t"
            f"indel_effective_coverage\tindel_number\tindel_burden\tindel_burden_upper_ci\tindel_burden_lower_ci\n"
        )
        for nn, sample in enumerate(samples):
            sample = sample.strip("/")
            stats_file = f"{sample}/{sample}_stats.txt"
            snv_burden_file = f"{sample}/{sample}_sbs_burden.txt"
            indel_burden_file = f"{sample}/{sample}_indel_burden.txt"

            if not os.path.exists(stats_file):
                raise FileNotFoundError(f"Stats file not found: {stats_file}")
            if not os.path.exists(snv_burden_file):
                raise FileNotFoundError(f"SNV burden file not found: {snv_burden_file}")
            if not os.path.exists(indel_burden_file):
                raise FileNotFoundError(
                    f"Indel burden file not found: {indel_burden_file}"
                )

            with open(stats_file) as stats:
                lines = stats.readlines()
                uniq_reads = int(lines[0].strip("\n").split("\t")[1])
                pf_reads = int(lines[1].strip("\n").split("\t")[1])
                pf_read_family = int(lines[2].strip("\n").split("\t")[1])
                eff_cov = int(lines[3].strip("\n").split("\t")[1])
                dup_rate = float(lines[8].strip("\n").split("\t")[1])
                read_efficiency = float(lines[9].strip("\n").split("\t")[1])
            with open(snv_burden_file) as f:
                lines = f.readlines()
                uncorrected_burden = float(lines[0].strip("\n").split("\t")[1])
                uncorrected_burden_lci = float(lines[1].strip("\n").split("\t")[1])
                uncorrected_burden_uci = float(lines[2].strip("\n").split("\t")[1])
                uncorrected_mutnum = int(lines[3].strip("\n").split("\t")[1])
                corrected_burden = float(lines[4].strip("\n").split("\t")[1])
                corrected_burden_lci = float(lines[5].strip("\n").split("\t")[1])
                corrected_burden_uci = float(lines[6].strip("\n").split("\t")[1])
                corrected_mutnum = float(lines[7].strip("\n").split("\t")[1])
                genome_cov = int(float(lines[8].strip("\n").split("\t")[1]))
            with open(indel_burden_file) as f:
                lines = f.readlines()
                indel_burden = float(lines[0].strip("\n").split("\t")[1])
                indel_burden_lci = float(lines[1].strip("\n").split("\t")[1])
                indel_burden_uci = float(lines[2].strip("\n").split("\t")[1])
                indel_num = int(lines[3].strip("\n").split("\t")[1])
            # indel effective coverage is in stats file line 5
            with open(stats_file) as stats:
                lines = stats.readlines()
                indel_eff_cov = int(lines[5].strip("\n").split("\t")[1])
            output.write(
                f"{sample}\t{pf_reads}\t{uniq_reads}\t{pf_read_family}\t{dup_rate}\t{read_efficiency}\t"
                f"{eff_cov}\t{uncorrected_mutnum}\t{uncorrected_burden}\t{uncorrected_burden_uci}\t{uncorrected_burden_lci}\t"
                f"{corrected_mutnum}\t{genome_cov}\t{corrected_burden}\t{corrected_burden_uci}\t{corrected_burden_lci}\t"
                f"{indel_eff_cov}\t{indel_num}\t{indel_burden}\t{indel_burden_uci}\t{indel_burden_lci}\n"
            )
            if nn == 0:
                sbs96_file = f"{sample}/{sample}_sbs_96_corrected.txt"
                if not os.path.exists(sbs96_file):
                    raise FileNotFoundError(
                        f"SBS96 corrected file not found: {sbs96_file}"
                    )
                sbs96_pd_now = pd.read_csv(sbs96_file, sep="\t", index_col=0)
                sbs96_pd = pd.DataFrame(index=sbs96_pd_now.index)
                sbs96_pd["MutationType"] = sbs96_pd.index
                sbs96_pd[sample] = sbs96_pd_now["number"]
            else:
                sbs96_file = f"{sample}/{sample}_sbs_96_corrected.txt"
                if not os.path.exists(sbs96_file):
                    raise FileNotFoundError(
                        f"SBS96 corrected file not found: {sbs96_file}"
                    )
                sbs96_pd_now = pd.read_csv(sbs96_file, sep="\t", index_col=0)
                sbs96_pd[sample] = sbs96_pd_now["number"]
    sbs96_pd.sort_index(inplace=True)
    sbs96_pd.to_csv(
        args.output.removesuffix(".txt") + "_corrected.SBS96.all", sep="\t", index=False
    )
