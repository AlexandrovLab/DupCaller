#!/usr/bin/env python3
# if __name__ == "__main__":
import pandas as pd
import os

from .funcs.misc import parse_stats_file


def do_summarize(args):
    samples = args.input
    with open(args.output, "w") as output:
        output.write(
            f"sample\tpass_filter_reads\tunique_reads\tread_families\tduplication_rate\tread_family_efficiency\t"
            f"sbs_base_coverage\tuncorrected_mutations\tuncorrected_burden\tuncorrected_burden_upper_ci\tuncorrected_burden_lower_ci\t"
            f"corrected_mutations\tmutations_per_genome\tgenome_length\tcorrected_burden\tcorrected_burden_upper_ci\tcorrected_burden_lower_ci\t"
            f"indel_base_coverage\tindel_number\tindel_burden\tindel_burden_upper_ci\tindel_burden_lower_ci\tdbs_base_coverage\n"
        )
        for nn, sample in enumerate(samples):
            sample = sample.strip("/")
            sample_name = os.path.basename(sample)
            stats_file = f"{sample}/{sample_name}_stats.txt"
            snv_burden_file = f"{sample}/SBS/{sample_name}_sbs_burden.txt"
            indel_burden_file = f"{sample}/INDEL/{sample_name}_indel_burden.txt"
            sbs96_file = f"{sample}/SBS/{sample_name}_sbs_96_corrected.txt"

            if not os.path.exists(stats_file):
                raise FileNotFoundError(f"Stats file not found: {stats_file}")
            if not os.path.exists(snv_burden_file):
                raise FileNotFoundError(f"SNV burden file not found: {snv_burden_file}")
            if not os.path.exists(indel_burden_file):
                raise FileNotFoundError(
                    f"Indel burden file not found: {indel_burden_file}"
                )
            if not os.path.exists(sbs96_file):
                raise FileNotFoundError(f"SBS96 corrected file not found: {sbs96_file}")

            stats = parse_stats_file(stats_file)
            uniq_reads = int(stats["Number of Read Families"])
            pf_reads = int(stats["Number of Pass-filter Reads"])
            pf_read_family = int(stats["Number of Effective Read Families"])
            sbs_base_cov = float(stats["SBS Base Coverage"])
            indel_base_cov = float(stats["Indel Base Coverage"])
            dbs_base_cov = float(stats["DBS Base Coverage"])
            dup_rate = float(stats["Pass-filter Duplication Rate"])
            read_efficiency = float(stats["Efficiency"])

            # _sbs_burden.txt line layout:
            # 0: Uncorrected burden
            # 1: Uncorrected burden 95% lower
            # 2: Uncorrected burden 95% upper
            # 3: Uncorrected mutation number
            # 4: Corrected burden
            # 5: Corrected burden 95% lower
            # 6: Corrected burden 95% upper
            # 7: Corrected mutation number
            # 8: Corrected mutation number (duplicate; no longer a distinct
            #    genome-extrapolated figure -- mutations_per_genome below is
            #    always numerically identical to corrected_mutations)
            # 9: Duplex coverage (effective coverage at min_group_size=1)
            # 10: Unmasked burden
            # 11: Unmasked burden 95% lower
            # 12: Unmasked burden 95% upper
            # 13: Reference base number
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
                mutations_per_genome = float(lines[8].strip("\n").split("\t")[1])
                genome_length = int(float(lines[13].strip("\n").split("\t")[1]))

            # _indel_burden.txt line layout:
            # 0: Uncorrected burden
            # 1: Uncorrected burden 95% lower
            # 2: Uncorrected burden 95% upper
            # 3: Uncorrected mutation number
            # 9: Duplex coverage (effective coverage at min_group_size=1)
            with open(indel_burden_file) as f:
                lines = f.readlines()
                indel_burden = float(lines[0].strip("\n").split("\t")[1])
                indel_burden_lci = float(lines[1].strip("\n").split("\t")[1])
                indel_burden_uci = float(lines[2].strip("\n").split("\t")[1])
                indel_num = int(lines[3].strip("\n").split("\t")[1])

            output.write(
                f"{sample_name}\t{pf_reads}\t{uniq_reads}\t{pf_read_family}\t{dup_rate}\t{read_efficiency}\t"
                f"{sbs_base_cov}\t{uncorrected_mutnum}\t{uncorrected_burden}\t{uncorrected_burden_uci}\t{uncorrected_burden_lci}\t"
                f"{corrected_mutnum}\t{mutations_per_genome}\t{genome_length}\t{corrected_burden}\t{corrected_burden_uci}\t{corrected_burden_lci}\t"
                f"{indel_base_cov}\t{indel_num}\t{indel_burden}\t{indel_burden_uci}\t{indel_burden_lci}\t{dbs_base_cov}\n"
            )

            sbs96_pd_now = pd.read_csv(sbs96_file, sep="\t", index_col=0)
            if nn == 0:
                sbs96_uncorrected = pd.DataFrame(
                    {"MutationType": sbs96_pd_now.index}, index=sbs96_pd_now.index
                )
                sbs96_corrected = pd.DataFrame(
                    {"MutationType": sbs96_pd_now.index}, index=sbs96_pd_now.index
                )
                sbs96_genome = pd.DataFrame(
                    {"MutationType": sbs96_pd_now.index}, index=sbs96_pd_now.index
                )
            sbs96_uncorrected[sample_name] = sbs96_pd_now["mutation_number_uncorrected"]
            sbs96_corrected[sample_name] = sbs96_pd_now["mutation_number_corrected"]
            sbs96_genome[sample_name] = sbs96_pd_now["mutation_number_genome"]

    prefix = args.output.removesuffix(".txt")
    for df, suffix in [
        (sbs96_uncorrected, "_SBS96_uncorrected.txt"),
        (sbs96_corrected, "_SBS96_corrected.txt"),
        (sbs96_genome, "_SBS96_genome.txt"),
    ]:
        df.sort_index(inplace=True)
        df.to_csv(prefix + suffix, sep="\t", index=False)
