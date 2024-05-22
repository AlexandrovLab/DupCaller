#!/usr/bin/env python3
#if __name__ == "__main__":
def do_summarize(args):
    samples = args.input
    with open(args.output, "w") as output:
        output.write(
            f"sample\tpass_filter_reads\tunique_reads\tread_families\tduplication_rate\tread_family_efficiency\tsnv_effective_coverage\tuncorrected_burden\tuncorrected_burden_upper_ci\tuncorrected_burden_lower_ci\tcorrected_burden\tcorrected_burden_upper_ci\tcorrected_burden_lower_ci\n"
        )
        for sample in samples:
            stats_file = f"{sample}/{sample}_stats.txt"
            snv_burden_file = f"{sample}/{sample}_burden.txt"
            #indel_burden_file = f"{sample}/{sample}_indel_burden.txt"
            with open(stats_file) as stats:
                lines = stats.readlines()
                uniq_reads = int(lines[0].strip("\n").split("\t")[1])
                pf_reads = int(lines[1].strip("\n").split("\t")[1])
                pf_read_family = int(lines[2].strip("\n").split("\t")[1])
                eff_cov = int(lines[3].strip("\n").split("\t")[1])
                dup_rate = float(lines[5].strip("\n").split("\t")[1])
                read_efficiency = float(lines[6].strip("\n").split("\t")[1])
            with open(snv_burden_file) as stats:
                lines = stats.readlines()
                uncorrected_burden = float(lines[0].strip("\n").split("\t")[1])
                uncorrected_burden_lci = float(lines[1].strip("\n").split("\t")[1])
                uncorrected_burden_uci = float(lines[2].strip("\n").split("\t")[1])
                corrected_burden = float(lines[3].strip("\n").split("\t")[1])
                corrected_burden_lci = float(lines[4].strip("\n").split("\t")[1])
                corrected_burden_uci = float(lines[5].strip("\n").split("\t")[1])
            """
            with open(indel_burden_file) as stats:
                lines = stats.readlines()
                indel_num = int(lines[0].strip("\n").split("\t")[1])
                indel_cov = int(lines[1].strip("\n").split("\t")[1])
                indel_naive_burden = float(lines[2].strip("\n").split("\t")[1])
            """
            output.write(
                f"{sample}\t{pf_reads}\t{uniq_reads}\t{pf_read_family}\t{dup_rate}\t{read_efficiency}\t{eff_cov}\t{uncorrected_burden}\t{uncorrected_burden_uci}\t{uncorrected_burden_lci}\t{corrected_burden}\t{corrected_burden_uci}\t{corrected_burden_lci}\n"
            )
