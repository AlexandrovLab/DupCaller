#!/usr/bin/env python3
# if __name__ == "__main__":
import pandas as pd
import os

from .funcs.misc import parse_stats_file
from .Estimate import INDEL83_TO_SIGPROFILER_LABELS


def do_summarize(args):
    samples = args.input
    with open(args.output, "w") as output:
        output.write(
            f"sample\tpass_filter_reads\tunique_reads\tread_families\tduplication_rate\tread_family_efficiency\t"
            f"sbs_base_coverage\tsbs_uncorrected_mutations\tsbs_uncorrected_burden\tsbs_uncorrected_burden_upper_ci\tsbs_uncorrected_burden_lower_ci\t"
            f"sbs_corrected_mutations\tsbs_mutations_per_genome\tsbs_mutations_per_genome_upper_ci\tsbs_mutations_per_genome_lower_ci\t"
            f"genome_length\tsbs_corrected_burden\tsbs_corrected_burden_upper_ci\tsbs_corrected_burden_lower_ci\t"
            f"indel_base_coverage\tindel_uncorrected_mutations\tindel_uncorrected_burden\tindel_uncorrected_burden_upper_ci\tindel_uncorrected_burden_lower_ci\t"
            f"indel_corrected_mutations\tindel_corrected_burden\tindel_corrected_burden_upper_ci\tindel_corrected_burden_lower_ci\t"
            f"indel_mutations_per_genome\tindel_mutations_per_genome_upper_ci\tindel_mutations_per_genome_lower_ci\t"
            f"dbs_base_coverage\tdbs_mutations_per_genome\tdbs_mutations_per_genome_upper_ci\tdbs_mutations_per_genome_lower_ci\n"
        )
        for nn, sample in enumerate(samples):
            sample = sample.strip("/")
            sample_name = os.path.basename(sample)
            stats_file = f"{sample}/{sample_name}_stats.txt"
            snv_burden_file = f"{sample}/SBS/{sample_name}_sbs_burden.txt"
            indel_burden_file = f"{sample}/INDEL/{sample_name}_indel_burden.txt"
            dbs_burden_file = f"{sample}/DBS/{sample_name}_dbs_burden.txt"
            sbs96_file = f"{sample}/SBS/{sample_name}_sbs_96_corrected.txt"
            indel83_file = f"{sample}/INDEL/{sample_name}_indel_83_corrected.txt"

            if not os.path.exists(stats_file):
                raise FileNotFoundError(f"Stats file not found: {stats_file}")
            if not os.path.exists(snv_burden_file):
                raise FileNotFoundError(f"SNV burden file not found: {snv_burden_file}")
            if not os.path.exists(indel_burden_file):
                raise FileNotFoundError(
                    f"Indel burden file not found: {indel_burden_file}"
                )
            if not os.path.exists(dbs_burden_file):
                raise FileNotFoundError(f"DBS burden file not found: {dbs_burden_file}")
            if not os.path.exists(sbs96_file):
                raise FileNotFoundError(f"SBS96 corrected file not found: {sbs96_file}")
            if not os.path.exists(indel83_file):
                raise FileNotFoundError(
                    f"ID83 corrected file not found: {indel83_file}"
                )

            stats = parse_stats_file(stats_file)
            uniq_reads = int(stats["Number of Read Families"])
            pf_reads = int(stats["Number of Pass-filter Reads"])
            pf_read_family = int(stats["Number of Effective Read Families"])
            sbs_base_cov = float(stats["SBS Base Coverage"])
            indel_base_cov = float(stats["Indel Base Coverage"])
            dbs_base_cov = float(stats["DBS Base Coverage"])
            dup_rate = float(stats["Pass-filter Duplication Rate"])
            read_efficiency = float(stats["Efficiency"])

            # parse_stats_file works on any tab-separated key\tvalue file,
            # not just _stats.txt -- used here to key-index into
            # _sbs_burden.txt/_indel_burden.txt/_dbs_burden.txt instead of
            # reading by fixed line number, which breaks silently if a line
            # is ever added/removed/reordered upstream in Estimate.py.
            snv_burden = parse_stats_file(snv_burden_file)
            uncorrected_burden = float(snv_burden["Uncorrected burden"])
            uncorrected_burden_lci = float(snv_burden["Uncorrected burden 95% lower"])
            uncorrected_burden_uci = float(snv_burden["Uncorrected burden 95% upper"])
            uncorrected_mutnum = int(snv_burden["Uncorrected mutation number"])
            corrected_burden = float(snv_burden["Corrected burden"])
            corrected_burden_lci = float(snv_burden["Corrected burden 95% lower"])
            corrected_burden_uci = float(snv_burden["Corrected burden 95% upper"])
            corrected_mutnum = float(snv_burden["Corrected mutation number"])
            mutations_per_genome = float(snv_burden["Mutation number per genome"])
            mutations_per_genome_lci = float(
                snv_burden["Mutation number per genome 95% lower"]
            )
            mutations_per_genome_uci = float(
                snv_burden["Mutation number per genome 95% upper"]
            )
            genome_length = int(float(snv_burden["Reference base number"]))

            # _indel_burden.txt mirrors _sbs_burden.txt's field names exactly
            # (see Estimate.py); the uncorrected block, corrected block, and
            # the per-genome figure are all consumed here.
            indel_burden = parse_stats_file(indel_burden_file)
            indel_uncorrected_burden = float(indel_burden["Uncorrected burden"])
            indel_uncorrected_burden_lci = float(
                indel_burden["Uncorrected burden 95% lower"]
            )
            indel_uncorrected_burden_uci = float(
                indel_burden["Uncorrected burden 95% upper"]
            )
            indel_uncorrected_mutnum = int(indel_burden["Uncorrected mutation number"])
            indel_corrected_burden = float(indel_burden["Corrected burden"])
            indel_corrected_burden_lci = float(
                indel_burden["Corrected burden 95% lower"]
            )
            indel_corrected_burden_uci = float(
                indel_burden["Corrected burden 95% upper"]
            )
            indel_corrected_mutnum = float(indel_burden["Corrected mutation number"])
            indel_mutations_per_genome = float(
                indel_burden["Mutation number per genome"]
            )
            indel_mutations_per_genome_lci = float(
                indel_burden["Mutation number per genome 95% lower"]
            )
            indel_mutations_per_genome_uci = float(
                indel_burden["Mutation number per genome 95% upper"]
            )

            # _dbs_burden.txt mirrors _sbs_burden.txt (no "Unmasked burden"
            # block -- see Estimate.py); only the per-genome figure is
            # consumed here.
            dbs_burden = parse_stats_file(dbs_burden_file)
            dbs_mutations_per_genome = float(dbs_burden["Mutation number per genome"])
            dbs_mutations_per_genome_lci = float(
                dbs_burden["Mutation number per genome 95% lower"]
            )
            dbs_mutations_per_genome_uci = float(
                dbs_burden["Mutation number per genome 95% upper"]
            )

            output.write(
                f"{sample_name}\t{pf_reads}\t{uniq_reads}\t{pf_read_family}\t{dup_rate}\t{read_efficiency}\t"
                f"{sbs_base_cov}\t{uncorrected_mutnum}\t{uncorrected_burden}\t{uncorrected_burden_uci}\t{uncorrected_burden_lci}\t"
                f"{corrected_mutnum}\t{mutations_per_genome}\t{mutations_per_genome_uci}\t{mutations_per_genome_lci}\t"
                f"{genome_length}\t{corrected_burden}\t{corrected_burden_uci}\t{corrected_burden_lci}\t"
                f"{indel_base_cov}\t{indel_uncorrected_mutnum}\t{indel_uncorrected_burden}\t"
                f"{indel_uncorrected_burden_uci}\t{indel_uncorrected_burden_lci}\t"
                f"{indel_corrected_mutnum}\t{indel_corrected_burden}\t"
                f"{indel_corrected_burden_uci}\t{indel_corrected_burden_lci}\t"
                f"{indel_mutations_per_genome}\t{indel_mutations_per_genome_uci}\t{indel_mutations_per_genome_lci}\t"
                f"{dbs_base_cov}\t{dbs_mutations_per_genome}\t{dbs_mutations_per_genome_uci}\t{dbs_mutations_per_genome_lci}\n"
            )

            sbs96_pd_now = pd.read_csv(sbs96_file, sep="\t", index_col=0)
            # _indel_83_corrected.txt is indexed by our internal ID83 syntax
            # (e.g. "Cdelhp1", "2delstr3") -- relabel to SigProfilerPlotting's
            # COSMIC-standard ID83 strings (e.g. "1:Del:C:0") via
            # INDEL83_TO_SIGPROFILER_LABELS, the same mapping Estimate.py's
            # own per-sample sigPlt.plotID call uses, so the aggregated
            # matrices below are directly consumable by SigProfilerPlotting/
            # SigProfilerExtractor without any further relabeling.
            indel83_pd_now = pd.read_csv(indel83_file, sep="\t", index_col=0).rename(
                index=INDEL83_TO_SIGPROFILER_LABELS
            )
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
                id83_uncorrected = pd.DataFrame(
                    {"MutationType": indel83_pd_now.index}, index=indel83_pd_now.index
                )
                id83_corrected = pd.DataFrame(
                    {"MutationType": indel83_pd_now.index}, index=indel83_pd_now.index
                )
                id83_genome = pd.DataFrame(
                    {"MutationType": indel83_pd_now.index}, index=indel83_pd_now.index
                )
            sbs96_uncorrected[sample_name] = sbs96_pd_now["mutation_number_uncorrected"]
            sbs96_corrected[sample_name] = sbs96_pd_now["mutation_number_corrected"]
            sbs96_genome[sample_name] = sbs96_pd_now["mutation_number_genome"]
            id83_uncorrected[sample_name] = indel83_pd_now[
                "mutation_number_uncorrected"
            ]
            id83_corrected[sample_name] = indel83_pd_now["mutation_number_corrected"]
            id83_genome[sample_name] = indel83_pd_now["mutation_number_genome"]

    prefix = args.output.removesuffix(".txt")
    for df, suffix in [
        (sbs96_uncorrected, "_SBS96_uncorrected.txt"),
        (sbs96_corrected, "_SBS96_corrected.txt"),
        (sbs96_genome, "_SBS96_genome.txt"),
        (id83_uncorrected, "_ID83_uncorrected.txt"),
        (id83_corrected, "_ID83_corrected.txt"),
        (id83_genome, "_ID83_genome.txt"),
    ]:
        df.sort_index(inplace=True)
        df.to_csv(prefix + suffix, sep="\t", index=False)
