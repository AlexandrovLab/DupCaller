import pandas as pd
import numpy as np


def sum_dataframes(dfs):
    """
    Takes multiple pandas DataFrame objects as input and returns the sum of them.

    Args:
    *dfs: Variable length argument containing pandas DataFrame objects

    Returns:
    pandas DataFrame: The sum of input DataFrames
    """
    # Concatenate the DataFrames along the index
    concatenated_df = pd.concat(dfs)

    # Group by index and sum the values
    summed_df = concatenated_df.groupby(level=concatenated_df.index.names).sum()

    return summed_df


def do_aggregate(args):
    dmg_id = list()
    dmg_tn = list()
    amp_id = list()
    amp_tn = list()
    for sample in args.input:
        dmg_id.append(pd.read_csv(sample + "/" + sample + ".dmg.id.txt", sep="\t"))
        dmg_tn.append(pd.read_csv(sample + "/" + sample + ".dmg.tn.txt", sep="\t"))
        amp_id.append(pd.read_csv(sample + "/" + sample + ".amp.id.txt", sep="\t"))
        amp_tn.append(pd.read_csv(sample + "/" + sample + ".amp.tn.txt", sep="\t"))
    # aggregate_dmg_id = sum_dataframes(trinuc_pds)
    # aggregate.to_csv(args.output, sep="\t", index=False)
