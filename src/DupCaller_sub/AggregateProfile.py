import os

import numpy as np
import pandas as pd

from .funcs.learn import estimate_sbs_srd_rates
from .funcs.misc import build_trinuc64_order

# File suffixes written by DupCaller.py call/learn under
# {folder}/ERROR/{basename(folder)}<suffix> -- see Caller.py's
# amp_tn_pd/dmg_tn_pd/amperri_hp_file etc. write-out for the source of
# truth on this naming.
_PROFILE_SUFFIXES = [
    ".amp.tn.txt",
    ".dmg.tn.txt",
    ".amp.hp.txt",
    ".dmg.hp.txt",
    ".amp.str.txt",
    ".dmg.str.txt",
]


def do_aggregate(args):
    if not args.input and not args.input_file:
        raise ValueError("Either -i/--input or -f/--input-file must be provided")
    if not args.output:
        raise ValueError("-o/--output must be provided")

    # -i/--input takes DupCallerCall *output folders*, whose error files
    # live at {folder}/ERROR/{basename(folder)}{suffix} (see Caller.py/
    # Learn.py's error_prefix -- named after the folder's own basename,
    # not its full path). -f/--input-file instead takes lines that are
    # already the full error-file *prefix* (folder + "/ERROR/" +
    # basename, ready to have suffix appended directly), per its own
    # help text -- so the two inputs are normalized into one list of
    # prefixes here rather than reusing a single (folder-shaped) path
    # formula for both.
    prefixes = []
    if args.input:
        for folder in args.input:
            folder = folder.rstrip("/")
            sample_name = os.path.basename(folder)
            prefixes.append(os.path.join(folder, "ERROR", sample_name))
    if args.input_file:
        with open(args.input_file) as fh:
            for line in fh:
                prefix = line.strip()
                if prefix:
                    prefixes.append(prefix)

    for suffix in _PROFILE_SUFFIXES:
        profiles = [
            pd.read_csv(f"{prefix}{suffix}", sep="\t", index_col=0)
            for prefix in prefixes
        ]
        aggregated = sum(profiles)
        aggregated.to_csv(args.output + suffix, sep="\t")

    # SBS SRD (single-read-damage) rate matrix: not a simple sum-of-tables
    # like the legacy profiles above, since it's an EM fit rather than a
    # raw count. Its sufficient statistic -- the per-trinuc-context,
    # per-base-quality BQ histogram of raw amp-error observations
    # (sbs_alt_bq_hist, see funcs/learn.py's profileTriNucMismatches/
    # estimate_sbs_srd_rates) -- IS additive across samples, so the
    # correct way to combine multiple samples' SRD estimates is to sum
    # their raw histograms first and re-run the EM fit once on the sum,
    # not to average/combine the samples' already-fit (64,4) rate
    # matrices directly.
    #
    # The histogram is saved by Learn.py/Caller.py as a QC diagnostic
    # under {folder}/tmp/{sample_name}.amp.tn.bqhist.npz (tmp/, a sibling
    # of ERROR/ -- not under the ERROR/-based ${prefix} used above), so
    # its path has to be re-derived from each ${prefix} rather than just
    # appending a new suffix.
    hist_sum = None
    for prefix in prefixes:
        error_dir, sample_name = os.path.split(prefix)
        folder = os.path.dirname(error_dir)
        bqhist_path = os.path.join(folder, "tmp", f"{sample_name}.amp.tn.bqhist.npz")
        if not os.path.exists(bqhist_path):
            raise FileNotFoundError(
                f"BQ histogram not found: {bqhist_path} (expected next to {prefix})"
            )
        with np.load(bqhist_path) as npz:
            hist = npz["hist"]
        hist_sum = hist if hist_sum is None else hist_sum + hist

    _, num2trinuc = build_trinuc64_order()
    srd_mat = estimate_sbs_srd_rates(hist_sum, args.pseudocount)
    srd_pd = pd.DataFrame(srd_mat, columns=["A", "T", "C", "G"], index=num2trinuc)
    srd_pd.to_csv(args.output + ".amp.tn.srd.txt", sep="\t")
