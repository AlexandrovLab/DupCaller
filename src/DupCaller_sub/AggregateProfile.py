import os

import numpy as np
import pandas as pd

from .funcs.learn import estimate_sbs_srd_rates
from .funcs.misc import build_trinuc64_order

# Count tables written by call/learn as {folder}/ERROR/{basename(folder)}<suffix>.
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

    # -i gives call output folders; -f gives error-file prefixes
    # ({folder}/ERROR/{sample}) directly. Normalize both to prefixes.
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

    # The SRD matrix is an EM fit, so it can't be summed directly; instead
    # sum the per-sample BQ histograms (additive) and refit once. The
    # histograms live in {folder}/tmp/, a sibling of ERROR/.
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
