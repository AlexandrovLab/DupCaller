import os

import pandas as pd

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
