#!/usr/bin/env python3
import sys
import glob
import pathlib
import pandas as pd

inputfolder = sys.argv[1]
outexcelfile = sys.argv[2]
files = glob.glob("**/ame_results.txt", recursive=True)
files.sort()
repname2file = dict()
with pd.ExcelWriter(outexcelfile) as writer:
    for f in files:
        x = pathlib.Path(f)
        y = list(
            filter(
                lambda x: "narrowPeak_motif_enrichment" in x,
                map(lambda x: str(x), x.parents),
            )
        )
        if len(y) != 1:
            continue
        repname = y[0].split(".")[0]
        repname2file[repname] = f
        df = pd.read_csv(f, sep="\t", header=0, skip_blank_lines=True)
        df.to_excel(writer, sheet_name=repname, index=False)
