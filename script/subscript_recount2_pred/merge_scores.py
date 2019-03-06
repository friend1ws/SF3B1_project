#! /usr/bin/env python

from __future__ import print_function
import sys, os, glob

score_dir = sys.argv[1]

all_files = glob.glob(score_dir + "/*.score.txt")

print('\t'.join(["Study", "Run", "Score"]))
for sfile in sorted(all_files):
    # project_name = os.path.basename(sfile).replace(".score.txt", "")
    if os.path.basename(sfile).startswith("TCGA"): continue
    with open(sfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "Study": continue
            if F[0] == "NA": continue
            print('\t'.join(F))


