#! /usr/bin/env python

import sys, os, glob

score_dir = sys.argv[1]

all_files = glob.glob(score_dir + "/*.score.txt")

print '\t'.join(["Project_Name", "Run_Name", "Score"])
for sfile in sorted(all_files):
    project_name = os.path.basename(sfile).replace(".score.txt", "")
    with open(sfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "run_name": continue
            print project_name + '\t' + F[0] + '\t' + F[1]


