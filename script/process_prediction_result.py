#! /usr/bin/env python

import os, glob

all_files = glob.glob("../output/recount2/prediction/*.score.txt")

project_info = {}
with open("../data/recount2/project_info_curated.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        project_info[F[0]] = F[1] + '\t' + F[2]


for sfile in sorted(all_files):
    project_name = os.path.basename(sfile).replace(".score.txt", "")
    with open(sfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[1] in ["score", "NA"]: continue
            if float(F[1]) >= 20.1:
                print project_name + '\t' + F[0] + '\t' + F[1] + '\t' + project_info[project_name]
   
