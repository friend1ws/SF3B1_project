#! /usr/bin/env python

import glob

SF3B1_keys = {}

with open("GRCh37_SF3B1_hotspot_kchiba.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = '\t'.join(F[:8])
        SF3B1_keys[key] = 1


all_mut_files = glob.glob("/home/yshira/mypaper/savnet_paper.bk170605/analysis/TCGA/data/mutation/*/*.filt.txt")

SF3B1_keys = {}

for mfile in sorted(all_mut_files):
    with open(mfile, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[5] != "exonic": continue
            if F[6] != "SF3B1": continue

            match_flag = 0
            for mkey in ["E622", "R625", "N626", "H662", "T663", "K666", "Q699", "K700", "G740", "K741", "G742", "D781", "E902", "Q903"]:
                if mkey in F[9]: match_flag = 1

            if match_flag == 1:
                FF = F[9].split(':')
                key = '\t'.join(F[:5]) + '\t' + "SF3B1" + '\t' + FF[3] + '\t' + FF[4]
                SF3B1_keys[key] = 1

for key in SF3B1_keys:
    print key

