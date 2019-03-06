#! /usr/bin/env python

from __future__ import print_function

bigwig2info = {}
with open("../data/recount2/TCGA.tsv", 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if len(F) != len(header2ind): continue
        G = {cname: F[i] for (cname, i) in header2ind.items()}
        bigwig2info[G["bigwig_file"].replace(".bw", "")] = \
            G["gdc_cases.project.project_id"] + '\t' + \
            G["gdc_cases.samples.portions.analytes.aliquots.submitter_id"]

with open("../output/recount2/prediction/TCGA.score.txt", 'r') as hin:
    header = hin.readline().rstrip('\n')
    print(header + '\t' + "Cancer_Type" + '\t' + "Barcode")
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[1] not in bigwig2info: continue
        print('\t'.join(F) + '\t' + bigwig2info[F[1]])

        

