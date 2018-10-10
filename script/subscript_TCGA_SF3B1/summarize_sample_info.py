#! /usr/bin/env python

import sys, os, glob

savnet_input_file_dir = sys.argv[1]
recount_TCGA_file = sys.argv[2]
sample_id_file = sys.argv[3]

sample_name2info = {}
all_savnet_input_file = glob.glob(savnet_input_file_dir + "/*.input.txt")
for savnet_input_file in sorted(all_savnet_input_file):

    cancer_type = os.path.basename(savnet_input_file).replace(".input.txt", '')
    with open(savnet_input_file) as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            sample_name = F[header2ind["Sample_Name"]]
            sj_file = F[header2ind["SJ_File"]]
            mutation_info = F[header2ind["Mutation_Info"]]
            weight = F[header2ind["Weight"]]

            sample_name2 = os.path.basename(sj_file).replace(".SJ.out.tab", '')
            sample_name2info[sample_name2] = cancer_type + '\t' + sample_name + '\t' + mutation_info


bw2info = {}
with open(recount_TCGA_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if len(F) < header2ind["bigwig_file"]: continue
        bw = F[header2ind["bigwig_file"]].replace(".bw", "")
        sample_name2 = F[header2ind["gdc_cases.samples.portions.analytes.aliquots.submitter_id"]]
        mapped_read_count = F[header2ind["mapped_read_count"]]
        avg_read_length = F[header2ind["avg_read_length"]]
 
        if sample_name2 not in sample_name2info: continue

        bw2info[bw] = sample_name2info[sample_name2] + '\t' + mapped_read_count + '\t' + avg_read_length


with open(sample_id_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[2] not in bw2info: continue
        print F[0] + '\t' + bw2info[F[2]]


