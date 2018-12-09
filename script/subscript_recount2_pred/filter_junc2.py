#! /usr/bin/env python

import sys, gzip

junc_count_file = sys.argv[1]
phenotype_file = sys.argv[2]
param_file = sys.argv[3]
sample_id_file = sys.argv[4]


junc_id2junc_info_alt = {}
junc_id2junc_info_ref = {}
with open(param_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[0].split(';')
        junc_id2junc_info_alt[FF[6]] = F[0]
        junc_id2junc_info_ref[FF[7]] = F[0]

sample_id2sample_name = {}
with open(sample_id_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        sample_id2sample_name[F[0]] = F[2]

sample2info = {}
with open(phenotype_file, 'r') as hin:

    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for (i, cname) in enumerate(header):
        header2ind[cname] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        run_name = F[header2ind["run"]]
        mapped_read_count = F[header2ind["mapped_read_count"]]
        avg_read_length = F[header2ind["avg_read_length"]]
        characteristics = F[header2ind["characteristics"]]
            
        sample2info[run_name] = mapped_read_count + '\t' + avg_read_length + '\t' + characteristics


key2counts = {}
with gzip.open(junc_count_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] not in junc_id2junc_info_alt and F[0] not in junc_id2junc_info_ref: continue

        sample_ids = F[1].split(',')
        counts = F[2].split(',')

        for i in range(len(sample_ids)):
            run_name = sample_id2sample_name[sample_ids[i]]
            mapped_read_count, avg_read_length, characteristics = sample2info[run_name].split('\t')

            if F[0] in junc_id2junc_info_alt:
                junc_info = junc_id2junc_info_alt[F[0]]
                key = '\t'.join([sample_ids[i], run_name, mapped_read_count, avg_read_length, junc_info])
                if key not in key2counts: key2counts[key] = ['0', '0']
                key2counts[key][0] = counts[i]

            if F[0] in junc_id2junc_info_ref:
                junc_info = junc_id2junc_info_ref[F[0]]
                key = '\t'.join([sample_ids[i], run_name, mapped_read_count, avg_read_length, junc_info])
                if key not in key2counts: key2counts[key] = ['0', '0']
                key2counts[key][1] = counts[i]


print '\t'.join(["Sample_ID", "Sample_Name", "Mapped_Read_Count", "Avg_Read_Length", "Splicing_Key", "Read_Count1", "Read_Count2"])      
for key in key2counts:
    print key + '\t' + '\t'.join(key2counts[key])

