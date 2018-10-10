#! /usr/bin/env python

import sys, gzip

sample_info_file = sys.argv[1]
junc_info_file = sys.argv[2]
coverage_file = sys.argv[3]

sample_id2info = {}
with open(sample_info_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        sample_id2info[F[0]] = '\t'.join(F)


junc_id2junc_key1 = {}
junc_id2junc_key2 = {}
with open(junc_info_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        junc_id2junc_key1[F[1]] = F[0] + ';' + F[1] + ';' + F[2]
        junc_id2junc_key2[F[2]] = F[0] + ';' + F[1] + ';' + F[2]


key2count = {}
with gzip.open(coverage_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if not (F[0] in junc_id2junc_key1) and not (F[0] in junc_id2junc_key2): continue

        sample_ids = F[1].split(',')
        counts = F[2].split(',')

        for i in range(len(sample_ids)):
            sample_id = sample_ids[i]
            count = counts[i]
            if sample_id not in sample_id2info: continue

            if F[0] in junc_id2junc_key1:
                key = sample_id2info[sample_id] + '\t' + junc_id2junc_key1[F[0]]
                if key not in key2count: key2count[key] = ["0", "0"]
                key2count[key][0] = count
            elif F[0] in junc_id2junc_key2:
                key = sample_id2info[sample_id] + '\t' + junc_id2junc_key2[F[0]]
                if key not in key2count: key2count[key] = ["0", "0"]
                key2count[key][1] = count
 

print "Sample_ID" + '\t' + "Cancer_Type" + '\t' + "Sample_Name" + '\t' + "Mutation_Info" + '\t' + "Mapped_Read_Count" + '\t' + "Avg_Read_Length" + '\t' + "Splicing_Key" + '\t' + "Read_Count1" + '\t' + "Read_Count2" 
for key in key2count:
    print key + '\t' + key2count[key][0] + '\t' + key2count[key][1]





