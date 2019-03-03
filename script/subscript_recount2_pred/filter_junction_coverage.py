#! /usr/bin/env python

import sys, gzip

junc_count_file = sys.argv[1]
param_file = sys.argv[2]
out_file = sys.argv[3]

junc_id2junc_info_alt = {}
junc_id2junc_info_ref = {}
with open(param_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[0].split(';')
        junc_id2junc_info_alt[FF[6]] = F[0]
        junc_id2junc_info_ref[FF[7]] = F[0]


hout = gzip.open(out_file, 'w')

with gzip.open(junc_count_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] not in junc_id2junc_info_alt and F[0] not in junc_id2junc_info_ref: continue
        print >> hout, '\t'.join(F)

