#! /usr/bin/env python

import sys

hotspot_check_file = sys.argv[1]
hotspot_info_file = sys.argv[2]
run_id_file = sys.argv[3]

mutkey2aa = {}
with open(hotspot_info_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mutkey = '\t'.join(F[:5])
        mutkey2aa[mutkey] = F[7]

run_id2info = {}
with open(hotspot_check_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        mutkey = '\t'.join(F[1:6])
        mut_info = mutkey2aa[mutkey]
        vaf = F[14]
        run_id = F[0]

        if run_id not in run_id2info: run_id2info[run_id] = []
        run_id2info[run_id].append(mut_info + ',' + vaf)

        
with open(run_id_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        info = ';'.join(run_id2info[F[1]]) if F[1] in run_id2info else "---"
        if float(F[2]) <= 5.0: continue
        print '\t'.join([F[0], F[1], F[2], info])

