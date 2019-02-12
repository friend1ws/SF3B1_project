#! /usr/bin/env python

import sys

junc_info_file = sys.argv[1]
param_file = sys.argv[2]

key2info1 = {}
with open(junc_info_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        key2info1[F[0]] = '\t'.join(F[1:])


key2info2 = {}
with open(param_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = ';'.join(F[0].split(';')[:-2])
        key2info2[key] = '\t'.join(F[1:])


print '\t'.join(["Junction_ID_Alt", "Junction_ID_Ref", "Junction_Key_Alt_GRCh38", "Junction_Key_Ref_GRCh38", "Junction_Key_Alt_GRCh37", "Junction_Key_Ref_GRCh37", \
                 "Alpha_0", "Beta_0", "Pi_0", "Alpha_1", "Beta_1", "Pi_1"])

for key in sorted(key2info1):
    if key in key2info2:
        print key2info1[key] + '\t' + key2info2[key]


