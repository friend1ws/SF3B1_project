#! /usr/bin/env python

import sys, pysam

input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]
id_file = sys.argv[4]

id_tb = pysam.TabixFile(id_file)

"""
key2id = {}
with open(id_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        FF = F[3].split('|')
        for i in [-2, 2, -1, 1, 0]:
            for j in [-2, 2, -1, 1, 0]:
                key = F[0] + '\t' + str(int(F[1]) + i) + '\t' + str(int(F[2]) + j)
                if key not in key_list: continue
                key2id[F[0] + '\t' + F[1] + '\t' + F[2]] = FF[0]
"""


def check_id(chr, start, end):

    tabixErrorFlag = 0
    try:
        records = id_tb.fetch(chr, start - 5, end + 5)
    except Exception as inst:
        tabixErrorFlag = 1
    
    id = None
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t') 
            if chr == record[0] and str(start) == record[1] and str(end) == record[2]:
               id = record[3].split('|')[0]

    return id


hout = open(output_file, 'w')
with open(input_file1, 'r') as hin1:
    with open(input_file2, 'r') as hin2:

        while True:
            line1 = hin1.readline()
            line2 = hin2.readline()
            if not line1 or not line2:
                break

            F1 = line1.rstrip('\n').split('\t')
            F2 = line2.rstrip('\n').split('\t')

            FF = F1[3].split(';')

            alt_38_start = str(int(F1[1]) + 1)
            alt_38_end = str(int(F1[2]) + 1)
            ref_38_start = str(int(F2[1]) + 1)
            ref_38_end = str(int(F2[2]) + 1)

            alt_38 = ','.join([F1[0], alt_38_start, alt_38_end]).replace('chr', '')
            ref_38 = ','.join([F2[0], ref_38_start, ref_38_end]).replace('chr', '')
            alt_37 = ','.join([FF[0], FF[1], FF[2]])
            if FF[5] == '+':
                ref_37 = ','.join([FF[0], FF[1], FF[4]])
            elif FF[5] == '-':
                ref_37 = ','.join([FF[0], FF[4], FF[2]])
            
            id1 = check_id(F1[0], int(F1[1]), int(F1[2]))
            id2 = check_id(F2[0], int(F2[1]), int(F2[2]))

            
            if id1 is not None and id2 is not None:
                # print >> hout, F1[3] + '\t' + id1 + '\t' + id2
                print >> hout, '\t'.join([F1[3], id1, id2, alt_38, ref_38, alt_37, ref_37])

# chr10	101600862	101601168	10;103360620;103360926;DPCD;103360959;+

hout.close()


