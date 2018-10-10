#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file1 = sys.argv[2]
output_file2 = sys.argv[3]

hout1 = open(output_file1, 'w')
hout2 = open(output_file2, 'w')
with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')

        print >> hout1, 'chr' + F[0] + '\t' + str(int(F[1]) - 1) + '\t' + str(int(F[2]) - 1) + '\t' + ';'.join(F)
        if F[5] == '+':
            print >> hout2, 'chr' + F[0] + '\t' + str(int(F[1]) - 1) + '\t' + str(int(F[4]) - 1) + '\t' + ';'.join(F)
        else:
            print >> hout2, 'chr' + F[0] + '\t' + str(int(F[4]) - 1) + '\t' + str(int(F[2]) - 1) + '\t' + ';'.join(F)

hout1.close()
hout2.close()

