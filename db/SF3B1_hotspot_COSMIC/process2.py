#! /usr/bin/env bash

import sys

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        F[0] = F[0].replace('chr', '')
        F[1] = str(int(F[1]) + 1)
        print '\t'.join(F)

