#! /usr/bin/env python

import sys

input_file = sys.argv[1]

print '\t'.join(["--env RUN_ID", "--output-recursive OUTPUT_DIR"])

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] == "Project_Name": continue
        if float(F[2]) <= 10.0: continue 
        run_id = F[1]
        output_dir = "s3://friend1ws-virginia-ecsub2/sra-toolkit/" + run_id
        print '\t'.join([run_id, output_dir])
        

