#! /usr/bin/env python

import sys

input_file = sys.argv[1]

print '\t'.join(["--env SAMPLE", "--input FASTQ1", "--input FASTQ2", "--output-recursive OUTPUT_DIR", "--input-recursive REFERENCE"])

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        sample = F[0]
        fastq1 = "s3://friend1ws-virginia-ecsub2/sra-toolkit/" + sample + "/" + sample + "_1.fastq"
        fastq2 = "s3://friend1ws-virginia-ecsub2/sra-toolkit/" + sample + "/" + sample + "_2.fastq"
        output_dir = "s3://friend1ws-virginia-ecsub2/star/" + sample
        reference = "s3://genomon-bucket/_GRCh37/reference/GRCh37.STAR-2.5.2a"
        
        print '\t'.join([sample, fastq1, fastq2, output_dir, reference])

 

