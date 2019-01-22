#!/bin/bash

set -o errexit
set -o nounset

fastq-dump -v --split-files ${RUN_ID} -O ${OUTPUT_DIR}
# SRR628583 -O /data/SRR628583 

# python ${SCRIPT}/wordcount.py ${INPUT_FILE} ${OUTPUT_FILE}
