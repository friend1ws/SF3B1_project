#! /usr/bin/env bash

PROJECT=$1
OUTPUT_DIR=$2

RECOUNT_PATH=http://duffel.rail.bio/recount

wget -P ${OUTPUT_DIR} ${RECOUNT_PATH}/${PROJECT}/${PROJECT}.junction_coverage.tsv.gz 

wget -P ${OUTPUT_DIR} ${RECOUNT_PATH}/${PROJECT}/${PROJECT}.tsv

python filter_junc.py \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz \
    ${OUTPUT_DIR}/${PROJECT}.tsv param_matrix.recount2.ninb.txt \
    sample_ids.tsv > \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.filt.txt


Rscript estimate.R \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.filt.txt \
    ${OUTPUT_DIR}/${PROJECT}.tsv \
    ${OUTPUT_DIR}/${PROJECT}.score.txt


rm -rf ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz
rm -rf ${OUTPUT_DIR}/${PROJECT}.junction_coverage.filt.txt
rm -rf ${OUTPUT_DIR}/${PROJECT}.tsv

