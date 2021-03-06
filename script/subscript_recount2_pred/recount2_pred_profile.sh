#! /usr/bin/env bash

PROJECT=$1
OUTPUT_DIR=$2

RECOUNT_PATH=http://duffel.rail.bio/recount

wget -P ${OUTPUT_DIR} ${RECOUNT_PATH}/${PROJECT}/${PROJECT}.junction_coverage.tsv.gz 

wget -P ${OUTPUT_DIR} ${RECOUNT_PATH}/${PROJECT}/${PROJECT}.tsv

python subscript_recount2_pred/filter_junc2.py \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz \
    ${OUTPUT_DIR}/${PROJECT}.tsv \
    ../output/recount2/TCGA/param_matrix.recount2.zibb.txt \
    ../data/recount2/sample_ids.tsv > \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.filt.txt

if [ ! -d ${OUTPUT_DIR}/figure/${PROJECT} ]
then
    mkdir -p ${OUTPUT_DIR}/figure/${PROJECT}
fi

Rscript subscript_recount2_pred/junc_profile.R \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.filt.txt \
    ../output/recount2/TCGA/param_matrix.recount2.zibb.txt \
    ${OUTPUT_DIR}/figure/${PROJECT}


rm -rf ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz
rm -rf ${OUTPUT_DIR}/${PROJECT}.junction_coverage.filt.txt
rm -rf ${OUTPUT_DIR}/${PROJECT}.tsv


