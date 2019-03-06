#! /usr/bin/env bash

PROJECT=$1
OUTPUT_DIR=$2

RECOUNT_PATH=http://duffel.rail.bio/recount

if [ ! -f ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.g ]
then
    wget -P ${OUTPUT_DIR} ${RECOUNT_PATH}/${PROJECT}/${PROJECT}.junction_coverage.tsv.gz 
fi

if [ ${PROJECT} = "SRP012682" ] || [ ${PROJECT} = "TCGA" ]
then
    echo "OK"
    mv ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz.tmp
    python subscript_recount2_pred/filter_junction_coverage.py \
        ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz.tmp \
        ../output/recount2/TCGA/param_matrix.recount2.zibb.txt \
        ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz
    rm -rf ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz.tmp
fi
    
Rscript subscript_recount2_pred/estimate_zibb.R \
    ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz \
    ${OUTPUT_DIR}/${PROJECT}.score.txt


rm -rf ${OUTPUT_DIR}/${PROJECT}.junction_coverage.tsv.gz

