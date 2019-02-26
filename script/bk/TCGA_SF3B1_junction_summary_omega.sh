#! /usr/bin/env bash

if [ ! -d ../output/omega/TCGA ]
then
    mkdir -p ../output/omega/TCGA
fi

python subscript_TCGA_SF3B1/gather_sj_count.py ../output/junc_list/sj_list.txt ../data/sample_list_files > ../output/omega/TCGA/TCGA_recount_SF3B1_junction_summary.txt


