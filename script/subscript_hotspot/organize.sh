#! /usr/bin/env bash

if [ ! -d ../output/sra/hotspot ]
then
    mkdir -p ../output/sra/hotspot
fi

tail +2 ../output/recount2/prediction/recount2.score.merged2.txt | awk '$3 > 0.0 {print}' - | sort -k3 -n -r > ../output/sra/hotspot/run_id_list2.txt

python subscript_hotspot/organize_result.py ../output/sra/hotspot/SF3B1.hotspot.result.txt ../db/SF3B1_hotspot_COSMIC/COSMIC_SF3B1_hotspot.txt ../output/sra/hotspot/run_id_list2.txt > ../output/sra/hotspot/SF3B1.hotspot.result.info.txt

Rscript subscript_hotspot/add_sra_meta.R

python subscript_hotspot/summarize_study.py ../output/sra/hotspot/SF3B1.hotspot.result.info2.txt ../data/recount2/project_info_curated.txt > ../output/sra/hotspot/study_summary.txt 
