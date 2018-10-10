#! /bin/bash

if [ ! -d ../output/recount2/TCGA ]
then
    mkdir -p ../output/recount2/TCGA
fi

if [ ! -f ../data/recount2/TCGA.tsv ]
then
    wget -P ../data/recount2/ http://duffel.rail.bio/recount/TCGA/TCGA.tsv
fi

if [ ! -f ../data/recount2/sample_ids.tsv ]
then
    wget -P ../data/recount2 https://jhubiostatistics.shinyapps.io/recount/sample_ids.tsv
fi

if [ ! -f ../data/recount2/TCGA.junction_coverage.tsv.gz ] 
then
    wget -P ../data/recount2/ http://duffel.rail.bio/recount/TCGA/TCGA.junction_coverage.tsv.gz
fi


python subscript_TCGA_SF3B1/summarize_sample_info.py ../data/sample_list_files/ ../data/recount2/TCGA.tsv ../data/recount2/sample_ids.tsv | sort -k2,2 -k3,3 > ../output/recount2/TCGA/my_sample_info.txt

python subscript_TCGA_SF3B1/summarize_count.py ../output/recount2/TCGA/my_sample_info.txt ../output/junc_list/sj_list_id.txt ../data/recount2/TCGA.junction_coverage.tsv.gz > ../output/recount2/TCGA/TCGA_recount_SF3B1_junction_summary.txt


