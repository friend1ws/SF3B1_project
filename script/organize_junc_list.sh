#! /bin/bash

if [ ! -d ../output/junc_list ]
then
    mkdir -p ../output/junc_list 
fi

if [ ! -d ../data/recount2 ]
then
    mkdir -p ../data/recount2
fi

if [ ! -f ../db/junction.bed.gz ]
then
    annot_utils junction --gene_model gencode --grc ../db/junction.bed.gz
fi

if [ ! -f ../db/hg19ToHg38.over.chain ]
then
    wget -P ../db http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    gunzip ../db/hg19ToHg38.over.chain.gz
fi

if [ ! -f ../data/recount2/TCGA.junction_id_with_transcripts.bed.gz ]
then
    wget -P ../data/recount2 http://duffel.rail.bio/recount/TCGA/TCGA.junction_id_with_transcripts.bed.gz
    zcat ../data/recount2/TCGA.junction_id_with_transcripts.bed.gz > ../data/recount2/TCGA.junction_id_with_transcripts.bed
    bgzip -f ../data/recount2/TCGA.junction_id_with_transcripts.bed
    tabix -p bed ../data/recount2/TCGA.junction_id_with_transcripts.bed.gz
fi


python subscript_organize_junc_list/read_sj_list.py > ../output/junc_list/sj_list.txt

python subscript_organize_junc_list/proc.py ../output/junc_list/sj_list.txt ../output/junc_list/sj_list_alt.bed.tmp ../output/junc_list/sj_list_ref.bed.tmp

liftOver ../output/junc_list/sj_list_ref.bed.tmp ../db/hg19ToHg38.over.chain ../output/junc_list/sj_list_ref.bed ../output/junc_list/unmap_ref.bed
liftOver ../output/junc_list/sj_list_alt.bed.tmp ../db/hg19ToHg38.over.chain ../output/junc_list/sj_list_alt.bed ../output/junc_list/unmap_ref.bed

python subscript_organize_junc_list/add_id.py ../output/junc_list/sj_list_alt.bed ../output/junc_list/sj_list_ref.bed ../output/junc_list/sj_list_id.txt ../data/recount2/TCGA.junction_id_with_transcripts.bed.gz 

rm -rf ../output/junc_list/unmap_ref.bed
rm -rf ../output/junc_list/sj_list_ref.bed.tmp
rm -rf ../output/junc_list/sj_list_ref.bed
rm -rf ../output/junc_list/sj_list_alt.bed.tmp
rm -rf ../output/junc_list/sj_list_alt.bed


