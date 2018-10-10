#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

if [ ! -d ../output/recount2/prediction ]
then
    mkdir ../output/recount2/prediction
fi

cat ../data/recount2/sample_ids.tsv | cut -f 2 | sort -u > ../data/recount2/project_list.txt

while read project
do
    echo "bash subscript_recount2_pred/recount2_pred.sh ${project} ../output/recount2/prediction"
    bash subscript_recount2_pred/recount2_pred.sh ${project} ../output/recount2/prediction
done < ../data/recount2/project_list.txt

