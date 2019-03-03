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
    if [ ! -f ../output/recount2/prediction/${project}.score.txt ]
    then
        echo "bash subscript_recount2_pred/recount2_pred_zibb.sh ${project} ../output/recount2/prediction"
        bash subscript_recount2_pred/recount2_pred_zibb.sh ${project} ../output/recount2/prediction
    fi
done < ../data/recount2/project_list.txt

python subscript_recount2_pred/merge_scores.py ../output/recount2/prediction > ../output/recount2/prediction/recount2.score.merged.txt

