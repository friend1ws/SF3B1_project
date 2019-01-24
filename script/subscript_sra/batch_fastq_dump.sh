#! /usr/bin/env bash

if [ ! -d ../output/sra ]
then
    mkdir -p ../output/sra/tasks
fi

tail +2 ../output/recount2/prediction/recount2.score.merged.txt | grep -v SRP012682 | awk '$3 > 5.0 {print}' - | sort -k3 -n -r | cut -f 2 > ../output/sra/tasks/run_id_list.txt


echo -e "--env RUN_ID\t--output-recursive OUTPUT_DIR" > ../output/sra/tasks/tasks-fastq-dump.tsv
while read line
do
    echo ${line} >&2
    test1=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key sra-toolkit/${line}/${line}_1.fastq 2> /dev/null`
    cid1=$?
    test2=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key sra-toolkit/${line}/${line}_2.fastq 2> /dev/null`
    cid2=$?

    if [ $cid1 == 255 -a $cid2 = 255 ]
    then
        echo -e "${line}\ts3://friend1ws-virginia-ecsub2/sra-toolkit/${line}" >> ../output/sra/tasks/tasks-fastq-dump.tsv
        # echo -e "${line}\t${cid1}\t${cid2}"
    fi
done < ../output/sra/tasks/run_id_list.txt

<<_
bash subscript_sra/check_buckets.sh > ../output/sra/tasks/run_id_list_check.txt 

awk '$2 == 255 && $3 == 255 {print}' ../output/sra/tasks/run_id_list_check.txt | cut -f 1 > ../output/sra/tasks/run_id_list_nonexist.txt

echo -e "--env RUN_ID\t--output-recursive OUTPUT_DIR" > ../output/sra/tasks/tasks-fastq-dump.tsv
while read line
do
    echo -e "${line}\ts3://friend1ws-virginia-ecsub2/sra-toolkit/${line}" >> ../output/sra/tasks/tasks-fastq-dump.tsv
done < ../output/sra/tasks/run_id_list_nonexist.txt
_

# ecsub submit --script subscript_sra/run-fastq-dump.sh --tasks ../output/sra/tasks/tasks-fastq-dump.tsv --aws-s3-bucket s3://friend1ws-virginia-ecsub2 --wdir /tmp/ecsub --image inutano/sra-toolkit --aws-ec2-instance-type m4.xlarge --disk-size 200 


