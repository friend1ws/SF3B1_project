#! /usr/bin/env bash

while read line
do
    echo ${line} >&2
    test1=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key sra-toolkit/${line}/${line}_1.fastq 2> /dev/null`
    cid1=$?
    test2=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key sra-toolkit/${line}/${line}_2.fastq 2> /dev/null`
    cid2=$?

    echo -e "${line}\t${cid1}\t${cid2}"

done < ../output/sra/tasks/run_id_list.txt

