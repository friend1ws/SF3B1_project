#! /usr/bin/env bash

bash subscript_sra/check_buckets.sh > ../output/sra/tasks/run_id_list_check.txt

awk '$2 == 0 && $3 == 0 {print}' ../output/sra/tasks/run_id_list_check.txt | cut -f 1 > ../output/sra/tasks/run_id_list_paired.txt

awk '$2 == 0 && $3 != 255 {print}' ../output/sra/tasks/run_id_list_check.txt | cut -f 1 > ../output/sra/tasks/run_id_list_single.txt
_

echo -e "--env SAMPLE\t--input FASTQ1\t--input FASTQ2\t--output-recursive OUTPUT_DIR\t--input-recursive REFERENCE" > ../output/sra/tasks/tasks-star-alignment-paired.tsv 
while read line
do
    echo ${line} >&2
    test1=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key star/${line}/${line}.Aligned.sortedByCoord.out.bam 2> /dev/null`
    cid1=$?
    test2=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key star/${line}/${line}.Aligned.sortedByCoord.out.bam.bai 2> /dev/null`
    cid2=$?

    if [ $cid1 != 0 -o $cid2 != 0 ]
    then
        fastq1="s3://friend1ws-virginia-ecsub2/sra-toolkit/${line}/${line}_1.fastq"
        fastq2="s3://friend1ws-virginia-ecsub2/sra-toolkit/${line}/${line}_2.fastq"
        output_dir="s3://friend1ws-virginia-ecsub2/star/${line}"
        reference="s3://genomon-bucket/_GRCh37/reference/GRCh37.STAR-2.5.2a"
        echo -e "${line}\t${fastq1}\t${fastq2}\t${output_dir}\t${reference}" >> ../output/sra/tasks/tasks-star-alignment-paired.tsv
    fi

done < ../output/sra/tasks/run_id_list_paired.txt

ecsub submit --script subscript_sra/star-alignment-paired.sh --tasks ../output/sra/tasks/tasks-star-alignment-paired.tsv --aws-s3-bucket s3://friend1ws-virginia-ecsub2 --wdir /tmp/ecsub --image genomon/star_alignment --aws-ec2-instance-type t2.2xlarge --disk-size 256 

echo -e "--env SAMPLE\t--input FASTQ1\t--output-recursive OUTPUT_DIR\t--input-recursive REFERENCE" > ../output/sra/tasks/tasks-star-alignment-single.tsv 
while read line
do
    echo ${line} >&2
    test1=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key star/${line}/${line}.Aligned.sortedByCoord.out.bam 2> /dev/null`
    cid1=$?
    test2=`aws s3api head-object --bucket friend1ws-virginia-ecsub2 --key star/${line}/${line}.Aligned.sortedByCoord.out.bam.bai 2> /dev/null`
    cid2=$?

    if [ $cid1 != 0 -o $cid2 != 0 ]
    then
        fastq1="s3://friend1ws-virginia-ecsub2/sra-toolkit/${line}/${line}_1.fastq"
        output_dir="s3://friend1ws-virginia-ecsub2/star/${line}"
        reference="s3://genomon-bucket/_GRCh37/reference/GRCh37.STAR-2.5.2a"
        echo -e "${line}\t${fastq1}\t${output_dir}\t${reference}" >> ../output/sra/tasks/tasks-star-alignment-single.tsv
    fi

done < ../output/sra/tasks/run_id_list_single.txt

ecsub submit --script subscript_sra/star-alignment-single.sh --tasks ../output/sra/tasks/tasks-star-alignment-single.tsv --aws-s3-bucket s3://friend1ws-virginia-ecsub2 --wdir /tmp/ecsub --image genomon/star_alignment --aws-ec2-instance-type t2.2xlarge --disk-size 256
_
