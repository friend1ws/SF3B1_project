#! /bin/bash
#$ -S /bin/bash
#$ -cwd

if [ ! -d output_s3 ]
then
    mkdir -p output_s3
fi

CONTROL_BAM=SF3B1_control.bam

echo -n > output_s3/bam_list.txt
for bfile in `aws s3 ls s3://friend1ws-virginia-ecsub2/star/ --recursive | grep Aligned.sortedByCoord.out.bam | grep -v bam.bai | tr -s ' ' | cut -f 4 -d ' '`
do
    echo "s3://friend1ws-virginia-ecsub2/${bfile}" >> output_s3/bam_list.txt
done


echo -n > output_s3/SF3B1.hotspot.result.txt
while read TUMOR_BAM 
do
    # BFILE=${TUMOR_BAM##*/}
    BFILE=`basename $TUMOR_BAM`
    SAMPLE_NAME=${BFILE%%.Aligned.sortedByCoord.out.bam}
    echo ${SAMPLE_NAME}

    echo "samtools view -bh ${TUMOR_BAM} 2:198256698-198299817 > output_s3/tmp.bam"
    samtools view -bh ${TUMOR_BAM} 2:198256698-198299817 > output_s3/tmp.bam
    samtools index output_s3/tmp.bam

    echo "hotspotCall -t 0.03 output_s3/tmp.bam ${CONTROL_BAM} output_s3/tmp1.txt ../db/SF3B1_hotspot_COSMIC/COSMIC_SF3B1_hotspot.txt"
    hotspotCall -t 0.03 output_s3/tmp.bam ${CONTROL_BAM} output_s3/tmp1.txt ../db/SF3B1_hotspot_COSMIC/COSMIC_SF3B1_hotspot.txt 
    sed -e '1d' output_s3/tmp1.txt > output_s3/tmp2.txt
    cat output_s3/tmp2.txt

    while read line
    do
        echo -e "${SAMPLE_NAME}\t${line}" >> output_s3/SF3B1.hotspot.result.txt
    done < output_s3/tmp2.txt

done < output_s3/bam_list.txt

