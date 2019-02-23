#! /bin/bash
#$ -S /bin/bash
#$ -cwd

if [ ! -d ../output/sra/hotspot/s3_hotspot ]
then
    mkdir -p ../output/sra/hotspot/s3_hotspot
fi

CONTROL_BAM=../data/SF3B1_control.bam

echo -n > ../output/sra/hotspot/bam_list.txt
for bfile in `aws s3 ls s3://friend1ws-virginia-ecsub2/star/ --recursive | grep Aligned.sortedByCoord.out.bam | grep -v bam.bai | tr -s ' ' | cut -f 4 -d ' '`
do
    echo "s3://friend1ws-virginia-ecsub2/${bfile}" >> ../output/sra/hotspot/bam_list.txt
done


echo -n > ../output/sra/hotspot/SF3B1.hotspot.result.txt
while read TUMOR_BAM 
do
    # BFILE=${TUMOR_BAM##*/}
    BFILE=`basename $TUMOR_BAM`
    SAMPLE_NAME=${BFILE%%.Aligned.sortedByCoord.out.bam}
    echo ${SAMPLE_NAME}

    echo "samtools view -bh ${TUMOR_BAM} 2:198256698-198299817 > ../output/sra/hotspot/tmp.bam"
    samtools view -bh ${TUMOR_BAM} 2:198256698-198299817 > ../output/sra/hotspot/tmp.bam
    samtools index ../output/sra/hotspot/tmp.bam

    echo "hotspotCall -t 0.03 ../output/sra/hotspot/tmp.bam ${CONTROL_BAM} ../output/sra/hotspot/tmp1.txt ../db/SF3B1_hotspot_COSMIC/COSMIC_SF3B1_hotspot.txt"
    hotspotCall -t 0.03 ../output/sra/hotspot/tmp.bam ${CONTROL_BAM} ../output/sra/hotspot/tmp1.txt ../db/SF3B1_hotspot_COSMIC/COSMIC_SF3B1_hotspot.txt 

    sed -e '1d' ../output/sra/hotspot/tmp1.txt > ../output/sra/hotspot/tmp2.txt
    cat ../output/sra/hotspot/tmp2.txt

    while read line
    do
        echo -e "${SAMPLE_NAME}\t${line}" >> ../output/sra/hotspot/SF3B1.hotspot.result.txt
    done < ../output/sra/hotspot/tmp2.txt

    rm -rf ${BFILE}.bai
    rm -rf ../output/sra/hotspot/tmp.bam
    rm -rf ../output/sra/hotspot/tmp.bam.bai
    rm -rf ../output/sra/hotspot/tmp1.txt
    rm -rf ../output/sra/hotspot/tmp2.txt

done < ../output/sra/hotspot/bam_list.txt

