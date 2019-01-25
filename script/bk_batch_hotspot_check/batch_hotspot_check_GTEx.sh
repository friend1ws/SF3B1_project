#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log/ -o log/

if [ ! -d ../output/hotspot/GTEx ]
then
    mkdir -p ../output/hotspot/GTEx
fi

CONTROL_BAM=/share/portal/GTEx/rna/Blood.Whole_Blood/hg19/alignment/star/GTEX-WFON-0005-SM-5S2RV/GTEX-WFON-0005-SM-5S2RV.Aligned.sortedByCoord.out.bam

ls /share/portal/GTEx/rna/*/hg19/alignment/star/*/*.Aligned.sortedByCoord.out.bam > ../output/hotspot/GTEx/bam_list.txt

echo -n > ../output/hotspot/GTEx/SF3B1.hotspot.result.txt
while read TUMOR_BAM 
do
    BFILE=`basename ${TUMOR_BAM}`
    SAMPLE_NAME=${BFILE%%.Aligned.sortedByCoord.out.bam}

    echo "hotspotCall -t 0.02 ${TUMOR_BAM} ${CONTROL_BAM} ../output/hotspot/GTEx/tmp1.txt ../db/GRCh37_SF3B1_hotspot_omega.txt"
    hotspotCall -t 0.02 ${TUMOR_BAM} ${CONTROL_BAM} ../output/hotspot/GTEx/tmp1.txt ../db/GRCh37_SF3B1_hotspot_omega.txt 
    sed -e '1d' ../output/hotspot/GTEx/tmp1.txt > ../output/hotspot/GTEx/tmp2.txt
    cat ../output/hotspot/GTEx/tmp2.txt

    while read line
    do
        echo -e "${SAMPLE_NAME}\t${line}" >> ../output/hotspot/GTEx/SF3B1.hotspot.result.txt
    done < ../output/hotspot/GTEx/tmp2.txt

done < ../output/hotspot/GTEx/bam_list.txt
 
