#!/bin/bash

set -o errexit
set -o nounset

OUTPUT_PREF=${OUTPUT_DIR}/${SAMPLE}
mkdir -p ${OUTPUT_DIR}

/usr/local/bin/STAR \
    --genomeDir ${REFERENCE} \
    --readFilesIn ${FASTQ1} ${FASTQ2} \
    --outFileNamePrefix ${OUTPUT_PREF}. \
    --runThreadN 6 --outSAMstrandField intronMotif --outSAMunmapped Within --alignMatesGapMax 500000 --alignIntronMax 500000 --alignSJstitchMismatchNmax -1 -1 -1 -1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --outSAMtype BAM Unsorted
 

rm $FASTQ1
rm $FASTQ2

/usr/local/bin/samtools sort \
    -T ${OUTPUT_PREF}.Aligned.sortedByCoord.out \
    -@ 6 -m 3G \
    ${OUTPUT_PREF}.Aligned.out.bam \
    -O bam > ${OUTPUT_PREF}.Aligned.sortedByCoord.out.bam

/usr/local/bin/samtools index \
    ${OUTPUT_PREF}.Aligned.sortedByCoord.out.bam

rm ${OUTPUT_PREF}.Aligned.out.bam

