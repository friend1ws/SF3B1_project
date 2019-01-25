#! /usr/bin/env bash

python process1.py V87_38_TARGETEDSCREENMUTANT.csv > COSMIC_SF3B1_hotspot.txt.tmp

liftOver -bedPlus=3 COSMIC_SF3B1_hotspot.txt.tmp ../hg38ToHg19.over.chain COSMIC_SF3B1_hotspot.txt.tmp2 unmap.txt

python process2.py COSMIC_SF3B1_hotspot.txt.tmp2 > COSMIC_SF3B1_hotspot.txt

rm -rf COSMIC_SF3B1_hotspot.txt.tmp
rm -rf unmap.txt
rm -rf COSMIC_SF3B1_hotspot.txt.tmp2

