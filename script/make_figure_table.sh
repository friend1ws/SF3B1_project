#! /bin/bash


echo -e "Sample_Name\tCancer_Type\tMutation_Status\tSF3B1ness_Score" > ../figure/TableS1.txt
sort -k2,2 -k1,1 ../output/recount2/TCGA/TCGA.pred.result.recount2.zibb.txt >> ../figure/TableS1.txt

