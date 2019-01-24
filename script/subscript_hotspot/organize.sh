#! /usr/bin/env bash

if [ ! -d ../output/sra/hotspot ]
then
    mkdir -p ../output/sra/hotspot
fi

tail +2 ../output/recount2/prediction/recount2.score.merged.txt | grep -v SRP012682 | awk '$3 > 5.0 {print}' - | sort -k3 -n -r > ../output/sra/hotspot/run_id_list2.txt


