#! /user/bin/env bash

python make_hotspot.py | sort -k1,1 -k2,2n -k3,3n > GRCh37_SF3B1_hotspot_omega.txt 

