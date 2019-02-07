#! /usr/bin/env python

import sys

input_file = sys.argv[1]
info_curated_file = sys.argv[2]

study2info = {}
with open(info_curated_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        study2info[F[0]] = F[1]

study2mutnum = {}
study2num = {}
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split('\t')
    for line in hin:
        F = {header[i]: elm for (i, elm) in enumerate(line.rstrip('\n').split('\t'))}
       
        if F["Mutation_Info"] == "---": continue
        if float(F["SF3B1ness_Score"]) < 100.0: continue
        if int(F["Bases"]) < 500000000: continue  
        if F["Study_Accession"] not in study2mutnum: study2mutnum[F["Study_Accession"]] = 0
        study2mutnum[F["Study_Accession"]] = study2mutnum[F["Study_Accession"]] + 1

        study2num[F["Study_Accession"]] = F["Study_Sample_Number"]

for study in sorted(study2mutnum):
    print study + '\t' + study2num[study] + '\t' + str(study2mutnum[study]) + '\t' + study2info[study]

"""
Study_accession Run_accession   SF3B1ness_score Mutation_Info   Bases   Study_Title     Study_Sample_number     Sample_Attribute
SRP063493       SRR2313066      4214.814        p.V701F,6,0.529 10777043880     Cancer associated SF3B1 hotspot mutations induce cryptic 3' splice site selection through use of a different branch point       72      source_name: Bone Marrow || biomaterial_provider: Weill Cornell Medicine || disease state: CLL || sf3b1 mutation: V701F
SRP063493       SRR2313063      3706.717        p.G742D,40,0.117;p.K700E,460,0.297      10717123860     Cancer associated SF3B1 hotspot mutations induce cryptic 3' splice site selection through use of a different branch point       72      source_name: Bone Marrow || biomaterial_provider: Weill Cornell Medicine || disease state: CLL || sf3b1 mutation: K700E|G742D
SRP050146       SRR1660313      3485.17 p.K700E,460,0.472       9826078800      RNA sequencing of bone marrow CD34+ cells from myelodysplastic syndrome patients with and without SF3B1 mutation and from healthy controls      17      source_name: CD34+ cells, MDS, SF3B1 mutated || cell type: bone marrow CD34+ cells || disease status: myelodysplastic syndromes (MDS) || mds subtype: RCMD-RS || sf3b1 status: mutated
SRP050146       SRR1660310      3345.079        p.R625L,25,0.488        16087041800     RNA sequencing of bone marrow CD34+ cells from myelodysplastic syndrome patients with and without SF3B1 mutation and from healthy controls      17      source_name: CD34+ cells, MDS, SF3B1 mutated || cell type: bone marrow CD34+ cells || disease status: myelodysplastic syndromes (MDS) || mds subtype: RCMD-RS || sf3b1 status: mutated
SRP050146       SRR1660312      3337.528        p.E622D,22,0.482        15371079400     RNA sequencing of bone marrow CD34+ cells from myelodysplastic syndrome patients with and without SF3B1 mutation and from healthy controls      17      source_name: CD34+ cells, MDS, SF3B1 mutated || cell type: bone marrow CD34+ cells || disease status: myelodysplastic syndromes (MDS) || mds subtype: RARS || sf3b1 status: mutate
"""
