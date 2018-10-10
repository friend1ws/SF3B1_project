#! /usr/bin/env python

import sys, os, glob

sj_list_file = sys.argv[1]
savnet_input_file_dir = sys.argv[2]

# target_cancer_type = ["BLCA", "BRCA", "LUAD", "SKCM", "UVM"]
target_cancer_type = ["BRCA", "SKCM", "UVM"]

sj_info = {}
sj_info2 = {}
with open(sj_list_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        key = F[0] + ':' + F[1] + '-' + F[2]
        sj_info[key] = ';'.join(F)
        
        if F[5] == '+':
            key2 = F[0] + ':' + F[1] + '-' + F[4]
            sj_info2[key2] = ';'.join(F)
        else:
            key2 = F[0] + ':' + F[4] + '-' + F[2]
            sj_info2[key2] = ';'.join(F)
   

print "Cancer_Type" + '\t' + "Sample_Name" + '\t' + "Weight" + '\t' + "Mutation_Info" + '\t' + "Splicing_Key" + '\t' + "Read_Count1" + '\t' + "Read_Count2"

all_savnet_input_file = glob.glob(savnet_input_file_dir + "/*.input.txt")
for savnet_input_file in sorted(all_savnet_input_file):

    cancer_type = os.path.basename(savnet_input_file).replace(".input.txt", '')
    # if cancer_type not in target_cancer_type: continue

    info2count = {}
    with open(savnet_input_file, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            sample_name = F[header2ind["Sample_Name"]]
            sj_file = F[header2ind["SJ_File"]]
            mutation_info = F[header2ind["Mutation_Info"]]
            weight = F[header2ind["Weight"]]

            info2count = {}
            sj_file = sj_file.replace("/sshare3", "")
            print >> sys.stderr, "Procesing: " + sj_file
            
            with open(sj_file, 'r') as hin2:
                for line2 in hin2:
                    F2 = line2.rstrip('\n').split('\t')
                    sj_key = F2[0] + ':' + F2[1] + '-' + F2[2]
                    if sj_key in sj_info:
                        info = cancer_type + '\t' + sample_name + '\t' + weight + '\t' + mutation_info + '\t' + sj_info[sj_key]
                        if info not in info2count: info2count[info] = [0, 0]
                        info2count[info][0] = F2[6]
                    if sj_key in sj_info2:
                        info = cancer_type + '\t' + sample_name + '\t' + weight + '\t' + mutation_info + '\t' + sj_info2[sj_key]
                        if info not in info2count: info2count[info] = [0, 0]
                        info2count[info][1] = F2[6]

            for info in info2count:
                count = info2count[info]
                print info + '\t' + str(count[0]) + '\t' + str(count[1])
            # print cancer_type + '\t' + sample_name + '\t' + weight + '\t' + mutation_info + '\t' + sj_info[sj_key] + '\t' + F2[6]



