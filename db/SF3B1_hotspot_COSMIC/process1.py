#! /usr/bin/env python

import sys, re

input_file = sys.argv[1]
mutation_aa_pos2count = {}
mutation_id2count = {}
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split(',')
    header = [x.strip() for x in header] 
    for line in hin:
        F = {header[i]: v for (i, v) in enumerate(line.rstrip('\n').split(','))}
        mutation_aa = F["MUTATION_AA"]
        match = re.search(r'^p.\w(\d+)\w$', mutation_aa)
        if match is None: continue
        mutation_aa_pos = match.group(1)
        if mutation_aa_pos not in mutation_aa_pos2count: mutation_aa_pos2count[mutation_aa_pos] = 0
        mutation_aa_pos2count[mutation_aa_pos] = mutation_aa_pos2count[mutation_aa_pos] + 1

        mutation_id = F["MUTATION_ID"]
        if mutation_id not in mutation_id2count: mutation_id2count[mutation_id] = 0
        mutation_id2count[mutation_id] = mutation_id2count[mutation_id] + 1


COMPLEMENT_DNA = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
mut_key2info = {}
with open(input_file, 'r') as hin:
    header = hin.readline().rstrip('\n').split(',')
    header = [x.strip() for x in header]
    for line in hin:
        F = {header[i]: v for (i, v) in enumerate(line.rstrip('\n').split(','))}
        mutation_aa = F["MUTATION_AA"]
        match = re.search(r'^p.\w(\d+)\w$', mutation_aa)
        if match is None: continue
        mutation_aa_pos = match.group(1)

        if int(mutation_aa_pos) > 801 or int(mutation_aa_pos) < 604: continue        
        # if mutation_aa_pos2count[mutation_aa_pos] < 5: continue

        mutation_genome_position = F["MUTATION_GENOME_POSITION"]
        match = re.search(r'^(\d+)\:(\d+)\-(\d+)$', mutation_genome_position)
        if match is None: continue
        genome_chr = 'chr' + match.group(1)
        genome_pos = match.group(2)

        mutation_cds = F["MUTATION_CDS"]
        match = re.search(r'^c.\d+(\w)>(\w)$', mutation_cds)
        if match is None: continue
        genome_ref = COMPLEMENT_DNA[match.group(1)]
        genome_alt = COMPLEMENT_DNA[match.group(2)]

        mutation_id = F["MUTATION_ID"]
        mutation_id_count = mutation_id2count[mutation_id]
        mut_key = '\t'.join([genome_chr, str(int(genome_pos) - 1), genome_pos, genome_ref, genome_alt])
        mut_key2info[mut_key] = '\t'.join([genome_chr, str(int(genome_pos) - 1), genome_pos, genome_ref, genome_alt, "SF3B1", mutation_cds, mutation_aa, mutation_id, str(mutation_id_count)])


for mut_key in sorted(mut_key2info):
    print mut_key2info[mut_key]

