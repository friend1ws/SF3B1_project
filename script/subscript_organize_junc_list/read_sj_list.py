#! /usr/bin/env python

import xlrd, pysam

junc_tb = pysam.TabixFile("../db/junction.bed.gz")
# key_exists = {}

def junc_check(tchr, tstart, tend, tgene):

    tabix_error_flag = 0
    try:
        junc_lines = junc_tb.fetch(tchr, tstart - 50, tend + 50)
    except Exception as inst:
        tabix_error_flag = 1

    is_annotated = 0
    return_key = None
    if tabix_error_flag == 0 and junc_lines is not None:
        for junc_line in junc_lines:

            key = '\t'.join([tchr, str(tstart), str(tend)])
            # if key in key_exists: continue

            junction = junc_line.split('\t')
            junction[1], junction[2] = int(junction[1]) + 1, int(junction[2]) 
            if tstrand == '+' and junction[5] == '+' and tchr == junction[0] and tstart == junction[1]:
                if tend > junction[2] - 50 and tend < junction[2]:
                    return_key = '\t'.join([tchr, str(tstart), str(tend), tgene, str(junction[2]), tstrand])
                    # key_exists[key] = 1

            if tstrand == '-' and junction[5] == '-' and tchr == junction[0] and tend == junction[2]:
                if tstart > junction[1] and tstart < junction[1] + 50:
                    return_key = '\t'.join([tchr, str(tstart), str(tend), tgene, str(junction[1]), tstrand])
                    # key_exists[key] = 1

            if tstart == junction[1] and tend == junction[2]: is_annotated = 1

    if is_annotated == 0:
        return return_key
    else:
        return None


junc_exists = {}
book = xlrd.open_workbook("../data/SF3B1_junction_list/ncomms10615-s2.xlsx")
sheet = book.sheet_by_index(0)
for row_index in range(2, sheet.nrows):
    row = sheet.row(row_index)
    if row[8].value == "Acceptor":
        tchr, tstart, tend, tstrand, tgene = [row[i].value for i in [0, 1, 2, 3, 4]]
        tchr, tstart, tend = tchr.replace('chr', ''), int(tstart) + 1, int(tend)
        junc_key = junc_check(tchr, tstart, tend, tgene)
        junc_key2 = '\t'.join([tchr, str(tstart), str(tend)])
        if junc_key is not None and junc_key2 not in junc_exists: 
            print junc_key
            junc_exists[junc_key2] = 1    


book = xlrd.open_workbook("../data/SF3B1_junction_list/mmc3.xls")
sheet = book.sheet_by_index(0)
for row_index in range(2, sheet.nrows):
    row = sheet.row(row_index)
    if row[7].value == "True":
        tchr, tstart, tend, tstrand, tgene = [row[i].value for i in [1, 2, 3, 4, 5]]
        tchr, tstart, tend = tchr.replace('chr', ''), int(tstart) + 1, int(tend)
        junc_check(tchr, tstart, tend, tgene) 
        junc_key = junc_check(tchr, tstart, tend, tgene)
        junc_key2 = '\t'.join([tchr, str(tstart), str(tend)])
        if junc_key is not None and junc_key2 not in junc_exists: 
            print junc_key
            junc_exists[junc_key2] = 1
    



