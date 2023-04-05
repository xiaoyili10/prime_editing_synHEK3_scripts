#!/usr/bin/env python3
import re
import pandas as pd
import sys

#This script computes editing status of all barcodes detected
#It further separates CTT and CTT with error
def count_edits(filename, output):
    barcode_table = pd.read_csv(filename, sep = "\t")
    barcode_ls = barcode_table["barcode"].tolist()
    barcode_ls = set(barcode_ls)
    
    dic_wt = {}
    dic_ctt = {}
    dic_ctt_err = {}
    dic_other = {}
    dic_count = {}
    count = 0

    for bc in barcode_ls:
        dic_wt[bc] = 0
        dic_ctt[bc] = 0
        dic_ctt_err[bc] = 0
        dic_other[bc] = 0
        dic_count[bc] = 0
    
    edits_table = open(filename, "r")
    line = edits_table.readline() #skip the header line
    line = edits_table.readline()
    count += 1
    
    while (line):
        line = line.rstrip()
        barcode = line.split("\t")[0]
        if barcode in barcode_ls:
            dic_count[barcode] += 1
            edits = line.split("\t")[2]
            if edits == "wt":
                dic_wt[barcode] += 1
            elif re.search(r"^ins-3C,ins-2T,ins-2T$", edits):
                dic_ctt[barcode] += 1
            elif re.search("ins-3C,ins-2T,ins-2T", edits):
                dic_ctt_err[barcode] += 1
            else:
                dic_other[barcode] += 1
    
        line = edits_table.readline()
        count += 1

    edits_table.close()
    
    
    with open(output, 'w') as outfile:
        outfile.write("barcode\ttotal\tunedited\tctt_ins\tctt_ins_error\tother\n")
        for key in barcode_ls:
            outfile.write(key + "\t" + str(dic_count[key]) + "\t" + str(dic_wt[key]) + "\t" + str(dic_ctt[key]) + "\t" + str(dic_ctt_err[key]) + "\t" + str(dic_other[key]) + "\n")
    outfile.close()
    
    return

if __name__ == "__main__":
    filename = sys.argv[1]
    output = sys.argv[2]
    count_edits(filename, output)