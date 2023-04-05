#!/usr/bin/env python3
import sys
import glob
import pandas as pd
import re
import numpy as np
from Levenshtein import distance
from scipy.spatial.distance import pdist, squareform
import operator

def edit_distance_cal(barcode_dic):
    barcode_ls = list(barcode_dic.keys())
    transformed_strings = np.array(barcode_ls).reshape(-1,1)
    distance_matrix = pdist(transformed_strings,lambda x,y: distance(x[0],y[0]))
    distance_matrix = squareform(distance_matrix)
    return(distance_matrix)

def collapse_barcode(barcode_dic, minimal_dist, minimal_counts):
    new_dic = {}
    temp_dic = barcode_dic
    while bool(temp_dic):
        bc_ls = list(temp_dic.keys())
        top_bc = bc_ls[0]
        select = edit_distance_cal(temp_dic)[0]
        to_collapse = list(np.array(bc_ls)[select <= minimal_dist])
        new_dic[top_bc] = sum({key: temp_dic[key] for key in to_collapse}.values())
        [temp_dic.pop(key) for key in to_collapse]
    
    # remove barcodes with a few reads
    new_dic = {key: value for key, value in new_dic.items() if value > minimal_counts}
    # remove barcodes with N
    new_dic = {key: value for key, value in new_dic.items() if not re.search("N", key)}
    return(new_dic)
    

def barcode_distance_correct(file):
    df = pd.read_csv(file, sep = "\t", header = None)
    df.columns = ["chr", "pos", "strand", "total", "barcode", "count"]
    shortname = file.split(".out.txt")[0]
    output = open(shortname + "_corrected.txt", "w")
    first_row = ["chr", "pos", "strand", "barcode", "count"]
    output.write('\t'.join(first_row) + "\n")

    prev_chrom, prev_pos, prev_str, ct, bc = df.loc[0, ["chr", "pos", "strand", "count", "barcode"]]
    dic_bc = {}
    dic_bc[bc] = ct
    for i in range(len(df)):
        chrom, pos, strand, ct, bc = df.loc[i, ["chr", "pos", "strand", "count", "barcode"]]
        if chrom == prev_chrom and pos == prev_pos and strand == prev_str:
            dic_bc[bc] = ct
        else:
            dic_bc_sort = dict(sorted(dic_bc.items(), key=operator.itemgetter(1), reverse=True))
            dic_bc_collapse = collapse_barcode(dic_bc_sort, minimal_dist = 2, minimal_counts = 3)
        
            current_entry = [prev_chrom] + [prev_pos] + [prev_str]
            for key in dic_bc_collapse.keys():
                barcode_info = [key] + [dic_bc_collapse[key]]
                to_append = current_entry + barcode_info
                output.write('\t'.join(str(x) for x in to_append) + "\n")
        
            prev_chrom = chrom
            prev_pos = pos
            prev_str = strand
        
            dic_bc = {}
            dic_bc[bc] = ct
            
    dic_bc_sort = dict(sorted(dic_bc.items(), key=operator.itemgetter(1), reverse=True))
    dic_bc_collapse = collapse_barcode(dic_bc_sort, minimal_dist = 2, minimal_counts = 3)
    current_entry = [prev_chrom] + [prev_pos] + [prev_str]
    for key in dic_bc_collapse.keys():
        barcode_info = [key] + [dic_bc_collapse[key]]
        to_append = current_entry + barcode_info
        output.write('\t'.join(str(x) for x in to_append) + "\n")
        
    output.close()

if __name__ == "__main__":
    sample = sys.argv[1]
    barcode_distance_correct(sample)
    
    
    
    
    