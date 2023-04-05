#!/usr/bin/env python3
import sys
import re
import pandas as pd
import operator
from collections import Counter

regular_chr = list(range(1,23)) + ["X", "Y", "M"]
chrom = ["chr" + str(s) for s in regular_chr]

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def clean_up(input):
    bedfile = open(input, "r")
    shortname = input.split(".")[0]
    temp = open(shortname + ".txt.temp", 'w')
    
    for line in bedfile:
        chr, start, end, name = line.split("\t")[0:4]
        strand = line.split("\t")[5]
        barcode = name.split(",")[0]
        cigar = line.split("\t")[7]
        seq = line.split("\t")[11]
        if strand == "+":
            pos = end
        elif strand == "-":
            pos = start
            seq = reverse_complement(seq)
        
        umi = seq[0:8]
        
        if chr in chrom: #Keep only reads mapped to standard chromosomes
            temp.write(barcode + "\t" + chr + "\t" + pos + "\t" + start + "\t" + end + "\t" + strand + "\t" + cigar + "\t" + umi + "\t" + seq + "\n")
    temp.close()
    bedfile.close()

def sort_data(input):
    shortname = input.split(".")[0]
    filename = shortname + ".txt.temp"
    df = pd.read_table(filename, header = None, sep = "\t")
    df.columns = ["barcode", "chr", "pos", "start", "end", "strand", "cigar", "UMI", "seq"]
    df = df.sort_values(["chr", "pos", "strand", "barcode", "UMI"], ascending = (True, True, True, True, True))
    
    df.to_csv(shortname + ".sort.txt.temp", sep = "\t", header = False, index = False)
    

def duplicate_removal(input):
    shortname = input.split(".")[0]
    filename = open(shortname + ".sort.txt.temp", 'r')
    output = open(shortname + ".out.txt", 'w')
    
    barcode_dic = {}

    # Read in the first line and store information in prev_ variables
    line = filename.readline()
    barcode = line.split("\t")[0]
    prev_chr, prev_pos = line.split("\t")[1:3]
    prev_str = line.split("\t")[5]
    barcode_dic[barcode] = 1
    
    integration_count = 1
    dup_count = 1
    
    line = filename.readline()
    
    while (line):
        barcode, chr, pos = line.split("\t")[0:3]
        strand = line.split("\t")[5]

        if chr == prev_chr and pos == prev_pos and strand == prev_str:
            dup_count += 1
            if barcode in barcode_dic.keys():
                barcode_dic[barcode] += 1
            else:
                barcode_dic[barcode] = 1
            line = filename.readline()
        else:
            barcode_dic_sort = dict(sorted(barcode_dic.items(), key=operator.itemgetter(1), reverse=True))
            
            for key in barcode_dic_sort.keys():
                output.write(prev_chr + "\t" + prev_pos + "\t" + prev_str + "\t" + str(dup_count) + "\t" + key + "\t" + str(barcode_dic_sort[key]) + "\n")
            barcode_dic = {}
            barcode_dic[barcode] = 1
            
            integration_count += 1
            dup_count = 1
            
            prev_chr = chr
            prev_pos = pos
            prev_str = strand
            
            line = filename.readline()
    
    barcode_dic_sort = dict(sorted(barcode_dic.items(), key=operator.itemgetter(1), reverse=True))
    for key in barcode_dic_sort.keys():
        output.write(prev_chr + "\t" + prev_pos + "\t" + prev_str + "\t" + str(dup_count) + "\t" + key + "\t" + str(barcode_dic_sort[key]) + "\n")
    
    filename.close()
    output.close()

def UMI_collapse(input):
    shortname = input.split(".")[0]
    filename = open(shortname + ".sort.txt.temp", 'r')
    output = open(shortname + ".out.txt", 'w')
    
    # Read in the first line and store information in prev_ variables
    line = filename.readline()
    prev_barcode, prev_chr, prev_pos = line.split("\t")[0:3]
    prev_str = line.split("\t")[5]
    umi = line.split("\t")[7]
    umi_ls = [umi]
    #barcode_dic[barcode] = 1
    #integration_count = 1
    dup_count = 1
    
    line = filename.readline()
    
    while (line):
        barcode, chr, pos = line.split("\t")[0:3]
        strand = line.split("\t")[5]
        umi = line.split("\t")[7]

        if chr == prev_chr and pos == prev_pos and strand == prev_str and barcode == prev_barcode:
            dup_count += 1
            umi_ls.append(umi)
            line = filename.readline()
        else:
            umi_dic = Counter(umi_ls)
            
            output.write(prev_chr + "\t" + prev_pos + "\t" + prev_str + "\t" + str(dup_count) + "\t" + prev_barcode + "\t" + str(len(umi_dic)) + "\n")
                
            # reset counter
            umi_ls = [umi]
            
            dup_count = 1
            
            prev_barcode = barcode
            prev_chr = chr
            prev_pos = pos
            prev_str = strand
            
            line = filename.readline()
    
    umi_dic = Counter(umi_ls)
    output.write(prev_chr + "\t" + prev_pos + "\t" + prev_str + "\t" + str(dup_count) + "\t" + prev_barcode + "\t" + str(len(umi_dic)) + "\n")
    
    filename.close()
    output.close()

if __name__ == "__main__":
    input = sys.argv[1]
    key = sys.argv[2]
    clean_up(input)
    sort_data(input)
    if key == "barcode":
        duplicate_removal(input)
    elif key == "UMI":
        UMI_collapse(input)
    else:
        print("Please enter \"barcode\" or \"UMI\" after file name.")