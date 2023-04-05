#!/usr/bin/env python3
import re
import pandas as pd
import operator
import edlib
import numpy as np
import sys

#This script calculate edits plus/minus 10 bp around the cut site

def get_alignment(query, subject, mode):
    path = edlib.align(query, subject, mode = mode, task = "path")
    aln = edlib.getNiceAlignment(path, query, subject)
    return(aln)

def edit_cal(query, subject):
    edit_outcome = []
    alignment = get_alignment(query, subject, mode = "NW")
    m_align = alignment["matched_aligned"]
    q_align = alignment["query_aligned"]
    t_align = alignment["target_aligned"]
    shift = 0
    len_PB = 11 #distance from ref start to the end of guide RNA
    
    if m_align == "|"*len(subject):
        edit_outcome.append("wt")
    else:   
        for pos in range(len(m_align)):
            if m_align[pos:pos+1] == ".":
                #re-calibrate mapping position in case of insertion before the cut site
                code = "sub"
                re_cal = len_PB + shift - pos #re-calibrated position
                edit_outcome.append("sub" + t_align[pos:pos+1] + '{0:+d}'.format(-re_cal) + q_align[pos:pos+1])
            elif m_align[pos:pos+1] == "-":
                if t_align[pos:pos+1] == "-":
                    shift += 1
                    re_cal = len_PB + shift - 1 - pos
                    edit_outcome.append("ins" + '{0:+d}'.format(-re_cal) + q_align[pos:pos+1])
                elif q_align[pos:pos+1] == "-":
                    re_cal = len_PB + shift - pos
                    edit_outcome.append("del" + '{0:+d}'.format(-re_cal) + t_align[pos:pos+1])
    
    return(",".join(edit_outcome))

def edit_process(input_file, output_file):
    ref = "TGAGCACGTGATGGCAGA"
    input = open(input_file)
    output = open(output_file, "w")
    output.write("barcode\tsequence\tedits\n")
    line = input.readline()
    
    while (line):
        line = line.rstrip("\n")
        bc = line.split(",")[0]
        seq = line.split(",")[1]
        edit = edit_cal(seq, ref)
        output.write(bc + "\t" + seq + "\t" + edit + "\n")
        
        line = input.readline()
    
    input.close()
    output.close()
    



if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    edit_process(input_file, output_file)