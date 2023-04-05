#!/usr/bin/env python3
import sys
import gzip
import re

def fastq_to_allele_table(fastq_file, output):
    f1 = gzip.open(fastq_file)
    f2 = open(output, "w")
    line1 = f1.readline()
    total_line = 0
    counter = 0
    
    while line1:
        total_line += 1
        name = line1.decode("utf-8")
        bc = (name.split(",")[0]).split("@")[1]
        
        line1 = f1.readline()
        seq = line1.decode("utf-8")
        if re.search(r'CTTCCTTTCC(.+?)GTCTG', seq): #edited, make sure there is something between the handles
            target = re.search(r'CTTCCTTTCC(.+?)GTCTG', seq).group(1)
            f2.write(reverse_complement(bc) + "," + reverse_complement(target) + "\n")
            counter += 1
        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()
    f1.close()
    f2.close()
    print("sample name: %s, total line: %f, line subjected to downstream edit analysis: %f" %(fastq_file, total_line, counter))

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])


if __name__ == "__main__":
    fastq_file = sys.argv[1]
    output = sys.argv[2]
    fastq_to_allele_table(fastq_file, output)