#!/usr/bin/env python3
import sys
import gzip
import re

#In this sequencing run, R1 >= 118 bp
#This script extracts barcodes from Read 1 and attach them to the read names in Read 1.

def barcode_attach(sample, input_folder, output_folder):
    Read1 = input_folder + "/" + sample + "_R1_001.fastq.gz"
    output_file = output_folder + "/" + sample + "_bc.fastq.gz"
    f0 = gzip.open(Read1)
    f1 = gzip.open(Read1)
    f2 = gzip.open(output_file, 'wb')
    
    line1 = f1.readline() #read name line
    total_line = 0
    bc_line = 0
    
    while (line1):
        total_line += 1
        line0 = f0.readline() #to the title line
        line1 = f1.readline() #to the sequence line

        seq1 = line1.decode("utf-8")
        
        if re.search(r'GCTTCTCGTC(.{16})GCCCTCTGGA', seq1):
            bc_line += 1
            bc = re.search(r'GCTTCTCGTC(.{16})GCCCTCTGGA', seq1).group(1)
            
            first_line = '@'.encode('utf-8') + bc.encode('utf-8') + ','.encode('utf-8') + line0[1:]
            f2.write(first_line)

            second_line = f0.readline()
            f2.write(second_line)

            third_line = f0.readline()
            f2.write(third_line)

            fourth_line = f0.readline()
            f2.write(fourth_line)
        else:
            line0 = f0.readline()
            line0 = f0.readline()
            line0 = f0.readline()


        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()


    f0.close()
    f1.close()
    f2.close()
    print("sample name: %s, total line: %f, line with barcode: %f" %(sample, total_line, bc_line))
    
    
if __name__ == "__main__":
    sample = sys.argv[1]
    input_folder = sys.argv[2]
    output_folder = sys.argv[3]
    barcode_attach(sample, input_folder, output_folder)