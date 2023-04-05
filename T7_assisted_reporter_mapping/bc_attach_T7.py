#!/usr/bin/env python3
import sys
import gzip
import re

#This script add barcode from R2 to R1 names (new libraries)

def UMI_attach(sample, input_folder, output_folder):
    #open the read1, read2, and output file
    Read1 = input_folder + "/" + sample + "_R1_001.fastq.gz"
    Read2 = input_folder + "/" + sample + "_R2_001.fastq.gz"
    output_file = output_folder + "/" + sample + "_R1_001.bc.fastq.gz"
    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    f3 = gzip.open(output_file, 'wb')
    
    line1 = f1.readline() #read name line
    line2 = f2.readline() #read name line
    total_line = 0
    barcode_line = 0
    
    while (line2):
        total_line += 1
        line2 = f2.readline() #to the sequence line

        seq = line2.decode("utf-8") 
        if re.search(r'AGGGC(.{16})GACGA', seq):
            barcode_line += 1
            bc = re.search(r'AGGGC(.{16})GACGA', seq).group(1)

            first_line = '@'.encode('utf-8') + bc.encode('utf-8') + ','.encode('utf-8') + line1[1:]
            f3.write(first_line)

            second_line = f1.readline()
            f3.write(second_line)

            third_line = f1.readline()
            f3.write(third_line)

            four_line = f1.readline()
            f3.write(four_line)
            
        else:
            line1 = f1.readline() #error here, should be line2 = f2.readline()
            line1 = f1.readline() #error
            line1 = f1.readline() #error

        line1 = f1.readline()

        line2 = f2.readline()
        line2 = f2.readline()
        line2 = f2.readline()

    f1.close()
    f2.close()
    f3.close()
    print("sample name: %s, total line: %f, lines with barcode: %f" %(sample, total_line, barcode_line))
    
if __name__ == "__main__":
    sample = sys.argv[1]
    input_folder = sys.argv[2]
    output_folder = sys.argv[3]
    UMI_attach(sample, input_folder, output_folder)