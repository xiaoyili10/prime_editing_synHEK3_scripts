#!/usr/bin/env python3
import sys
import re

#This script remove reads that do not align well at the TTAA sites or align to several locations.

def read_filtering(samfile):
    
    f1 = open(samfile)
    shortname = samfile.split('.')[0]
    f2 = open(shortname + '.sam', 'w')
    pass_filter = 0
    total = 0
    
    for line in f1:
        if (line[0] == '@'):
            f2.write(line)
        else:
            total += 1
            #name = (((line.split('\t'))[0]).split(','))
            flag = (line.split('\t'))[1]
            cigar = (line.split('\t'))[5]
            read = (line.split('\t'))[9]
            tags = (line.split('\t'))[11:]
            
            
            if not re.search("XA:Z:", " ".join(tags)):
                if int(flag) == 0:
                    if re.search(r"TTAA$", read):
                        if re.search(r"[0-9]+M$", cigar):
                            match_len = re.split('M|I|D|N|S|H|P|=|X', cigar)[-2]
                            if int(match_len) >=4:
                                pass_filter += 1
                                f2.write(line)
                elif int(flag) == 16:
                    if re.search(r"^TTAA", read):
                        if re.search(r"^[0-9]+M", cigar):
                            match_len = cigar.split('M')[0]
                            if int(match_len) >=4:
                                pass_filter += 1
                                f2.write(line)
    print("sample name: %s, total line: %f, line pass filter: %f" %(shortname, total, pass_filter))
    f1.close()
    f2.close()

if __name__ == "__main__":
    samfile = sys.argv[1]
    read_filtering(samfile)