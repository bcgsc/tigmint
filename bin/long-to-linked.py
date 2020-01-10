#!/usr/bin/env python3 

'''
Simulate linked reads from long ONT reads
usage: python3 long-to-linked.py long_reads.fa

'''

import sys

# Split long reads and output a list containing these reads
def split_long_read(read_string, n=500):
    return [read_string.strip()[0+i:n+i] for i in range(0, len(read_string), n)]

#print(split_long_read("abcdefghijklmnopqrs"))

# Load long ONT reads ** add here
with open(sys.argv[1],'r') as long_reads:
    for i, line in enumerate(long_reads):
        if (i + 1) % 2 != 0:
            header = line.strip()
        if (i + 1) % 2 == 0: # Line containing read in fasta file
            split_reads = split_long_read(line) # Split long read into a list of shorter reads
            for read in split_reads:
                print(header.replace(" ", "") + "_" + str(split_reads.index(read)) + " BX:Z:" + str(int((i+1)/2)))
                print(read)
