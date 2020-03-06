#!/usr/bin/env python3 

'''
Simulate linked reads from long ONT reads
usage: gunzip -c [long_read_file.fq.gz] | python3 long-to-linked.py -l [cut length] | gzip > file_name.fq.gz

'''

import sys
import argparse

# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser( \
            description="Split long reads into pseudo-linked reads with a length of l.")
    parser.add_argument( "-l", "--length", \
            metavar="cut length", dest="cut_length", default=500, type=int, \
            help="Length for pseudo-linked reads to be cut")
    return parser.parse_args()

# Split long reads and output a list containing these reads
def split_long_read(read_string, length):
    return [read_string.strip()[0+i:l+i] for i in range(0, len(read_string), l)]

def cut_reads(length):
    
    # Read first line to determine file type
    line = sys.stdin.readline()
#    print(line, file=sys.stdout)

    # Fasta file type
    # Load header of first fasta entry
    if line[0] == ">":
        header = line.strip().replace(" ", "_")
        # Load remainder of reads
        for i, line in enumerate(sys.stdin):
            if (i + 1) % 2 == 0:
                header = line.strip().replace(" ", "_")
            if (i + 1) % 2 == 1: # Line containing read in fasta file
                split_reads = split_long_read(line, length) # Split long read into a list of shorter reads
                for read in split_reads:
                    print(header + "_" + str(split_reads.index(read)) + " BX:Z:" + str(int((i+1)/2)), file=sys.stdout)
                    print(read, file=sys.stdout)

    # Fastq file type
    # Load header of first fastq entry
    elif line[0] == "@":
        header = line.strip().replace(" ", "_")
        fa_header = ">" + header[1:]
#       print(header, file=sys.stdout)
    # Load remainder of reads
        for i, line in enumerate(sys.stdin):
            if (i % 4) == 3:
                header = line.strip().replace(" ", "_")
                fa_header = ">" + header[1:]
#               print(line, file=sys.stdout)
            if (i % 4) == 0:
                split_reads = split_long_read(line, length)
                for read in split_reads:
                    print(fa_header + "_" + str(split_reads.index(read)) + " BX:Z:" + str(int((i + 1) / 4)), file=sys.stdout)
                    print(read, file=sys.stdout)


if __name__ == "__main__":
    arguments = parse_arguments()
    l = arguments.cut_length
    cut_reads(l)

