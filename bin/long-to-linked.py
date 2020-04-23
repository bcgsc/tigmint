#!/usr/bin/env python3 
"""
Simulate linked reads from long ONT reads
Usage: gunzip -c [long_read_file.fa.gz] | python3 long-to-linked.py -l [cut length] | gzip > file_name.fa.gz
"""

import sys
import gzip
import argparse
from read_fasta import read_fasta

# Parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser( \
            description="Split long reads into pseudo-linked reads with a length of l.")
    parser.add_argument( "-l", "--length", \
            metavar="cut length", dest="cut_length", default=500, type=int, \
            help="Length for pseudo-linked reads to be cut")
    parser.add_argument("-r", "--reads", \
            metavar="long read file", dest="long_reads", required=True, \
            help="Long read file")
    parser.add_argument("-m", "--min_size",\
            metavar="minimum read size", dest="min_size", default=2000, type=int, \
            help="Minimum read length to be considered a molecule")
    return parser.parse_args()

# Split long reads and output a list containing these reads
def split_long_read(read_string, length):
    """
    Take an input read string and length l. Output a list containing substrings of read of length l.
    """
    return [read_string.strip()[0+i:l+i] for i in range(0, len(read_string), l)]

def cut_reads(length, long_reads, min_size):
    """ 
    Open gzipped long read file. Print cut reads to stdout.
    """
    with gzip.open(long_reads, 'rt') as long_reads:
        i = 0
        for header, seq, _, _ in read_fasta(long_reads):
            i += 1
            if len(seq) >= min_size:
                split_reads = split_long_read(seq, length)
                for read in split_reads:
                    print(">" + header + "_" + str(split_reads.index(read)) + " BX:Z:" + str(i), file=sys.stdout)
                    print(read, file=sys.stdout)

if __name__ == "__main__":
    arguments = parse_arguments()
    l = arguments.cut_length
    ont = arguments.long_reads
    min_size = arguments.min_size
    cut_reads(l, ont, min_size)

