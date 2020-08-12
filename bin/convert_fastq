#!/usr/bin/env python3
"""
Convert long reads in fastq to fasta.
Usage: python3 convert-fastq.py infile.fq.gz | gzip > outfile.fa.gz
"""

import sys
import gzip
from read_fasta import read_fasta

def convert_fq_fa(fastq):
    for header, seq, _, _ in read_fasta(fastq):
        print(">", header, "\n", seq, sep="", file=sys.stdout)
    return

if __name__ == "__main__":
    with gzip.open(sys.argv[1], "rt") as in_fastq:
        convert_fq_fa(in_fastq)
