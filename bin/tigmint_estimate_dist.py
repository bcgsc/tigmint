#!/usr/bin/env python3

'''
Sample the long reads, and estimate the optimal 'dist' parameter setting
@author: Lauren Coombe (lcoombe)
'''

import argparse
import numpy as np
from read_fasta import read_fasta

def get_n_read_lengths(reads_filename, num_reads, lower_bound):
    "Collect the read lengths for n reads"
    read_lengths = []
    read_count = 0
    with open(reads_filename, 'r') as fin:
        for _, seq, _, _ in read_fasta(fin):
            seq_len = len(seq)
            if seq_len > lower_bound:
                read_lengths.append(seq_len)
                read_count += 1
            if read_count >= num_reads:
                break
    return read_lengths


def main():
    "Read the first n reads, and estimate the optimal 'dist' setting based on the percentile"
    parser = argparse.ArgumentParser(description="Sample the reads from the given fasta/fastq file "
                                                 "and estimate the Tigmint-long dist parameter "
                                                 "based on the chosen length percentile")
    parser.add_argument("READS", help="Input reads file (fa or fq). Set to '-' to read from stdin",
                        type=str)
    parser.add_argument("-n", help="Number of reads to sample [10000]", type=int, required=False,
                        default=10000)
    parser.add_argument("-p", help="Read length percentile to use for dist "
                                   "parameter estimation [50]",
                        type=int, default=50, required=False)
    parser.add_argument("-d", help="Lower bound read length for percentile calculation [1000]",
                        type=int, default=1000,
                        required=False)
    parser.add_argument("-o", help="Output file name [tigmint-long.params.tsv]", type=str,
                        default="tigmint-long.params.tsv",
                        required=False)
    parser.add_argument("-v", "--version",
                        action="version",
                        version="tigmint_estimate_dist.py 1.2.5")

    args = parser.parse_args()

    args.READS = "/dev/stdin" if args.READS == "-" else args.READS

    read_lengths = get_n_read_lengths(args.READS, args.n, args.d)

    dist = int(np.percentile(read_lengths, args.p))

    param_file = open(args.o, 'a')
    print("read_p{p}\t{dist}".format(p=args.p, dist=dist), file=param_file)
    param_file.close()


if __name__ == "__main__":
    main()
