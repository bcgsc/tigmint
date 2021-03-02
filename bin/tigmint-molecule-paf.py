#!/usr/bin/env python3
"""
Group linked reads from paf format into molecules.
@author: Lauren Coombe
Based on tigmint-molecule written by Justin Chu, Janet Li and Shaun Jackman
"""
import argparse
import re
import sys
from tigmint_molecule import FileFormat, Molecule, MolecIdentifier


class ReadMapping:
    "Class representing a linked read mapping"

    def __init__(self, rname, start, end, barcode):
        self.rname = rname
        self.start = start
        self.end = end
        self.barcode = barcode

    def __str__(self):
        return "{0}: {1}-{2} BX:{3}".format(self.rname, self.start, self.end, self.barcode)

    def __hash__(self):
        return hash((self.rname, self.start, self.end, self.barcode))

    def __eq__(self, other):
        return str(self) == str(other)


class MolecIdentifierPaf:
    """Group molecules into barcodes"""

    def print_current_molecule(self, ref, start, end, barcode, num_reads, new_molec_file):
        if num_reads >= self.opt.min_reads and end - start >= self.opt.min_size:
            print(ref, start, end, barcode, num_reads, sep="\t", file=new_molec_file)


    def print_new_molecule(self, barcode, intervals, new_molec_file):
        "Given a collection of intervals from a barcode, find the molecule extents and print them"
        barcode_re = re.compile(r'^BX:Z:(\S+)')

        barcode_match = re.search(barcode_re, barcode)
        assert barcode_match
        barcode_short = barcode_match.group(1)

        for ref in intervals:
            current_start = None
            current_end = None
            num_reads = 0
            for interval in sorted(intervals[ref], key=lambda x: (x.start, x.end)):
                if current_start is None:
                    current_start = interval.start
                    current_end = interval.end
                elif interval.start - current_end <= self.opt.max_dist:
                    current_end = interval.end
                else:
                    # Start new molecule
                    self.print_current_molecule(ref, current_start, current_end, barcode_short,
                                                num_reads, new_molec_file)
                    current_start = interval.start
                    current_end = interval.end
                    num_reads = 0

                num_reads += 1
            self.print_current_molecule(ref, current_start, current_end, barcode_short,
                                        num_reads, new_molec_file)

    def run(self):
        "Run molecule identification"

        if self.opt.out_molecules_filename:
            out_molecules_file = open(self.opt.out_molecules_filename, 'w')
        else:
            out_molecules_file = sys.stdout

        if self.opt.output_format == FileFormat.TSV:
            print("Rname\tStart\tEnd\tSize\tBX\tMI\tReads",
                  file=out_molecules_file)

        prev_barcode = None
        cur_intervals = {}
        #new_molec_id = 0

        with open(self.opt.PAF, 'r') as in_paf_file:
            for paf_entry in in_paf_file:
                paf_entry = paf_entry.strip().split("\t")
                rname, start, end, barcode = paf_entry[5], int(paf_entry[7]), int(paf_entry[8]), paf_entry[18]
                if prev_barcode != barcode:
                    if prev_barcode is not None:
                        self.print_new_molecule(prev_barcode, cur_intervals, out_molecules_file)
                        cur_intervals = {}
                    prev_barcode = barcode

                if rname not in cur_intervals:
                    cur_intervals[rname] = set()
                cur_intervals[rname].add(ReadMapping(rname, start, end, barcode))
            self.print_new_molecule(prev_barcode, cur_intervals, out_molecules_file)

    def parse_arguments(self):
        "Parse input arguments"
        parser = argparse.ArgumentParser(
            description="Group linked reads simulated from long reads into molecules. "
                        "Read a PAF file and output a BED file. ")
        parser.add_argument(
            '--version', action='version', version='tigmint-molecule-paf 1.2.2')
        parser.add_argument(
            metavar="PAF", dest="PAF",
            help="Input PAF file, - for stdin")
        parser.add_argument(
            "-o", "--output", dest="out_molecules_filename",
            help="Output molecule file [stdout]",
            metavar="FILE")
        parser.add_argument(
            "--bed", action="store_const", dest="output_format", const=FileFormat.BED,
            default=FileFormat.BED,
            help="Output in BED format [default]")
        parser.add_argument(
            "--tsv", action="store_const", dest="output_format", const=FileFormat.TSV,
            help="Output in TSV format")
        parser.add_argument(
            "-d", "--dist", dest="max_dist", type=int, default=50000,
            help="Maximum distance between reads in the same molecule [50000]",
            metavar="N")
        parser.add_argument(
            "-m", "--reads", dest="min_reads", type=int, default=4,
            help="Minimum number of reads per molecule (duplicates are filtered out) [4]",
            metavar="N")
        parser.add_argument(
            "-s", "--size", dest="min_size", type=int, default=2000,
            help="Minimum molecule size [2000]",
            metavar="N")
        parser.add_argument(
            "-p", "--params", dest="param_file", type=str, default=None)

        self.opt = parser.parse_args()

        # Use calculated max dist if parameter file is provided
        if self.opt.param_file:
            self.opt.max_dist = MolecIdentifier.get_dist(self)

        self.opt.PAF = "/dev/stdin" if self.opt.PAF == "-" else self.opt.PAF


    def __init__(self):
        self.parse_arguments()

def main():
    "Create a MolecIdentifier instance"
    MolecIdentifierPaf().run()

if __name__ == "__main__":
    main()