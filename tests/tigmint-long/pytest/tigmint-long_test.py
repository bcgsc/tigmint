#!/usr/bin/env python3
"""Pytest tests for Tigmint with linked reads."""

import shlex
import subprocess
import os
import gzip

def run_calculate_span(reads, gsize):
    """Calculate the sequence coverage and recommended span value for long reads."""
    open_reads = subprocess.Popen(shlex.split("gunzip -c {0}.fa.gz".format(reads)), stdout=subprocess.PIPE)
    input_reads = open_reads.stdout
    calculate = subprocess.Popen(shlex.split("../../../bin/calculate-span -g {0}".format(gsize)), stdin=input_reads,
        stdout=subprocess.PIPE, universal_newlines=False)
    open_reads.stdout.close()
    span = calculate.communicate()[0].strip()
    return_code = calculate.returncode
    assert return_code == 0
    return span

def run_long_to_linked(reads, cut, minsize):
    """Cut long reads to shorter reads of size cut."""
    open_reads = subprocess.Popen(shlex.split("gunzip -c {0}.fa.gz".format(reads)), stdout=subprocess.PIPE)
    input_reads = open_reads.stdout
    cut = subprocess.Popen(shlex.split("../../../bin/long-to-linked -r - -l {0} -m {1}".format(cut, minsize)), stdin=input_reads,
        stdout=subprocess.PIPE, universal_newlines=False)
    outs = cut.communicate()[0].splitlines()
    return_code = cut.returncode
    assert return_code == 0
    return outs

def run_tigmint_molecule(draft, reads):
    """Test tigmint-molecule with a sample of linked reads and draft assembly."""
    command = "../../../bin/tigmint-make {0}.{1}.cut500.as0.65.nm500.molecule.size2000.bed draft={0} reads={1}".format(draft, reads)
    command_shlex = shlex.split(command)
    return_code = subprocess.call(command_shlex)
    assert return_code == 0
    with open("../expected_outputs/test_contig.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed".format(draft, reads), "r") as expected_molecules:
        with open("./{0}.{1}.cut500.as0.65.nm500.molecule.size2000.bed".format(draft, reads), "r") as outfile:
            output_molecules = outfile.readlines()
            for i, actual in enumerate(expected_molecules):
                assert output_molecules[i] == actual
    dir_files = os.listdir()
    bamfile = "{0}.{1}.cut500.sortbx.bam".format(draft, reads)
    if bamfile in dir_files:
        os.remove(bamfile)

def run_tigmint_cut(draft, reads):
    """Test tigmint-cut with output from run_tigmint_molecule."""
    command = "../../../bin/tigmint-make {0}.{1}.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa draft={0} reads={1} span=2".format(draft, reads)
    command_shlex = shlex.split(command)
    return_code = subprocess.call(command_shlex)
    assert return_code == 0
    with open("../expected_outputs/test_contig.test_longreads.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa.bed", "r") as expected_breakpoints:
        with open("./{0}.{1}.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa.bed".format(draft, reads), "r") as outfile:
            output_breakpoints = outfile.readlines()
            for i, actual in enumerate(expected_breakpoints):
                assert output_breakpoints[i] == actual
    with open("../expected_outputs/test_contig.test_longreads.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa", "r") as expected_breaktigs:
        with open("./{0}.{1}.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa".format(draft, reads), "r") as outfile:
            output_breaktigs = outfile.readlines()
            for i, actual in enumerate(expected_breaktigs):
                if i % 2 == 1:
                    assert output_breaktigs[i] == actual
    
    created_files = [".fa.fai", ".{0}.cut500.as0.65.nm500.molecule.size2000.bed".format(reads), 
                    ".{0}.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa".format(reads), 
                    ".{0}.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa.bed".format(reads)]
    dir_files = os.listdir()
    for suffix in created_files:
        output_file = draft + suffix
        if output_file in dir_files:
            os.remove(output_file)

def test_calculate_span_g100000():
    """Test calculate-span output with gsize = 100000."""
    output = run_calculate_span("pytest_longreads", 100000)
    assert output == b"20"

def test_calculate_span_g8072592():
    """Test calculate-span output with gsize = 8072592."""
    output = run_calculate_span("pytest_longreads", 8072592)
    assert output == b"0"

def test_long_to_linked_default():
    """Test long-to-linked output with cut length = 500, minimum molecule size = 2000."""
    cut_reads = run_long_to_linked("pytest_longreads", 500, 2000)
    assert len(cut_reads) == 32528
    cut_reads.pop()
    assert cut_reads.pop() == b">66928_ch187_17 BX:Z:489"

def test_long_to_linked_all_filtered():
    """
    Test long-to-linked output with cut length = 500, minimum molecule size = 47000.
    This should result in all reads being filtered, and no resulting cut reads.
    """
    cut_reads = run_long_to_linked("pytest_longreads", 500, 47000)
    assert len(cut_reads) == 0

def test_tigmint_molecule():
    run_tigmint_molecule("pytest_contig", "pytest_longreads")

def test_cut_reads():
    if "pytest_longreads.cut500.fa.gz" in os.listdir():
        with gzip.open("../expected_outputs/test_longreads.cut500.fa.gz", "rt") as expected_cuts:
            with gzip.open("pytest_longreads.cut500.fa.gz", "rt") as outfile:
                output_cuts = outfile.readlines()
                for i, actual in enumerate(expected_cuts):
                    assert output_cuts[i] == actual
        os.remove("pytest_longreads.cut500.fa.gz")

def test_tigmint_cut():
    run_tigmint_cut("pytest_contig", "pytest_longreads")