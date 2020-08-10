#!/usr/bin/env python3
"""Pytest tests for Tigmint with linked reads."""

import shlex
import subprocess
import os

def run_tigmint_molecule(draft, reads):
    """Test tigmint-molecule with a sample of linked reads and draft assembly."""
    command = "../../../bin/tigmint-make {0}.{1}.as0.65.nm5.molecule.size2000.bed draft={0} reads={1}".format(draft, reads)
    command_shlex = shlex.split(command)
    return_code = subprocess.call(command_shlex)
    assert return_code == 0
    with open("../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", "r") as expected_molecules:
        with open("./{0}.{1}.as0.65.nm5.molecule.size2000.bed".format(draft, reads), "r") as outfile:
            output_molecules = outfile.readlines()
            for i, actual in enumerate(expected_molecules):
                assert output_molecules[i] == actual
    dir_files = os.listdir()
    created_files = [".fa.amb", ".fa.ann", ".fa.bwt", ".fa.pac", ".fa.sa", ".{0}.sortbx.bam".format(reads)]
    for suffix in created_files:
        output_file = draft + suffix
        if output_file in dir_files:
            os.remove(output_file)

def run_tigmint_cut(draft, reads):
    """Test tigmint-cut with output from run_tigmint_molecule."""
    command = "../../../bin/tigmint-make {0}.{1}.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa draft={0} reads={1}".format(draft, reads)
    command_shlex = shlex.split(command)
    return_code = subprocess.call(command_shlex)
    assert return_code == 0
    with open("../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed", "r") as expected_breakpoints:
        with open("./{0}.{1}.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed".format(draft, reads), "r") as outfile:
            output_breakpoints = outfile.readlines()
            for i, actual in enumerate(expected_breakpoints):
                assert output_breakpoints[i] == actual
    with open("../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa", "r") as expected_breaktigs:
        with open("./{0}.{1}.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa".format(draft, reads), "r") as outfile:
            output_breaktigs = outfile.readlines()
            for i, actual in enumerate(expected_breaktigs):
                if i % 2 == 1:
                    assert output_breaktigs[i] == actual
    
    created_files = [".fa.fai", ".{0}.as0.65.nm5.molecule.size2000.bed".format(reads), 
                    ".{0}.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa".format(reads), 
                    ".{0}.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed".format(reads)]
    dir_files = os.listdir()
    for suffix in created_files:
        output_file = draft + suffix
        if output_file in dir_files:
            os.remove(output_file)
            
def test_tigmint_molecule():
    run_tigmint_molecule("pytest_contig", "pytest_linkedreads")

def test_tigmint_cut():
    run_tigmint_cut("pytest_contig", "pytest_linkedreads")