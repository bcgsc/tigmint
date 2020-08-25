#!/usr/bin/env python3
"""Pytests for Tigmint."""

import pytest
import shlex
import subprocess
import os
import gzip

def long_to_linked(length=500, minsize=2000):
    """Test long-to-linked."""
    open_reads = subprocess.Popen(shlex.split("gunzip -c test_longreads.fa.gz"), stdout=subprocess.PIPE, universal_newlines=True)
    input_reads = open_reads.stdout
    long_to_linked = subprocess.Popen(shlex.split("../bin/long-to-linked -l%i -r - -m%i" % (length, minsize)), 
        stdin=input_reads, stdout=subprocess.PIPE, universal_newlines=True)
    cut_reads = long_to_linked.communicate()[0].splitlines()
    assert long_to_linked.returncode == 0
    for read in cut_reads:
        yield read

def tigmint_molecule(bamfile):
    """Test tigmint-molecule with a given alignment bam file."""
    tigmint_molecule = subprocess.Popen(shlex.split("../bin/tigmint-molecule -a0.65 -n5 -q0 -d50000 -s2000 %s" % bamfile), stdout=subprocess.PIPE)
    sorted_molecules = subprocess.Popen(shlex.split("sort -k1,1 -k2,2n -k3,3n"), stdin=tigmint_molecule.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    tigmint_molecule.wait()
    molecules = sorted_molecules.communicate()[0].splitlines()
    assert tigmint_molecule.returncode == 0
    for molecule in molecules:
        yield molecule

@pytest.fixture
def tigmint_cut():
    """Test tigmint-cut."""
    outfiles = []
    def _run(draft, reads, molecule_bed, out_fasta, span):
        outfiles.append(out_fasta)
        outfiles.append(out_fasta + ".bed")
        tigmint_cut = subprocess.call(shlex.split("../bin/tigmint-cut -p8 -w1000 -n %s -t0 -o %s -r %s -g 100000 %s %s" % (span, out_fasta, reads, draft, molecule_bed)))
        assert tigmint_cut == 0
        return out_fasta
    yield _run
    # Test teardown; remove created files
    for outfile in outfiles:
        if os.path.exists(outfile):
            os.remove(outfile)

# Tests

def test_long_to_linked_default():
    """Test long-to-linked script with  default parameters."""
    with gzip.open("test_longreads.cut500.fa.gz", "rt") as exp:
        for obs in long_to_linked():
            assert exp.readline().strip() == obs

def test_long_to_linked_all_filtered():
    """Test long-to-linked script with all molecules filtered."""
    for obs in long_to_linked(minsize=47000):
        assert obs == None

def test_tigmint_molecule_linked_default():
    """Test tigmint-molecule with linked reads and default parameters."""
    with open("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", "r") as exp:
        for obs in tigmint_molecule("test_contig.test_linkedreads.sortbx.bam"):
            assert exp.readline().strip() == obs

def test_tigmint_molecule_long_default():
    """Test tigmint-molecule with long read and default parameters."""
    with open("test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed", "r") as exp:
        for obs in tigmint_molecule("test_contig_long.test_longreads.cut500.sortbx.bam"):
            assert exp.readline().strip() == obs

def test_tigmint_cut_linked_span20(tigmint_cut):
    """Test tigmint-cut with linked reads and default parameters."""
    real_breaktigs = "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
    test_breaktigs = "pytest_" + real_breaktigs
    test_breaktigs_bed = test_breaktigs + ".bed"
    tigmint_cut("test_contig.fa", "test_linkedreads.fq.gz",
                "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", test_breaktigs, 20)
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
                for i, exp in enumerate(exp_breaktigs):
                    obs = obs_breaktigs.readline()
                    if i % 2 == 1:
                        assert exp.strip() == obs.strip()
    exp_bed = ["test\t0\t3051\ttest-1", "test\t3051\t3196\ttest-2",
                "test\t3196\t6182\ttest-3"]
    with open(test_breaktigs_bed) as obs_bed:
        for i, obs in enumerate(obs_bed):
            assert obs.strip() == exp_bed[i]

def test_tigmint_cut_long_span2(tigmint_cut):
    """Test tigmint-cut with long reads and default parameters."""
    real_breaktigs = "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa"
    test_breaktigs = "pytest_" + real_breaktigs
    test_breaktigs_bed = test_breaktigs + ".bed"
    tigmint_cut("test_contig_long.fa", "test_longreads.fq.gz",
                "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed", test_breaktigs, 2)
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
            for i, exp in enumerate(exp_breaktigs):
                obs = obs_breaktigs.readline()
                if i % 2 == 1:
                    assert exp.strip() == obs.strip()
    exp_bed = ["test\t0\t2585\ttest-1", "test\t2585\t2685\ttest-2",
                "test\t2685\t4998\ttest-3"]
    with open(test_breaktigs_bed) as obs_bed:
        for i, obs in enumerate(obs_bed):
            assert obs.strip() == exp_bed[i]

def test_tigmint_cut_long_spanauto(tigmint_cut):
    """Test tigmint-cut with long reads and automatic span calculation."""
    real_breaktigs = "test_contig_long.fa"  # No cuts should be made with gsize=100000
    test_breaktigs = "pytest_test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa"
    test_breaktigs_bed = test_breaktigs + ".bed"
    tigmint_cut("test_contig_long.fa", "test_longreads.fa.gz",
                "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed", 
                test_breaktigs, "auto")
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
            for i, exp in enumerate(exp_breaktigs):
                obs = obs_breaktigs.readline()
                if i % 2 == 1:
                    assert exp.strip() == obs.strip()
    exp_bed = "test\t0\t4998\ttest"
    with open(test_breaktigs_bed) as obs_bed:
        assert obs_bed.readline().strip() == exp_bed

