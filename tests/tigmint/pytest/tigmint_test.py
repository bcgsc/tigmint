#!/usr/bin/env python3
"""Pytest tests for Tigmint with linked reads."""

import pytest
import shlex
import subprocess
import os

@pytest.fixture
def setup_molecule(tmpdir):
    """Setup code to create required alignment and index files for tigmint-molecule."""
    draft = tmpdir.join("test_contig.fa")
    reads = tmpdir.join("test_linkedreads.fq.gz")
    bam = tmpdir.join("test_contig.test_linkedreads.sortbx.bam")
    copy_draft = subprocess.call(shlex.split("cp ../test_contig.fa %s" % draft))
    copy_reads = subprocess.call(shlex.split("cp ../test_linkedreads.fq.gz %s" % reads))
    assert copy_draft == 0 and copy_reads == 0
    bwa_index = subprocess.call(shlex.split("bwa index %s" % draft))
    assert bwa_index == 0
    align_reads = subprocess.Popen(shlex.split("bwa mem -t8 -pC %s %s" % (draft, reads)), stdout=subprocess.PIPE)
    view_reads = subprocess.Popen(shlex.split("samtools view -u -F4"), stdin=align_reads.stdout, stdout=subprocess.PIPE)
    sort_bam = subprocess.Popen(shlex.split("samtools sort -@8 -tBX -o %s" % bam), stdin=view_reads.stdout)
    align_reads.wait()
    view_reads.wait()
    sort_bam.wait()
    assert align_reads.returncode == 0
    assert view_reads.returncode == 0
    assert sort_bam.returncode == 0

@pytest.fixture
def tigmint_molecule(tmpdir):
    """Test tigmint-molecule."""
    bam = tmpdir.join("test_contig.test_linkedreads.sortbx.bam")
    bed = tmpdir.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed")
    tigmint_molecule = subprocess.Popen(shlex.split("../../../bin/tigmint-molecule -a0.65 -n5 -q0 -d50000 -s2000 %s" % bam), stdout=subprocess.PIPE)
    sorted_molecules = subprocess.Popen(shlex.split("sort -k1,1 -k2,2n -k3,3n"), stdin=tigmint_molecule.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    tigmint_molecule.wait()
    sorted_molecules.wait()
    molecules = sorted_molecules.stdout
    assert tigmint_molecule.returncode == 0
    assert sorted_molecules.returncode == 0
    yield molecules

@pytest.fixture
def setup_cut(tmpdir):
    """Setup code to create required files for tigmint-cut."""
    draft = tmpdir.join("test_contig.fa")
    reads = tmpdir.join("test_linkedreads.fq.gz")
    copy_draft = subprocess.call(shlex.split("cp ../test_contig.fa %s" % draft))
    copy_reads = subprocess.call(shlex.split("cp ../test_linkedreads.fq.gz %s" % reads))
    assert copy_draft == 0 and copy_reads == 0
    molecules_bed = tmpdir.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed")
    copy_bed = subprocess.call(shlex.split("cp ../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed %s" % molecules_bed))
    assert copy_bed == 0
    samtools_index = subprocess.call(shlex.split("samtools faidx %s" % draft))
    assert samtools_index == 0

@pytest.fixture
def tigmint_cut_fasta(tmpdir):
    """Test tigmint-cut."""
    reads = tmpdir.join("test_linkedreads.fq.gz")
    draft = tmpdir.join("test_contig.fa")
    molecules_bed = tmpdir.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed")
    breaktigs_fasta = tmpdir.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa")
    tigmint_cut = subprocess.call(shlex.split("../../../bin/tigmint-cut -p8 -w1000 -n20 -t0 -o %s -r %s -g 100000 %s %s" % (breaktigs_fasta, reads, draft, molecules_bed)))
    assert tigmint_cut == 0
    yield breaktigs_fasta.open()
            
def test_tigmint_molecule_default(tmpdir, setup_molecule, tigmint_molecule):
    with open("../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", "r") as expected:
        for observed in tigmint_molecule:
            assert expected.readline() == observed

def test_tigmint_cut_default(tmpdir, setup_cut, tigmint_cut_fasta):
    with open("../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa", 
        "r") as expected_fasta:
        for i, observed in enumerate(tigmint_cut_fasta):
            expected = expected_fasta.readline()
            if i % 2 == 1:
                assert observed == expected
    observed_bed = tmpdir.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed").open()
    with open("../expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed",
        "r") as expected_bed:
        for expected_breakpoint in expected_bed:
            observed_breakpoint = observed_bed.readline()
            assert expected_breakpoint == observed_breakpoint

    