#!/usr/bin/env python3
"""Pytest tests for Tigmint with linked reads."""

import pytest
import shlex
import subprocess
import os
import gzip

@pytest.fixture(scope="session")
def setup_all(tmpdir_factory):
    path = tmpdir_factory.getbasetemp()
    draft_path = path.join("test_contig.fa")
    linked_reads_path = path.join("test_linkedreads.fq.gz")
    long_reads_path = path.join("test_longreads.fa.gz")
    cut_long_reads_path = path.join("test_longreads.cut500.fa.gz")
    copy_draft = subprocess.call(shlex.split("cp test_contig.fa %s" % draft_path))
    copy_linked_reads = subprocess.call(shlex.split("cp test_linkedreads.fq.gz %s" % linked_reads_path))
    copy_long_reads = subprocess.call(shlex.split("cp test_longreads.fa.gz %s" % long_reads_path))
    copy_cut_long_reads = subprocess.call(shlex.split("cp expected_outputs/test_longreads.cut500.fa.gz %s" % cut_long_reads_path))
    assert copy_draft == 0
    assert copy_linked_reads == 0
    assert copy_long_reads == 0
    assert copy_cut_long_reads == 0
    linked_bam = path.join("test_contig.test_linkedreads.sortbx.bam")
    long_bam = path.join("test_contig.test_longreads.cut500.sortbx.bam")
    bwa_index = subprocess.call(shlex.split("bwa index %s" % draft_path))
    assert bwa_index == 0
    samtools_index = subprocess.call(shlex.split("samtools faidx %s" % draft_path))
    assert samtools_index == 0
    linked_reads_bed_path = path.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed")
    copy_linked_reads_bed = subprocess.call(shlex.split("cp expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed %s" % linked_reads_bed_path))
    assert copy_linked_reads_bed == 0
    long_reads_bed_path = path.join("test_contig.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed")
    copy_long_reads_bed = subprocess.call(shlex.split("cp expected_outputs/test_contig.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed %s" % long_reads_bed_path))
    assert copy_long_reads_bed == 0

    # Align linked reads
    align_linked_reads = subprocess.Popen(shlex.split("bwa mem -t8 -pC %s %s" % (draft_path, linked_reads_path)), stdout=subprocess.PIPE)
    view_linked_reads = subprocess.Popen(shlex.split("samtools view -u -F4"), stdin=align_linked_reads.stdout, stdout=subprocess.PIPE)
    sort_linked_bam = subprocess.Popen(shlex.split("samtools sort -@8 -tBX -o %s" % linked_bam), stdin=view_linked_reads.stdout)
    align_linked_reads.wait()
    view_linked_reads.wait()
    sort_linked_bam.wait()
    assert align_linked_reads.returncode == 0
    assert view_linked_reads.returncode == 0
    assert sort_linked_bam.returncode == 0

    # Align long reads
    align_long_reads = subprocess.Popen(shlex.split("minimap2 -y -t8 -ax map-ont --secondary=no %s %s" % (draft_path, cut_long_reads_path)), stdout=subprocess.PIPE)
    view_long_reads = subprocess.Popen(shlex.split("samtools view -b -u -F4"), stdin=align_long_reads.stdout, stdout=subprocess.PIPE)
    sort_long_bam = subprocess.Popen(shlex.split("samtools sort -@8 -tBX -o %s" % long_bam), stdin=view_long_reads.stdout)
    align_long_reads.wait()
    view_long_reads.wait()
    sort_long_bam.wait()
    assert align_long_reads.returncode == 0
    assert view_long_reads.returncode == 0
    assert sort_long_bam.returncode == 0

def tigmint_molecule(tmpdir_factory, readtype):
    """Test tigmint-molecule for a given read type."""
    path = tmpdir_factory.getbasetemp()
    bam_path = path.join("test_contig.test_%sreads.sortbx.bam" % readtype)
    bed_path = path.join("test_contig.test_%sreads.as0.65.nm5.molecule.size2000.bed" % readtype)
    tigmint_molecule = subprocess.Popen(shlex.split("../bin/tigmint-molecule -a0.65 -n5 -q0 -d50000 -s2000 %s" % bam_path), stdout=subprocess.PIPE)
    sorted_molecules = subprocess.Popen(shlex.split("sort -k1,1 -k2,2n -k3,3n"), stdin=tigmint_molecule.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    tigmint_molecule.wait()
    sorted_molecules.wait()
    molecules = sorted_molecules.stdout
    assert tigmint_molecule.returncode == 0
    assert sorted_molecules.returncode == 0
    return molecules

def tigmint_cut(tmpdir_factory, readtype):
    """Test tigmint-cut."""
    path = tmpdir_factory.getbasetemp()
    reads = path.join("test_%sreads.fq.gz" % readtype)
    draft = path.join("test_contig.fa")
    molecules_bed = path.join("test_contig.test_%sreads.as0.65.nm5.molecule.size2000.bed" % readtype)
    breaktigs_fasta = path.join("test_contig.test_%sreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa" % readtype)
    tigmint_cut = subprocess.call(shlex.split("../bin/tigmint-cut -p8 -w1000 -n20 -t0 -o %s -r %s -g 100000 %s %s" % (breaktigs_fasta, reads, draft, molecules_bed)))
    assert tigmint_cut == 0
    return breaktigs_fasta.open()

def long_to_linked(tmpdir_factory, length=500, minsize=2000):
    path = tmpdir_factory.getbasetemp()
    long_reads_path = path.join("test_longreads.fa.gz")
    open_reads = subprocess.Popen(shlex.split("gunzip -c %s" % long_reads_path), stdout=subprocess.PIPE, universal_newlines=True)
    input_reads = open_reads.stdout
    long_to_linked = subprocess.Popen(shlex.split("../bin/long-to-linked -l%i -r - -m%i" % (length, minsize)), 
        stdin=input_reads, stdout=subprocess.PIPE, universal_newlines=True)
    cut_reads = long_to_linked.communicate()[0].splitlines()
    assert long_to_linked.returncode == 0
    return cut_reads

def test_tigmint_molecule_linked_default(tmpdir_factory, setup_all):
    with open("expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", "r") as expected:
        for observed in tigmint_molecule(tmpdir_factory, "linked"):
            assert observed == expected.readline()

def test_tigmint_molecule_long_default():
    pass

def test_long_to_linked_default(tmpdir_factory, setup_all):
    """Test long-to-linked script with the default parameters."""
    assert len(long_to_linked(tmpdir_factory)) == 32528

def test_long_to_linked_all_filtered(tmpdir_factory, setup_all):
    """Test long-to-linked script with all molecules filtered."""
    assert len(long_to_linked(tmpdir_factory, minsize=47000)) == 0

def test_tigmint_cut_linked_default(tmpdir_factory, setup_all):
    path = tmpdir_factory.getbasetemp()
    with open("expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa", 
        "r") as expected_fasta:
        for i, observed in enumerate(tigmint_cut(tmpdir_factory, "linked")):
            expected = expected_fasta.readline()
            if i % 2 == 1:
                assert observed == expected
    observed_bed = path.join("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed").open()
    with open("expected_outputs/test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed",
        "r") as expected_bed:
        for expected_breakpoint in expected_bed:
            observed_breakpoint = observed_bed.readline()
            assert expected_breakpoint == observed_breakpoint

# ***ADD TIGMINT-LONG TESTS***