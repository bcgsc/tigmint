#!/usr/bin/env python3
"""Pytests for Tigmint."""

import pytest
import shlex
import subprocess
import os
import gzip

def long_to_linked(length=500, minsize=2000, span="auto", G=100000, dist="default"):
    """Run long-to-linked."""
    open_reads = subprocess.Popen(shlex.split("gunzip -c test_longreads.fa.gz"), stdout=subprocess.PIPE, universal_newlines=True)
    input_reads = open_reads.stdout
    output_param_file = "tigmint-long.span_G_%s.span_%s.dist_%s.tsv" % (G, str(span), dist)
    params = []
    if span == "auto":
        params.append('-s')
    if dist == "auto":
        params.append('-d')
    long_to_linked = subprocess.Popen(shlex.split("../bin/long-to-linked -l%i -m%i -g %i -o %s %s" % (length, minsize, G, output_param_file, " ".join(params))), 
            stdin=input_reads, stdout=subprocess.PIPE, universal_newlines=True)

    cut_reads = long_to_linked.communicate()[0].splitlines()
    assert long_to_linked.returncode == 0
    for read in cut_reads:
        yield read

def tigmint_molecule(bamfile, dist="default"):
    """Run tigmint_molecule.py with a given alignment bam file."""
    if dist == "default":
        tigmint_molecule = subprocess.Popen(shlex.split("../bin/tigmint_molecule.py -a0.65 -n5 -q0 -d50000 -s2000 %s" % bamfile), stdout=subprocess.PIPE)
    elif dist == "auto":
        tigmint_molecule = subprocess.Popen(shlex.split("../bin/tigmint_molecule.py -a0.65 -n5 -q0 -d50000 -s2000 -p tigmint-long.span_G_1000000.span_20.dist_auto.tsv %s" % bamfile), stdout=subprocess.PIPE)
    sorted_molecules = subprocess.Popen(shlex.split("sort -k1,1 -k2,2n -k3,3n"), stdin=tigmint_molecule.stdout, stdout=subprocess.PIPE, universal_newlines=True)
    tigmint_molecule.wait()
    molecules = sorted_molecules.communicate()[0].splitlines()
    assert tigmint_molecule.returncode == 0
    for molecule in molecules:
        yield molecule

@pytest.fixture
def tigmint_cut():
    """Run tigmint-cut."""
    outfiles = []
    def _run(draft, reads, molecule_bed, out_fasta, span=2, auto_span=False, G=100000):
        outfiles.append(out_fasta)
        outfiles.append(out_fasta + ".bed")
        if auto_span:
            param_file = "tigmint-long.span_G_%s.span_auto.dist_default.tsv" % G
            tigmint_cut = subprocess.call(shlex.split("../bin/tigmint-cut -p8 -w1000 -n %s -t0 -o %s -f %s \
                %s %s" % (span, out_fasta, param_file, draft, molecule_bed)))
        else:
            tigmint_cut = subprocess.call(shlex.split("../bin/tigmint-cut -p8 -w1000 -n %s -t0 -o %s \
                %s %s" % (span, out_fasta, draft, molecule_bed)))
        assert tigmint_cut == 0
        return out_fasta
    yield _run
    # Test teardown; remove created files
    for outfile in outfiles:
        if os.path.exists(outfile):
            os.remove(outfile)

@pytest.fixture()
def tigmint_pipeline():
    """Run entire Tigmint pipelines."""
    os.chdir("test_installation/")
    run_pipelines = subprocess.call(shlex.split("./run_tigmint_demo.sh"))
    yield run_pipelines
    outfiles = ["test_contig.fa.amb", "test_contig.fa.ann", "test_contig.fa.bwt",
                "test_contig.fa.fai", "test_contig.fa.pac", "test_contig.fa.sa",
                "test_contig_long.cut500.tigmint.fa", "test_contig.tigmint.fa", "test_contig_long.fa.fai",
                "test_contig_long.test_longreads.cut500.molecule.size2000.bed",
                "test_contig_long.test_longreads.cut500.molecule.size2000.trim0.window1000.span10.breaktigs.fa",
                "test_contig_long.test_longreads.cut500.molecule.size2000.trim0.window1000.span10.breaktigs.fa.bed",
                "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed",
                "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa",
                "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed",
                "test_contig.test_linkedreads.sortbx.bam",
                "test_longreads.cut500.fq.gz", "test_longreads.tigmint-long.span.txt",
                "test_longreads.tigmint-long.params.tsv"]
    for outfile in outfiles:
        if os.path.exists(outfile):
            os.remove(outfile)
    

# Tests

def test_long_to_linked_default():
    """Test long-to-linked script with default parameters."""
    with gzip.open("test_longreads.cut500.fq.gz", "rt") as exp:
        for obs in long_to_linked():
            assert exp.readline().strip() == obs
    with open("tigmint-long.span_G_100000.span_auto.dist_default.tsv") as span:
        assert span.readline().strip() == "span\t20"

def test_long_to_linked_all_filtered():
    """Test long-to-linked script with all molecules filtered."""
    for obs in long_to_linked(minsize=47000):
        assert obs == None

def test_long_to_linked_large_genome():
    """Test long-to-linked script with a large genome size."""
    with gzip.open("test_longreads.cut500.fq.gz", "rt") as exp:
        for obs in long_to_linked(G=1000000):
            assert exp.readline().strip() == obs
    with open("tigmint-long.span_G_1000000.span_auto.dist_default.tsv") as span:
        assert span.readline().strip() == "span\t2"

def test_long_to_linked_no_auto_params():
    """Test long-to-linked script with a given span."""
    with gzip.open("test_longreads.cut500.fq.gz", "rt") as exp:
        for obs in long_to_linked(span=20, G=1000):
            assert exp.readline().strip() == obs
    assert not os.access("tigmint-long.span_G_1000.span_20.dist_default.tsv", os.F_OK)

def test_long_to_linked_auto_dist():
    """Test long-to-linked script with dist=auto."""
    with gzip.open("test_longreads.cut500.fq.gz", "rt") as exp:
        for obs in long_to_linked(span=20, G=1000000, dist="auto"):
            assert exp.readline().strip() == obs
    with open("tigmint-long.span_G_1000000.span_20.dist_auto.tsv", "rt") as param_file:
        assert param_file.readline().strip() == "read_p50\t17154"

def test_tigmint_molecule_linked_default():
    """Test tigmint_molecule.py with linked reads and default parameters."""
    with open("test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", "r") as exp:
        for obs in tigmint_molecule("test_contig.test_linkedreads.sortbx.bam"):
            assert exp.readline().strip() == obs

def test_tigmint_molecule_long_default():
    """Test tigmint_molecule.py with long read and default parameters."""
    with open("test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed", "r") as exp:
        for obs in tigmint_molecule("test_contig_long.test_longreads.cut500.sortbx.bam"):
            assert exp.readline().strip() == obs

def test_tigmint_molecule_long_dist_auto():
    """Test tigmint_molecule.py with dist=read length p5 from calculated file."""
    with open("test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.dist_auto.bed", "r") as exp:
        for obs in tigmint_molecule("test_contig_long.test_longreads.cut500.sortbx.bam", dist="auto"):
            assert exp.readline().strip() == obs

def test_tigmint_cut_linked_span20(tigmint_cut):
    """Test tigmint-cut with linked reads and default parameters."""
    real_breaktigs = "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
    test_breaktigs = "pytest_" + real_breaktigs
    test_breaktigs_bed = test_breaktigs + ".bed"
    tigmint_cut("test_contig.fa", "test_linkedreads.fq.gz",
                "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed", test_breaktigs, span=20)
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
                for i, exp in enumerate(exp_breaktigs):
                    obs = obs_breaktigs.readline()
                    if i % 2 == 1:
                        assert exp == obs
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
                "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed", test_breaktigs)
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
            for i, exp in enumerate(exp_breaktigs):
                obs = obs_breaktigs.readline()
                if i % 2 == 1:
                    assert exp == obs
    exp_bed = ["test\t0\t2585\ttest-1", "test\t2585\t2685\ttest-2",
                "test\t2685\t4998\ttest-3"]
    with open(test_breaktigs_bed) as obs_bed:
        for i, obs in enumerate(obs_bed):
            assert obs.strip() == exp_bed[i]

def test_tigmint_cut_long_spanauto(tigmint_cut):
    """Test tigmint-cut with long reads, automatic span calculation and G=100000."""
    # No cuts should be made with G=100000
    real_breaktigs = "test_contig_long.fa"
    test_breaktigs = "pytest_test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa"
    test_breaktigs_bed = test_breaktigs + ".bed"
    tigmint_cut("test_contig_long.fa", "test_longreads.fa.gz",
                "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed",
                test_breaktigs, auto_span=True)
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
            for i, exp in enumerate(exp_breaktigs):
                obs = obs_breaktigs.readline()
                if i % 2 == 1:
                    assert exp == obs
    exp_bed = "test\t0\t4998\ttest"
    with open(test_breaktigs_bed) as obs_bed:
        assert obs_bed.readline().strip() == exp_bed

def test_tigmint_cut_long_spanauto_largeG(tigmint_cut):
    """Test tigmint-cut with long reads, automatic span calculation and G=1000000."""
    # Auto calculation should result in span=2
    real_breaktigs = "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.trim0.window1000.span2.breaktigs.fa"
    test_breaktigs = "pytest_" + real_breaktigs
    test_breaktigs_bed = test_breaktigs + ".bed"
    tigmint_cut("test_contig_long.fa", "test_longreads.fa.gz",
                "test_contig_long.test_longreads.cut500.as0.65.nm500.molecule.size2000.bed",
                test_breaktigs, auto_span=True, G=1000000)
    with open(test_breaktigs) as obs_breaktigs:
        with open(real_breaktigs) as exp_breaktigs:
            for i, exp in enumerate(exp_breaktigs):
                obs = obs_breaktigs.readline()
                if i % 2 == 1:
                    assert exp == obs
    exp_bed = ["test\t0\t2585\ttest-1", "test\t2585\t2685\ttest-2",
                "test\t2685\t4998\ttest-3"]
    with open(test_breaktigs_bed) as obs_bed:
        for i, obs in enumerate(obs_bed):
            assert obs.strip() == exp_bed[i]

def test_pipeline(tigmint_pipeline):
    """Test entire pipeline with long and linked reads."""
    assert tigmint_pipeline == 0
    
    # Compare tigmint outputs
    # Alignments
    bam = "test_contig.test_linkedreads.sortbx.bam"
    exp = subprocess.Popen(shlex.split("samtools view %s" % ("expected_outputs/" + bam)),
                            stdout=subprocess.PIPE, universal_newlines=True)
    obs = subprocess.Popen(shlex.split("samtools view %s" % bam), stdout=subprocess.PIPE,
                            universal_newlines=True)
    exp_bam = exp.communicate()[0].splitlines()
    obs_bam = obs.communicate()[0].splitlines()
    for i, exp_alignment in enumerate(exp_bam):
        assert exp_alignment == obs_bam[i]
    
    # Other readable files
    tigmint_outputs = ["test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.bed",
        "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa",
        "test_contig.test_linkedreads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed"]
    for output in tigmint_outputs:
        expected_output = "expected_outputs/" + output
        with open(expected_output) as exp:
            with open(output) as obs:
                obs_content = obs.readlines()
                if output.endswith(".fa"):
                    for i, exp_line in enumerate(exp):
                        if i % 2 == 1:
                            assert exp_line == obs_content[i]
                else:
                    for i, exp_line in enumerate(exp):
                        assert exp_line == obs_content[i]
    
    # Compare tigmint-long outputs
    # Cut reads
    cut_reads = "test_longreads.cut500.fq.gz"
    exp = subprocess.Popen(shlex.split("gunzip -c %s" % ("expected_outputs/" + cut_reads)),
                            stdout=subprocess.PIPE, universal_newlines=True)
    obs = subprocess.Popen(shlex.split("gunzip -c %s" % cut_reads), stdout=subprocess.PIPE,
                            universal_newlines=True)
    exp_cut_reads = exp.communicate()[0].splitlines()
    obs_cut_reads = obs.communicate()[0].splitlines()
    for i, exp_read in enumerate(exp_cut_reads):
        assert exp_read == obs_cut_reads[i]

    # Other output files
    tigmint_long_outputs = ["test_longreads.tigmint-long.params.tsv",
        "test_contig_long.test_longreads.cut500.molecule.size2000.bed",
        "test_contig_long.test_longreads.cut500.molecule.size2000.trim0.window1000.span10.breaktigs.fa",
        "test_contig_long.test_longreads.cut500.molecule.size2000.trim0.window1000.span10.breaktigs.fa.bed"]
    for output in tigmint_long_outputs:
        expected_output = "expected_outputs/" + output
        with open(expected_output) as exp:
            with open(output) as obs:
                obs_content = obs.readlines()
                if output.endswith(".fa"):
                    for i, exp_line in enumerate(exp):
                        if i % 2 == 1:
                            assert exp_line == obs_content[i]
                else:
                    for i, exp_line in enumerate(exp):
                        assert exp_line == obs_content[i]
