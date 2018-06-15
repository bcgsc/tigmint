<img src="http://sjackman.ca/img/tigmint.png" style="width:4in">

# Correct misassemblies using linked reads

Cut sequences at positions with few spanning molecules.

Written by [Shaun Jackman](http://sjackman.ca), Lauren Coombe, and Justin Chu.

[bioRxiv doi:10.1101/304253](https://www.biorxiv.org/content/early/2018/04/20/304253) &middot; [Slides](http://sjackman.ca/tigmint-recomb-slides) &middot; [Poster](https://f1000research.com/posters/6-1406)

# Description

Tigmint identifies and corrects misassemblies using linked reads from 10x Genomics Chromium. The reads are first aligned to the assembly, and the extents of the large DNA molecules are inferred from the alignments of the reads. The physical coverage of the large molecules is more consistent and less prone to coverage dropouts than that of the short read sequencing data. The sequences are cut at positions that have insufficient spanning molecules. Tigmint outputs a BED file of these cut points, and a FASTA file of the cut sequences.

Each window of a specified fixed size is checked for a minimum number of spanning molecules. Sequences are cut at those positions where a window with sufficient coverage is followed by some number of windows with insufficient coverage is then followed again by a window with sufficient coverage.

# Installation

## Install Tigmint using Brew

Install [Linuxbrew](http://linuxbrew.sh/) on Linux or Windows Subsystem for Linux (WSL), or 
install [Homebrew](https://brew.sh/) on macOS, and then run the command

	brew install tigmint

## Install Tigmint from the source code

Download and extract the source code. Compiling is not needed.

```
git clone https://github.com/bcgsc/tigmint && cd tigmint
```
or
```
curl -L https://github.com/bcgsc/tigmint/archive/master.tar.gz | tar xz && mv tigmint-master tigmint && cd tigmint
```

# Dependencies

## Install Python package dependencies
```sh
pip3 install intervaltree pybedtools pysam statistics
```

Tigmint uses Bedtools, BWA and Samtools. These dependencies may be installed using [Homebrew](https://brew.sh) on macOS or [Linuxbrew](http://linuxbrew.sh) on Linux.

## Install the dependencies of Tigmint
```sh
brew install bedtools bwa samtools
```

## Install the dependencies of ARCS (optional)
```sh
brew tap brewsci/bio
brew install arcs links-scaffolder
```

## Install the dependencies for calculating assembly metrics (optional)
```sh
brew install abyss seqtk
```

# Usage

To run Tigmint on the draft assembly `draft.fa` with the reads `reads.fq.gz`, which have been run through `longranger basic`:

```sh
samtools faidx draft.fa
bwa index draft.fa
bwa mem -t8 -p -C draft.fa reads.fq.gz | samtools sort -@8 -tBX -o draft.reads.sortbx.bam
tigmint-molecule draft.reads.sortbx.bam | sort -k1,1 -k2,2n -k3,3n >draft.reads.molecule.bed
tigmint-cut -p8 -o draft.tigmint.fa draft.fa draft.reads.molecule.bed
```

- `bwa mem -C` is used to copy the BX tag from the FASTQ header to the SAM tags.
- `samtools sort -tBX` is used to sort first by barcode and then position.

Alternatively, you can run the Tigmint pipeline using the Makefile driver script `tigmint-make`. To run Tigmint on the draft assembly `myassembly.fa` with the reads `myreads.fq.gz`, which have been run through `longranger basic`:

```sh
tigmint-make tigmint draft=myassembly reads=myreads
```

To run both Tigmint and scaffold the corrected assembly with [ARCS](https://github.com/bcgsc/arcs):

```sh
tigmint-make arcs draft=myassembly reads=myreads
```

To run Tigmint, ARCS, and calculate assembly metrics using the reference genome `GRCh38.fa`:

```sh
tigmint-make metrics draft=myassembly reads=myreads ref=GRCh38 G=3088269832
```

# Note

+ `tigmint-make` is a Makefile script, and so any `make` options may also be used with `tigmint-make`, such as `-n` (`--dry-run`).
+ The file extension of the assembly must be `.fa` and the reads `.fq.gz`, and the extension is not included in the parameters `draft` and `reads`. These specific file name requirements result from implementing the pipeline in GNU Make.

# tigmint-make commands

+ `tigmint`: Run Tigmint, and produce a file named `$draft.tigmint.fa`
+ `arcs`: Run Tigmint and ARCS, and produce a file name `$draft.tigmint.arcs.fa`
+ `metrics`: Run, Tigmint, ARCS, and calculate assembly metrics using `abyss-fac` and `abyss-samtobreak`, and produce TSV files.

# Parameters of Tigmint

+ `draft`: Name of the draft assembly, `draft.fa`
+ `reads`: Name of the reads, `reads.fq.gz`
+ `span=20`: Number of spanning molecules threshold
+ `window=1000`: Window size (bp) for checking spanning molecules
+ `minsize=2000`: Minimum molecule size
+ `as=0.65`: Minimum AS/read length ratio
+ `nm=5`: Maximum number of mismatches
+ `dist=50000`: Maximum distance (bp) between reads to be considered the same molecule
+ `mapq=0`: Mapping quality threshold
+ `trim=0`: Number of bases to trim off contigs following cuts
+ `t=8`: Number of threads

# Parameters of ARCS
+ `c=5`
+ `e=30000`
+ `r=0.05`

# Parameters of LINKS
+ `a=0.1`
+ `l=10`

# Parameters for calculating assembly metrics

+ `ref`: Reference genome, `ref.fa`, for calculating assembly contiguity metrics
+ `G`: Size of the reference genome, for calculating NG50 and NGA50

# Tips

- If your barcoded reads are in multiple FASTQ files, the initial alignments of the barcoded reads to the draft assembly can be done in parallel and merged prior to running Tigmint.
- When aligning with BWA-MEM, use the `-C` option to include the barcode in the BX tag of the alignments.
- Sort by BX tag using `samtools sort -tBX`.
- Merge multiple BAM files using `samtools merge -tBX`.

# Support

After first looking for existing issue at <https://github.com/bcgsc/tigmint/issues>, please report a new issue at <https://github.com/bcgsc/tigmint/issues/new>. Please report the names of your input files, the exact command line that you are using, and the entire output of Tigmint.

# Pipeline

[![Tigmint pipeline illustration](pipeline.gv.png)](pipeline.gv.svg)
