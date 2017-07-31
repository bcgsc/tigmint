# Correct misassemblies using Linked Reads

Split sequences at positions with low depth of coverage and high number of molecule starts.

Written by [Shaun Jackman](http://sjackman.ca)

# Usage

To run Tigmint on the draft assembly `assembly.fa` with the reads `reads.fq.gz`, which have been run through `longranger basic`:

```sh
tigmint-make tigmint draft=myassembly reads=myreads
```

To run Tigmint and ARCS on the draft assembly `assembly.fa` with the reads `reads.fq.gz`, which have been run through `longranger basic`:

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

# Commands

+ `tigmint`: Run Tigmint, and produce a file named `$draft.tigmint.fa`
+ `arcs`: Run Tigmint and ARCS, and produce a file name `$draft.tigmint.arcs.fa`
+ `metrics`: Run, Tigmint, ARCS, and calculate assembly metrics using `abyss-fac` and `abyss-samtobreak`, and produce TSV files.

# Parameters of Tigmint

+ `draft`: Name of the draft assembly, `draft.fa`
+ `reads`: Name of the reads, `reads.fq.gz`
+ `depth_threshold=100`: Depth of coverage threshold
+ `starts_threshold=4`: Number of molecule starts threshold
+ `minsize=2000`: Minimum molecule size
+ `as=100`: Minimum alignment score
+ `nm=5`: Maximum number of mismatches
+ `t=8`: Number of threads
+ `gzip=gzip`: gzip compression program, use `pigz -p8` for parallelized compression

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

# Support

After first looking for existing issue at <https://github.com/bcgsc/tigmint/issues>, please report a new issue at <https://github.com/bcgsc/tigmint/issues/new>. Please report the names of your input files, the exact command line that you are using, and the entire output of `tigmint-make`.
