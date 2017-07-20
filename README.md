# Correct misassemblies using Tigmint

Split sequences at positions with low depth of coverage and high number of molecule starts.

Written by [Shaun Jackman](http://sjackman.ca)

# Usage

To run Tigmint on the draft assembly `assembly.fa` with the reads `reads.fq.gz`, which have been run through `longranger basic`:

```sh
tigmint-make draft=myassembly reads=myreads
```

To run Tigmint and calculate assembly metrics using the reference genome `GRCh38.fa`:

```sh
tigmint-make draft=myassembly reads=myreads ref=GRCh38 G=3088269832
```

Note: `tigmint-make` is a Makefile script, and so any `make` options may also be used with `tigmint-make`, such as `-n` (`--dry-run`).

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
