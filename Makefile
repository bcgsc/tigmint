.DELETE_ON_ERROR:
.SECONDARY:

SHELL=bash -eu -o pipefail

all: tigmint-make.cwl tigmint-make.gv.svg

check: mt.tigmint.fa

# GraphViz

# Create phony input files.
draft.fa reads.fq.gz:
	touch $@

# Convert the pipepline graph to GraphViz using makefile2graph.
tigmint-make.gv: tigmint-make draft.fa reads.fq.gz
	makefile2graph -f $< \
		| gsed -r \
			-e 's/label="(all|arcs|tigmint)".*]/label="\1", shape=ellipse, style=filled]/' \
			-e 's/label="(draft.tigmint.fa|draft.tigmint.arcs.fa)".*]/label="\1", shape=parallelogram, style=filled]/' \
			-e 's/color="green"/shape=parallelogram, style=filled/;s/color="red"/shape=rectangle/' \
		| tred >$@

# Render a GraphViz file to PNG.
%.gv.png: %.gv
	dot -Tpng -Gsize=16,16 -o $@ $<

# Render a GraphViz file to SVG.
%.gv.svg: %.gv
	dot -Tsvg -o $@ $<

# Common Workflow Language (CWL)

# Fetch lindenb/xml-patch-make.
xml-patch-make/stylesheets/graph2cwl.xsl:
	git clone https://github.com/lindenb/xml-patch-make

# Compile xml-patch-make
xml-patch-make/make-4.1/make-4.1/make: xml-patch-make/stylesheets/graph2cwl.xsl
	make -C xml-patch-make

# Generate a generic XML pipeline.
tigmint-make.xml: tigmint-make xml-patch-make/make-4.1/make-4.1/make
	xml-patch-make/make-4.1/make-4.1/make --xml $@ -f $<

# Generate a XML pipeline.
%.xml: %.fa %.lrsim.fq.gz xml-patch-make/make-4.1/make-4.1/make
	xml-patch-make/make-4.1/make-4.1/make --xml $@ -f tigmint-make tigmint \
		draft=$* reads=$*.lrsim ref=$* G=16569

# Convert Makefile XML to CWL.
%.cwl: %.xml xml-patch-make/stylesheets/graph2cwl.xsl
	xsltproc -o $@ --stringparam shellpath $(PWD)/$*.cwl.sh xml-patch-make/stylesheets/graph2cwl.xsl $<
	chmod +x $@.sh

# Create a CWL JSON driver script.
%.cwl.json: %.cwl
	printf '#!/usr/bin/env cwl-runner\n{ "cwl:tool": "$<#main" }\n' >$@

# Generate test data.

# Download the human mitochondrial genome.
mt.fa:
	curl ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz \
		| seqtk seq >$@

# Simulate linked reads using LRSIM.
%.lrsim_S1_L001_R1_001.fastq.gz %.lrsim_S1_L001_R2_001.fastq.gz: %.fa
	simulateLinkedReads -g $< -p $*.lrsim -o -x0.005 -f4 -t1 -m1 -z1

# Convert paired FASTQ to the interleaved FASTQ format of longranger basic.
%.fq.gz: %_S1_L001_R1_001.fastq.gz %_S1_L001_R2_001.fastq.gz
	seqtk mergepe $^ \
		| paste - - - - - - - - \
		| awk -v OFS='\n' '{ print $$1 " BX:Z:" substr($$3, 1, 15) "-1", substr($$3, 16), "+", substr($$5, 16), $$6 " BX:Z:" substr($$3, 1, 15) "-1", $$8, $$9, $$10 }' \
		| gzip >$@

# Run Tigmint on the mitochondrial test data.
mt.tigmint.fa: %.tigmint.fa: %.fa %.lrsim.fq.gz
	./tigmint-make tigmint draft=$* reads=$*.lrsim depth_threshold=150 starts_threshold=2 ref=$* G=16569
