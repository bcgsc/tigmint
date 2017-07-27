.DELETE_ON_ERROR:
.SECONDARY:

SHELL=sh -eu -o pipefail

all: tigmint-make.cwl

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

# Run a CWL pipeline.
%.tigmint.fa: %.cwl.json
	cwl-runner $<
