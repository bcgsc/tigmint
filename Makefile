.DELETE_ON_ERROR:
.SECONDARY:

all: tigmint-make.cwl

# Fetch lindenb/xml-patch-make.
xml-patch-make/stylesheets/graph2cwl.xsl:
	git clone https://github.com/lindenb/xml-patch-make

# Compile xml-patch-make
xml-patch-make/make-4.1/make-4.1/make: xml-patch-make/stylesheets/graph2cwl.xsl
	make -C xml-patch-make

# Convert a Makefile to XML.
tigmint-make.xml: %.xml: % xml-patch-make/make-4.1/make-4.1/make
	xml-patch-make/make-4.1/make-4.1/make --xml $@ -f $<

# Convert Makefile XML to CWL.
%.cwl: %.xml xml-patch-make/stylesheets/graph2cwl.xsl
	xsltproc -o $@ --stringparam shellpath $(PWD)/$*.cwl.sh xml-patch-make/stylesheets/graph2cwl.xsl $<
