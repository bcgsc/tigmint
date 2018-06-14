.DELETE_ON_ERROR:
.SECONDARY:

# Number of threads.
t=16

SHELL=bash -eu -o pipefail

all: tigmint-make.cwl tigmint-make.gv.svg

check: mt.tigmint.fa
	diff -b mt.tigmint.fa.wc <(wc $<)

clean:
	rm -f \
		mt.mt.halved.lrsim.as0.65.nm5.molecule.size2000.bed \
		mt.mt.halved.lrsim.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa \
		mt.mt.halved.lrsim.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bed \
		mt.tigmint.fa

# Install Tigmint in this directory.
prefix=/usr/local
bindir=$(prefix)/bin

# Install Tigmint.
install:
	install -d $(DESTDIR)$(bindir)
	install bin/* $(DESTDIR)$(bindir)/

# GraphViz

# Create phony input files.
draft.fa reads.fq.gz:
	touch $@

# Convert the pipepline graph to GraphViz using makefile2graph.
tigmint-make.gv: bin/tigmint-make draft.fa reads.fq.gz
	makefile2graph -f $< all \
		| gsed -r \
			-e 's/label="(all|arcs|tigmint)".*]/label="\1", shape=ellipse, style=filled]/' \
			-e 's/label="(draft.tigmint.fa|draft.tigmint.arcs.fa)".*]/label="\1", shape=parallelogram, style=filled]/' \
			-e 's/color="green"/shape=parallelogram, style=filled/;s/color="red"/shape=rectangle/' \
		| tred >$@

# Render a GraphViz file to PNG.
%.gv.png: %.gv
	dot -Tpng -Gsize=4 -Gdpi=300 -o $@ $<

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
tigmint-make.xml: bin/tigmint-make xml-patch-make/make-4.1/make-4.1/make
	xml-patch-make/make-4.1/make-4.1/make --xml $@ -f $< all

# Generate a XML pipeline.
%.xml: %.fa %.lrsim.fq.gz xml-patch-make/make-4.1/make-4.1/make
	xml-patch-make/make-4.1/make-4.1/make --xml $@ -f bin/tigmint-make tigmint \
		draft=$* reads=$*.lrsim ref=$* G=16569

# Convert Makefile XML to CWL.
%.cwl: %.xml xml-patch-make/stylesheets/graph2cwl.xsl
	xsltproc -o $@ --stringparam shellpath $(PWD)/$*.cwl.sh xml-patch-make/stylesheets/graph2cwl.xsl $<

# Create a CWL JSON driver script.
%.cwl.json: %.cwl
	printf '#!/usr/bin/env cwl-runner\n{ "cwl:tool": "$<#main" }\n' >$@

# Generate human mitochondrial test data.

# Download the human mitochondrial genome.
mt.fa:
	curl ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz \
		| seqtk seq >$@

# Cut the mitochondrial genome in two.
mt1.fa: mt.fa
	samtools faidx $< MT:1-8284 | seqtk seq | sed 's/>.*/>MT1/' >$@

# Cut the mitochondrial genome in two.
mt2.fa: mt.fa
	samtools faidx $< MT:8285-16569 | seqtk seq | sed 's/>.*/>MT2/' >$@

# Concatenate reads from the two halves of the genome.
mt.halved.lrsim.fq.gz: mt1.lrsim.fq.gz mt2.lrsim.fq.gz
	cat $^ >$@

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
mt.tigmint.fa: %.tigmint.fa: %.fa %.halved.lrsim.fq.gz
	bin/tigmint-make tigmint draft=$* reads=$*.halved.lrsim depth_threshold=150 starts_threshold=2 ref=$* G=16569

# Generate yeast test data.

# Download the yeast genome.
scerevisiae/scerevisiae.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ensembl.org/pub/release-92/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
		| seqtk seq >$@

# Simulate linked reads using LRSIM.
scerevisiae/%.lrsim_S1_L001_R1_001.fastq.gz scerevisiae/%.lrsim_S1_L001_R2_001.fastq.gz: scerevisiae/%.fa
	cd $(@D) && simulateLinkedReads -g $(<F) -p $*.lrsim -o -x2 -f50 -t5 -m10 -z$t

# Add assembly errors.
%.fuse.fa: %.fa
	gsed '3~4d' $< | seqtk seq >$@

# Run Tigmint on the yeast test data.
scerevisiae/%.tigmint.fa: scerevisiae/%.fa scerevisiae/scerevisiae.lrsim.fq.gz
	$(PWD)/bin/tigmint-make -C $(@D) tigmint draft=$* reads=scerevisiae.lrsim ref=scerevisiae G=12157105

# Bedtools

# Create a bedgraph coverage track from a BED file.
%.bedgraph: %.bed mt.fa.fai
	bedtools genomecov -bg -g mt.fa.fai -i $< >$@

# Convert BED to BAM.
scerevisiae/scerevisiae.fuse.%.bed.bam: scerevisiae/scerevisiae.fuse.%.bed scerevisiae/scerevisiae.fuse.fa.fai
	awk '$$2 != $$3' $< | bedtools bedtobam -i - -g scerevisiae/scerevisiae.fuse.fa.fai | samtools sort -@$t -Obam -o $@

# Samtools

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# ChromeQC

# Report summary statistics of a Chromium library using ChromeQC.
%.chromeqc.html: %.bed
	Rscript -e 'rmarkdown::render("chromeqc/report/summary.rmd", "html_document", "$(PWD)/$@", params = list(molecules_bed="$(PWD)/$<"))'
