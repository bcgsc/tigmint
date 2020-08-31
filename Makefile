.DELETE_ON_ERROR:
.SECONDARY:

# Number of threads.
t=16

# Parallel gzip.
gzip=pigz -p$t

SHELL=bash -eu -o pipefail

all: tigmint-make.cwl tigmint-make.gv.svg

check: mt/mt.tigmint.fa
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
mt/mt.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz \
		| seqtk seq >$@

# Cut the mitochondrial genome in two.
mt/mt1.fa: mt/mt.fa
	samtools faidx $< MT:1-8284 | seqtk seq | sed 's/>.*/>MT1/' >$@

# Cut the mitochondrial genome in two.
mt/mt2.fa: mt/mt.fa
	samtools faidx $< MT:8285-16569 | seqtk seq | sed 's/>.*/>MT2/' >$@

# Concatenate reads from the two halves of the genome.
mt/mt.halved.lrsim.fq.gz: mt/mt1.lrsim.fq.gz mt/mt2.lrsim.fq.gz
	cat $^ >$@

# Simulate linked reads using LRSIM.
mt/%.lrsim_S1_L001_R1_001.fastq.gz mt/%.lrsim_S1_L001_R2_001.fastq.gz: mt/%.fa
	cd $(@D) && simulateLinkedReads -g $(<F) -p $*.lrsim -o -x0.005 -f4 -t1 -m1 -z1

# Convert paired FASTQ to the interleaved FASTQ format of longranger basic.
%.fq.gz: %_S1_L001_R1_001.fastq.gz %_S1_L001_R2_001.fastq.gz
	seqtk mergepe $^ \
		| paste - - - - - - - - \
		| awk -v OFS='\n' '{ print $$1 " BX:Z:" substr($$3, 1, 15) "-1", substr($$3, 16), "+", substr($$5, 16), $$6 " BX:Z:" substr($$3, 1, 15) "-1", $$8, $$9, $$10 }' \
		| $(gzip) >$@

# Run Tigmint on the mitochondrial test data.
mt/%.tigmint.fa: mt/%.fa mt/%.halved.lrsim.fq.gz
	$(PWD)/bin/tigmint-make -C $(@D) tigmint draft=$* reads=$*.halved.lrsim ref=$* G=16569

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
scerevisiae/%.fuse.fa: scerevisiae/%.fa
	gsed '3~4d' $< | seqtk seq >$@

# Run Tigmint on the yeast test data.
scerevisiae/%.tigmint.fa: scerevisiae/%.fa scerevisiae/scerevisiae.lrsim.fq.gz
	$(PWD)/bin/tigmint-make -C $(@D) tigmint draft=$* reads=scerevisiae.lrsim ref=scerevisiae G=12157105

# Convert BED to BAM.
scerevisiae/scerevisiae.fuse.%.bed.bam: scerevisiae/scerevisiae.fuse.%.bed scerevisiae/scerevisiae.fuse.fa.fai
	bedtools bedtobam -i $< -g scerevisiae/scerevisiae.fuse.fa.fai | samtools sort -@$t -Obam -o $@

# Generate Caenorhabditis elegans test data.

# Download the C. elegans genome.
celegans/celegans.orig.fa:
	mkdir -p $(@D)
	curl ftp://ftp.ensembl.org/pub/release-92/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz \
		| seqtk seq >$@

# Remove the mitochondrial chromosome.
celegans/celegans.fa: celegans/celegans.orig.fa
	sed '/MtDNA/,+1d' $< >$@

# Simulate linked reads using LRSIM.
celegans/%.lrsim_S1_L001_R1_001.fastq.gz celegans/%.lrsim_S1_L001_R2_001.fastq.gz: celegans/%.fa
	cd $(@D) && simulateLinkedReads -g $(<F) -p $*.lrsim -o -x20 -f100 -t20 -m10 -z$t

# Add assembly errors.
celegans/%.fuse.fa: celegans/%.fa
	gsed '3~6d;5~6d' $< | seqtk seq >$@

# Run Tigmint on the C. elegans test data.
celegans/%.tigmint.fa: celegans/%.fa celegans/celegans.lrsim.fq.gz
	$(PWD)/bin/tigmint-make -C $(@D) tigmint draft=$* reads=celegans.lrsim ref=celegans G=100286401

# Convert BED to BAM.
celegans/celegans.fuse.%.bed.bam: celegans/celegans.fuse.%.bed celegans/celegans.fuse.fa.fai
	bedtools bedtobam -i $< -g celegans/celegans.fuse.fa.fai | samtools sort -@$t -Obam -o $@

# Bedtools

# Create a bedgraph coverage track from a BED file.
%.bedgraph: %.bed mt.fa.fai
	bedtools genomecov -bg -g mt.fa.fai -i $< >$@

# Samtools

# Sort a BX-sorted BAM file by position.
%.sort.bam: %.sortbx.bam
	samtools sort -@$t -T$$(mktemp -u -t $(@F).XXXXXX) -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# ChromeQC

# Report summary statistics of a Chromium library using ChromeQC.
%.chromeqc.html: %.bed
	Rscript -e 'rmarkdown::render("chromeqc/report/summary.rmd", "html_document", "$(PWD)/$@", params = list(molecules_bed="$(PWD)/$<"))'
