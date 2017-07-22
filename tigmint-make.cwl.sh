#!/bin/bash
function die () {
    echo 1>&2 "ERROR: $0 : $@"
    exit 1
}
if [ "$#" -ne 1 ]; then
    die "Illegal number of parameters"
fi

oldpwd=${PWD}




function run2() {
	cd '/Users/sjackman/work/tigmint'
	
	
	##touch  -c "reads.fq.gz""
	cd "${oldpwd}"
	touch "ok2.txt"
	}




function run1() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "reads.fq.gz" ]; then
		die "File reads.fq.gz missing"
	fi
	
		gunzip -c reads.fq.gz | gawk '  		{ bx = "NA" }  		match($0, "BX:Z:([ACGT]*)-1", x) { bx = x[1] }  		bx == "NA" { getline; getline; getline; next }  		{ print $1 "_" bx " " $2; getline; print; getline; print; getline; print }'  		| gzip >reads.bx.fq.gz
	
	##touch  -c "reads.bx.fq.gz""
	cd "${oldpwd}"
	touch "ok1.txt"
	}




function run4() {
	cd '/Users/sjackman/work/tigmint'
	
	
	##touch  -c "draft.fa""
	cd "${oldpwd}"
	touch "ok4.txt"
	}




function run3() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.fa" ]; then
		die "File draft.fa missing"
	fi
	
		bwa index draft.fa
	
	##touch  -c "draft.fa.bwt""
	cd "${oldpwd}"
	touch "ok3.txt"
	}




function run5() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "reads.bx.fq.gz" ]; then
		die "File reads.bx.fq.gz missing"
	fi
	
	if [ ! -f "draft.fa.bwt" ]; then
		die "File draft.fa.bwt missing"
	fi
	
		bwa mem -t8 -pC draft.fa reads.bx.fq.gz | samtools view -h -F4 | samtools sort -@8 -o draft.reads.bx.bam
	
	##touch  -c "draft.reads.bx.bam""
	cd "${oldpwd}"
	touch "ok5.txt"
	}




function run6() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.bam" ]; then
		die "File draft.reads.bx.bam missing"
	fi
	
		samtools index draft.reads.bx.bam
	
	##touch  -c "draft.reads.bx.bam.bai""
	cd "${oldpwd}"
	touch "ok6.txt"
	}




function run7() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.bam" ]; then
		die "File draft.reads.bx.bam missing"
	fi
	
		samtools view -h -F4 draft.reads.bx.bam | gawk -F'\t' '  			/^@/ { print; next }  			{ as = 0 }  			match($0, "AS:.:([^\t]*)", x) { as = x[1] }  			as >= 100'  		| samtools view -@8 -b -o draft.reads.bx.as100.bam
	
	##touch  -c "draft.reads.bx.as100.bam""
	cd "${oldpwd}"
	touch "ok7.txt"
	}




function run8() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.bam" ]; then
		die "File draft.reads.bx.as100.bam missing"
	fi
	
		samtools view -h -F4 draft.reads.bx.as100.bam  		| gawk -F'\t' '  			/^@/ { print; next }  			{ nm = 999999999 }  			match($0, "NM:i:([^\t]*)", x) { nm = x[1] }  			nm < 5'  		| samtools view -@8 -b -o draft.reads.bx.as100.nm5.bam
	
	##touch  -c "draft.reads.bx.as100.nm5.bam""
	cd "${oldpwd}"
	touch "ok8.txt"
	}




function run9() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam" ]; then
		die "File draft.reads.bx.as100.nm5.bam missing"
	fi
	
		samtools index draft.reads.bx.as100.nm5.bam
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.bai""
	cd "${oldpwd}"
	touch "ok9.txt"
	}




function run10() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam" ]; then
		die "File draft.reads.bx.as100.nm5.bam missing"
	fi
	
		samtools view -F4 draft.reads.bx.as100.nm5.bam | gawk -F'\t' '  		BEGIN { print "Flags\tRname\tPos\tMapq\tAS\tNM\tBX\tMI" }  		{ as = bx = mi = nm = "NA" }  		match($0, "AS:.:([^\t]*)", x) { as = x[1] }  		match($0, "NM:.:([^\t]*)", x) { nm = x[1] }  		match($0, "BX:Z:([^\t]*)", x) { bx = x[1] }  		match($0, "MI:i:([^\t]*)", x) { mi = x[1] }  		{ print $2 "\t" $3 "\t" $4 "\t" $5 "\t" as "\t" nm "\t" bx "\t" mi }' >draft.reads.bx.as100.nm5.bam.bx.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.bx.tsv""
	cd "${oldpwd}"
	touch "ok10.txt"
	}




function run11() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.bx.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.bx.tsv missing"
	fi
	
		./mi.r draft.reads.bx.as100.nm5.bam.bx.tsv draft.reads.bx.as100.nm5.bam.mi.bx.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.tsv""
	cd "${oldpwd}"
	touch "ok11.txt"
	}




function run12() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.tsv missing"
	fi
	
		mlr --tsvlite  		then stats1 -g BX,MI,Rname -a count,min,p50,max -f Pos,Mapq,AS,NM  		then rename Pos_min,Start,Pos_max,End,Mapq_p50,Mapq_median,AS_p50,AS_median,NM_p50,NM_median,Pos_count,Reads  		then put '$Size = $End - $Start'  		then cut -o -f Rname,Start,End,Size,BX,MI,Reads,Mapq_median,AS_median,NM_median  		then filter '$Reads >= 4'  		draft.reads.bx.as100.nm5.bam.mi.bx.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv""
	cd "${oldpwd}"
	touch "ok12.txt"
	}




function run13() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv missing"
	fi
	
		Rscript -e 'rmarkdown::render("summary.rmd", "html_document", "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html", params = list(input_tsv="draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv", output_tsv="draft.reads.bx.as100.nm5.bam.mi.summary.tsv"))'
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html""
	cd "${oldpwd}"
	touch "ok13.txt"
	}




function run14() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv missing"
	fi
	
		mlr --tsvlite --headerless-csv-output  		put '$Start = $Start - 1; $End = $End - 1'  		then put '$Name = "Reads=" . $Reads . ",Size=" . $Size . ",Mapq=" . $Mapq_median . ",AS=" . $AS_median . ",NM=" . $NM_median . ",BX=" . $BX . ",MI=" . $MI'  		then cut -o -f Rname,Start,End,Name,Reads draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed""
	cd "${oldpwd}"
	touch "ok14.txt"
	}




function run15() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed missing"
	fi
	
		awk '$3 - $2 >= 2000' draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed""
	cd "${oldpwd}"
	touch "ok15.txt"
	}




function run16() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.fa" ]; then
		die "File draft.fa missing"
	fi
	
		samtools faidx draft.fa
	
	##touch  -c "draft.fa.fai""
	cd "${oldpwd}"
	touch "ok16.txt"
	}




function run17() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed missing"
	fi
	
	if [ ! -f "draft.fa.fai" ]; then
		die "File draft.fa.fai missing"
	fi
	
		(printf "Rname\tDepth\tCount\tRsize\tFraction\n"; awk '$2 != $3' draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed | bedtools genomecov -g draft.fa.fai -i -) >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv""
	cd "${oldpwd}"
	touch "ok17.txt"
	}




function run18() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv missing"
	fi
	
		mlr --tsvlite  		then filter '$Rname == "genome" && $Depth > 0'  		then step -a rsum -f Fraction  		then put -q '@Depth_count += $Count; if (is_null(@p25) && $Fraction_rsum >= 0.25) { @p25 = $Depth }; if (is_null(@p50) && $Fraction_rsum >= 0.50) { @p50 = $Depth }; if (is_null(@p75) && $Fraction_rsum >= 0.75) { @p75 = $Depth } end { emitf @Depth_count, @p25, @p50, @p75 }'  		then rename p25,Depth_p25,p50,Depth_p50,p75,Depth_p75  		then put '$Depth_IQR = $Depth_p75 - $Depth_p25'  		draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv""
	cd "${oldpwd}"
	touch "ok18.txt"
	}




function run19() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed missing"
	fi
	
	if [ ! -f "draft.fa.fai" ]; then
		die "File draft.fa.fai missing"
	fi
	
		(printf "Rname\tPos\tDepth\n"; awk '$2 != $3' draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed | bedtools genomecov -d -g draft.fa.fai -i -) >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv""
	cd "${oldpwd}"
	touch "ok19.txt"
	}




function run20() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv missing"
	fi
	
		mlr --tsvlite filter '$Depth < 100' draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv""
	cd "${oldpwd}"
	touch "ok20.txt"
	}




function run21() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv missing"
	fi
	
		mlr --tsvlite  		then filter '$Size >= 2000'  		then count-distinct -f Rname,Start  		then rename Start,Pos,count,Starts  		then sort -f Rname -n Pos  		draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv""
	cd "${oldpwd}"
	touch "ok21.txt"
	}




function run22() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv missing"
	fi
	
		mlr --tsvlite join -u -j Rname,Pos -f draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv""
	cd "${oldpwd}"
	touch "ok22.txt"
	}




function run23() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv missing"
	fi
	
		mlr --tsvlite filter '$Depth < 100 && $Starts >= 4 && $Pos >= 1000' draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv""
	cd "${oldpwd}"
	touch "ok23.txt"
	}




function run24() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv missing"
	fi
	
	if [ ! -f "draft.fa.fai" ]; then
		die "File draft.fa.fai missing"
	fi
	
		Rscript -e 'rmarkdown::render("breaktigs.rmd", "html_notebook", "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.nb.html", params = list(input_tsv="draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv", input_fai="draft.fa.fai", output_bed="draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed"))'
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed""
	cd "${oldpwd}"
	touch "ok24.txt"
	}




function run25() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed missing"
	fi
	
	if [ ! -f "draft.fa" ]; then
		die "File draft.fa missing"
	fi
	
		bedtools getfasta -name -fi draft.fa -bed draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed | sed 's/::/ /;s/^NN*//;s/NN*$//' >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa""
	cd "${oldpwd}"
	touch "ok25.txt"
	}




function run26() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa missing"
	fi
	
		ln -sf draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa draft.tigmint.fa
	
	##touch  -c "draft.tigmint.fa""
	cd "${oldpwd}"
	touch "ok26.txt"
	}




function run27() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.bam.bai" ]; then
		die "File draft.reads.bx.bam.bai missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.bai" ]; then
		die "File draft.reads.bx.as100.nm5.bam.bai missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa missing"
	fi
	
	if [ ! -f "draft.tigmint.fa" ]; then
		die "File draft.tigmint.fa missing"
	fi
	
	
	##touch "__27_phony.flag""
	cd "${oldpwd}"
	touch "ok27.txt"
	}




function run28() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa missing"
	fi
	
		bwa index draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt""
	cd "${oldpwd}"
	touch "ok28.txt"
	}




function run29() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt missing"
	fi
	
	if [ ! -f "reads.bx.fq.gz" ]; then
		die "File reads.bx.fq.gz missing"
	fi
	
		bwa mem -t8 -pC draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa reads.bx.fq.gz | samtools view -@8 -h -F4 -o draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam""
	cd "${oldpwd}"
	touch "ok29.txt"
	}




function run30() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa missing"
	fi
	
		arcs -s98 -c5 -l0 -z500 -m4-20000 -d0 -e30000 -r0.05 -v  		-f draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa  		-b draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs  		-g draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.dist.gv  		--tsv=draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.dist.tsv  		--barcode-counts=draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam.barcode-counts.tsv  		draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv""
	cd "${oldpwd}"
	touch "ok30.txt"
	}




function run31() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa missing"
	fi
	
		bin/arcs-makeTSVfile draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv""
	cd "${oldpwd}"
	touch "ok31.txt"
	}




function run32() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv missing"
	fi
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa missing"
	fi
	
		cp draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.tigpair_checkpoint.tsv
		LINKS -k20 -l10 -t2 -a0.1 -x1 -s /dev/null -f draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa -b draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa""
	cd "${oldpwd}"
	touch "ok32.txt"
	}




function run33() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa missing"
	fi
	
		gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa >draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa
	
	##touch  -c "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa""
	cd "${oldpwd}"
	touch "ok33.txt"
	}




function run34() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa missing"
	fi
	
		ln -sf draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa draft.tigmint.arcs.fa
	
	##touch  -c "draft.tigmint.arcs.fa""
	cd "${oldpwd}"
	touch "ok34.txt"
	}




function run35() {
	cd '/Users/sjackman/work/tigmint'
	
	if [ ! -f "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa" ]; then
		die "File draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa missing"
	fi
	
	if [ ! -f "draft.tigmint.arcs.fa" ]; then
		die "File draft.tigmint.arcs.fa missing"
	fi
	
	
	##touch "__35_phony.flag""
	cd "${oldpwd}"
	touch "ok35.txt"
	}




function run36() {
	cd '/Users/sjackman/work/tigmint'
	
	
	##touch "__36_phony.flag""
	cd "${oldpwd}"
	touch "ok36.txt"
	}




case "$1" in
        2)
    run2
    ;;
        1)
    run1
    ;;
        4)
    run4
    ;;
        3)
    run3
    ;;
        5)
    run5
    ;;
        6)
    run6
    ;;
        7)
    run7
    ;;
        8)
    run8
    ;;
        9)
    run9
    ;;
        10)
    run10
    ;;
        11)
    run11
    ;;
        12)
    run12
    ;;
        13)
    run13
    ;;
        14)
    run14
    ;;
        15)
    run15
    ;;
        16)
    run16
    ;;
        17)
    run17
    ;;
        18)
    run18
    ;;
        19)
    run19
    ;;
        20)
    run20
    ;;
        21)
    run21
    ;;
        22)
    run22
    ;;
        23)
    run23
    ;;
        24)
    run24
    ;;
        25)
    run25
    ;;
        26)
    run26
    ;;
        27)
    run27
    ;;
        28)
    run28
    ;;
        29)
    run29
    ;;
        30)
    run30
    ;;
        31)
    run31
    ;;
        32)
    run32
    ;;
        33)
    run33
    ;;
        34)
    run34
    ;;
        35)
    run35
    ;;
        36)
    run36
    ;;
	*)
	die "Undefined target id=$1"
	;;
esac

