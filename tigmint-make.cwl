#!/usr/bin/env cwl-runner
cwlVersion: v1.0
$graph:


- id: tool2
  class: CommandLineTool
  label: "draft.fa"
  doc: "draft.fa"
  inputs:
    2_target:
      label: "draft.fa"
      doc: "draft.fa"
      type: string
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok2.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - touch ok2.txt


- id: tool1
  class: CommandLineTool
  label: "draft.fa.bwt"
  doc: "draft.fa.bwt"
  inputs:
    1_target:
      label: "draft.fa.bwt"
      doc: "draft.fa.bwt"
      type: string
    1_dep2:
      label: "draft.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok1.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        bwa index draft.fa;
        ); touch ok1.txt


- id: tool4
  class: CommandLineTool
  label: "reads.fq.gz"
  doc: "reads.fq.gz"
  inputs:
    4_target:
      label: "reads.fq.gz"
      doc: "reads.fq.gz"
      type: string
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok4.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - touch ok4.txt


- id: tool3
  class: CommandLineTool
  label: "draft.reads.sortbx.bam"
  doc: "draft.reads.sortbx.bam"
  inputs:
    3_target:
      label: "draft.reads.sortbx.bam"
      doc: "draft.reads.sortbx.bam"
      type: string
    3_dep4:
      label: "reads.fq.gz"
      type: File
    3_dep1:
      label: "draft.fa.bwt"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok3.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
         bwa mem -t8 -pC draft.fa reads.fq.gz | samtools view -u -F4 | samtools sort -@8 -tBX -T$(mktemp -u -t draft.reads.sortbx.bam.XXXXXX) -o draft.reads.sortbx.bam;
        ); touch ok3.txt


- id: tool5
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.bed"
  doc: "draft.reads.as0.65.nm5.molecule.bed"
  inputs:
    5_target:
      label: "draft.reads.as0.65.nm5.molecule.bed"
      doc: "draft.reads.as0.65.nm5.molecule.bed"
      type: string
    5_dep3:
      label: "draft.reads.sortbx.bam"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok5.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
         bin/tigmint-molecule -a0.65 -n5 -q0 -d50000 -o draft.reads.as0.65.nm5.molecule.bed draft.reads.sortbx.bam;
        ); touch ok5.txt


- id: tool6
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.bed"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.bed"
  inputs:
    6_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.bed"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.bed"
      type: string
    6_dep5:
      label: "draft.reads.as0.65.nm5.molecule.bed"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok6.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        awk '$3 - $2 >= 2000' draft.reads.as0.65.nm5.molecule.bed >draft.reads.as0.65.nm5.molecule.size2000.bed;
        ); touch ok6.txt


- id: tool7
  class: CommandLineTool
  label: "draft.fa.fai"
  doc: "draft.fa.fai"
  inputs:
    7_target:
      label: "draft.fa.fai"
      doc: "draft.fa.fai"
      type: string
    7_dep2:
      label: "draft.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok7.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        samtools faidx draft.fa;
        ); touch ok7.txt


- id: tool8
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
  inputs:
    8_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: string
    8_dep6:
      label: "draft.reads.as0.65.nm5.molecule.size2000.bed"
      type: File
    8_dep2:
      label: "draft.fa"
      type: File
    8_dep7:
      label: "draft.fa.fai"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok8.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
         bin/tigmint-cut -p8 -w1000 -n20 -t0 -o draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa draft.fa draft.reads.as0.65.nm5.molecule.size2000.bed;
        ); touch ok8.txt


- id: tool9
  class: CommandLineTool
  label: "draft.tigmint.fa"
  doc: "draft.tigmint.fa"
  inputs:
    9_target:
      label: "draft.tigmint.fa"
      doc: "draft.tigmint.fa"
      type: string
    9_dep8:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok9.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        ln -sf draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa draft.tigmint.fa;
        ); touch ok9.txt


- id: tool10
  class: CommandLineTool
  label: "tigmint"
  doc: "tigmint"
  inputs:
    10_target:
      label: "tigmint"
      doc: "tigmint"
      type: string
    10_dep8:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: File
    10_dep9:
      label: "draft.tigmint.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok10.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - touch ok10.txt


- id: tool11
  class: CommandLineTool
  label: "reads.bx.fq.gz"
  doc: "reads.bx.fq.gz"
  inputs:
    11_target:
      label: "reads.bx.fq.gz"
      doc: "reads.bx.fq.gz"
      type: string
    11_dep4:
      label: "reads.fq.gz"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok11.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        gunzip -c reads.fq.gz | gawk '  		{ bx = "NA" }  		match($0, "BX:Z:([ACGT]*)-1", x) { bx = x[1] }  		bx == "NA" { getline; getline; getline; next }  		{ print $1 "_" bx " " $2; getline; print; getline; print; getline; print }'  		| pigz -p8 >reads.bx.fq.gz;
        ); touch ok11.txt


- id: tool12
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bwt"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bwt"
  inputs:
    12_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bwt"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bwt"
      type: string
    12_dep8:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok12.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        bwa index draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa;
        ); touch ok12.txt


- id: tool13
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam"
  inputs:
    13_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam"
      type: string
    13_dep12:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa.bwt"
      type: File
    13_dep11:
      label: "reads.bx.fq.gz"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok13.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        bwa mem -t8 -pC draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa reads.bx.fq.gz | samtools view -@8 -h -F4 -o draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam;
        ); touch ok13.txt


- id: tool14
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs_original.gv"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs_original.gv"
  inputs:
    14_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs_original.gv"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs_original.gv"
      type: string
    14_dep13:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam"
      type: File
    14_dep8:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok14.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        arcs -s98 -c5 -l0 -z500 -m4-20000 -d0 -e30000 -r0.05 -v  		-f draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa  		-b draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs  		-g draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.dist.gv  		--tsv=draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.dist.tsv  		--barcode-counts=draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam.barcode-counts.tsv  		draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.bx.sortn.bam;
        ); touch ok14.txt


- id: tool15
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv"
  inputs:
    15_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv"
      type: string
    15_dep14:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs_original.gv"
      type: File
    15_dep8:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok15.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        bin/tigmint-arcs-tsv draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs_original.gv draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa;
        ); touch ok15.txt


- id: tool16
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
  inputs:
    16_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
      type: string
    16_dep15:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv"
      type: File
    16_dep8:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok16.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        cp draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.links.tsv draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.tigpair_checkpoint.tsv;
        LINKS -k20 -l10 -t2 -a0.1 -x1 -s /dev/null -f draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.fa -b draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links;
        ); touch ok16.txt


- id: tool17
  class: CommandLineTool
  label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
  doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
  inputs:
    17_target:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      doc: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      type: string
    17_dep16:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok17.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa >draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa;
        ); touch ok17.txt


- id: tool18
  class: CommandLineTool
  label: "draft.tigmint.arcs.fa"
  doc: "draft.tigmint.arcs.fa"
  inputs:
    18_target:
      label: "draft.tigmint.arcs.fa"
      doc: "draft.tigmint.arcs.fa"
      type: string
    18_dep17:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok18.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - >
        (cd '/Users/sjackman/work/tigmint';
        ln -sf draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa draft.tigmint.arcs.fa;
        ); touch ok18.txt


- id: tool19
  class: CommandLineTool
  label: "arcs"
  doc: "arcs"
  inputs:
    19_target:
      label: "arcs"
      doc: "arcs"
      type: string
    19_dep17:
      label: "draft.reads.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      type: File
    19_dep18:
      label: "draft.tigmint.arcs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok19.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - touch ok19.txt


- id: tool20
  class: CommandLineTool
  label: "all"
  doc: "all"
  inputs:
    20_target:
      label: "all"
      doc: "all"
      type: string
    20_dep10:
      label: "tigmint"
      type: File
    20_dep19:
      label: "arcs"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok20.txt"
  baseCommand:
    - bash
    - -eu
    - -o
    - pipefail
    - -c
    - touch ok20.txt



- id: main
  class: Workflow
  inputs: []
  outputs:
    outfile:
      type: File
      outputSource: step20/output
  steps:

    step2:
      in:
        2_target: { default: __2.ok.flag }
      out: [ output ]
      run: "#tool2"

    step1:
      in:
        1_target: { default: __1.ok.flag }
        1_dep2: step2/output
      out: [ output ]
      run: "#tool1"

    step4:
      in:
        4_target: { default: __4.ok.flag }
      out: [ output ]
      run: "#tool4"

    step3:
      in:
        3_target: { default: __3.ok.flag }
        3_dep4: step4/output
        3_dep1: step1/output
      out: [ output ]
      run: "#tool3"

    step5:
      in:
        5_target: { default: __5.ok.flag }
        5_dep3: step3/output
      out: [ output ]
      run: "#tool5"

    step6:
      in:
        6_target: { default: __6.ok.flag }
        6_dep5: step5/output
      out: [ output ]
      run: "#tool6"

    step7:
      in:
        7_target: { default: __7.ok.flag }
        7_dep2: step2/output
      out: [ output ]
      run: "#tool7"

    step8:
      in:
        8_target: { default: __8.ok.flag }
        8_dep6: step6/output
        8_dep2: step2/output
        8_dep7: step7/output
      out: [ output ]
      run: "#tool8"

    step9:
      in:
        9_target: { default: __9.ok.flag }
        9_dep8: step8/output
      out: [ output ]
      run: "#tool9"

    step10:
      in:
        10_target: { default: __10.ok.flag }
        10_dep8: step8/output
        10_dep9: step9/output
      out: [ output ]
      run: "#tool10"

    step11:
      in:
        11_target: { default: __11.ok.flag }
        11_dep4: step4/output
      out: [ output ]
      run: "#tool11"

    step12:
      in:
        12_target: { default: __12.ok.flag }
        12_dep8: step8/output
      out: [ output ]
      run: "#tool12"

    step13:
      in:
        13_target: { default: __13.ok.flag }
        13_dep12: step12/output
        13_dep11: step11/output
      out: [ output ]
      run: "#tool13"

    step14:
      in:
        14_target: { default: __14.ok.flag }
        14_dep13: step13/output
        14_dep8: step8/output
      out: [ output ]
      run: "#tool14"

    step15:
      in:
        15_target: { default: __15.ok.flag }
        15_dep14: step14/output
        15_dep8: step8/output
      out: [ output ]
      run: "#tool15"

    step16:
      in:
        16_target: { default: __16.ok.flag }
        16_dep15: step15/output
        16_dep8: step8/output
      out: [ output ]
      run: "#tool16"

    step17:
      in:
        17_target: { default: __17.ok.flag }
        17_dep16: step16/output
      out: [ output ]
      run: "#tool17"

    step18:
      in:
        18_target: { default: __18.ok.flag }
        18_dep17: step17/output
      out: [ output ]
      run: "#tool18"

    step19:
      in:
        19_target: { default: __19.ok.flag }
        19_dep17: step17/output
        19_dep18: step18/output
      out: [ output ]
      run: "#tool19"

    step20:
      in:
        20_target: { default: __20.ok.flag }
        20_dep10: step10/output
        20_dep19: step19/output
      out: [ output ]
      run: "#tool20"

