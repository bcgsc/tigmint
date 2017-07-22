#!/usr/bin/env cwl-runner
cwlVersion: v1.0
$graph:


- id: tool2
  class: CommandLineTool
  label: "reads.fq.gz"
  doc: "reads.fq.gz"
  inputs:
    2_target:
      label: "reads.fq.gz"
      doc: "reads.fq.gz"
      type: string
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok2.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "2" ]



- id: tool1
  class: CommandLineTool
  label: "reads.bx.fq.gz"
  doc: "reads.bx.fq.gz"
  inputs:
    1_target:
      label: "reads.bx.fq.gz"
      doc: "reads.bx.fq.gz"
      type: string
    1_dep2:
      label: "reads.fq.gz"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok1.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "1" ]



- id: tool4
  class: CommandLineTool
  label: "draft.fa"
  doc: "draft.fa"
  inputs:
    4_target:
      label: "draft.fa"
      doc: "draft.fa"
      type: string
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok4.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "4" ]



- id: tool3
  class: CommandLineTool
  label: "draft.fa.bwt"
  doc: "draft.fa.bwt"
  inputs:
    3_target:
      label: "draft.fa.bwt"
      doc: "draft.fa.bwt"
      type: string
    3_dep4:
      label: "draft.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok3.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "3" ]



- id: tool5
  class: CommandLineTool
  label: "draft.reads.bx.bam"
  doc: "draft.reads.bx.bam"
  inputs:
    5_target:
      label: "draft.reads.bx.bam"
      doc: "draft.reads.bx.bam"
      type: string
    5_dep1:
      label: "reads.bx.fq.gz"
      type: File
    5_dep3:
      label: "draft.fa.bwt"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok5.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "5" ]



- id: tool6
  class: CommandLineTool
  label: "draft.reads.bx.bam.bai"
  doc: "draft.reads.bx.bam.bai"
  inputs:
    6_target:
      label: "draft.reads.bx.bam.bai"
      doc: "draft.reads.bx.bam.bai"
      type: string
    6_dep5:
      label: "draft.reads.bx.bam"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok6.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "6" ]



- id: tool7
  class: CommandLineTool
  label: "draft.reads.bx.as100.bam"
  doc: "draft.reads.bx.as100.bam"
  inputs:
    7_target:
      label: "draft.reads.bx.as100.bam"
      doc: "draft.reads.bx.as100.bam"
      type: string
    7_dep5:
      label: "draft.reads.bx.bam"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok7.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "7" ]



- id: tool8
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam"
  doc: "draft.reads.bx.as100.nm5.bam"
  inputs:
    8_target:
      label: "draft.reads.bx.as100.nm5.bam"
      doc: "draft.reads.bx.as100.nm5.bam"
      type: string
    8_dep7:
      label: "draft.reads.bx.as100.bam"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok8.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "8" ]



- id: tool9
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.bai"
  doc: "draft.reads.bx.as100.nm5.bam.bai"
  inputs:
    9_target:
      label: "draft.reads.bx.as100.nm5.bam.bai"
      doc: "draft.reads.bx.as100.nm5.bam.bai"
      type: string
    9_dep8:
      label: "draft.reads.bx.as100.nm5.bam"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok9.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "9" ]



- id: tool10
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.bx.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.bx.tsv"
  inputs:
    10_target:
      label: "draft.reads.bx.as100.nm5.bam.bx.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.bx.tsv"
      type: string
    10_dep8:
      label: "draft.reads.bx.as100.nm5.bam"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok10.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "10" ]



- id: tool11
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.tsv"
  inputs:
    11_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.tsv"
      type: string
    11_dep10:
      label: "draft.reads.bx.as100.nm5.bam.bx.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok11.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "11" ]



- id: tool12
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
  inputs:
    12_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
      type: string
    12_dep11:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok12.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "12" ]



- id: tool13
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html"
  inputs:
    13_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html"
      type: string
    13_dep12:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok13.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "13" ]



- id: tool14
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed"
  inputs:
    14_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed"
      type: string
    14_dep12:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok14.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "14" ]



- id: tool15
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed"
  inputs:
    15_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed"
      type: string
    15_dep14:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.bed"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok15.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "15" ]



- id: tool16
  class: CommandLineTool
  label: "draft.fa.fai"
  doc: "draft.fa.fai"
  inputs:
    16_target:
      label: "draft.fa.fai"
      doc: "draft.fa.fai"
      type: string
    16_dep4:
      label: "draft.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok16.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "16" ]



- id: tool17
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv"
  inputs:
    17_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv"
      type: string
    17_dep15:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed"
      type: File
    17_dep16:
      label: "draft.fa.fai"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok17.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "17" ]



- id: tool18
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv"
  inputs:
    18_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv"
      type: string
    18_dep17:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok18.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "18" ]



- id: tool19
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv"
  inputs:
    19_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv"
      type: string
    19_dep15:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed"
      type: File
    19_dep16:
      label: "draft.fa.fai"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok19.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "19" ]



- id: tool20
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv"
  inputs:
    20_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv"
      type: string
    20_dep19:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok20.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "20" ]



- id: tool21
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv"
  inputs:
    21_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv"
      type: string
    21_dep12:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok21.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "21" ]



- id: tool22
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv"
  inputs:
    22_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv"
      type: string
    22_dep20:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth100.tsv"
      type: File
    22_dep21:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.starts.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok22.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "22" ]



- id: tool23
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv"
  inputs:
    23_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv"
      type: string
    23_dep22:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.tsv"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok23.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "23" ]



- id: tool24
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed"
  inputs:
    24_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed"
      type: string
    24_dep23:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tsv"
      type: File
    24_dep16:
      label: "draft.fa.fai"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok24.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "24" ]



- id: tool25
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
  inputs:
    25_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: string
    25_dep24:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.bed"
      type: File
    25_dep4:
      label: "draft.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok25.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "25" ]



- id: tool26
  class: CommandLineTool
  label: "draft.tigmint.fa"
  doc: "draft.tigmint.fa"
  inputs:
    26_target:
      label: "draft.tigmint.fa"
      doc: "draft.tigmint.fa"
      type: string
    26_dep25:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok26.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "26" ]



- id: tool27
  class: CommandLineTool
  label: "tigmint"
  doc: "tigmint"
  inputs:
    27_target:
      label: "tigmint"
      doc: "tigmint"
      type: string
    27_dep6:
      label: "draft.reads.bx.bam.bai"
      type: File
    27_dep9:
      label: "draft.reads.bx.as100.nm5.bam.bai"
      type: File
    27_dep13:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.summary.html"
      type: File
    27_dep18:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv"
      type: File
    27_dep25:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: File
    27_dep26:
      label: "draft.tigmint.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok27.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "27" ]



- id: tool28
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt"
  inputs:
    28_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt"
      type: string
    28_dep25:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok28.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "28" ]



- id: tool29
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam"
  inputs:
    29_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam"
      type: string
    29_dep28:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa.bwt"
      type: File
    29_dep1:
      label: "reads.bx.fq.gz"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok29.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "29" ]



- id: tool30
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv"
  inputs:
    30_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv"
      type: string
    30_dep29:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.bx.sortn.bam"
      type: File
    30_dep25:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok30.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "30" ]



- id: tool31
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv"
  inputs:
    31_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv"
      type: string
    31_dep30:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs_original.gv"
      type: File
    31_dep25:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok31.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "31" ]



- id: tool32
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
  inputs:
    32_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
      type: string
    32_dep31:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.links.tsv"
      type: File
    32_dep25:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok32.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "32" ]



- id: tool33
  class: CommandLineTool
  label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
  doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
  inputs:
    33_target:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      doc: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      type: string
    33_dep32:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.scaffolds.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok33.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "33" ]



- id: tool34
  class: CommandLineTool
  label: "draft.tigmint.arcs.fa"
  doc: "draft.tigmint.arcs.fa"
  inputs:
    34_target:
      label: "draft.tigmint.arcs.fa"
      doc: "draft.tigmint.arcs.fa"
      type: string
    34_dep33:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok34.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "34" ]



- id: tool35
  class: CommandLineTool
  label: "arcs"
  doc: "arcs"
  inputs:
    35_target:
      label: "arcs"
      doc: "arcs"
      type: string
    35_dep33:
      label: "draft.reads.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth100.starts.breakpoints.tigs.reads.c5_e30000_r0.05.arcs.a0.1_l10.links.fa"
      type: File
    35_dep34:
      label: "draft.tigmint.arcs.fa"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok35.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "35" ]



- id: tool36
  class: CommandLineTool
  label: "all"
  doc: "all"
  inputs:
    36_target:
      label: "all"
      doc: "all"
      type: string
    36_dep27:
      label: "tigmint"
      type: File
    36_dep35:
      label: "arcs"
      type: File
  outputs:
    output:
      type: File
      outputBinding:
          glob: "ok36.txt"
  baseCommand: [ /Users/sjackman/work/tigmint/tigmint-make.cwl.sh, "36" ]




- id: main
  class: Workflow
  inputs: []
  outputs:
    outfile:
      type: File
      outputSource: step36/output
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
      out: [ output ]
      run: "#tool3"

    step5:
      in:
        5_target: { default: __5.ok.flag }
        5_dep1: step1/output
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
        7_dep5: step5/output
      out: [ output ]
      run: "#tool7"

    step8:
      in:
        8_target: { default: __8.ok.flag }
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
      out: [ output ]
      run: "#tool10"

    step11:
      in:
        11_target: { default: __11.ok.flag }
        11_dep10: step10/output
      out: [ output ]
      run: "#tool11"

    step12:
      in:
        12_target: { default: __12.ok.flag }
        12_dep11: step11/output
      out: [ output ]
      run: "#tool12"

    step13:
      in:
        13_target: { default: __13.ok.flag }
        13_dep12: step12/output
      out: [ output ]
      run: "#tool13"

    step14:
      in:
        14_target: { default: __14.ok.flag }
        14_dep12: step12/output
      out: [ output ]
      run: "#tool14"

    step15:
      in:
        15_target: { default: __15.ok.flag }
        15_dep14: step14/output
      out: [ output ]
      run: "#tool15"

    step16:
      in:
        16_target: { default: __16.ok.flag }
        16_dep4: step4/output
      out: [ output ]
      run: "#tool16"

    step17:
      in:
        17_target: { default: __17.ok.flag }
        17_dep15: step15/output
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
        19_dep15: step15/output
        19_dep16: step16/output
      out: [ output ]
      run: "#tool19"

    step20:
      in:
        20_target: { default: __20.ok.flag }
        20_dep19: step19/output
      out: [ output ]
      run: "#tool20"

    step21:
      in:
        21_target: { default: __21.ok.flag }
        21_dep12: step12/output
      out: [ output ]
      run: "#tool21"

    step22:
      in:
        22_target: { default: __22.ok.flag }
        22_dep20: step20/output
        22_dep21: step21/output
      out: [ output ]
      run: "#tool22"

    step23:
      in:
        23_target: { default: __23.ok.flag }
        23_dep22: step22/output
      out: [ output ]
      run: "#tool23"

    step24:
      in:
        24_target: { default: __24.ok.flag }
        24_dep23: step23/output
        24_dep16: step16/output
      out: [ output ]
      run: "#tool24"

    step25:
      in:
        25_target: { default: __25.ok.flag }
        25_dep24: step24/output
        25_dep4: step4/output
      out: [ output ]
      run: "#tool25"

    step26:
      in:
        26_target: { default: __26.ok.flag }
        26_dep25: step25/output
      out: [ output ]
      run: "#tool26"

    step27:
      in:
        27_target: { default: __27.ok.flag }
        27_dep6: step6/output
        27_dep9: step9/output
        27_dep13: step13/output
        27_dep18: step18/output
        27_dep25: step25/output
        27_dep26: step26/output
      out: [ output ]
      run: "#tool27"

    step28:
      in:
        28_target: { default: __28.ok.flag }
        28_dep25: step25/output
      out: [ output ]
      run: "#tool28"

    step29:
      in:
        29_target: { default: __29.ok.flag }
        29_dep28: step28/output
        29_dep1: step1/output
      out: [ output ]
      run: "#tool29"

    step30:
      in:
        30_target: { default: __30.ok.flag }
        30_dep29: step29/output
        30_dep25: step25/output
      out: [ output ]
      run: "#tool30"

    step31:
      in:
        31_target: { default: __31.ok.flag }
        31_dep30: step30/output
        31_dep25: step25/output
      out: [ output ]
      run: "#tool31"

    step32:
      in:
        32_target: { default: __32.ok.flag }
        32_dep31: step31/output
        32_dep25: step25/output
      out: [ output ]
      run: "#tool32"

    step33:
      in:
        33_target: { default: __33.ok.flag }
        33_dep32: step32/output
      out: [ output ]
      run: "#tool33"

    step34:
      in:
        34_target: { default: __34.ok.flag }
        34_dep33: step33/output
      out: [ output ]
      run: "#tool34"

    step35:
      in:
        35_target: { default: __35.ok.flag }
        35_dep33: step33/output
        35_dep34: step34/output
      out: [ output ]
      run: "#tool35"

    step36:
      in:
        36_target: { default: __36.ok.flag }
        36_dep27: step27/output
        36_dep35: step35/output
      out: [ output ]
      run: "#tool36"

