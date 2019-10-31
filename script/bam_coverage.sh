#!/bin/bash

# G3 WT D1
workDir=/home/gandri/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results_old/fastq_aln/
bam=${workDir}G3-WT-d1_S1_L001_AlignmentSort.bam
strandP=${workDir}G3-WT-d1_S1_L001_strand_P.bedgraph
strandN=${workDir}G3-WT-d1_S1_L001_strand_N.bedgraph
unstrand=${workDir}G3-WT-d1_S1_L001.bedgraph

bedtools genomecov -ibam ${bam} -bg -strand + > ${strandP}
bedtools genomecov -ibam ${bam} -bg -strand - > ${strandN}
bedtools genomecov -ibam ${bam} -bg > ${unstrand}
