#!/bin/bash
source /Users/geoffroyandrieux/.bash_profile
cd /Users/geoffroyandrieux/Research/CASTSeq/pipelineGit/script/
Rscript CAST-Seq.R --pipeline "crispr"\
				   --sampleDname "G3_TOY"\
				   --sampleName "G3_treated"\
				   --controlName "G3_UNtreated"\
				   --homeD "../../"\
				   --distCutoff 1500\
				   --species "hg"\
				   --random 1000\
				   --saveReads "no"\
				   --cpu 4
