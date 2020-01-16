#!/bin/sh
# filter out reads with length > 300
#gunzip -c G3-WT-d1_S1_L001_merged.fastq.gz | awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>=300{print a”\n”b”\n”c”\n”$0;}' - | gzip -c - > G3-WT-d1_S1_L001_merged_filt.fastq.gz

#awk '{y= i++ % 4 ; L[y]=$0; if(y==3 && length(L[1])<=27) {printf("%s\n%s\n%s\n%s\n",L[0],L[1],L[2],L[3]);}}'



gunzip -c G3-WT-d1_S1_L001_merged.fastq.gz | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 0 && length(seq) <= 300) {print header, seq, qheader, qseq}}' - | gzip -c - > G3-WT-d1_S1_L001_merged_filt.fastq.gz

