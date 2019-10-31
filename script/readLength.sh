#!/bin/sh
# compute reads length distribution from a fastq file
zcat G3-WT-d1_S1_L001_trim2.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > G3-WT-d1_S1_L001_trim2.readLength.txt