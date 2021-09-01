#!/bin/bash

# sh FANCF_decoy_rep1_fastq_aln.sh /home/gandri/offTargets/Giando/pipelineGit/ FANCF_decoy_rep1 FANCF-withDecoyGel_S1_L001 UT-FANCF-withDecoyGel_S2_L001

############# DEFINE INPUTS ##################
CPU=12

# $1
#homeDir=/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/
homeDir=$1

# $2
#mainDir=${homeDir}FANCF_decoy/
mainDir=${homeDir}/samples/$2/

# $3 $4
#Samples="FANCF-decoy_L001 FANCF-UT-decoy_L001"
Samples="$3 $4"

# Set up the directories and annotation files
annotDir=${homeDir}annotations/
dataDir=${mainDir}data/
outDir=${mainDir}results/fastq_aln/
outBedDir=${mainDir}results/guide_aln/
mkdir -p ${outDir}
mkdir -p ${outBedDir}

pos=${dataDir}pos.fa
neg=${dataDir}neg.fa
mispriming=${dataDir}mispriming.fa
#linker=${dataDir}linker.fa
linker=${dataDir}linker_RC.fa
hth=${dataDir}headTOhead.fa

adapters=${annotDir}TruSeq4-PE.fa

# Crispresso parameters
#crispGuide=GTGAGTAGAGCGGAGGCAGGAGG
#crispSeq=TGCTCTTCAGCCTTTTGCAGTTTATCAGGATGAGGATGACCAGCATGTTGCCCACAAAACCAAAGATGAACACCAGTGAGTAGAGCGGAGGCAGGAGGCGGGCTGCGATTTGCTTCACATTGATTTTTTGGCAGGGCTCCGATGTATAATAATTGATGTCATAGATTGGACTTGACACTTGATAATCCATCTTGTTCCACCctgtgcataaataaaaagtga

# Repeat masker parameters
chrSize=${annotDir}hg38.chrom.sizes
genomeFasta=${annotDir}bowtie2Index/genome.fa

############# ALIGNMENT ##################
for fname in ${Samples}
do

	logFile=${outDir}${fname}_pipeline.log

	pair1=${dataDir}fastq/${fname}_R1_001.fastq.gz
	pair2=${dataDir}fastq/${fname}_R2_001.fastq.gz

	# Un tar
	#gunzip ${dataDir}fastq/${fname}_R1_001.fastq.gz
	#gunzip ${dataDir}fastq/${fname}_R2_001.fastq.gz


	############## Quality control fastq files #####################
	fastqc -o ${outDir} --noextract ${pair1}
	fastqc -o ${outDir} --noextract ${pair2}

	################# FLASH PAIRING ###################
	echo "START FLASH PAIRING" > ${logFile}

	# FLASH
	flash -t ${CPU} -o ${fname} -d ${outDir} -m 15 -M 250 ${pair2} ${pair1} > ${logFile}# double

	# MERGE PAIRED AND UNPAIRED R2
	# we used notCombined_1.fastq because we used R2 as R1 with flash
	#cat ${outDir}${fname}.extendedFrags.fastq ${outDir}${fname}.notCombined_1.fastq > ${outDir}${fname}_merged.fastq
	cat ${outDir}${fname}.extendedFrags.fastq ${outDir}${fname}.notCombined_1.fastq | gzip -c > ${outDir}${fname}_merged.fastq.gz
	gzip ${outDir}${fname}.extendedFrags.fastq
	gzip ${outDir}${fname}.notCombined_1.fastq
	gzip ${outDir}${fname}.notCombined_2.fastq
	
	echo "_merged.fastq.gz: " >> ${logFile}
	echo $(zcat ${outDir}${fname}_merged.fastq.gz|wc -l)/4|bc >> ${logFile}
	
	
	################# FILTER OUT READS > 300bp ################### 
	#echo "FILTER OUT READS > 300bp" >> ${logFile}
	#gunzip -c ${outDir}${fname}_merged.fastq.gz | awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)<=300{print a”\n”b”\n”c”\n”$0;}' - | gzip -c - > ${outDir}${fname}_merged_filt.fastq.gz
	#gunzip -c ${outDir}${fname}_merged.fastq.gz | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 0 && length(seq) <= 300) {print header, seq, qheader, qseq}}' - | gzip -c - > ${outDir}${fname}_merged_filt.fastq.gz


	############# FILTERING ###################
	echo "START FILTERING" >> ${logFile}
	# filter 1: find reads that contain the CCR5 sequence before the cut site
	bbduk.sh -Xmx4g in=${outDir}${fname}_merged.fastq.gz ref=${pos} outm=${outDir}${fname}_pos.fastq.gz outu=${outDir}${fname}_pos_NOmatch.fastq.gz k=25 mm=f edist=2 ow=t rcomp=f

	echo "_pos.fastq.gz: " >> ${logFile}
	echo $(zcat ${outDir}${fname}_pos.fastq.gz|wc -l)/4|bc >> ${logFile}

	# filter 2: eliminate CCR2 aspecifics
	bbduk.sh -Xmx4g in=${outDir}${fname}_pos.fastq.gz ref=${mispriming} out=${outDir}${fname}_Filt2.fastq.gz k=50 mm=f edist=0 ow=t rcomp=f

	echo "_Filt2.fastq.gz: " >> ${logFile}
	echo $(zcat ${outDir}${fname}_Filt2.fastq.gz|wc -l)/4|bc >> ${logFile}

	### check filtering
	#bbduk.sh in=${outDir}${fname}_Filt2.fastq ref=${dataDir}SeqCCR5pos.fa k=25 mm=f edist=2 ow=t    

	############## TRIMMING #########################
	echo "START TRIMMING" >> ${logFile}
	#trimming of the linker and eliminate short reads <20+ (ml=20)
	#rcomp = t
	bbduk.sh -Xmx4g in=${outDir}${fname}_Filt2.fastq.gz ref=${linker} out=${outDir}${fname}_trim1.fastq.gz k=20 mm=f edist=2 ow=t ktrim=r ml=20 rcomp=f

	echo "_trim1.fastq.gz: " >> ${logFile}
	echo $(zcat ${outDir}${fname}_trim1.fastq.gz|wc -l)/4|bc >> ${logFile}

	# trimming of the adapters
	bbduk.sh -Xmx4g in=${outDir}${fname}_trim1.fastq.gz ref=${adapters} out=${outDir}${fname}_trim2.fastq.gz k=20 mm=f edist=2 ow=t ktrim=r ml=20

	echo "_trim2.fastq.gz: " >> ${logFile}
	echo $(zcat ${outDir}${fname}_trim2.fastq.gz|wc -l)/4|bc >> ${logFile}

	############
	# CRISPRESSO
	echo "START CRISPRESSO" >> ${logFile}

	# filter 3: find reads that not contain the CCR5 sequence after the cut site
	bbduk.sh -Xmx4g in=${outDir}${fname}_trim2.fastq.gz ref=${neg} outm=${outDir}${fname}_Filt3_pos.fastq.gz outu=${outDir}${fname}_Filt3_neg.fastq.gz k=25 mm=f hdist=3 edist=2 ow=t rcomp=f

	# run crispresso with fastq.gz
	#source activate py27
	#CRISPResso -r1 ${outDir}${fname}_Filt3_pos.fastq.gz -a ${crispSeq} -g ${crispGuide} -o ${outDir} --min_identity_score 60 -p ${CPU}
	#source deactivate


	#trimming of the CCR5 sequence before the cut site
	bbduk.sh -Xmx4g in=${outDir}${fname}_trim2.fastq.gz ref=${pos} out=${outDir}${fname}_trim3.fastq.gz k=25 mm=f edist=2 ow=t ktrim=l rcomp=f ml=30

	echo "_trim3.fastq.gz: " >> ${logFile}
	echo $(zcat ${outDir}${fname}_trim3.fastq.gz|wc -l)/4|bc >> ${logFile}

	######### FINAL filter ############ 
	echo "START FINAL FILTER" >> ${logFile}
	# create the separate file for the head to head CCR5 sequences
	bbduk.sh -Xmx4g in=${outDir}${fname}_trim3.fastq.gz ref=${hth} outm=${outDir}${fname}_Headtoheadtrimmed.fastq.gz outu=${outDir}${fname}_final.fastq.gz k=20 mm=f edist=2 ow=t rcomp=f


	############## check quality #########################
	fastqc -o ${outDir} --noextract ${outDir}${fname}_merged.fastq.gz
	#fastqc -o ${outDir} --noextract ${outDir}${fname}_merged_filt.fastq.gz
	#fastqc -o ${outDir} --noextract ${outDir}${fname}_pos.fastq
	#fastqc -o ${outDir} --noextract ${outDir}${fname}_Filt2.fastq
	#fastqc -o ${outDir} --noextract ${outDir}${fname}_trim1.fastq
	#fastqc -o ${outDir} --noextract ${outDir}${fname}_trim2.fastq
	fastqc -o ${outDir} --noextract ${outDir}${fname}_trim3.fastq.gz
	#fastqc -o ${outDir} --noextract ${outDir}${fname}_final.fastq


	############## ALIGNMENT ###########################
	##bowtie alignment of trimmed but not collapsed seq##
	bowtie2 -x ${annotDir}bowtie2Index/genome -U ${outDir}${fname}_trim3.fastq.gz --very-sensitive -p ${CPU} | samtools view -bS -> ${outDir}${fname}_Alignment.bam

	echo "_Alignment.bam: " >> ${logFile}
	echo $(samtools view -c -F 260 ${outDir}${fname}_Alignment.bam) >> ${logFile}
	
	### eliminate mapq < 15 AND select mapq < 15 for repeat masker
	samtools view -h -b -q 15 -U ${outDir}${fname}_lowQ.bam ${outDir}${fname}_Alignment.bam > ${outDir}${fname}_mapqfiltered.bam
	samtools bam2fq -F 260 ${outDir}${fname}_lowQ.bam | seqtk seq -A > ${outDir}${fname}_lowQ.fa
	rm ${outDir}${fname}_Alignment.bam

	#### sort
	samtools sort ${outDir}${fname}_mapqfiltered.bam > ${outDir}${fname}_AlignmentSort.bam
	rm ${outDir}${fname}_mapqfiltered.bam
	
	echo "_AlignmentSort.bam: " >> ${logFile}
	echo $(samtools view -c -F 260 ${outDir}${fname}_AlignmentSort.bam) >> ${logFile}
	
	
	#### convert from bam to bed
	bedtools bamtobed -i ${outDir}${fname}_AlignmentSort.bam > ${outBedDir}${fname}_Alignment.bed
	#convert2bed --input=bam --output=bed < ${outDir}${fname}_AlignmentSort.bam > ${outDir}${fname}_Alignment_V2.bed

	echo "_Alignment.bed: " >> ${logFile}
	echo $(wc -l ${outBedDir}${fname}_Alignment.bed) >> ${logFile}


	###############
	# REPEAT MASKER

	# generate random reads bedrolls, allow overlap, do not keep chromosomes distribution.
	bedTMP=$(mktemp)
	shuffleBedTMP=$(mktemp)

	bedtools bamtobed -i ${outDir}${fname}_lowQ.bam > ${bedTMP}
	bedtools shuffle -i ${bedTMP} -g ${chrSize} > ${shuffleBedTMP}
	bedtools getfasta -fi ${genomeFasta} -bed ${shuffleBedTMP} -fo ${outDir}${fname}_lowQ.shuffle.fa

done






