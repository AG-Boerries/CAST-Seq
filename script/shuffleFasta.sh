#!/bin/bash

# Commands
bam2bedCMD="/home/gandri/Programs/bedtools2-master/bin/bedtools bamtobed"
shuffleCMD="/home/gandri/Programs/bedtools2-master/bin/bedtools shuffle"
getFastaCMD="/home/gandri/Programs/bedtools2-master/bin/bedtools getfasta"

# Parameters
bam=/home/gandri/offTargets/Giando/pipeline/FANCF_decoy/results/fastq_aln/FANCF-decoy_S1_L001_lowQ.bam
bedTMP=$(mktemp)
shuffleBedTMP=$(mktemp)
#bedTMP=/home/gandri/offTargets/Giando/pipeline/FANCF_decoy/results/fastq_aln/FANCF-decoy_S1_L001_lowQ.bed
#shuffleBedTMP=/home/gandri/offTargets/Giando/pipeline/FANCF_decoy/results/fastq_aln/FANCF-decoy_S1_L001_lowQ.shuffle.bed
fastaOut=/home/gandri/offTargets/Giando/pipeline/FANCF_decoy/results/fastq_aln/FANCF-decoy_S1_L001_lowQ.shuffle.fa

chrSize=/home/gandri/offTargets/Giando/pipeline/annotations/hg38.chrom.sizes
genomeFasta=/storage01/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa

# bam -> bed
${bam2bedCMD} -i ${bam} > ${bedTMP}

# bed -> shuffle bed
${shuffleCMD} -i ${bedTMP} -g ${chrSize} > ${shuffleBedTMP}

# Shuffle Bed -> fasta
${getFastaCMD} -fi ${genomeFasta} -bed ${shuffleBedTMP} -fo ${fastaOut}