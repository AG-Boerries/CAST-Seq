#!/bin/bash
#PBS -S /bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=12
#PBS -l walltime=200:00:00
#PBS -M geoffroy.andrieux@mol-med.uni-freiburg.de
#PBS -m ae
#PBS -N Talen
#PBS -e /home/gandri/offTargets/Giando/pipeline/script/err.log
#PBS -o /home/gandri/offTargets/Giando/pipeline/script/out.log
#PBS -d /home/gandri/offTargets/Giando/pipeline/script
Rscript /home/gandri/offTargets/Giando/pipeline/script/guide_alignmentDual.R

