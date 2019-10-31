
set.seed(1234)

NBCPU = 12

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################




############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

##########################################################################################
################                   REQUIRED LIBRARIES                 ####################
##########################################################################################
 
require(openxlsx)
require(GenomicRanges)
require(Biostrings)
require(data.table)
require(ggplot2)
require(ggseqlogo)
#require(textreadr)
require(parallel)
require(ChIPseeker)
require(org.Hs.eg.db)
require(clusterProfiler)
require(rtracklayer)
require(biomaRt)
require(tools)
require(BSgenome.Hsapiens.UCSC.hg38)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(karyoploteR)


##########################################################################################
################                     PARAMETERS                         ################## 
##########################################################################################

#####################################################
# SAMPLE SPECIFIC PARAMETERS (NEED TO BE CHANGED !!!)

# SET SAMPLE DIRECTORY NAME
sampleDname <- "G3_TOY"

# SET INPUT TEST FILE
sampleName <- "G3_treated"

# SET INPUT CONTROL FILE
controlName <- "G3_untreated"

# SET REFERENCE FOLDER
homeD <- "./"

##################
# OTHER PARAMETERS

scriptD <- file.path(homeD, "script")
annotD <- file.path(homeD, "annotations")

sampleD <- file.path(homeD, "samples", sampleDname)

dataD <- file.path(sampleD, "data")
resultD <- file.path(sampleD, "results", "guide_aln")
dir.create(resultD, showWarnings = FALSE)

# SET GUIDE SEQ
refSeq <- toupper(as.character(readDNAStringSet(file.path(dataD, "gRNA.fa"))))

# SET ON-TARGET SITE
otsBed <- file.path(dataD, "ots.bed")
otsDistance <- 50
surrounding_size <- 20000

# SET FLANKING REGIONS (+/- flanking size)
flankingSize <- 2500

# SET NUMBER OF RANDOM SEQUENCES
nb.rd <- 10000

# SET WIDTH
w = 250

# SET DISTANCE CUTOFF (FOR CLUSTERS)
distance.cutoff <- 1500

# SET PVALUE CUTOFF
pv.cutoff <- 0.05

# SET GENOME VERSION
myGenome.size <- file.path(annotD, "hg38.chrom.sizes")

# HG38 TSS TES
geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)

# SET ONCO ENTREZ LIST
oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
onco.width <- 3000

# SET HISTONE BROAD PEAKS FILES
histoneFiles <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = TRUE)
names(histoneFiles) <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = FALSE)
names(histoneFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histoneFiles))

##########################################################################################
################                     SOURCE SCRIPTS                  #####################
##########################################################################################

source(file.path(scriptD, "fastq_aln.R"))
source(file.path(scriptD, "bed2sequence.R"))
source(file.path(scriptD, "bedTools_fct.R"))
source(file.path(scriptD, "delta.R"))
source(file.path(scriptD, "cluster.R"))
source(file.path(scriptD, "checkMAPQ.R"))
source(file.path(scriptD, "hitsEnrichment.R"))
source(file.path(scriptD, "guide_alignment.R"))
source(file.path(scriptD, "flankingRegions.R"))
source(file.path(scriptD, "defineGroups.R"))
source(file.path(scriptD, "pieChart.R"))
source(file.path(scriptD, "annotateGenes.R"))
source(file.path(scriptD, "vs_allGenes.R"))
source(file.path(scriptD, "vs_onco.R"))
source(file.path(scriptD, "histone_forestPlot.R"))
source(file.path(scriptD, "scoring_system.R"))
source(file.path(scriptD, "chrPlot.R"))
source(file.path(scriptD, "ogi.R"))
source(file.path(scriptD, "pipeline.R"))

################     SET PYTHON FOR LCS    ################ 
require(reticulate)
use_python("/home/gandri/miniconda3/bin/python")
#use_python("/usr/bin/python")
source_python(file.path(scriptD, "lcs.py"))

##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################

runPipeline()

################     RUN PIPELINE OVERLAP   #############################################
#source(file.path(scriptD, "pipeline_overlap.R"))

#ovlD <- file.path(sampleD, "results", "overlap_aln")
#sampleName <- "G3_WT_ovl"

#runPipelineOverlap()