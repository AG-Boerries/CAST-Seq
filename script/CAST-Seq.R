
set.seed(1234)

start_time <- Sys.time()


#####################################################
###                                               ###
###               CAST-SEQ PIPELINE               ###
###                                               ###
#####################################################




##########################################################################################
################                     PARSE                              ################## 
##########################################################################################
require(optparse)

option_list = list(
	make_option(c("--pipeline"), type="character", default=NULL, 
              help="name of the pipeline you want to use. Choose between crispr, crispr_2ot or talen", metavar="character"),
	make_option(c("--sampleDname"), type="character", default=NULL, 
              help="name of sample directory", metavar="character"),   
	make_option(c("--sampleName"), type="character", default=NULL, 
              help="XXX name of test file. XXX_R1_001.fastq.gz AND XXX_R2_001.fastq.gz should exist [default= %default]", metavar="character"),
    make_option(c("--controlName"), type="character", default=NULL, 
              help="XXX name of control file. XXX_R1_001.fastq.gz AND XXX_R2_001.fastq.gz should exist [default= %default]", metavar="character"),
    make_option(c("--homeD"), type="character", default=NULL, 
              help="name of home directory [default= %default]", metavar="character"),
	make_option(c("--fastqD"), type="character", default=NULL, 
	            help="name of fastq directory [default= %default]", metavar="character"),
    make_option(c("--grna"), type="character", default="gRNA.fa", 
              help="name of gRNA fasta file [default= %default]", metavar="character"),
    make_option(c("--grnaR"), type="character", default="gRNA_Right.fa", 
              help="name of gRNA (RIGHT) fasta file (for talen pipeline only)[default= %default]", metavar="character"),
    make_option(c("--grnaL"), type="character", default="gRNA_Left.fa", 
              help="name of gRNA (LEFT) fasta file (for talen pipeline only)[default= %default]", metavar="character"),
    make_option(c("--onTarget"), type="character", default="ots.bed", 
              help="name of ON-target bed file [default= %default]", metavar="character"),              
	make_option(c("--otsDistance"), type="integer", default=50, 
              help="distance (bp) from the ON-target. Reads +/- this distance will be removed [default= %default]", metavar="integer"),     
	make_option(c("--surrounding_size"), type="integer", default=20000, 
              help="distance (bp) from the ON-target. Use for the scoring system [default= %default]", metavar="integer"),                 
    make_option(c("--flank1"), type="character", default="flank1.fa", 
              help="name of first flanking sequence [default= %default]", metavar="character"),         
    make_option(c("--flank2"), type="character", default="flank2.fa", 
              help="name of second flanking sequence [default= %default]", metavar="character"), 
    make_option(c("--flankingSize"), type="integer", default=2500, 
              help="distance to consider for HMT [default= %default]", metavar="integer"),         
    make_option(c("--random"), type="integer", default=10000, 
              help="number of random sequences to generate [default= %default]", metavar="integer"),        
    make_option(c("--width"), type="integer", default=250, 
              help="distance to extend the putative sites [default= %default]", metavar="integer"),            
    make_option(c("--distCutoff"), type="integer", default=1500, 
              help="distance to merge hits together [default= %default]", metavar="integer"),          
    make_option(c("--pvCutoff"), type="numeric", default=0.05, 
              help="pvalue threshold [default= %default]", metavar="numeric"),   
	make_option(c("--scoreCutoff"), type="numeric", default=NULL, 
              help="gRNA alignment score threshold [default= %default]", metavar="numeric"),          
    make_option(c("--hitsCutoff"), type="integer", default=1, 
              help="minimum number of hits per site [default= %default]", metavar="integer"),    
    make_option(c("--saveReads"), type="character", default="no", 
              help="should reads fastq sequences be saved [default= %default]", metavar="character"),         
    make_option(c("--species"), type="character", default="hg", 
              help="name of sample species [default= %default]", metavar="character"),   
    make_option(c("--ovlDname"), type="character", default=NULL, 
              help="name of overlap directory within sample directory (only for crispr_overlap or talen_overlap pipelines)", metavar="character"),       
    make_option(c("--ovlName"), type="character", default=NULL, 
              help="name of overlap sample within overlap directory (only for crispr_overlap or talen_overlap pipelines) Must be XXX if the files is called XXX.xlsx", metavar="character"),   
    make_option(c("--replicates"), type="character", default=NULL, 
              help="name of sample to be used in the overlap analysis. [default= %default]", metavar="character"),            
    make_option(c("--repNames"), type="character", default=NULL, 
              help="labels of the replicates to be used in the overlap analysis. [default= %default]", metavar="character"),          
    make_option(c("--repDname"), type="character", default=NULL, 
              help="name of a representative replicate (used to find the appropriate replicate files). [default= %default]", metavar="character"),       
    make_option(c("--ovl"), type="integer", default=1, 
              help="number of significant samples to be considered in the overlap analysis. [default= %default]", metavar="character"),          
	make_option(c("--cpu"), type="integer", default=2, 
              help="number of CPUs [default= %default]", metavar="integer"),            
    make_option(c("--pythonPath"), type="character", default="/Users/geoffroyandrieux/miniconda3/bin/python", 
              help="python path [default= %default]", metavar="character")             
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$pipeline)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (pipeline).n", call.=FALSE)
}

if(is.null(opt$sampleDname)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (sampleDname).n", call.=FALSE)
}

if(is.null(opt$sampleName) & !(opt$pipeline %in% c("crispr_overlap", "talen_overlap", "double_nickase_overlap"))){
  print_help(opt_parser)
  stop("At least one argument must be supplied (sampleName).n", call.=FALSE)
}

if(is.null(opt$controlName) & !(opt$pipeline %in% c("crispr_overlap", "talen_overlap", "double_nickase_overlap"))){
  print_help(opt_parser)
  stop("At least one argument must be supplied (controlName).n", call.=FALSE)
}

if(is.null(opt$homeD)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (homeD).n", call.=FALSE)
}

# Overlap pipeline parameters
if(is.null(opt$ovlDname) & opt$pipeline %in% c("crispr_overlap", "talen_overlap", "double_nickase_overlap")){
  print_help(opt_parser)
  stop("At least one argument must be supplied when using overlap or talen overlap pipelines (ovlDname).n", call.=FALSE)
}

if(is.null(opt$ovlName) & opt$pipeline %in% c("crispr_overlap", "talen_overlap", "double_nickase_overlap")){
  print_help(opt_parser)
  stop("At least one argument must be supplied when using overlap or talen overlap pipelines (ovlName).n", call.=FALSE)
}


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
require(karyoploteR)
require(UpSetR)
require(circlize)
library(tidyr)

##########################################################################################
################                     PARAMETERS                         ################## 
##########################################################################################

NBCPU <- opt$cpu

#####################################################
# SAMPLE SPECIFIC PARAMETERS (NEED TO BE CHANGED !!!)

# SET SAMPLE DIRECTORY NAME
sampleDname <- opt$sampleDname

# SET INPUT TEST FILE
sampleName <- opt$sampleName

# SET INPUT CONTROL FILE
controlName <- opt$controlName

# SET REFERENCE FOLDER
homeD <- opt$homeD


##################
# OTHER PARAMETERS

scriptD <- file.path(homeD, "script")

sampleD <- file.path(homeD, "samples", sampleDname)

dataD <- file.path(sampleD, "data")
resultD <- file.path(sampleD, "results", "guide_aln")
dir.create(resultD, showWarnings = FALSE)

# SET FASTQ FOLDER
if(is.null(opt$fastqD)){
  fastqD <- file.path(dataD, "fastq")
}else{
  fastqD <- opt$fastqD
}



# SET GUIDE SEQ

if(opt$pipeline == "talen" | opt$pipeline == "talen_overlap" | opt$pipeline == "double_nickase" | opt$pipeline == "double_nickase_overlap"){
	refSeq.left <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grnaL))))
	refSeq.right <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grnaR))))
}else{
	refSeq <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grna))))
}

# SET ON-TARGET SITE
otsBed <- file.path(dataD, opt$onTarget)
otsDistance <- opt$otsDistance
surrounding_size <- opt$surrounding_size

# SET FLANKING REGIONS (+/- flanking size)
if(file.exists(file.path(dataD, opt$flank1))){
	flank1.sq <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$flank1))))
}else{
	flank1.sq <- NULL
}
	
if(file.exists(file.path(dataD, opt$flank2))){
	flank2.sq <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$flank2))))
}else{
	flank2.sq <- NULL
}

flankingSize <- opt$flankingSize

# SET NUMBER OF RANDOM SEQUENCES
nb.rd <- opt$random

# SET WIDTH TO ADD AROUND PUTATIVE SITES
w = opt$width

# SET DISTANCE CUTOFF (FOR HITS)
distance.cutoff <- opt$distCutoff

# SET PVALUE CUTOFF
#filtName <- ""
pv.cutoff <- opt$pvCutoff
score.cutoff <- opt$scoreCutoff
hits.cutoff <- opt$hitsCutoff

# SAVE READS IN FASTA FORMAT
if(opt$saveReads == "yes"){
	saveReads <- TRUE
}else saveReads <- FALSE


# OVERLAP PIPELINE PARAMETERS
ovlD <- file.path(homeD, "samples_overlap", opt$ovlDname)
ovlName <- opt$ovlName
replicates <- opt$replicates
repNames <- opt$repNames
repD <- file.path(homeD, "samples", opt$repDname)
nb.ovl <- opt$ovl

# REMOVE TMP FILES
rmTMP <- TRUE

# COMPRESS RESULTS
tozip <- TRUE


##########################################################################################
################                  SPECIES PARAMETERS                    ################## 
##########################################################################################


if(opt$species == "hg"){
	# Human
	require(org.Hs.eg.db)
	require(BSgenome.Hsapiens.UCSC.hg38)
	require(TxDb.Hsapiens.UCSC.hg38.knownGene)
	TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene
	ORG <- org.Hs.eg.db
	ORG.STR <- "org.Hs.eg.db"
	GNM <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
	
	# ANNOTATION DIRECTORY
	annotD <- file.path(homeD, "annotations/human")

	# SET GENOME VERSION
	myGenome.size <- file.path(annotD, "chrom.sizes")

	# HG38 TSS TES
	geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)

	# SET ONCO ENTREZ LIST
	oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
	onco.width <- 3000
	
	# SET CIRCOS PLOT SPECIES
	circos.sp <- "hg38"

	# SET HISTONE BROAD PEAKS FILES
	histoneFiles <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = TRUE)
	names(histoneFiles) <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = FALSE)
	names(histoneFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histoneFiles))
}else if(opt$species == "mmul"){
	#print_help(opt_parser)
	#stop("mmul not yet compatible", call.=FALSE)
	
	# Macaca
	require(org.Mmu.eg.db)
	require(BSgenome.Mmulatta.UCSC.rheMac8)
	require(TxDb.Mmulatta.UCSC.rheMac8.refGene)
	TXDB <- TxDb.Mmulatta.UCSC.rheMac8.refGene
	ORG <- org.Mmu.eg.db
	ORG.STR <- "org.Mmu.eg.db"
	GNM <- BSgenome.Mmulatta.UCSC.rheMac8::Mmulatta
	
	# ANNOTATION DIRECTORY
	annotD <- file.path(homeD, "annotations/Macaca_mulatta")

	# SET GENOME VERSION
	myGenome.size <- file.path(annotD, "chrom.sizes")

	# HG38 TSS TES
	geneMat <- read.delim(file.path(annotD, "Mmul8.0.1_TSS_TES.txt"), header = FALSE)

	# SET ONCO ENTREZ LIST
	oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
	onco.width <- 3000
	
	# SET CIRCOS PLOT SPECIES
	circos.sp <- ""

	# SET HISTONE BROAD PEAKS FILES
	#histoneFiles <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = TRUE)
	#names(histoneFiles) <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = FALSE)
	#names(histoneFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histoneFiles))
}else if(opt$species == "mm"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Macaca
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  GNM <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus")
  
  # SET GENOME VERSION
  myGenome.size <- file.path(annotD, "chrom.sizes")
  
  # HG38 TSS TES
  geneMat <- read.delim(file.path(annotD, "Mm10_TSS_TES.txt"), header = FALSE)
  
  # SET ONCO ENTREZ LIST
  oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
  onco.width <- 3000
  
  # SET CIRCOS PLOT SPECIES
  circos.sp <- "mm10"
  
  # SET HISTONE BROAD PEAKS FILES
  #histoneFiles <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = TRUE)
  #names(histoneFiles) <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = FALSE)
  #names(histoneFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histoneFiles))
}else{
	print_help(opt_parser)
	stop("wrong species (must be hg).n", call.=FALSE)
}


##########################################################################################
################                     SOURCE SCRIPTS                  #####################
##########################################################################################

source(file.path(scriptD, "fastq_aln.R"))
source(file.path(scriptD, "bed2sequence.R"))
source(file.path(scriptD, "bedTools_fct.R"))
source(file.path(scriptD, "countReads.R"))
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
#source(file.path(scriptD, "vs_onco.R"))
#source(file.path(scriptD, "histone_forestPlot.R"))
source(file.path(scriptD, "scoring_system.R"))
source(file.path(scriptD, "chrPlot.R"))
source(file.path(scriptD, "ogi.R"))
source(file.path(scriptD, "getReads.R"))
source(file.path(scriptD, "pipeline.R"))
#source(file.path(scriptD, "pipeline_2OT.R"))
source(file.path(scriptD, "pipeline_overlap.R"))
source(file.path(scriptD, "circlize.R"))
source(file.path(scriptD, "site_overlap.R"))
source(file.path(scriptD, "ONreadout.R"))

if(opt$pipeline == "talen"){
	source(file.path(scriptD, "pipeline_TALEN.R"))
	source(file.path(scriptD, "guide_alignmentDual.R"))
}

if(opt$pipeline == "talen_overlap"){
  source(file.path(scriptD, "pipeline_TALENoverlap.R"))
  source(file.path(scriptD, "guide_alignmentDual.R"))
}

if(opt$pipeline == "double_nickase"){
  source(file.path(scriptD, "pipelineDoubleNickase.R"))
  source(file.path(scriptD, "guide_alignmentDual.R"))
}

if(opt$pipeline == "double_nickase_overlap"){
  source(file.path(scriptD, "pipelineDoubleNickase_overlap.R"))
  source(file.path(scriptD, "guide_alignmentDual.R"))
}

################     SET PYTHON FOR LCS    ################ 
require(reticulate)
#use_python("/home/gandri/miniconda3/bin/python")
#use_python("/usr/bin/python")
print(paste0("pythonPath: ", opt$pythonPath))
use_python(opt$pythonPath)
source_python(file.path(scriptD, "lcs.py"))

##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################

if(opt$pipeline == "crispr"){
	print("START CRISPR PIPELINE")
	runPipeline()
}else if(opt$pipeline == "crispr_2ot"){
	print("START CRISPR 2 ON-TARGETS PIPELINE")
	runPipeline_2OT()
}else if(opt$pipeline == "crispr_overlap"){
	print("START CRISPR OVERLAP PIPELINE")
	runPipelineOverlap()
}else if(opt$pipeline == "talen"){
	print("START TALEN PIPELINE")
	runPipelineTALEN()
}else if(opt$pipeline == "talen_overlap"){
	print("START TALEN OVERLAP PIPELINE")
	runPipelineTALENoverlap()
}else if(opt$pipeline == "double_nickase"){
  print("START DOUBLE NICKASE PIPELINE")
  runPipelineDoubleNickase()
}else if(opt$pipeline == "double_nickase_overlap"){
  print("START DOUBLE NICKASE OVERLAP PIPELINE")
  runPipelineDoubleNickaseOverlap()
}else{
  print_help(opt_parser)
  stop("wrong pipeline name. Must be crispr, crispr_2ot, crispr_overlap, talen or talen_overlap, double_nickase or double_nickase_overlap", call.=FALSE)
}


end_time <- Sys.time()
end_time - start_time
