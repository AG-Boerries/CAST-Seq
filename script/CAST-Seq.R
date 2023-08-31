
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
	make_option(c("--pname"), type="character", default=NULL, 
	            help="name of current project sample", metavar="character"),   
	make_option(c("--sampleDname"), type="character", default=NULL, 
              help="name of sample directory", metavar="character"),   
	make_option(c("--tsamp"), type="character", default=NULL, 
              help="XXX name of TREATED files. XXX_R1_001.fastq.gz AND XXX_R2_001.fastq.gz should exist [default= %default]", metavar="character"),
  make_option(c("--usamp"), type="character", default=NULL, 
              help="XXX name of UNTREATED file. XXX_R1_001.fastq.gz AND XXX_R2_001.fastq.gz should exist [default= %default]", metavar="character"),
	make_option(c("--tnames"), type="character", default=NULL, 
	            help="Label of TREATED files. If NULL, will be t1, t2, etc. [default= %default]", metavar="character"),
	make_option(c("--unames"), type="character", default=NULL, 
	            help="Label of UNTREATED files. If NULL, will be u1, u2, etc. [default= %default]", metavar="character"),
	make_option(c("--homeD"), type="character", default=NULL, 
              help="name of home directory [default= %default]", metavar="character"),
	make_option(c("--scriptD"), type="character", default="script", 
	            help="name of script directory [default= %default]", metavar="character"),
	make_option(c("--tfastqD"), type="character", default=NULL, 
	            help="name of fastq directory [default= %default]", metavar="character"),
	make_option(c("--ufastqD"), type="character", default=NULL, 
	            help="name of fastq directory [default= %default]", metavar="character"),
	make_option(c("--fastqExt1"), type="character", default="_R1_001.fastq.gz", 
	            help="extension of R1 fastq file [default= %default]", metavar="character"),
	make_option(c("--fastqExt2"), type="character", default="_R2_001.fastq.gz", 
	            help="extension of R1 fastq file [default= %default]", metavar="character"),
  make_option(c("--grna"), type="character", default="gRNA.fa", 
              help="name of gRNA fasta file [default= %default]", metavar="character"),
  make_option(c("--grnaR"), type="character", default="gRNA_Right.fa", 
              help="name of gRNA (RIGHT) fasta file (for talen pipeline only)[default= %default]", metavar="character"),
  make_option(c("--grnaL"), type="character", default="gRNA_Left.fa", 
              help="name of gRNA (LEFT) fasta file (for talen pipeline only)[default= %default]", metavar="character"),
  make_option(c("--onTarget"), type="character", default="ots.bed", 
              help="name of ON-target bed file [default= %default]", metavar="character"), 
	make_option(c("--distCov"), type="integer", default=100, 
	            help="distance from the maximum covered bin from where the gRNA will be aligned [default= %default]", metavar="integer"), 
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
    make_option(c("--random"), type="integer", default=1000, 
              help="number of random sequences to generate [default= %default]", metavar="integer"),        
    make_option(c("--width"), type="integer", default=250, 
              help="distance to extend the putative sites [default= %default]", metavar="integer"),       
	make_option(c("--bpOrient"), type="character", default="forward", 
              help="Orientation of the bait primer (forward or reverse) [default= %default]", metavar="character"),            
    make_option(c("--distCutoff"), type="integer", default=1500, 
              help="distance to merge hits together [default= %default]", metavar="integer"),          
    make_option(c("--pvCutoff"), type="numeric", default=0.005, 
              help="pvalue threshold [default= %default]", metavar="numeric"),   
	make_option(c("--scoreCutoff"), type="numeric", default=NULL, 
              help="gRNA alignment score threshold [default= %default]", metavar="numeric"),          
    make_option(c("--hitsCutoff"), type="integer", default=1, 
              help="minimum number of hits per site [default= %default]", metavar="integer"),    
    make_option(c("--saveReads"), type="character", default="no", 
              help="should reads fastq sequences be saved [default= %default]", metavar="character"),         
    make_option(c("--species"), type="character", default="hg", 
              help="name of sample species [default= %default]", metavar="character"),   
   make_option(c("--ovl"), type="integer", default=1, 
              help="number of samples to be considered in the overlap analysis. [default= %default]", metavar="integer"),    
	make_option(c("--signif"), type="integer", default=1, 
	            help="number of significant samples to be considered in the overlap analysis. [default= %default]", metavar="integer"),          
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

if(is.null(opt$pname)){
  print_help(opt_parser)
  stop("At least one argument must be supplied when using overlap or talen overlap pipelines (pname).n", call.=FALSE)
}

if(is.null(opt$sampleDname)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (sampleDname).n", call.=FALSE)
}

if(is.null(opt$tsamp)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (tsamp).n", call.=FALSE)
}

if(is.null(opt$usamp)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (usamp).n", call.=FALSE)
}

if(is.null(opt$homeD)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (homeD).n", call.=FALSE)
}


##########################################################################################
################                   REQUIRED LIBRARIES                 ####################
##########################################################################################

#BiocManager::install(c("openxlsx", "GenomicRanges", "Biostrings", "data.table", "ggplot2", "ggseqlogo", "parallel", "ChIPseeker", "org.Hs.eg.db", "clusterProfiler", "rtracklayer", "biomaRt", "UpSetR", "circlize", "tidyr", "pheatmap", "viridis"))

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
#require(karyoploteR)
require(UpSetR)
require(circlize)
require(tidyr)
require(pheatmap)
require(viridis)

##########################################################################################
################                     PARAMETERS                         ################## 
##########################################################################################

NBCPU <- opt$cpu

#####################################################
# SAMPLE SPECIFIC PARAMETERS (NEED TO BE CHANGED !!!)

# SET REFERENCE FOLDER
homeD <- opt$homeD

# SET SAMPLE PROJECT NAME
projectName <- opt$pname

# SET SAMPLE DIRECTORY NAME
sampleDname <- opt$sampleDname

scriptD <- file.path(homeD, opt$scriptD)

sampleD <- file.path(homeD, "samples_V2", sampleDname)

dataD <- file.path(sampleD, "data")

# SET FASTQ FOLDER
if(is.null(opt$tfastqD)){
  fastqD.treated <- file.path(dataD, "fastq")
}else{
  fastqD.treated <- opt$tfastqD
}

if(is.null(opt$ufastqD)){
  fastqD.untreated <- fastqD.treated
}else{
  fastqD.untreated <- opt$ufastqD
}


# SET INPUT TEST FILE
samples.treated <- opt$tsamp
samples.treated <- strsplit(samples.treated, split = ",|, ")[[1]]

labels.treated <- opt$tnames
if(!is.null(labels.treated)){
  labels.treated <- strsplit(labels.treated, split = ",|, ")[[1]]
  if(length(labels.treated) != length(samples.treated)){
    print("tsamp and tnames have different length; tnames will be modify")
    labels.treated <- paste0("T", seq(1, length(samples.treated)))
  }
}else{
  labels.treated <- paste0("T", seq(1, length(samples.treated)))
}

print(paste0("fastqD.treated: ", fastqD.treated))
fastqD.treated <- strsplit(fastqD.treated, split = ",|, ")[[1]]
if(length(fastqD.treated) == 1){
  print("Use the same fastq directory for all treated")
  fastqD.treated <- rep(fastqD.treated, length(samples.treated))
}else if(length(fastqD.treated) != length(samples.treated)){
  stop("lengths of samples.treated and fastqD.treated do NOT  match")
}


# SET INPUT CONTROL FILE
samples.untreated <- opt$usamp
samples.untreated <- strsplit(samples.untreated, split = ",|, ")[[1]]

# CHECK LENGTH OF samples.treated and samples.untreated
if(length(samples.untreated) != length(samples.treated)){
  if(length(samples.untreated) == 1){
    print("Use the same untreated sample for all")
    samples.untreated <- rep(samples.untreated, length(samples.treated))
  }else{
    stop("lengths of samples.treated and samples.untreated do NOT  match")
  }
}

labels.untreated <- opt$unames
if(!is.null(labels.untreated)){
  labels.untreated <- strsplit(labels.untreated, split = ",|, ")[[1]]
  if(length(labels.untreated) != length(samples.untreated)){
    print("usamp and unames have different length; unames will be modify")
    if(length(labels.untreated) == 1){
      labels.untreated <- rep(labels.untreated, length(samples.untreated))
    }else{
      labels.untreated <- paste0("U", seq(1, length(samples.untreated)))
    }

  }
}else{
  labels.untreated <- paste0("U", seq(1, length(samples.untreated)))
}

print(paste0("fastqD.untreated: ", fastqD.untreated))
fastqD.untreated <- strsplit(fastqD.untreated, split = ",|, ")[[1]]
if(length(fastqD.untreated) == 1){
  print("Use the same fastq directory for all UNtreated")
  fastqD.untreated <- rep(fastqD.untreated, length(samples.treated))
}else if(length(fastqD.untreated) != length(samples.untreated)){
  stop("lengths of samples.treated and fastqD.untreated do NOT  match")
}



##################
# OTHER PARAMETERS


# FASTQ EXTENSION
fastqExt1 <- opt$fastqExt1
fastqExt2 <- opt$fastqExt2


# SET GUIDE SEQ
if(opt$pipeline == "talen"){
	refSeq.left <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grnaL))))
	refSeq.right <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grnaR))))
}else if(opt$pipeline %in% c("crispr2", "double_nickase")){
  refSeq1 <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grnaL))))
  refSeq2 <- toupper(as.character(readDNAStringSet(file.path(dataD, opt$grnaR))))
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
w <- opt$width

# SET DISTANCE CUTOFF (FOR HITS)
distance.cutoff <- opt$distCutoff

# BAIT PRIMER ORIENTATION
bpOrient <- opt$bpOrient

# SET COVERAGE DISTANCE CUTOFF
distance.cov <- opt$distCov

# SET PVALUE CUTOFF
#filtName <- ""
pv.cutoff <- opt$pvCutoff
score.cutoff <- opt$scoreCutoff
hits.cutoff <- opt$hitsCutoff

# SAVE READS IN FASTA FORMAT
if(opt$saveReads == "yes"){
	saveReads <- TRUE
}else saveReads <- FALSE


# OVERLAP PARAMETERS
#ovlName <- opt$ovlName

#repD <- opt$repDname
#if(!is.null(repD)){
#	repD.split <- strsplit(repD, split = ",")[[1]]
#	repD.split <- file.path(homeD, "samples_V2", repD.split)
#}

nb.ovl <- opt$ovl
nb.signif <- opt$signif

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
	#histoneFiles <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = TRUE)
	#names(histoneFiles) <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = FALSE)
	#names(histoneFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histoneFiles))
}else if(opt$species == "hg_PSMAeTRuCKI"){
  # Human
  require(org.Hs.eg.db)
  require(BSgenome.hg38.PSMAeTRuCKI)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ORG <- org.Hs.eg.db
  ORG.STR <- "org.Hs.eg.db"
  GNM <- BSgenome.hg38.PSMAeTRuCKI::hg38PSMAeTRuCKI
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/human_PSMAeTRuCKI")
  
  # SET GENOME VERSION
  myGenome.size <- file.path(annotD, "chrom.sizes")
  
  # HG38 TSS TES
  geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)
  
  # SET ONCO ENTREZ LIST
  oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
  onco.width <- 3000
  
  # SET CIRCOS PLOT SPECIES
  circos.sp <- "hg38"
  
}else if(opt$species == "hg_CD247KI"){
  # Human
  require(org.Hs.eg.db)
  require(BSgenome.hg38.CD247KI)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ORG <- org.Hs.eg.db
  ORG.STR <- "org.Hs.eg.db"
  GNM <- BSgenome.hg38.CD247KI::hg38CD247KI
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/human_CD247KI")
  
  # SET GENOME VERSION
  myGenome.size <- file.path(annotD, "chrom.sizes")
  
  # HG38 TSS TES
  geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)
  
  # SET ONCO ENTREZ LIST
  oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
  onco.width <- 3000
  
  # SET CIRCOS PLOT SPECIES
  circos.sp <- "hg38"
  
}else if(opt$species == "hg_chr11N5254785to5254984"){
  # Human
  require(org.Hs.eg.db)
  require(BSgenome.hg38.chr11N5254785to5254984)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ORG <- org.Hs.eg.db
  ORG.STR <- "org.Hs.eg.db"
  GNM <- BSgenome.hg38.chr11N5254785to5254984::hg38chr11N5254785to5254984
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/human_chr11N5254785to5254984")
  
  # SET GENOME VERSION
  myGenome.size <- file.path(annotD, "chrom.sizes")
  
  # HG38 TSS TES
  geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)
  
  # SET ONCO ENTREZ LIST
  oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
  onco.width <- 3000
  
  # SET CIRCOS PLOT SPECIES
  circos.sp <- "hg38"
  
}else if(opt$species == "hg_NCF1"){
  # Human
  require(org.Hs.eg.db)
  require(BSgenome.hg38.NCF1)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene
  ORG <- org.Hs.eg.db
  ORG.STR <- "org.Hs.eg.db"
  GNM <- BSgenome.hg38.NCF1::hg38NCF1
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/hg_NCF1")
  
  # SET GENOME VERSION
  myGenome.size <- file.path(annotD, "chrom.sizes")
  
  # HG38 TSS TES
  geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)
  
  # SET ONCO ENTREZ LIST
  oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
  onco.width <- 3000
  
  # SET CIRCOS PLOT SPECIES
  circos.sp <- "hg38"
  
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
}else if(opt$species == "mm_PCSK9"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10PCSK9)
  GNM <- BSgenome.mm10PCSK9::mm10PCSK9
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus_PCSK9")
  
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
}else if(opt$species == "mm_AdV"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10AdV)
  GNM <- BSgenome.mm10AdV::mm10AdV
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus_AdV")
  
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
}else if(opt$species == "mm_AAVg12"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10AAVg12)
  GNM <- BSgenome.mm10AAVg12::mm10AAVg12
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus_AAVg12")
  
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
}else if(opt$species == "mm_AAVAAT"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10AAVAAT)
  GNM <- BSgenome.mm10AAVAAT::mm10AAVAAT
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus_AAV_AAT")
  
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
}else if(opt$species == "mm_AAVTBG"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10AAVTBG)
  GNM <- BSgenome.mm10AAVTBG::mm10AAVTBG
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus_AAV_TBG")
  
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
}else if(opt$species == "mm_PCSK9gMH"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10PCSK9gMH)
  GNM <- BSgenome.mm10PCSK9gMH::mm10PCSK9gMH
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/Mus_musculus_PCSK9_gMH")
  
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
}else if(opt$species == "mm10HITIAAV"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10.HITIAAV)
  GNM <- BSgenome.mm10.HITIAAV::mm10HITIAAV
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/mm10HITIAAV")
  
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
}else if(opt$species == "mm10HITIAAVwogRNA"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10.HITIAAVwogRNA)
  GNM <- BSgenome.mm10.HITIAAVwogRNA::mm10HITIAAVwogRNA
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/mm10HITIAAVwogRNA")
  
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
}else if(opt$species == "mm10MK1"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10.MK1)
  GNM <- BSgenome.mm10.MK1::mm10MK1
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/mm10MK1")
  
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
}else if(opt$species == "mm10MK2"){
  #print_help(opt_parser)
  #stop("mmul not yet compatible", call.=FALSE)
  
  # Mouse
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(TxDb.Mmusculus.UCSC.mm10.knownGene)
  TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ORG <- org.Mm.eg.db
  ORG.STR <- "org.Mm.eg.db"
  
  require(BSgenome.mm10.MK2)
  GNM <- BSgenome.mm10.MK2::mm10MK2
  
  # ANNOTATION DIRECTORY
  annotD <- file.path(homeD, "annotations/mm10MK2")
  
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
source(file.path(scriptD, "pipeline_crispr_V2.R"))
#source(file.path(scriptD, "pipeline_2OT.R"))
#source(file.path(scriptD, "pipeline_overlap.R"))
source(file.path(scriptD, "circlize.R"))
source(file.path(scriptD, "site_overlap.R"))
source(file.path(scriptD, "ONreadout.R"))
source(file.path(scriptD, "barcodeHoppingFilter.R"))
source(file.path(scriptD, "read_length_distribution.R"))

if(opt$pipeline == "talen"){
	source(file.path(scriptD, "pipeline_TALEN_V2.R"))
	#source(file.path(scriptD, "guide_alignmentDual.R"))
  source(file.path(scriptD, "guide_alignmentDual_SPLIT.R"))
}

if(opt$pipeline == "double_nickase"){
  source(file.path(scriptD, "pipeline_double_nickase_V2.R"))
  #source(file.path(scriptD, "guide_alignmentDual.R"))
}

if(opt$pipeline == "crispr2"){
  source(file.path(scriptD, "pipeline_crispr2_V2.R"))
}


################     SET PYTHON FOR LCS    ################ 
require(reticulate)
#use_python("/home/gandri/miniconda3/bin/python")
#use_python("/usr/bin/python")
print(paste0("pythonPath: ", opt$pythonPath))
#use_python(opt$pythonPath)
use_condaenv("pyCAST", required = TRUE)

source_python(file.path(scriptD, "lcs.py"))

##########################################################################################
############                        CHECK INPUT FILES                         ############
##########################################################################################

print("CHECK INPUT FILES")

if(file.exists(otsBed)){
	print("otsBed: OK")
}else{
	stop(paste0("otsBed: ", otsBed, " does NOT exist."))
}

if(file.exists(myGenome.size)){
	print("myGenome.size: OK")
}else{
	stop(paste0("myGenome.size: ", myGenome.size, " does NOT exist."))
}

if(file.exists(oncoEntrez)){
	print("oncoEntrez: OK")
}else{
	stop(paste0("oncoEntrez: ", oncoEntrez, " does NOT exist."))
}

posFile=file.path(dataD, "pos.fa")
if(file.exists(posFile)){
	print("posFile: OK")
}else{
	stop(paste0("posFile: ", posFile, " does NOT exist."))
}

negFile=file.path(dataD, "neg.fa")
if(file.exists(negFile)){
	print("negFile: OK")
}else{
	print("negFile: does NOT exist. XXXXX neg.fa has been created.")
	neg <- c("neg" = "XXXXX")
	writeXStringSet(BStringSet(neg), negFile)
}

misprimingFile=file.path(dataD, "mispriming.fa")
if(file.exists(misprimingFile)){
	print("misprimingFile: OK")
}else{
	print("misprimingFile: does NOT exist. XXXXX mispriming.fa has been created.")
	msp <- c("msp" = "XXXXX")
	writeXStringSet(BStringSet(msp), misprimingFile)
}

hthFile=file.path(dataD, "headTOhead.fa")
if(file.exists(hthFile)){
	print("hthFile: OK")
}else{
	print("hthFile: does NOT exist. XXXXX headTOhead.fa has been created.")
	hth <- c("hth" = "XXXXX")
	writeXStringSet(BStringSet(hth), hthFile)
}

linkerFile=file.path(dataD, "linker_RC.fa")
if(file.exists(linkerFile)){
	print("linkerFile: OK")
}else{
	stop(paste0("linkerFile: ", linkerFile, " does NOT exist."))
}

adaptersFile=file.path(annotD, "TruSeq4-PE.fa")
if(file.exists(adaptersFile)){
	print("adaptersFile: OK")
}else{
	stop(paste0("adaptersFile: ", adaptersFile, " does NOT exist."))
}


# CHECK FASTQ FILES

print(samples.treated)
lapply(1:length(samples.treated), function(i){
  
  pair1=file.path(fastqD.treated[i], paste0(samples.treated[i], fastqExt1))
  if(file.exists(pair1)){
    print("pair1: OK")
  }else{
    stop(paste0("pair1: ", pair1, " does NOT exist."))
  }
  
  pair2=file.path(fastqD.treated[i], paste0(samples.treated[i], fastqExt2))
  if(file.exists(pair2)){
    print("pair2: OK")
  }else{
    stop(paste0("pair2: ", pair2, " does NOT exist."))
  }
  
  return(NA)
})

print(samples.untreated)
lapply(1:length(samples.untreated), function(i){
  
  pair1=file.path(fastqD.untreated[i], paste0(samples.untreated[i], fastqExt1))
  if(file.exists(pair1)){
    print("pair1: OK")
  }else{
    stop(paste0("pair1: ", pair1, " does NOT exist."))
  }
  
  pair2=file.path(fastqD.untreated[i], paste0(samples.untreated[i], fastqExt2))
  if(file.exists(pair2)){
    print("pair2: OK")
  }else{
    stop(paste0("pair2: ", pair2, " does NOT exist."))
  }
  
  return(NA)
})


##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################

if(opt$pipeline == "crispr"){
  print("START CRISPR PIPELINE")
  runPipelineCRISPR()
}else if(opt$pipeline == "crispr_2ot"){
  print("START CRISPR 2 ON-TARGETS PIPELINE")
  runPipeline_2OT()
}else if(opt$pipeline == "talen"){
  print("START TALEN PIPELINE")
  runPipelineTALEN()
}else if(opt$pipeline == "double_nickase"){
  print("START DOUBLE NICKASE PIPELINE")
  runPipelineDOUBLENICKASE()
}else if(opt$pipeline == "crispr2"){
  print("START CRISPR 2 gRNAs PIPELINE")
  runPipeline_CRISPR2()
}else{
  print_help(opt_parser)
  stop("wrong pipeline name. Must be crispr, crispr2, crispr_2ot, talen, double_nickase", call.=FALSE)
}



end_time <- Sys.time()
end_time - start_time
