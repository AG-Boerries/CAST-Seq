
runPipelineTALEN <- function()
{
##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################


if(!file.exists(file.path(resultD, paste0(sampleName, "_Alignment.bed")))){
################     FASTQ ALIGNMENT    ################ 
print("################     FASTQ ALIGNMENT    ################")

fastqAln(functionstring=paste0("sh ", file.path(scriptD,"fastq_aln.sh")),
	homeFolder = homeD,
	annotFolder = annotD,
	sampleFolder = sampleDname,
	fastqFolder = fastqD,
	testSample = sampleName,
	utSample = controlName,
	cpu = NBCPU
	)	
	
print("################     COUNT READS    ################")
countReads(fastqFolder = fastqD,
           sampleFolder = sampleD,
           testSample = sampleName,
           utSample = controlName)

}

################     DELTA AND HITS    ################ 
print("################     DELTA AND HITS    ################ ")
# bed files (test and control) should be available in resultD

# test sample
sampleBed <- file.path(resultD, paste0(sampleName, "_Alignment.bed"))
		   
# control sample		   
controlBed <- file.path(resultD, paste0(controlName, "_Alignment.bed"))	
	   
# run delta analysis
getDelta(sampleBed, otsBed, otsDistance, distance = distance.cutoff)		   
getDeltaShuffle(gsub(".bed$", "_delta.bed", sampleBed), nb = 10, distance = distance.cutoff, genome.size = myGenome.size)				   
getDelta(controlBed, otsBed, otsDistance, distance = distance.cutoff)		   
getDeltaShuffle(gsub(".bed$", "_delta.bed", controlBed), nb = 10, distance = distance.cutoff, genome.size = myGenome.size)		

# Density without on-target sites
#deltaDensity(gsub(".bed$", "_delta.bed", sampleBed), otsBed, surrounding_size)
#deltaDensity(gsub(".bed$", "_delta.bed", controlBed), otsBed, surrounding_size)
	  	  
# run hits analysis
getHits(gsub(".bed$", "_delta.bed", sampleBed), distance.cutoff)
getHits(gsub(".bed$", "_delta.bed", controlBed), distance.cutoff)

################     TEST VS. CONTROL ENRICHMENT    ################
print("################     TEST VS. CONTROL ENRICHMENT    ################")
# USE RAW FASTQ FILE !!!!!!!!!
#nbReads.sample <- as.numeric(nbReadFastqgz(file.path(dataD, "fastq", paste0(sampleName, "_R2_001.fastq.gz"))))
#nbReads.control <- as.numeric(nbReadFastqgz(file.path(dataD, "fastq", paste0(controlName, "_R2_001.fastq.gz"))))

nbReads.sample <- as.numeric(nbReadFastqgz(file.path(fastqD, paste0(sampleName, "_R2_001.fastq.gz"))))
nbReads.control <- as.numeric(nbReadFastqgz(file.path(fastqD, paste0(controlName, "_R2_001.fastq.gz"))))

doEnrichment(gsub(".bed$", "_hits.bed", sampleBed), gsub(".bed$", "_hits.bed", controlBed),
			 nbReads.sample, nbReads.control, w, myGenome.size)
	
doEnrichment(gsub(".bed$", "_hits.bed", controlBed), gsub(".bed$", "_hits.bed", sampleBed),
			 nbReads.control, nbReads.sample, w, myGenome.size)
	
# ADD AVERAGE MAPQ
addMAPQ(file.path(resultD, paste0(sampleName, "_w", w, ".xlsx")),
	file.path(sampleD, "results", "fastq_aln", paste0(sampleName, "_AlignmentSort.bam")))
addMAPQ(file.path(resultD, paste0(controlName, "_w", w, ".xlsx")),
	file.path(sampleD, "results", "fastq_aln", paste0(controlName, "_AlignmentSort.bam")))

##########################################################################################
############                DESIGNER NUCLEASE TREATED SAMPLE                  ############
##########################################################################################

################     GUIDE SEQ ALIGNMENT    ################ 
print("################     GUIDE SEQ ALIGNMENT    ################")

# DO GUIDE ALIGNMENT ON REAL SEQUENCES
guideD <- resultD
getGuideAlignmentDual(inputF = file.path(resultD, paste0(sampleName, "_w", w, ".xlsx")),
				  guideLeft = refSeq.left, guideRight = refSeq.right,
				  alnFolder = guideD,
				  gnm = GNM
				  )
#file.remove(list.files(guideD, pattern = "_TMP.txt", full.names = TRUE))					  
				  
# GENERATE RANDOM SEQUENCE BED
randomD <- file.path(sampleD, "results", "random")
dir.create(randomD, showWarnings = FALSE)
randomName <- paste0("random_w",w)

# Check if bed file already exists, if yes skip the analysis
if(!file.exists(file.path(randomD, paste0(randomName, ".bed")))){
	getRandomBed(l = w*2, n = nb.rd,
		outFile = file.path(randomD, paste0(randomName, ".bed")),
		opt.string = paste("-g", myGenome.size, sep = " ")
		)

	# DO GUIDE ALIGNMENT ON RANDOM SEQUENCES
	getGuideAlignmentDual(inputF = file.path(randomD, paste0(randomName, ".bed")),
					  guideLeft = refSeq.left, guideRight = refSeq.right,
					  alnFolder = randomD,
					  gnm = GNM
					  )	
	#file.remove(list.files(randomD, pattern = "_TMP.txt", full.names = TRUE))
	}

# ASSIGN PVALUE FOR EVERY SINGLE COMBINATION
assignPV(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
	file.path(randomD, paste0(randomName, "_aln_stat.xlsx")))
	
# ASSIGN THE BEST COMBINATION
assignBestCombination(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")))	

# filt name
filtName <- ""
if(!is.null(hits.cutoff)) filtName <- c(filtName, "hits", hits.cutoff)
if(!is.null(pv.cutoff)) filtName <- c(filtName, "pv", pv.cutoff)
if(!is.null(score.cutoff)) filtName <- c(filtName, "score", score.cutoff)
filtName <- paste0(filtName, collapse = "_")	


################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
addFlanking(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")), otsBed, flankingSize)

# RANDOM SEQUENCES
# Check if output file already exists, if yes skip the analysis
if(!file.exists(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")))){
	addFlanking(file.path(randomD, paste0(randomName, "_aln_stat.xlsx")), otsBed, flankingSize)
	}

################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")
assignGroupsTALEN(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_group_summary.xlsx")),
	hits = NULL,
	score = NULL, pv = NULL
	)	
	
if(filtName != ""){		
groupSummary(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_group_summary_", filtName, ".xlsx")),
	hits = hits.cutoff,
	score = score.cutoff, pv = pv.cutoff
	)	
}		 

################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGeneTALEN(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")))
annotateGeneTALEN(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")))

# Barplot per group
geneBarplot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenesTALEN(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	geneMat = geneMat,
	genes.width = 0, site.width = 100000)


################     RETURN FINAL XLSX FILE    ################ 
print("################     RETURN FINAL XLSX FILE    ################")
finalize(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))	


################     HISTONE MARKS    ################ 
#print("################     HISTONE MARKS    ################ ")

#histoneForestPlot(inputF = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	randomF = file.path(randomD, paste0(randomName, ".bed")),
#	histFiles = histoneFiles)

################     CHR PLOT    ################ 
if(filtName != ""){	
chrPlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_aln_", filtName, "_chrPlot.pdf")),
	hits = hits.cutoff, score = score.cutoff, pv = pv.cutoff)
chrPlotAside(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_aln_", filtName, "_chrPlot")),
	hits = hits.cutoff, score = score.cutoff, pv = pv.cutoff)
}




#if(filtName != ""){	
################     CHR PLOT    ################ 
#circlizePipeline(siteFile = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#				 zoom.size = 50000, label = FALSE, 
#                bestScore.cutoff = 9, bestFlank.cutoff = 25,
#                 gene.bed = NULL, ots.bed = FALSE, 
#                outFile = file.path(guideD, paste0(sampleName, "_w", w, "_circlize.pdf")),
#                species = "hg38")
#}



##########################################################################################
############                        UNTREATED SAMPLE                          ############
##########################################################################################



################     GUIDE SEQ ALIGNMENT    ################ 
print("################     GUIDE SEQ ALIGNMENT    ################")

# DO GUIDE ALIGNMENT ON REAL SEQUENCES
guideD <- resultD
getGuideAlignmentDual(inputF = file.path(resultD, paste0(controlName, "_w", w, ".xlsx")),
				  guideLeft = refSeq.left, guideRight = refSeq.right,
				  alnFolder = guideD,
				  gnm = GNM
				  )
#file.remove(list.files(guideD, pattern = "_TMP.txt", full.names = TRUE))					  
				  

# ASSIGN PVALUE FOR EVERY SINGLE COMBINATION
assignPV(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")),
	file.path(randomD, paste0(randomName, "_aln_stat.xlsx")))
	
# ASSIGN THE BEST COMBINATION
assignBestCombination(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")))	


################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
addFlanking(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")), otsBed, flankingSize)


################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")
assignGroupsTALEN(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_group_summary.xlsx")),
	hits = NULL,
	score = NULL, pv = NULL
	)	
	
if(filtName != ""){		
groupSummary(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_group_summary_", filtName, ".xlsx")),
	hits = hits.cutoff,
	score = score.cutoff, pv = pv.cutoff
	)	
}		 

################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGeneTALEN(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")))
annotateGeneTALEN(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")))

# Barplot per group
geneBarplot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenesTALEN(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	geneMat = geneMat,
	genes.width = 0, site.width = 100000)


################     RETURN FINAL XLSX FILE    ################ 
print("################     RETURN FINAL XLSX FILE    ################")
finalize(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))	


################     HISTONE MARKS    ################ 
#print("################     HISTONE MARKS    ################ ")

#histoneForestPlot(inputF = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	randomF = file.path(randomD, paste0(randomName, ".bed")),
#	histFiles = histoneFiles)

################     CHR PLOT    ################ 
if(filtName != ""){	
chrPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_aln_", filtName, "_chrPlot.pdf")),
	hits = hits.cutoff, score = score.cutoff, pv = pv.cutoff)
chrPlotAside(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_aln_", filtName, "_chrPlot")),
	hits = hits.cutoff, score = score.cutoff, pv = pv.cutoff)
}




#if(filtName != ""){	
################     CHR PLOT    ################ 
#circlizePipeline(siteFile = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#				 zoom.size = 50000, label = FALSE, 
#                 bestScore.cutoff = 9, bestFlank.cutoff = 25,
#                 gene.bed = NULL, ots.bed = FALSE, 
#                 outFile = file.path(guideD, paste0(controlName, "_w", w, "_circlize.pdf")),
#                 species = "hg38")
#}




##########################################################################################
############                  OVERALL GENOMIC INSTABILITY                     ############
##########################################################################################

#getOGI(rawFastq.de = file.path(dataD, "fastq", paste0(sampleName, "_R2_001.fastq.gz")),
#	filtFastq.de = file.path(sampleD, "results/fastq_aln", paste0(sampleName, "_pos.fastq.gz")),
#	groupF.de = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx")),
#	rawFastq.u = file.path(dataD, "fastq", paste0(controlName, "_R2_001.fastq.gz")),
#	filtFastq.u = file.path(sampleD, "results/fastq_aln", paste0(controlName, "_pos.fastq.gz")),
#	groupF.u = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx")),
#	otsF = otsBed,
#	outF = file.path(guideD, paste0(sampleName, "_w", w, "_OGI.txt"))
#	)

##########################################################################################
############                          SAVE READS                              ############
##########################################################################################


if(saveReads){
	print("################     SAVE READS    ################")
	
	# Sample
	dir.create(file.path(guideD, paste0(sampleName, "_reads")), showWarnings = FALSE)
	getReads(inputFile = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")),
			 bamFile = file.path(sampleD, "results/fastq_aln", paste0(sampleName, "_AlignmentSort.bam")),
			 outputDir = file.path(guideD, paste0(sampleName, "_reads"))
			 )
	addNbReads(inputFile = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")),
			 bamFile = file.path(sampleD, "results/fastq_aln", paste0(sampleName, "_AlignmentSort.bam")))		 
	
	# Control
	dir.create(file.path(guideD, paste0(controlName, "_reads")), showWarnings = FALSE)
	getReads(inputFile = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")),
			 bamFile = file.path(sampleD, "results/fastq_aln", paste0(controlName, "_AlignmentSort.bam")),
			 outputDir = file.path(guideD, paste0(controlName, "_reads"))
			 )	
 	addNbReads(inputFile = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")),
			 bamFile = file.path(sampleD, "results/fastq_aln", paste0(controlName, "_AlignmentSort.bam")))	
}


##########################################################################################
############                        CHR PLOT (CIRCLIZE)                       ############
##########################################################################################
if(filtName != ""){
  print("############    CHR PLOT (CIRCLIZE)    ############")
  ################     CHR PLOT (CIRCLIZE)    ################ 
  circlizePipelineTALEN(siteFile = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
                   zoom.size = 25000, label = FALSE, 
                   PV.cutoff = pv.cutoff,
                   bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                   showNBS = TRUE,
                   gene.bed = NULL, ots.bed = otsBed, 
                   outFile = file.path(guideD, paste0(sampleName, "_w", w, "_circlize_25k.pdf")),
                   species = circos.sp)
  
  circlizePipelineTALEN(siteFile = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
                        zoom.size = 25000, label = FALSE, 
                        PV.cutoff = pv.cutoff,
                        bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                        showNBS = FALSE,
                        gene.bed = NULL, ots.bed = otsBed, 
                        outFile = file.path(guideD, paste0(sampleName, "_w", w, "_circlize_25k_woNBS.pdf")),
                        species = circos.sp)
  
}

##########################################################################################
############                           HITS BARPLOT                           ############
##########################################################################################
print("############    HITS BARPLOT    ############")

hitsBarplot(file.path(guideD, paste0(sampleName, "_w", w, "_FINAL.xlsx")),
            pv = pv.cutoff, top = 50, showNBS = TRUE, log = TRUE,
            outName = file.path(guideD, paste0(sampleName, "_w", w, "_hits_barplot.pdf")))
hitsBarplot(file.path(guideD, paste0(sampleName, "_w", w, "_FINAL.xlsx")),
            pv = pv.cutoff, top = 50, showNBS = FALSE, log = TRUE,
            outName = file.path(guideD, paste0(sampleName, "_w", w, "_hits_barplot_woNBS.pdf")))

##########################################################################################
############                        REMOVE TMP FILES                          ############
##########################################################################################

if(rmTMP){
  print("############    REMOVE TMP FILES    ############")
  mypattern <- c("_final.fastq.gz$", "_Filt3_neg.fastq.gz$", "_trim2.fastq.gz$",
                 "_trim1.fastq.gz$", "_Filt2.fastq.gz$", "_pos.fastq.gz$", "_pos_NOmatch.fastq.gz$",
                 "_merged.fastq.gz$", ".unassembled.reverse.fastq.gz$", ".unassembled.forward.fastq.gz$",
                 "_assembled.fastq.gz$", "_unassembled.R1.fastq.gz$", "_unassembled.R2.fastq.gz", 
                 ".assembled.fastq.gz$", ".discarded.fastq$")
  torm <- unlist(lapply(mypattern, function(mp) list.files(sampleD, pattern = mp, full.names = TRUE, recursive = TRUE)))
  file.remove(torm)
  
  mypattern <- c("_aln_stat_FLANK_GROUP_GENES.xlsx$", "_aln_stat_FLANK_GROUP.xlsx$",
                 "_aln_stat_FLANK.xlsx$", "_aln_stat.xlsx$")
  torm <- unlist(lapply(mypattern, function(mp) list.files(guideD, pattern = mp, full.names = TRUE, recursive = TRUE)))
  file.remove(torm)
}


##########################################################################################
############                        COMPRESS RESULTS                          ############
##########################################################################################

if(tozip){
  print("############    COMPRESS RESULTS    ############")
  setwd(file.path(sampleD, "results"))
  zip.files <- c("guide_aln",
                 file.path("fastq_aln", paste0(sampleName, "_QC.xlsx"))
                 #file.path("fastq_aln", paste0(sampleName, "_pipeline.log")),
                 #file.path("fastq_aln", paste0(controlName, "_pipeline.log"))
  )
  utils::zip(zipfile = file.path(paste0(sampleName, "_CASTseq.zip")), files = zip.files)
}


}



