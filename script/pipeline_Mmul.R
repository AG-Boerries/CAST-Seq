
runPipeline <- function()
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
	testSample = sampleName,
	utSample = controlName,
	cpu = NBCPU
	)
}


################     DELTA AND CLUSTERS    ################ 
print("################     DELTA AND CLUSTERS    ################ ")
# bed files (test and control) should be available in resultD

# test sample
sampleBed <- file.path(resultD, paste0(sampleName, "_Alignment.bed"))
		   
# control sample		   
controlBed <- file.path(resultD, paste0(controlName, "_Alignment.bed"))	
	   
# run delta analysis
getDelta(sampleBed, otsBed, otsDistance)		   
getDeltaShuffle(gsub(".bed$", "_delta.bed", sampleBed), nb = 10, distance = distance.cutoff, genome.size = myGenome.size)				   
getDelta(controlBed, otsBed, otsDistance)		   
getDeltaShuffle(gsub(".bed$", "_delta.bed", controlBed), nb = 10, distance = distance.cutoff, genome.size = myGenome.size)		

# Density without on-target sites
#deltaDensity(gsub(".bed$", "_delta.bed", sampleBed), otsBed, surrounding_size)
#deltaDensity(gsub(".bed$", "_delta.bed", controlBed), otsBed, surrounding_size)
	  	  
# run cluster analysis
getCluster(gsub(".bed$", "_delta.bed", sampleBed), distance.cutoff)
getCluster(gsub(".bed$", "_delta.bed", controlBed), distance.cutoff)

################     TEST VS. CONTROL ENRICHMENT    ################
print("################     TEST VS. CONTROL ENRICHMENT    ################")
# USE RAW FASTQ FILE !!!!!!!!!
nbReads.sample <- as.numeric(nbReadFastqgz(file.path(dataD, "fastq", paste0(sampleName, "_R2_001.fastq.gz"))))
nbReads.control <- as.numeric(nbReadFastqgz(file.path(dataD, "fastq", paste0(controlName, "_R2_001.fastq.gz"))))

doEnrichment(gsub(".bed$", "_cluster.bed", sampleBed), gsub(".bed$", "_cluster.bed", controlBed),
			 nbReads.sample, nbReads.control, w, myGenome.size)
	
doEnrichment(gsub(".bed$", "_cluster.bed", controlBed), gsub(".bed$", "_cluster.bed", sampleBed),
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
getGuideAlignment(inputF = file.path(resultD, paste0(sampleName, "_w", w, ".xlsx")),
				  guide = refSeq,
				  alnFolder = guideD,
				  gnm = BSgenome.Mmulatta.UCSC.rheMac8::Mmulatta
				  )
file.remove(list.files(guideD, pattern = "_TMP.txt", full.names = TRUE))					  
				  
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
	getGuideAlignment(inputF = file.path(randomD, paste0(randomName, ".bed")),
					  guide = refSeq,
					  alnFolder = randomD,
					  gnm = BSgenome.Mmulatta.UCSC.rheMac8::Mmulatta
					  )	
	file.remove(list.files(randomD, pattern = "_TMP.txt", full.names = TRUE))
	}
	
# filt name
filtName <- ""
if(!is.null(clusters.cutoff)) filtName <- c(filtName, "clust", clusters.cutoff)
if(!is.null(pv.cutoff)) filtName <- c(filtName, "pv", clusters.cutoff)
if(!is.null(score.cutoff)) filtName <- c(filtName, "score", score.cutoff)
filtName <- paste0(filtName, collapse = "_")		

# PLOT GUIDE ALIGNMENT
guidePlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
		  file.path(guideD, paste0(sampleName, "_w", w, "_aln_heatmap.pdf")),
		  score = NULL, pv = NULL, ref = refSeq)# ALL
		  
if(filtName != ""){		  
	guidePlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
		  file.path(guideD, paste0(sampleName, "_w", w, "_aln_heatmap_", filtName, ".pdf")),
		  clusters = clusters.cutoff,
		  score = score.cutoff, pv = pv.cutoff, ref = refSeq)# Significant		  
}	  
	  
# LOGO PLOT
logoPlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
		 file.path(guideD, paste0(sampleName, "_w", w, "_aln_logo.pdf")),
		 score = NULL, pv = NULL, ref = refSeq)# ALL

if(filtName != ""){		
logoPlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
		 file.path(guideD, paste0(sampleName, "_w", w, "_aln_logo_", filtName, ".pdf")),
		 clusters = clusters.cutoff,
		 score = score.cutoff, pv = pv.cutoff, ref = refSeq)# Significant
}

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
assignGroups(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_group_summary.xlsx")),
	clusters = NULL,
	score = NULL, pv = NULL
	)	
	
if(filtName != ""){			
groupSummary(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_group_summary_", filtName, ".xlsx")),
	clusters = clusters.cutoff,
	score = score.cutoff, pv = pv.cutoff
	)	
}			 
			 
################     PIE CHARTS    ###############
print("################     PIE CHARTS    ###############")
piePlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	    file.path(guideD, paste0(sampleName, "_w", w, "_aln_piechart.pdf")),
	    score = NULL, pv = NULL)
	  
if(filtName != ""){			    
piePlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	    file.path(guideD, paste0(sampleName, "_w", w, "_aln_piechart_", filtName, ".pdf")),
		clusters = clusters.cutoff,
		score = score.cutoff, pv = pv.cutoff)	    
}

################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############
print("################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############")
pcBarplot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_heatmap.xlsx")))# ALL
if(file.exists(paste0(sampleName, "_w", w, "_aln_heatmap_clust1_pv0.05.xlsx"))) pcBarplot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_heatmap_clust1_pv0.05.xlsx")))# Significant

################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGene(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")))
annotateGene(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")))

# Barplot per group
geneBarplot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenes(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	geneMat = geneMat,
	genes.width = 0, site.width = 100000)


################     HISTONE MARKS    ################ 
#print("################     HISTONE MARKS    ################ ")

#histoneForestPlot(inputF = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	randomF = file.path(randomD, paste0(randomName, ".bed")),
#	histFiles = histoneFiles)


################     SCORING SYSTEM    ################ 
#print("################     SCORING SYSTEM    ################")

#addScore(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
#	pv.cutoff,
#	otsBed,
#	surrounding_size)
	
#scoreDensity(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")))	


################     CHR PLOT    ################ 
chrPlot(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_aln_clust1_pv0.05_chrPlot.pdf")), clusters = 1, score = NULL, pv = pv.cutoff)
chrPlotAside(file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(sampleName, "_w", w, "_aln_clust1_pv0.05_chrPlot")), clusters = 1, score = NULL, pv = pv.cutoff)



##########################################################################################
############                        UNTREATED SAMPLE                          ############
##########################################################################################

################     GUIDE SEQ ALIGNMENT    ################ 
print("################     GUIDE SEQ ALIGNMENT    ################")

# DO GUIDE ALIGNMENT ON REAL SEQUENCES
guideD <- resultD
getGuideAlignment(inputF = file.path(resultD, paste0(controlName, "_w", w, ".xlsx")),
				  guide = refSeq,
				  alnFolder = guideD,
				  gnm = BSgenome.Mmulatta.UCSC.rheMac8::Mmulatta
				  )
file.remove(list.files(guideD, pattern = "_TMP.txt", full.names = TRUE))					  
	
				  
# PLOT GUIDE ALIGNMENT
guidePlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")),
		  file.path(guideD, paste0(controlName, "_w", w, "_aln_heatmap.pdf")),
		  score = NULL, pv = NULL, ref = refSeq)# ALL
		
if(filtName != ""){			  
guidePlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")),
		  file.path(guideD, paste0(controlName, "_w", w, "_aln_heatmap_", filtName, ".pdf")),
		  clusters = clusters.cutoff,
		  score = score.cutoff, pv = pv.cutoff, ref = refSeq)# Significant		  
}	
	  
# LOGO PLOT
logoPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")),
		 file.path(guideD, paste0(controlName, "_w", w, "_aln_logo.pdf")),
		 score = NULL, pv = NULL, ref = refSeq)# ALL

if(filtName != ""){	
logoPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")),
		 file.path(guideD, paste0(controlName, "_w", w, "_aln_logo_", filtName, ".pdf")),
		 clusters = clusters.cutoff,
		 score = score.cutoff, pv = pv.cutoff, ref = refSeq)# Significant
}

################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
addFlanking(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat.xlsx")), otsBed, flankingSize)


################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")

assignGroups(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_group_summary.xlsx")),
	clusters = clusters.cutoff,
	score = score.cutoff, pv = pv.cutoff
	)	
	
if(filtName != ""){		
groupSummary(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_group_summary_", filtName, ".xlsx")),
	clusters = clusters.cutoff,
	score = score.cutoff, pv = pv.cutoff
	)	
}			 
			 
################     PIE CHARTS    ###############
print("################     PIE CHARTS    ###############")
piePlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	    file.path(guideD, paste0(controlName, "_w", w, "_aln_piechart.pdf")),
	    score = NULL, pv = NULL)
	 
if(filtName != ""){		    
piePlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	    file.path(guideD, paste0(controlName, "_w", w, "_aln_piechart_", filtName, ".pdf")),
	    clusters = clusters.cutoff,
		score = score.cutoff, pv = pv.cutoff)	    
}

################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############
print("################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############")
pcBarplot(file.path(guideD, paste0(controlName, "_w", w, "_aln_heatmap.xlsx")))# ALL
if(file.exists(paste0(controlName, "_w", w, "_aln_heatmap_clust1_pv0.05.xlsx"))) pcBarplot(file.path(guideD, paste0(controlName, "_w", w, "_aln_heatmap_clust1_pv0.05.xlsx")))# Significant

################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGene(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")))

# Barplot per group
geneBarplot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenes(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	geneMat = geneMat,
	genes.width = 0, site.width = 100000)


################     HISTONE MARKS    ################ 
#print("################     HISTONE MARKS    ################ ")

#histoneForestPlot(inputF = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	randomF = file.path(randomD, paste0(randomName, ".bed")),
#	histFiles = histoneFiles)


################     SCORING SYSTEM    ################ 
#print("################     SCORING SYSTEM    ################")

#addScore(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
#	pv.cutoff,
#	otsBed,
#	surrounding_size)
	
#scoreDensity(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")))	


################     CHR PLOT    ################ 
chrPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_aln_clust1_pv0.05_chrPlot.pdf")), clusters = 1, score = NULL, pv = pv.cutoff)
chrPlotAside(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_aln_clust1_pv0.05_chrPlot")), clusters = 1, score = NULL, pv = pv.cutoff)




##########################################################################################
############                  OVERALL GENOMIC INSTABILITY                     ############
##########################################################################################

#getOGI(rawFastq.de = file.path(dataD, "fastq", paste0(sampleName, "_R2_001.fastq.gz")),
#	filtFastq.de = file.path(sampleD, "results/fastq_aln", paste0(sampleName, "_pos.fastq.gz")),
#	groupF.de = file.path(guideD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")),
#	rawFastq.u = file.path(dataD, "fastq", paste0(controlName, "_R2_001.fastq.gz")),
#	filtFastq.u = file.path(sampleD, "results/fastq_aln", paste0(controlName, "_pos.fastq.gz")),
#	groupF.u = file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")),
#	otsF = otsBed,
#	outF = file.path(guideD, paste0(sampleName, "_w", w, "_OGI.txt"))
#	)


}

