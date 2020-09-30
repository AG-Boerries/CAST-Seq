
runPipelineOverlap <- function()
{

##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################

################     CHECK INPUT    ################ 
print("################     CHECK INPUT    ################")


print(ovlD)
print(ovlName)


##########################################################################################
############                DESIGNER NUCLEASE TREATED SAMPLE                  ############
##########################################################################################

################     GUIDE SEQ ALIGNMENT    ################ 
print("################     GUIDE SEQ ALIGNMENT    ################")

# DO GUIDE ALIGNMENT ON REAL SEQUENCES
getGuideAlignmentDual(inputF = file.path(ovlD, paste0(ovlName, ".xlsx")),
					  guideLeft = refSeq.left, guideRight = refSeq.right,
					  alnFolder = ovlD,
					  gnm = GNM
					  )			  
file.remove(list.files(ovlD, pattern = "_TMP.txt", full.names = TRUE))					  
				  
# GENERATE RANDOM SEQUENCE BED
randomD <- file.path(sampleD, "results", "random")
dir.create(randomD, showWarnings = FALSE)
randomName <- paste0("random_w", w)


# ASSIGN PVALUE FOR EVERY SINGLE COMBINATION
assignPV(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")),
	file.path(randomD, paste0(randomName, "_aln_stat.xlsx")))
	
# ASSIGN THE BEST COMBINATION
assignBestCombination(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")))	


################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
if(is.null(flank1.sq)){
	addFlanking(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")), otsBed, flankingSize)
}else{
	addFlankingFromSq(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")), flank1.sq, flank2.sq)
}

################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")
assignGroupsTALEN(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)			 
			 
groupSummary(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(ovlD, paste0(ovlName, "_group_summary.xlsx")),
	hits = hits.cutoff,
	score = score.cutoff,
	pv = pv.cutoff
	)	
		 
			 
################     PIE CHARTS    ###############
print("################     PIE CHARTS    ###############")
piePlot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")),
	    file.path(ovlD, paste0(ovlName, "_aln_piechart.pdf")),
	    score = NULL, pv = NULL)
	    

################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############
print("################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############")
pcBarplot(file.path(ovlD, paste0(ovlName, "_aln_heatmap.xlsx")))# ALL


################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGeneTALEN(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")))

# Barplot per group
geneBarplot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenesTALEN(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	geneMat = geneMat,
	genes.width = 0, site.width = 100000)


################     HISTONE MARKS    ################ 
print("################     HISTONE MARKS    ################ ")

histoneForestPlot(inputF = file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	randomF = file.path(randomD, paste0(randomName, ".bed")),
	histFiles = histoneFiles)

################     SCORING SYSTEM    ################ 
#print("################     SCORING SYSTEM    ################")

#addScore(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")),
#	file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
#	pv.cutoff,
#	otsBed,
#	surrounding_size)
	
#scoreDensity(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx")))	


################     CHR PLOT    ################ 
chrPlot(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_aln_", filtName, "_chrPlot.pdf")),
	hits = hits.cutoff, score = score.cutoff, pv = pv.cutoff)
chrPlotAside(file.path(guideD, paste0(controlName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(guideD, paste0(controlName, "_w", w, "_aln_", filtName, "_chrPlot")),
	hits = hits.cutoff, score = score.cutoff, pv = pv.cutoff)




}





runPipelineTALENoverlap <- function()
{

##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################


##########################################################################################
############                DESIGNER NUCLEASE TREATED SAMPLE                  ############
##########################################################################################

################     GUIDE SEQ ALIGNMENT    ################ 
print("################     GUIDE SEQ ALIGNMENT    ################")

# DO GUIDE ALIGNMENT ON REAL SEQUENCES
guideD <- resultD
getGuideAlignmentDual(inputF = file.path(ovlD, paste0(sampleName, "_w", w, ".xlsx")),
				  guideLeft = refSeq.left, guideRight = refSeq.right,
				  alnFolder = ovlD
				  )			  
				  
# GENERATE RANDOM SEQUENCE BED
randomD <- file.path(sampleD, "results", "random")
dir.create(randomD, showWarnings = FALSE)
randomName <- paste0("random_w",w)

# ASSIGN PVALUE FOR EVERY SINGLE COMBINATION
assignPV(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
	file.path(randomD, paste0(randomName, "_aln_stat.xlsx")))
	
# ASSIGN THE BEST COMBINATION
assignBestCombination(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")))	


################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
addFlanking(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")), otsBed, flankingSize)

################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")
assignGroupsTALEN(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(ovlD, paste0(sampleName, "_w", w, "_group_summary.xlsx")),
	clusters = NULL,
	score = NULL, pv = NULL
	)	 


################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGeneTALEN(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")))
annotateGeneTALEN(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")))

# Barplot per group
geneBarplot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenesTALEN(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	genes.width = 0, site.width = 100000)

################     ONCOGENE ANNOTATION    ################ 
#print("################     ONCOGENE ANNOTATION    ################")
#addOnco(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")), oncoEntrez, onco.width)


################     HISTONE MARKS    ################ 
print("################     HISTONE MARKS    ################ ")

histoneForestPlot(inputF = file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	randomF = file.path(randomD, paste0(randomName, ".bed")),
	histFiles = histoneFiles)

################     SCORING SYSTEM    ################ 
#print("################     SCORING SYSTEM    ################")

#addScore(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")),
#	file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
#	pv.cutoff,
#	otsBed,
#	surrounding_size)
	
#scoreDensity(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx")))	


################     CHR PLOT    ################ 
chrPlot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(ovlD, paste0(sampleName, "_w", w, "_aln_chrPlot.pdf")), clusters = NULL, score = NULL, pv = NULL)
chrPlotAside(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(ovlD, paste0(sampleName, "_w", w, "_aln_chrPlot")), clusters = NULL, score = NULL, pv = NULL)




}



