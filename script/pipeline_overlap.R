




runPipelineOverlap <- function()
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
getGuideAlignment(inputF = file.path(ovlD, paste0(sampleName, "_w", w, ".xlsx")),
				  guide = refSeq,
				  alnFolder = ovlD
				  )
file.remove(list.files(ovlD, pattern = "_TMP.txt", full.names = TRUE))					  
				  
# GENERATE RANDOM SEQUENCE BED
randomD <- file.path(sampleD, "results", "random")
dir.create(randomD, showWarnings = FALSE)
randomName <- paste0("random_w",w)

# PLOT GUIDE ALIGNMENT
guidePlot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
		  file.path(ovlD, paste0(sampleName, "_w", w, "_aln_heatmap.pdf")),
		  score = NULL, pv = NULL, ref = refSeq)# ALL
		  
# LOGO PLOT
logoPlot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")),
		 file.path(ovlD, paste0(sampleName, "_w", w, "_aln_logo.pdf")),
		 score = NULL, pv = NULL, ref = refSeq)# ALL



################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
addFlanking(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat.xlsx")), otsBed, flankingSize)

################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")
assignGroups(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(ovlD, paste0(sampleName, "_w", w, "_group_summary.xlsx")),
	clusters = NULL,
	score = NULL, pv = NULL
	)	
		 
			 
################     PIE CHARTS    ###############
print("################     PIE CHARTS    ###############")
piePlot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")),
	    file.path(ovlD, paste0(sampleName, "_w", w, "_aln_piechart.pdf")),
	    score = NULL, pv = NULL)
	    

################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############
print("################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############")
pcBarplot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_heatmap.xlsx")))# ALL


################     GENE ANNOTATION    ################ 
print("################     GENE ANNOTATION    ################ ")
# Annotation
annotateGene(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP.xlsx")))
annotateGene(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")))

# Barplot per group
geneBarplot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenes(file.path(ovlD, paste0(sampleName, "_w", w, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
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

