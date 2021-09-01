
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
getGuideAlignment(inputF = file.path(ovlD, paste0(ovlName, ".xlsx")),
				  guide = refSeq,
				  alnFolder = ovlD,
				  gnm = GNM
				  )
file.remove(list.files(ovlD, pattern = "_TMP.txt", full.names = TRUE))					  
				  
# GENERATE RANDOM SEQUENCE BED
randomD <- file.path(sampleD, "results", "random")
dir.create(randomD, showWarnings = FALSE)
randomName <- paste0("random_w", w)

# PLOT GUIDE ALIGNMENT
guidePlot(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")),
		  file.path(ovlD, paste0(ovlName, "_aln_heatmap.pdf")),
		  score = NULL, pv = NULL, ref = refSeq)# ALL
		  
# LOGO PLOT
logoPlot(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")),
		 file.path(ovlD, paste0(ovlName, "_aln_logo.pdf")),
		 score = NULL, pv = NULL, ref = refSeq)# ALL



################     FLANKING REGIONS    ################ 
print("################     FLANKING REGIONS    ################ ")
# REAL SEQUENCES
if(is.null(flank1.sq)){
	addFlanking(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")), otsBed, flankingSize)
}else{
	addFlankingFromSq(file.path(ovlD, paste0(ovlName, "_aln_stat.xlsx")), flank1.sq, flank2.sq)
}

print(file.exists(file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx"))))

################     DEFINE GROUPS    ################ 
print("################     DEFINE GROUPS    ################ ")
assignGroups(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK.xlsx")),
			 file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
			 otsBed,
			 pv.cutoff)
			 
groupSummary(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")),
	file.path(ovlD, paste0(ovlName, "_group_summary.xlsx")),
	hits = NULL,
	score = NULL, pv = NULL
	)	
	
guidePlot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")),
          file.path(ovlD, paste0(ovlName, "_aln_heatmap_OMT.pdf")),
          score = NULL, pv = NULL, ref = refSeq, OMTonly = TRUE)# ALL	 
			 
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
annotateGene(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")))

# Barplot per group
geneBarplot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")))

# Forest plot per group
geneForestPlot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
			   file.path(randomD, paste0(randomName, "_aln_stat_FLANK_GENES.xlsx")))

################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
addGenes(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	oncoFile = oncoEntrez,
	geneMat = geneMat,
	genes.width = 0, site.width = 100000)
	
################     RETURN FINAL XLSX FILE    ################ 
print("################     RETURN FINAL XLSX FILE    ################")
finalizeOverlap(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")))	
		


################     HISTONE MARKS    ################ 
#print("################     HISTONE MARKS    ################ ")

#histoneForestPlot(inputF = file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
#	randomF = file.path(randomD, paste0(randomName, ".bed")),
#	histFiles = histoneFiles)

################     SCORING SYSTEM    ################ 
#print("################     SCORING SYSTEM    ################")

#addScore(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")),
#	file.path(randomD, paste0(randomName, "_aln_stat_FLANK.xlsx")),
#	pv.cutoff,
#	otsBed,
#	surrounding_size)
	
#scoreDensity(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx")))	


################     CHR PLOT    ################ 
chrPlot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(ovlD, paste0(ovlName, "_aln_chrPlot.pdf")), hits = NULL, score = NULL, pv = NULL)
chrPlotAside(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
	file.path(ovlD, paste0(ovlName, "_aln_chrPlot")), hits = NULL, score = NULL, pv = NULL)

##########################################################################################
############                        CHR PLOT (CIRCLIZE)                       ############
##########################################################################################
print("############    CHR PLOT (CIRCLIZE)    ############")
################     CHR PLOT (CIRCLIZE)    ################ 
circlizePipeline(siteFile = file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
                   zoom.size = 25000, label = FALSE, 
                   PV.cutoff = NULL,
                   bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                   showNBS = TRUE,
                   gene.bed = NULL, ots.bed = otsBed, 
                   realigned = TRUE,
                   outFile = file.path(ovlD, paste0(ovlName, "_circlize_25k.pdf")),
                   species = circos.sp)
  

circlizePipeline(siteFile = file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP_GENES.xlsx")),
                 zoom.size = 25000, label = FALSE, 
                 PV.cutoff = NULL,
                 bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                 showNBS = FALSE,
                 gene.bed = NULL, ots.bed = otsBed, 
                 realigned = TRUE,
                 outFile = file.path(ovlD, paste0(ovlName, "_circlize_25k_woNBS.pdf")),
                 species = circos.sp)

##########################################################################################
############                           HITS BARPLOT                           ############
##########################################################################################
print("############    HITS BARPLOT    ############")

hitsBarplot(file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
            pv = NULL, top = 50, showNBS = TRUE, log = TRUE,
            outName = file.path(ovlD, paste0(ovlName, "_hits_barplot.pdf")))
hitsBarplot(file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
            pv = NULL, top = 50, showNBS = FALSE, log = TRUE,
            outName = file.path(ovlD, paste0(ovlName, "_hits_barplot_woNBS.pdf")))

##########################################################################################
############                        REMOVE TMP FILES                          ############
##########################################################################################

if(rmTMP){
  print("############    REMOVE TMP FILES    ############")
  mypattern <- c("_aln_stat_FLANK_GROUP_GENES.xlsx$", "_aln_stat_FLANK_GROUP.xlsx$",
                 "_aln_stat_FLANK.xlsx$", "_aln_stat.xlsx$")
  torm <- unlist(lapply(mypattern, function(mp) list.files(ovlD, pattern = mp, full.names = TRUE, recursive = TRUE)))
  file.remove(torm)
}



}

