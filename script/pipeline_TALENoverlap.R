
runPipelineTALENoverlap <- function()
{

##########################################################################################
############                           RUN PIPELINE                           ############
##########################################################################################

################     CHECK INPUT    ################ 
print("################     CHECK INPUT    ################")

print(ovlD)
print(ovlName)

##########################################################################################
############                         OVERLAP ANALYSIS                         ############
##########################################################################################
  
print("################    OVERLAP ANALYSIS    ################")

print(replicates)
print(repNames)
print(repD)
  
replicates.split <- strsplit(replicates, split = ",")[[1]]
repNames.split <- strsplit(repNames, split = ",")[[1]]
  
siteFiles <- paste0(replicates.split, "_w", w, "_FINAL.xlsx")
#siteFiles <- sapply(siteFiles, function(i) list.files(repD, pattern = i, recursive=TRUE, full.names = TRUE))
siteFiles <- sapply(1:length(siteFiles), function(i) list.files(repD.split[i], pattern = siteFiles[i], recursive=TRUE, full.names = TRUE))
siteFiles <- lapply(siteFiles, function(i) i[grepl("\\/results\\/", i)])
siteFiles <- unlist(siteFiles)
names(siteFiles) <- repNames.split
  
print(siteFiles)
if(length(siteFiles)<2) stop("at least two files are needed to perform overlap analysis")

compList <- dfComparisonList(siteFiles, names(siteFiles), width = distance.cutoff, nb.signif = nb.ovl, NBS = TRUE)
write.xlsx(compList, file.path(ovlD, paste0(ovlName, ".xlsx")), overwrite = TRUE)
makeUpset(file.path(ovlD, paste0(ovlName, ".xlsx")))


##########################################################################################
############                DESIGNER NUCLEASE TREATED SAMPLE                  ############
##########################################################################################

################     GUIDE SEQ ALIGNMENT    ################ 
print("################     GUIDE SEQ ALIGNMENT    ################")

print(file.path(ovlD, paste0(ovlName, ".xlsx")))
print(file.exists(file.path(ovlD, paste0(ovlName, ".xlsx"))))

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
#print("################     PIE CHARTS    ###############")
#piePlot(file.path(ovlD, paste0(ovlName, "_aln_stat_FLANK_GROUP.xlsx")),
#	    file.path(ovlD, paste0(ovlName, "_aln_piechart.pdf")),
#	    score = NULL, pv = NULL)
	    

################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############
#print("################     MISMATCHES / INDEL PERCENTAGE BARPLOT    ###############")
#pcBarplot(file.path(ovlD, paste0(ovlName, "_aln_heatmap.xlsx")))# ALL


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
chrPlot(file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
        file.path(ovlD, paste0(ovlName, "_aln_chrPlot.pdf")), hits = NULL, score = NULL, pv = NULL)
chrPlotAside(file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
             file.path(ovlD, paste0(ovlName, "_aln_chrPlot")), hits = NULL, score = NULL, pv = NULL)



##########################################################################################
############                        CHR PLOT (CIRCLIZE)                       ############
##########################################################################################
print("############    CHR PLOT (CIRCLIZE)    ############")
################     CHR PLOT (CIRCLIZE)    ################ 


	tryCatch(
    	{
	circlizePipelineTALEN(siteFile = file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
                        zoom.size = 25000, label = FALSE, 
                        PV.cutoff = NULL,
                        bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                        showNBS = TRUE,
                        gene.bed = NULL, ots.bed = otsBed, 
                        realigned = TRUE,
                        outFile = file.path(ovlD, paste0(ovlName, "_circlize_25k.pdf")),
                        species = circos.sp)
    },
    error = function(e){
    print(read.delim(otsBed, header = FALSE))
	print("no sites on defined otsBed, use max gRNA score")
		circlizePipelineTALEN(siteFile = file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
                        zoom.size = 25000, label = FALSE, 
                        PV.cutoff = NULL,
                        bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                        showNBS = TRUE,
                        gene.bed = NULL, ots.bed = NULL, 
                        realigned = FALSE,
                        outFile = file.path(ovlD, paste0(ovlName, "_circlize_25k.pdf")),
                        species = circos.sp)

		    }
	)
	
	
	tryCatch(
    	{
		circlizePipelineTALEN(siteFile = file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
                      zoom.size = 25000, label = FALSE, 
                      PV.cutoff = NULL,
                      bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                      showNBS = FALSE,
                      gene.bed = NULL, ots.bed = otsBed, 
                      realigned = TRUE,
                      outFile = file.path(ovlD, paste0(ovlName, "_circlize_25k_woNBS.pdf")),
                      species = circos.sp)
    },
    error = function(e){
    print(read.delim(otsBed, header = FALSE))
	print("no sites on defined otsBed, use max gRNA score")
		circlizePipelineTALEN(siteFile = file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
                      zoom.size = 25000, label = FALSE, 
                      PV.cutoff = NULL,
                      bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                      showNBS = FALSE,
                      gene.bed = NULL, ots.bed = NULL, 
                      realigned = FALSE,
                      outFile = file.path(ovlD, paste0(ovlName, "_circlize_25k_woNBS.pdf")),
                      species = circos.sp)

		    }
	)




  
##########################################################################################
############                           HITS BARPLOT                           ############
##########################################################################################
print("############    HITS BARPLOT    ############")

hitsBarplot(file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
            pv = NULL, top = 50, showNBS = TRUE, log = TRUE,
            outName = file.path(ovlD, paste0(ovlName, "_w", w, "_hits_barplot.pdf")))
hitsBarplot(file.path(ovlD, paste0(ovlName, "_FINAL.xlsx")),
            pv = NULL, top = 50, showNBS = FALSE, log = TRUE,
            outName = file.path(ovlD, paste0(ovlName, "_w", w, "_hits_barplot_woNBS.pdf")))


##########################################################################################
############                        REMOVE TMP FILES                          ############
##########################################################################################

if(rmTMP){
  print("############    REMOVE TMP FILES    ############")
  mypattern <- c("_aln_stat_FLANK_GROUP_GENES.xlsx$", "_aln_stat_FLANK_GROUP.xlsx$",
                 "_aln_stat_FLANK.xlsx$", "_aln_stat.xlsx$")
  torm <- unlist(lapply(mypattern, function(mp) list.files(ovlD, pattern = mp, full.names = TRUE, recursive = TRUE)))
  torm <- torm[file.exists(torm)]
  file.remove(torm)
}

}




