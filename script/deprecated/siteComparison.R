
library(GenomicRanges)
library(openxlsx)
library(rtracklayer)
library(ChIPpeakAnno)
library(UpSetR)

source("~/Scripts/work/tools/venn_script.r")

hg19to38 <- function(m)
{
	chain <- import.chain("~/Research/database/ucsc/hg19ToHg38.over.chain")
	
	m$strand <- "*"
	m.gr <- makeGRangesFromDataFrame(m, seqnames.field = "V1",
		start.field = "V2", end.field = "V3", strand.field = "strand",
		keep.extra.columns = TRUE)

	m38.gr = liftOver(m.gr, chain)

	m38.gr = unlist(m38.gr)
	genome(m38.gr) = "hg38"
	return(data.frame(m38.gr))
}


dfComparison <- function(f1, f2, name1, name2)
{
	df1 <- read.xlsx(f1, sheet = 1)
	df2 <- read.xlsx(f2, sheet = 1)
	
	df1 <- df1[, -4]# no strand
	df2 <- df2[, -4]# no strand
	
	#df1 <- df1[,  c("chromosome", "start", "end", "read", "collapseCluster", "read.ctl",
	#	"collapseCluster.ctl", "width.raw", "width", "OddRatio", "pvalue",
	#	"adj.pvalue", "artificial", "avg.MAPQ", "alignment.id", "pattern", "subject",
	#	"aln.sum", "score", "group")]
		
	#df2 <- df2[, c("chromosome", "start", "end", "read", "collapseCluster", "read.ctl",
	#	"collapseCluster.ctl", "width.raw", "width", "OddRatio", "pvalue",
	#	"adj.pvalue", "artificial", "avg.MAPQ", "alignment.id", "pattern", "subject",
	#	"aln.sum", "score", "group")]	
		

	gr1 <- makeGRangesFromDataFrame(df1, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	gr2 <- makeGRangesFromDataFrame(df2, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	
	gr.ovl <- findOverlaps(query = gr1, subject = gr2, type = "any", maxgap = width)
	if(length(queryHits(gr.ovl)) == 0)
		{
		df.ovl <- matrix(NA, nrow = 1, ncol = ncol(df1) + ncol(df2))
		df1.spe <- df1
		df2.spe <- df2
		}else{
		df.ovl <- data.frame(df1[queryHits(gr.ovl),], df2[subjectHits(gr.ovl),])
		df1.spe <- df1[-queryHits(gr.ovl),]
		df2.spe <- df2[-subjectHits(gr.ovl),]
		}
	colnames(df.ovl) <- c(paste(colnames(df1), name1, sep = "."),
		paste(colnames(df2), name2, sep = "."))
		
	if(length(queryHits(gr.ovl)) != 0){
		pv.name1 <- paste0("adj.pvalue.", name1)
		pv.name2 <- paste0("adj.pvalue.", name2)
		toKeep <- df.ovl[, pv.name1] < 0.05 | df.ovl[, pv.name2] < 0.05
		df.ovl.filt <- df.ovl[toKeep, ]
		
		idx.name1 <- grep(name1, colnames(df.ovl.filt))
		idx.name2 <- grep(name2, colnames(df.ovl.filt))
		
		df.ovl.filt1 <- df.ovl.filt[, idx.name1]
		colnames(df.ovl.filt1) <- colnames(df1)
		df.ovl.filt2 <- df.ovl.filt[, idx.name2]
		colnames(df.ovl.filt2) <- colnames(df2)
		
		df.ovl.filt <- rbind(df.ovl.filt1, df.ovl.filt2)
		df.ovl.filt <- df.ovl.filt[order(df.ovl.filt$chromosome, df.ovl.filt$start, df.ovl.filt$end), ]
		
		df.ovl.filt <- df.ovl.filt[!duplicated(df.ovl.filt), ]# remove duplicates
		
		# merge overlap (+/- 2000bp)
		df.ovl.filt <- df.ovl.filt[, 1:7]
		
		idx <- 0
		df.final <- c()
		for(chr in unique(df.ovl.filt[,1]))
			{
			tempdf <- subset(df.ovl.filt,  chromosome==chr)
			df.final <- rbind(df.final, tempdf[1,])
			idx <- idx + 1
			if(nrow(tempdf) > 1){
				for(i in 2:nrow(tempdf))
					{
					if(tempdf[i, "start"] <= (df.final[idx, "end"] + width))
						{
						df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] <-
							df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] + tempdf[i, c("read", "hits", "read.ctl", "hits.ctl")]# add read, hits...
						if(tempdf[i, "end"] > df.final[idx, "end"]) df.final[idx, "end"] <- tempdf[i, "end"]# change end if new end is higher
						
						}else{
							df.final <- rbind(df.final, tempdf[i, ])
							idx <- idx + 1
							}
					}
				}
			}
		
		compL <- list(common.filtered = df.final, common = df.ovl, df1.spe = df1.spe, df2.spe = df2.spe)
		names(compL) <- c("common.filtered", "common", paste0(name1, ".spe"), paste0(name2, ".spe"))
		}else{
		compL <- list(common = df.ovl, df1.spe = df1.spe, df2.spe = df2.spe)
		names(compL) <- c("common", paste0(name1, ".spe"), paste0(name2, ".spe"))
		}
	
	return(compL)
}


getReads <- function(dfA, dfB, width)
{
	reads <- rep(0, nrow(dfA))
	
	grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)

	gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
	
	idxA <- queryHits(gr.ovl)
	if(length(idxA) == 0) return(NA)
	
	idxA.unique <- unique(idxA)
	idxB <- subjectHits(gr.ovl)
	
	reads[idxA.unique] <- sapply(idxA.unique, function(i){
		currentIndex <- idxB[idxA == i]
		return(sum(dfB$read[currentIndex]))
	})

	return(reads)
}


getHits <- function(dfA, dfB, width)
{
	hits <- rep(0, nrow(dfA))
	
	grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)

	gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
	
	idxA <- queryHits(gr.ovl)
	if(length(idxA) == 0) return(NA)
	
	idxA.unique <- unique(idxA)
	idxB <- subjectHits(gr.ovl)
	
	hits[idxA.unique] <- sapply(idxA.unique, function(i){
		currentIndex <- idxB[idxA == i]
		return(sum(dfB$hits[currentIndex]))
	})

	return(hits)
}


getCommonIndices <- function(dfA, dfB, width)
{
	grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
		
	gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
	if(length(queryHits(gr.ovl)) == 0) return(NA)
	
	return(sort(unique(queryHits(gr.ovl))))
}



getCommon <- function(dfA, dfB, width)
{
	grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
		
	gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
	if(length(queryHits(gr.ovl)) == 0) return(NA)

	df.ovl <- data.frame(dfA[queryHits(gr.ovl),], dfB[subjectHits(gr.ovl),])

	pv.nameA <- grep("^adj.pvalue", colnames(df.ovl))[1]
	pv.nameB <- grep("^adj.pvalue", colnames(df.ovl))[2]
	toKeep <- df.ovl[, pv.nameA] < 0.05 | df.ovl[, pv.nameB] < 0.05# AND or OR ???
	df.ovl.filt <- df.ovl[toKeep, ]
	
	df.ovl.filt1 <- df.ovl.filt[, 1:ncol(dfA)]
	colnames(df.ovl.filt1) <- colnames(dfA)
	df.ovl.filt2 <- df.ovl.filt[, (ncol(dfA)+1):ncol(df.ovl.filt)]
	colnames(df.ovl.filt2) <- colnames(dfB)
	
	df.ovl.filt <- rbind(df.ovl.filt1, df.ovl.filt2)
	df.ovl.filt <- df.ovl.filt[order(df.ovl.filt$chromosome, df.ovl.filt$start, df.ovl.filt$end), ]
	
	df.ovl.filt <- df.ovl.filt[!duplicated(df.ovl.filt), ]# remove duplicates

	return(df.ovl.filt)
}


getCommon_fromList <- function(dfList, width)
{
	combMat <- combn(1:length(dfList), 2)
	df.ovl.filt <- lapply(1:ncol(combMat), function(i)
		getCommon(dfList[[combMat[1,i]]], dfList[[combMat[2,i]]], width)
		)
	df.ovl.filt <- do.call(rbind, df.ovl.filt)	
	df.ovl.filt <- df.ovl.filt[order(df.ovl.filt$chromosome, df.ovl.filt$start, df.ovl.filt$end), ]
	df.ovl.filt <- df.ovl.filt[!duplicated(df.ovl.filt), ]# remove duplicates	

	return(df.ovl.filt)
}



dfComparisonList <- function(fList, nList, width = 2000, NBS = TRUE)
{
	dfList <- lapply(fList, read.xlsx, sheet = 1)
	dfList <- lapply(dfList, function(i) i[,-4])
	#grList <- lapply(dfList, function(i) makeGRangesFromDataFrame(i, seqnames.field = "chromosome",
	#	start.field = "start", end.field = "end",
	#	keep.extra.columns = TRUE, ignore.strand = TRUE))
	
	# Remove NBS
	if(!NBS) dfList <- lapply(dfList, function(i) i[i$group != "NBS", ])
		
	# Overlaping sites	
	df.ovl.filt.raw <- getCommon_fromList(dfList, width)	

	# merge overlap (+/- 2000bp)
	df.ovl.filt <- df.ovl.filt.raw[, 1:7]
	
	idx <- 0
	df.final <- c()
	for(chr in unique(df.ovl.filt[,1]))
		{
		tempdf <- subset(df.ovl.filt,  chromosome==chr)
		df.final <- rbind(df.final, tempdf[1,])
		idx <- idx + 1
		if(nrow(tempdf) > 1){
			for(i in 2:nrow(tempdf))
				{
				if(tempdf[i, "start"] <= (df.final[idx, "end"] + width))
					{
					df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] <-
						df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] + tempdf[i, c("read", "hits", "read.ctl", "hits.ctl")]# add read, hits...
					if(tempdf[i, "end"] > df.final[idx, "end"]) df.final[idx, "end"] <- tempdf[i, "end"]# change end if new end is higher
					
					}else{
						df.final <- rbind(df.final, tempdf[i, ])
						idx <- idx + 1
						}
				}
			}
		}
	gr.final <- makeGRangesFromDataFrame(df.final, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)	
		
	# Specific sites
	speList <- lapply(dfList, function(i){
		gr <- makeGRangesFromDataFrame(i, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)
		gr.ovl <- findOverlaps(query = gr, subject = gr.final, type = "any", maxgap = width)
		df.spe <- i[-queryHits(gr.ovl),]
		return(df.spe)
		})
	names(speList) <- nList	
	
	# Intersect common sites vs. individual sites
	readList <- lapply(dfList, getReads, dfA = df.final, width = 0)
	readM <- do.call(cbind, readList)
	colnames(readM) <- paste0(nList, "_read")
	
	hitList <- lapply(dfList, getHits, dfA = df.final, width = 0)
	hitM <- do.call(cbind, hitList)
	colnames(hitM) <- paste0(nList, "_hits")

	df.final <- cbind(df.final, readM, hitM)	
	
	# Set the final number of "read" and "hits"	
	df.final$read <- rowSums(df.final[, grep("_read", colnames(df.final))])	
	df.final$hits <- rowSums(df.final[, grep("_hits", colnames(df.final))])	
	
	# Re-order rows
	df.final <- df.final[order(-df.final$hits, -df.final$read), ]
	
	# return
	compL <- c(list(common.filtered = df.final), list(common = df.ovl.filt.raw), speList)
	return(compL)
	
}


makeUpset <- function(ovlFile){
	mysheets <- getSheetNames(ovlFile)
  	mList <- lapply(mysheets, read.xlsx, xlsxFile = ovlFile)
  	names(mList) <- mysheets
  	
  	ovlMat <- mList[["common.filtered"]]
  	
  	# get nb of reads
  	nbMat <- ovlMat[, grepl("*_read$", colnames(ovlMat))]

	siteNames <- paste(ovlMat$chromosome, ovlMat$start, ovlMat$end, sep = "_")
	sampleNames <- gsub("_read$", "", colnames(nbMat))
	rownames(nbMat) <- siteNames
	colnames(nbMat) <- sampleNames
	
	# GET siteList from nbMat
	siteList <- lapply(1:ncol(nbMat), function(i){
		return(rownames(nbMat)[as.numeric(nbMat[, i]) != 0])	
	})
	names(siteList) <- colnames(nbMat)
	
	# ADD specific sites
	for(i in names(siteList)){
		speMat <- mList[[i]]
		speMat <- speMat[speMat$adj.pvalue < 0.05, ]
		siteList[[i]] <- c(siteList[[i]], paste(speMat$chromosome, speMat$start, speMat$end, sep = "_"))
	}
		
	# UpSet plot

	pdf(gsub(".xlsx", "_UpSetR.pdf", ovlFile), width = 10, onefile = FALSE)
  	print(upset(fromList(siteList), order.by = "freq", nsets = length(siteList),
              sets = names(siteList), keep.order = TRUE,
              mainbar.y.label = "Intersection", sets.x.label = "sites per sample"))
  	dev.off()	
  	
  	
	# Save Venn
  vennOut <- doVenn(siteList)
  	
  		
}

dfComparisonV2 <- function(f1, f2, name1, name2)
{
	df1 <- read.xlsx(f1, sheet = 1)
	df2 <- read.xlsx(f2, sheet = 1)
	
	df1 <- df1[, -4]# no strand
	#df1 <- df1[,  c("chromosome", "start", "end", "read", "collapseCluster", "read.ctl",
	#	"collapseCluster.ctl", "width.raw", "width", "OddRatio", "pvalue",
	#	"adj.pvalue", "artificial", "avg.MAPQ", "alignment.id", "pattern", "subject",
	#	"aln.sum", "score", "group")]
		
	df2 <- df2[, c("seqnames", "start", "end", "Read")]	
		

	gr1 <- makeGRangesFromDataFrame(df1, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	gr2 <- makeGRangesFromDataFrame(df2, seqnames.field = "seqnames",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	
	gr.ovl <- findOverlaps(query = gr1, subject = gr2, type = "any", maxgap = 2000)
	if(length(queryHits(gr.ovl)) == 0)
		{
		df.ovl <- matrix(NA, nrow = 1, ncol = ncol(df1) + ncol(df2))
		df1.spe <- df1
		df2.spe <- df2
		}
	else{
		df.ovl <- data.frame(df1[queryHits(gr.ovl),], df2[subjectHits(gr.ovl),])
		df1.spe <- df1[-queryHits(gr.ovl),]
		df2.spe <- df2[-subjectHits(gr.ovl),]
		}
	colnames(df.ovl) <- c(paste(colnames(df1), name1, sep = "."),
		paste(colnames(df2), name2, sep = "."))
	
	compL <- list(common = df.ovl, df1.spe = df1.spe, df2.spe = df2.spe)
	names(compL) <- c("common", paste0(name1, ".spe"), paste0(name2, ".spe"))
	
	return(compL)
}


filtering <- function(f, width = 2000)
{
	m <- read.xlsx(f, sheet = 1)
	m <- m[, 1:7]
	
	m <- m[order(m$chromosome, m$start, m$end), ]
	
	idx <- 0
	df.final <- c()
	for(chr in unique(m[,1]))
		{
		tempdf <- subset(m,  chromosome==chr)
		df.final <- rbind(df.final, tempdf[1,])
		idx <- idx + 1
		if(nrow(tempdf) > 1){
			for(i in 2:nrow(tempdf))
				{
				if(tempdf[i, "start"] <= (df.final[idx, "end"] + width))
					{
					df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] <-
						df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] + tempdf[i, c("read", "hits", "read.ctl", "hits.ctl")]# add read, collapseCluster...
					if(tempdf[i, "end"] > df.final[idx, "end"]) df.final[idx, "end"] <- tempdf[i, "end"]# change end if new end is higher
					
					}else{
						df.final <- rbind(df.final, tempdf[i, ])
						idx <- idx + 1
						}
				}
			}
		}
	return(df.final)
}





if(FALSE)
{


#######
# HBG1 196
f1 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Fwd_196/results/guide_aln",
	"Fwd-196-1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Fwd_196_rep2/results/guide_aln",
	"Fwd-196-2_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f3 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Rev_196/results/guide_aln",
	"Rev-196-1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f4 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Rev_196_rep2/results/guide_aln",
	"Rev-196-2_S4_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")	
	
# COMPARISON
outDir <- "/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/overlap_HBG1"
dir.create(outDir, showWarnings = FALSE)
setwd(file.path(outDir))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("Fwd196", "Fwd196_rep2", "Rev196", "Rev196_rep2"), width = 1000)
write.xlsx(compList, "HBG1_196_w250.xlsx")


#######
# HBG1 197
f1 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Fwd_197/results/guide_aln",
	"Fwd-197-1_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Fwd_197_rep2/results/guide_aln",
	"Fwd-197-2_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f3 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Rev_197/results/guide_aln",
	"Rev-197-1_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f4 <- file.path("/Users/gandrieux/cluster/cluster/offTargets/Giando/pipelineGit/samples/HBG1Rev_197_rep2/results/guide_aln",
	"Rev-197-2_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")	
	
# COMPARISON
outDir <- "/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/overlap_HBG1"
dir.create(outDir, showWarnings = FALSE)
setwd(file.path(outDir))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("Fwd197", "Fwd197_rep2", "Rev197", "Rev197_rep2"), width = 1000)
write.xlsx(compList, "HBG1_197_w250.xlsx")



##########################################################################################
# FINAL LIST (06.05.19)

#######
# G3 WT
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D4/results/guide_aln/",
	"G3-WT-d4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D4_rep2/results/guide_aln/",
	"G3-d4-WT_S6_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D14/results/guide_aln/",
	"G3-WT-d14_S3_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")	
	
# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("D1", "D4", "D4_rep2", "D14"), width = 1000)
write.xlsx(compList, "G3_WT_ovl_w250.xlsx")

#########
# G3 HiFi
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D1/results/guide_aln/",
	"G3-Hifi-d1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D4/results/guide_aln/",
	"G3-Hifi-d4_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D4_rep2/results/guide_aln/",
	"G3-d4-HF_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D14/results/guide_aln/",
	"G3-Hifi-d14_S6_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		

# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("D1", "D4", "D4_rep2", "D14"), width = 1000)
write.xlsx(compList, "G3_HiFi_ovl_w250.xlsx")

##########
# HiFi 399
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results/guide_aln/",
	"399HifiD1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results/guide_aln/",
	"399HifiD4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4_rep2/results/guide_aln/",
	"399-d4-HF_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D14/results/guide_aln/",
	"399HifiD14_S3_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")

# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("D1", "D4", "D4_rep2", "D14"), width = 2000)
write.xlsx(compList, "HiFi399_ovl_w250.xlsx")


########
# WT 399

f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D1/results/guide_aln/",
	"399WtD1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D4/results/guide_aln/",
	"399WtD4_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D4_rep2/results/guide_aln/",
	"399-d4-WT_S3_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D14/results/guide_aln/",
	"399WtD14_S6_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")

# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("D1", "D4", "D4_rep2", "D14"), width = 2000)
write.xlsx(compList, "WT399_ovl_w250.xlsx")

#############
# FANCF decoy
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/FANCF_decoy_rep1/results/guide_aln/",
	"FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/FANCF_decoy_rep2/results/guide_aln/",
	"FANCF-decoy_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2), c("rep1", "rep2"), width = 2000)
write.xlsx(compList, "FANCF_decoy_ovl_w250.xlsx")

#######
# VEGFA
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/VEGFA_decoy_rep1/results/guide_aln/",
	"VEGFA-withDecoyGel_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/VEGFA_decoy_rep2/results/guide_aln/",
	"VEGFA-decoy_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_SCORE.xlsx")
	
# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2), c("rep1", "rep2"), width = 2000)
write.xlsx(compList, "VEGFA_decoy_ovl_w250.xlsx")

########
# HBBTal

f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D1/results/guide_aln/",
	"HBBTalD1_S7_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D4/results/guide_aln/",
	"HBBTalD4_S8_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D4_rep2/results/guide_aln/",
	"HBB-D4-TAL_S8_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")	
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D14/results/guide_aln/",
	"HBB-Tal-d14_S8_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
		
# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2, f3, f4), c("D1", "D4", "D4_rep2", "D14"), width = 1000)
write.xlsx(compList, "HBBTal_ovl_w250.xlsx")


###################
# NOT FOR THE PAPER


# OT
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/OT/results/guide_aln/",
	"OT2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/OT6/results/guide_aln/",
	"OT6_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2), c("OT2", "OT6"), width = 2000)
write.xlsx(compList, "OT_ovl_w250.xlsx")


# TAL
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/TAL_rep1/results/guide_aln/",
	"Tal1_S3_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/TAL_rep2/results/guide_aln/",
	"Tal2_S4_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
	
# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
compList <- dfComparisonList(c(f1, f2), c("rep1", "rep2"), width = 2000)
write.xlsx(compList, "TAL_ovl_w250.xlsx")



##########################################################################################
# CASTSeq vs. CircleSeq

######################################
# FANCF
f1 <- file.path("~/cluster/master/offTargets/Giando/pipeline/FANCF_decoy_rep1/results/overlap_aln/",
	"FANCF_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")
	
f2 <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"circle_FANCF.U2OS.hg38.xlsx")

# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison"))
compList <- dfComparisonV2(f1, f2, "CASTseq", "CIRCLEseq")
write.xlsx(compList, "FANCF_decoy_CASTseq_vs_CIRCLEseq.xlsx")


######################################
# FANCF UNION
f1 <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando/",
	"FANCF_decoy_UNION.xlsx")
	
f2 <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"circle_FANCF.U2OS.hg38.xlsx")
	
# UNION
f1u <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando/",
	"FANCF_decoy_UNION_merged.xlsx")
write.xlsx(filtering(f1, 2000), f1u)

# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison"))
compList <- dfComparisonV2(f1u, f2, "CASTseq", "CIRCLEseq")
write.xlsx(compList, "FANCF_decoy_unfiltered_CASTseq_vs_CIRCLEseq.xlsx")


######################################
# VEGFA
f1 <- file.path("~/cluster/master/offTargets/Giando/pipeline/VEGFA_decoy_rep1/results/overlap_aln/",
	"VEGFA_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")
	
f2 <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"circle_VEGFA.hg38.xlsx")

# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison"))
compList <- dfComparisonV2(f1, f2, "CASTseq", "CIRCLEseq")
write.xlsx(compList, "VEGFA_decoy_CASTseq_vs_CIRCLEseq.xlsx")


######################################
# VEGFA
f1 <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando/",
	"VEGFA_decoy_UNION.xlsx")
	
f2 <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"circle_VEGFA.hg38.xlsx")
	
# UNION
f1u <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando/",
	"VEGFA_decoy_UNION_merged.xlsx")
write.xlsx(filtering(f1, 2000), f1u)


# COMPARISON
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison"))
compList <- dfComparisonV2(f1u, f2, "CASTseq", "CIRCLEseq")
write.xlsx(compList, "VEGFA_decoy_unfiltered_CASTseq_vs_CIRCLEseq.xlsx")



##########################################################################################
# CASTSeq vs. CircleSeq vs. GUIDESeq

# GUIDESeq hg19 to hg38
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/"))
guideM <- read.xlsx("GUIDE-seq tables.xlsx", sheet = 1)

vegfa <- guideM[guideM$Cells == "U2OS" & guideM$Targetsite == "VEGFA_site3", ]
colnames(vegfa)[c(1:3, 6)] <- c("V1", "V2", "V3", "strand")
vegfa.hg38 <- hg19to38(vegfa)

fancf <- guideM[guideM$Cells == "U2OS" & guideM$Targetsite == "FANCF", ]
colnames(fancf)[c(1:3, 6)] <- c("V1", "V2", "V3", "strand")
fancf.hg38 <- hg19to38(fancf)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/"))
write.xlsx(vegfa.hg38, "GUIDE-seq_VEGFA.hg38.xlsx")
write.xlsx(fancf.hg38, "GUIDE-seq_FANCF.hg38.xlsx")


#######
# FANCF
f.cast.filt <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/FANCF_decoy_rep1/results/overlap_aln/",
	"FANCF_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

#f.cast <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando/",
#	"FANCF_decoy_UNION_merged.xlsx")

f.cast <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519",
	"FANCF_decoy_UNION.xlsx")
	
f.circle <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"circle_FANCF.U2OS.hg38.xlsx")

f.guide <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"GUIDE-seq_FANCF.hg38.xlsx")


df.cast.filt <- read.xlsx(f.cast.filt)
gr.cast.filt <- makeGRangesFromDataFrame(df.cast.filt, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)

df.cast <- read.xlsx(f.cast)
gr.cast <- makeGRangesFromDataFrame(df.cast[, 1:3], seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)

df.circle <- read.xlsx(f.circle)
df.circle <- df.circle[, !grepl("strand", colnames(df.circle), ignore.case = TRUE)]
gr.circle <- makeGRangesFromDataFrame(df.circle, seqnames.field = "seqnames",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)
			
df.guide <- read.xlsx(f.guide)
df.guide <- df.guide[, !grepl("strand", colnames(df.guide), ignore.case = TRUE)]
gr.guide <- makeGRangesFromDataFrame(df.guide, seqnames.field = "seqnames",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)
			
ovlObj <- findOverlapsOfPeaks(gr.cast.filt, gr.cast, gr.circle, gr.guide, maxgap=2000)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/vs_GUIDE_vs_CIRCLE"))
ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))

openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "FANCF_CASTSeq_CIRCLESeq_GUIDESeq.xlsx")

pdf("FANCF_CASTSeq_CIRCLESeq_GUIDESeq_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

#######
# VEGFA
f.cast.filt <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/VEGFA_decoy_rep1/results/overlap_aln/",
	"VEGFA_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

#f.cast <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando/",
#	"VEGFA_decoy_UNION_merged.xlsx")

f.cast <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519/",
	"VEGFA_decoy_UNION.xlsx")
	
f.circle <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"circle_VEGFA.hg38.xlsx")

f.guide <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc/",
	"GUIDE-seq_VEGFA.hg38.xlsx")



df.cast.filt <- read.xlsx(f.cast.filt)
gr.cast.filt <- makeGRangesFromDataFrame(df.cast.filt, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)

df.cast <- read.xlsx(f.cast)
gr.cast <- makeGRangesFromDataFrame(df.cast[, 1:3], seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)
gr.cast <- unique(gr.cast)

df.circle <- read.xlsx(f.circle)
df.circle <- df.circle[, !grepl("strand", colnames(df.circle), ignore.case = TRUE)]
gr.circle <- makeGRangesFromDataFrame(df.circle, seqnames.field = "seqnames",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)
			
df.guide <- read.xlsx(f.guide)
df.guide <- df.guide[, !grepl("strand", colnames(df.guide), ignore.case = TRUE)]
gr.guide <- makeGRangesFromDataFrame(df.guide, seqnames.field = "seqnames",
			start.field = "start", end.field = "end",
			keep.extra.columns = TRUE, ignore.strand = TRUE)
			
ovlObj <- findOverlapsOfPeaks(gr.cast.filt, gr.cast, gr.circle, gr.guide, maxgap=2000)
makeVennDiagram(ovlObj)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/vs_GUIDE_vs_CIRCLE"))
ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))

openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "VEGFA_CASTSeq_CIRCLESeq_GUIDESeq.xlsx")

pdf("VEGFA_CASTSeq_CIRCLESeq_GUIDESeq_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

# Get peaks
toString(a$peakNames[[1]])

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))

openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "VEGFA_CASTSeq_CIRCLESeq_GUIDESeq.xlsx")

}

###########################################################
# MIN LENGTH 20, 28, 30, 32


# G3 WT D1
f.ml20 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results20/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml28 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results28/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml30 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results30/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml32 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results32/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

df.ml20 <- read.xlsx(f.ml20)
gr.ml20 <- makeGRangesFromDataFrame(df.ml20, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.ml28 <- read.xlsx(f.ml28)
gr.ml28 <- makeGRangesFromDataFrame(df.ml28, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.ml30 <- read.xlsx(f.ml30)
gr.ml30 <- makeGRangesFromDataFrame(df.ml30, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.ml32 <- read.xlsx(f.ml32)
gr.ml32 <- makeGRangesFromDataFrame(df.ml32, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)						


ovlObj <- findOverlapsOfPeaks(gr.ml20, gr.ml28, gr.ml30, gr.ml32, maxgap=0)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/ml"))
pdf("G3_WT_D1_minLength_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "G3_WT_D1_minLength.xlsx")



# HiFi399 D1
f.ml20 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results20/guide_aln/",
	"399HifiD1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml28 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results28/guide_aln/",
	"399HifiD1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml30 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results30/guide_aln/",
	"399HifiD1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml32 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results32/guide_aln/",
	"399HifiD1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

df.ml20 <- read.xlsx(f.ml20)
gr.ml20 <- makeGRangesFromDataFrame(df.ml20, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.ml28 <- read.xlsx(f.ml28)
gr.ml28 <- makeGRangesFromDataFrame(df.ml28, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.ml30 <- read.xlsx(f.ml30)
gr.ml30 <- makeGRangesFromDataFrame(df.ml30, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.ml32 <- read.xlsx(f.ml32)
gr.ml32 <- makeGRangesFromDataFrame(df.ml32, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)						


ovlObj <- findOverlapsOfPeaks(gr.ml20, gr.ml28, gr.ml30, gr.ml32, maxgap=0)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/ml"))
pdf("HiFi399_D1_minLength_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "HiFi399_D1_minLength.xlsx")




# HiFi399 D4
f.ml20 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results20/guide_aln/",
	"399HifiD4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml28 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results28/guide_aln/",
	"399HifiD4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml30 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results30/guide_aln/",
	"399HifiD4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.ml32 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results32/guide_aln/",
	"399HifiD4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

df.ml20 <- read.xlsx(f.ml20)
gr.ml20 <- makeGRangesFromDataFrame(df.ml20, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.ml28 <- read.xlsx(f.ml28)
gr.ml28 <- makeGRangesFromDataFrame(df.ml28, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.ml30 <- read.xlsx(f.ml30)
gr.ml30 <- makeGRangesFromDataFrame(df.ml30, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.ml32 <- read.xlsx(f.ml32)
gr.ml32 <- makeGRangesFromDataFrame(df.ml32, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)						


ovlObj <- findOverlapsOfPeaks(gr.ml20, gr.ml28, gr.ml30, gr.ml32, maxgap=0)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/ml"))
pdf("HiFi399_D4_minLength_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "HiFi399_D4_minLength.xlsx")




##################################################################
# FANCF DECOY REP1 FOR, REV, BOTH

f.for <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/FANCF_decoy_rep1/results_FOR/guide_aln/",
	"FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.rev <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/FANCF_decoy_rep1/results_REV/guide_aln/",
	"FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.both <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/FANCF_decoy_rep1/results_both/guide_aln/",
	"FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

df.for <- read.xlsx(f.for)
gr.for <- makeGRangesFromDataFrame(df.for, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.rev <- read.xlsx(f.rev)
gr.rev <- makeGRangesFromDataFrame(df.rev, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.both <- read.xlsx(f.both)
gr.both <- makeGRangesFromDataFrame(df.both, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			

ovlObj <- findOverlapsOfPeaks(gr.for, gr.rev, gr.both, maxgap=0)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/linker"))
pdf("FANCF_decoy_rep1_LINKER_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "FANCF_decoy_rep1_LINKER.xlsx")


# G3_WT_D1 FOR, REV, BOTH

f.for <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results_FOR/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.rev <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results_REV/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.both <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results_both/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

df.for <- read.xlsx(f.for)
gr.for <- makeGRangesFromDataFrame(df.for, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.rev <- read.xlsx(f.rev)
gr.rev <- makeGRangesFromDataFrame(df.rev, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
df.both <- read.xlsx(f.both)
gr.both <- makeGRangesFromDataFrame(df.both, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			

ovlObj <- findOverlapsOfPeaks(gr.for, gr.rev, gr.both, maxgap=0)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/linker"))
pdf("G3_WT_D1_LINKER_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "G3_WT_D1_LINKER.xlsx")


##################################################################
# G3 Rev 1 vs. G3 For 1

g3r <- file.path("~/Research/CASTSeq/revision/",
                   "G3Rev_WT_D1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
g3r1 <- file.path("~/Research/CASTSeq/revision/REVd1WT1_S19_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
g3f <- file.path("~/Research/CASTSeq/revision/G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

df.g3r <- read.xlsx(g3r)
df.g3r <- df.g3r[df.g3r$adj.pvalue < 0.05, ]
df.g3r <- df.g3r[df.g3r$group != "CBS", ]
gr.g3r <- makeGRangesFromDataFrame(df.g3r, seqnames.field = "chromosome",
                                   start.field = "start", end.field = "end",
                                   keep.extra.columns = FALSE, ignore.strand = TRUE)

df.g3r1 <- read.xlsx(g3r1)
df.g3r1 <- df.g3r1[df.g3r1$adj.pvalue < 0.05, ]
df.g3r1 <- df.g3r1[df.g3r1$group != "CBS", ]
gr.g3r1 <- makeGRangesFromDataFrame(df.g3r1, seqnames.field = "chromosome",
                                   start.field = "start", end.field = "end",
                                   keep.extra.columns = FALSE, ignore.strand = TRUE)

df.g3f <- read.xlsx(g3f)
df.g3f <- df.g3f[df.g3f$adj.pvalue < 0.05, ]
df.g3f <- df.g3f[df.g3f$group != "CBS", ]
gr.g3f <- makeGRangesFromDataFrame(df.g3f, seqnames.field = "chromosome",
                                    start.field = "start", end.field = "end",
                                    keep.extra.columns = FALSE, ignore.strand = TRUE)


ovlObj <- findOverlapsOfPeaks(gr.g3r, gr.g3r1, gr.g3f, maxgap=0)

setwd(file.path("~/Research/CASTSeq/revision/overlap"))
pdf("G3_F_vs_R_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "G3_F_vs_R.xlsx")


##################################################################
# G3 Rev vs. For NEW PIPELINE 20.05.20

siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D1/results/guide_aln/G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D4/results/guide_aln/G3-WT-d4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D4_rep2/results/guide_aln/G3-d4-WT_S6_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D14/results/guide_aln/G3-WT-d14_S3_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D26/results/guide_aln/G3_UT_D26_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3Rev_WT_D1/results/guide_aln/G3Rev_WT_D1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3Rev_WT_D14/results/guide_aln/G3Rev_WT_D14_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3Rev_WT_D26/results/guide_aln/G3Rev_WT_D26_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
               )
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "G3_For_Rev_280520_NEW.xlsx")

makeUpset("G3_For_Rev_280520_NEW.xlsx")

# FOR ONLY
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D1/results/guide_aln/G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D4/results/guide_aln/G3-WT-d4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D4_rep2/results/guide_aln/G3-d4-WT_S6_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D14/results/guide_aln/G3-WT-d14_S3_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
               #file.path("~/Research/CASTSeq/pipelineGit/samples/G3_WT_D26/results/guide_aln/G3_UT_D26_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "G3_For_280520.xlsx")

makeUpset("G3_For_280520.xlsx")

# REV ONLY
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/G3Rev_WT_D1/results/guide_aln/G3Rev_WT_D1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/G3Rev_WT_D14/results/guide_aln/G3Rev_WT_D14_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
               #file.path("~/Research/CASTSeq/pipelineGit/samples/G3Rev_WT_D26/results/guide_aln/G3Rev_WT_D26_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "G3_Rev_280520.xlsx")

makeUpset("G3_Rev_280520.xlsx")

################
# CCR5 SITE2 FOR
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D1/results_1500/guide_aln/399WtD1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D4/results_1500/guide_aln/399WtD4_S5_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D4_rep2/results_1500/guide_aln/399-d4-WT_S3_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D14/results_1500/guide_aln/399WtD14_S6_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               #file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D26/results_1500/guide_aln/399WTd26_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/CCR5_2907_TN_UTN/results1500/guide_aln/CCR5-2907-TN_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
               )
#names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])
names(siteFiles) <- c("G3_399_D1", "G3_399_D4", "G3_399_D4rep2", "G3_399_D14", "CCR5_2907_TN_UTN")


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "CCR5_2_FOR_180920.xlsx")

makeUpset("CCR5_2_FOR_180920.xlsx")

################
# CCR5 SITE2 REV
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/CCR5_0808_TO_UTO/results1500/guide_aln/CCR5-08-08-TO_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/CCR5_2907_TO_UTO/results1500/guide_aln/CCR5-2907-TO_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
#names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])
names(siteFiles) <- c("CCR5_0808_TO_UTO", "CCR5_2907_TO_UTO")

setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "CCR5_2_REV_180920.xlsx")

makeUpset("CCR5_2_REV_180920.xlsx")


########################
# CCR5 SITE2 FOR AND REV
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D1/results_1500/guide_aln/399WtD1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D4/results_1500/guide_aln/399WtD4_S5_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D4_rep2/results_1500/guide_aln/399-d4-WT_S3_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D14/results_1500/guide_aln/399WtD14_S6_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               #file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/WT399_D26/results_1500/guide_aln/399WTd26_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/CCR5_2907_TN_UTN/results1500/guide_aln/CCR5-2907-TN_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/CCR5_0808_TO_UTO/results1500/guide_aln/CCR5-08-08-TO_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/CCR5_2907_TO_UTO/results1500/guide_aln/CCR5-2907-TO_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
#names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])
names(siteFiles) <- c("G3_399_D1", "G3_399_D4", "G3_399_D4rep2", "G3_399_D14",
                      "CCR5_2907_TN_UTN", "CCR5_0808_TO_UTO", "CCR5_2907_TO_UTO")


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "CCR5_2_FOR_REV_180920.xlsx")

makeUpset("CCR5_2_FOR_REV_180920.xlsx")







# RAG1A SET1 T1 to T4
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T1_UT1/results/guide_aln/RAG1-A-set-1-T1_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T2_UT2/results/guide_aln/RAG1-A-set-1-T2_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T3_UT3/results/guide_aln/RAG1-A-set-1-T3_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T4_UT4/results/guide_aln/RAG1-A-set-1-T4_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "RAG1A_set1_T1T4.xlsx")

makeUpset("RAG1A_set1_T1T4.xlsx")


# RAG1A SET1 T1 to T4 AND set2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T1_UT1/results/guide_aln/RAG1-A-set-1-T1_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T2_UT2/results/guide_aln/RAG1-A-set-1-T2_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T3_UT3/results/guide_aln/RAG1-A-set-1-T3_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set1_T4_UT4/results/guide_aln/RAG1-A-set-1-T4_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_set2_T1_UT1/results/guide_aln/RAG1-A-set-2-T1_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])


setwd(file.path("~/Research/CASTSeq/revision/overlap/"))
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0)
write.xlsx(compList, "RAG1A_set1_T1T4_set2_T1.xlsx")

makeUpset("RAG1A_set1_T1T4_set2_T1.xlsx")


siteList <- lapply(siteFiles, read.xlsx)
siteList <- lapply(siteList, function(i) i[i$adj.pvalue < 0.05, ])
siteList <- lapply(siteList, function(i) i[i$group != "CBS", ])

siteList.gr <- lapply(siteList, makeGRangesFromDataFrame, seqnames.field = "chromosome",
                                   start.field = "start", end.field = "end",
                                   keep.extra.columns = FALSE, ignore.strand = TRUE)

F1 <-  siteList.gr[["G3_WT_D1"]]
F4 <-  siteList.gr[["G3_WT_D4"]]
F14 <-  siteList.gr[["G3_WT_D14"]]
F26 <-  siteList.gr[["G3_WT_D26"]]

R1 <-  siteList.gr[["G3Rev_WT_D1"]]
R14 <-  siteList.gr[["G3Rev_WT_D14"]]
R26 <-  siteList.gr[["G3Rev_WT_D26"]]

# G3 For
ovlObj <- findOverlapsOfPeaks(F1, F4, F14, F26,
                              maxgap=250)

setwd(file.path("~/Research/CASTSeq/revision/overlap"))
pdf("G3_FOR_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "G3_FOR_Venn.xlsx")

# G3 REV
ovlObj <- findOverlapsOfPeaks(R1, R14, R26,
                              maxgap=250)

setwd(file.path("~/Research/CASTSeq/revision/overlap"))
pdf("G3_REV_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "G3_REV_Venn.xlsx")


# G3 FOR and REV
ovlObj <- findOverlapsOfPeaks(F1, R1, F14, R14,
                              maxgap=250)

setwd(file.path("~/Research/CASTSeq/revision/overlap"))
pdf("G3_FOR_REV_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "G3_FOR_REV_Venn.xlsx")




#############################################################
# HBG1Rev mm0, mm2

f.mm0 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBG1Rev/results_0/guide_aln/",
	"HBG1Rev_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.mm2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBG1Rev/results_2/guide_aln/",
	"HBG1Rev_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")


df.mm0 <- read.xlsx(f.mm0)
gr.mm0 <- makeGRangesFromDataFrame(df.mm0, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.mm2 <- read.xlsx(f.mm2)
gr.mm2 <- makeGRangesFromDataFrame(df.mm2, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
ovlObj <- findOverlapsOfPeaks(gr.mm0, gr.mm2, maxgap=0)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/mismatches"))
pdf("HBG1Rev_Mismatches_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "HBG1Rev_Mismatches.xlsx")


#######################################################################
# RAG1A and RAG1B 05.04.20

f.RAG1A.r1 <- file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A/results_k25/guide_aln/UT-RAG1A_S1_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.RAG1A.r2 <- file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1A_040520/results/guide_aln/Treated-RAG1A_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.known <- file.path("~/Research/CASTSeq/RAG1A and B off-targets HTGTS paper.xlsx")


df.RAG1A.r1 <- read.xlsx(f.RAG1A.r1, sheet = 1)
df.RAG1A.r1 <- df.RAG1A.r1[df.RAG1A.r1$adj.pvalue < 0.05, ]
df.RAG1A.r1 <- df.RAG1A.r1[df.RAG1A.r1$group != "CBS", ]

gr.RAG1A.r1 <- makeGRangesFromDataFrame(df.RAG1A.r1, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)

df.RAG1A.r2 <- read.xlsx(f.RAG1A.r2, sheet = 1)
df.RAG1A.r2 <- df.RAG1A.r2[df.RAG1A.r2$adj.pvalue < 0.05, ]
df.RAG1A.r2 <- df.RAG1A.r2[df.RAG1A.r2$group != "CBS", ]

gr.RAG1A.r2 <- makeGRangesFromDataFrame(df.RAG1A.r2, seqnames.field = "chromosome",
                                        start.field = "start", end.field = "end",
                                        keep.extra.columns = FALSE, ignore.strand = TRUE)

df.known <- read.xlsx(f.known, sheet = "RAG1A")
df.known <- hg19to38(data.frame(V1 = df.known$Chr, V2 = df.known$Start, V3 = df.known$End))
gr.known <- makeGRangesFromDataFrame(df.known, seqnames.field = "seqnames",
                                        start.field = "start", end.field = "end",
                                        keep.extra.columns = FALSE, ignore.strand = TRUE)
				
ovlObj <- findOverlapsOfPeaks(gr.RAG1A.r1, gr.RAG1A.r2, gr.known, maxgap=0)

setwd(file.path("~/Research/CASTSeq/overlap/"))
pdf("RAG1A_woCBS_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "RAG1A_woCBS_overlap.xlsx")


#########

f.RAG1B <- file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1B_040520/results/guide_aln/Treated-RAG1B_S4_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.RAG1B_rep2 <- file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1B_210520/results/guide_aln/T-jk_S2_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

f.known <- file.path("~/Research/CASTSeq/RAG1A and B off-targets HTGTS paper.xlsx")

df.RAG1B <- read.xlsx(f.RAG1B, sheet = 1)
df.RAG1B  <- df.RAG1B[df.RAG1B$adj.pvalue < 0.05, ]
df.RAG1B <- df.RAG1B[df.RAG1B$group != "CBS", ]

gr.RAG1B <- makeGRangesFromDataFrame(df.RAG1B, seqnames.field = "chromosome",
                                        start.field = "start", end.field = "end",
                                        keep.extra.columns = FALSE, ignore.strand = TRUE)

df.RAG1B_rep2 <- read.xlsx(f.RAG1B_rep2, sheet = 1)
df.RAG1B_rep2  <- df.RAG1B_rep2[df.RAG1B_rep2$adj.pvalue < 0.05, ]
df.RAG1B_rep2 <- df.RAG1B_rep2[df.RAG1B_rep2$group != "CBS", ]

gr.RAG1B_rep2 <- makeGRangesFromDataFrame(df.RAG1B_rep2, seqnames.field = "chromosome",
                                     start.field = "start", end.field = "end",
                                     keep.extra.columns = FALSE, ignore.strand = TRUE)


df.known <- read.xlsx(f.known, sheet = "RAG1B")
df.known <- hg19to38(data.frame(V1 = df.known$Chr, V2 = df.known$Start, V3 = df.known$End))
gr.known <- makeGRangesFromDataFrame(df.known, seqnames.field = "seqnames",
                                     start.field = "start", end.field = "end",
                                     keep.extra.columns = FALSE, ignore.strand = TRUE)

ovlObj <- findOverlapsOfPeaks(gr.RAG1B, gr.RAG1B_rep2, gr.known, maxgap=0)

setwd(file.path("~/Research/CASTSeq/overlap/"))
pdf("RAG1B_woCBS_Venn_V2.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "RAG1B_woCBS_overlap_V2.xlsx")



#####################################################
# UNC13D

f.UNC13D <- file.path("~/Research/CASTSeq/pipelineGit/samples/UNC13D/results/guide_aln/Treated-UNC13D_S4_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
f.UNC13D_rep2 <- file.path("~/Research/CASTSeq/pipelineGit/samples/UNC13D_260520/results/guide_aln/T-kay_S4_L001_w250_aln_stat_FLANK_GROUP_GENES.xlsx")

siteFiles <- c(UNC13D = f.UNC13D,
               UNC13D_rep2 = f.UNC13D_rep2)

# read files
#names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])

siteList <- lapply(siteFiles, read.xlsx)
siteList <- lapply(siteList, function(i) i[i$adj.pvalue < 0.05, ])
siteList <- lapply(siteList, function(i) i[i$group != "CBS", ])

siteList.gr <- lapply(siteList, makeGRangesFromDataFrame, seqnames.field = "chromosome",
                      start.field = "start", end.field = "end",
                      keep.extra.columns = FALSE, ignore.strand = TRUE)

UNC13D  <-  siteList.gr[["UNC13D"]]
UNC13D_rep2 <-  siteList.gr[["UNC13D_rep2"]]


ovlObj <- findOverlapsOfPeaks(UNC13D, UNC13D_rep2, maxgap=1000)

setwd(file.path("~/Research/CASTSeq/overlap/"))
pdf("UNC13D_woCBS_qv0.05_Venn.pdf")
makeVennDiagram(ovlObj)
dev.off()

ovlPeaks <- ovlObj$peaklist
ovlPeaks.df <- lapply(ovlPeaks, as.data.frame)
ovlPeaks.df <- lapply(ovlPeaks.df, function(i) i[, 1:3])
names(ovlPeaks.df) <- gsub("///", "_", names(ovlPeaks.df))
names(ovlPeaks.df) <- gsub("gr.", "", names(ovlPeaks.df))

nbMat <- data.frame(VennGroup = names(ovlPeaks.df), NB = unlist(lapply(ovlPeaks.df, nrow)))
openxlsx::write.xlsx(c(list(Overview = nbMat), ovlPeaks.df), "UNC13D_woCBS_qv0.05_overlap.xlsx")








