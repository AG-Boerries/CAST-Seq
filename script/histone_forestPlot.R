#library(ggplot2)
#library(openxlsx)
#library(ChIPseeker)
#library(clusterProfiler)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

annotateExonIntron <- function(m)
{
	rawAnn <- m$annotation
	newAnn <- rawAnn
	
	firstExon <- grep("exon 1 of", rawAnn)    
	firstIntron <- grep("intron 1 of", rawAnn)
	exon <- grep("^Exon ", rawAnn) 
	intron <- grep("^Intron ", rawAnn)
	
	newAnn[exon] <- "other Exon"
	newAnn[intron] <- "other Intron"
	
	newAnn[firstExon] <- "first Exon"
	newAnn[firstIntron] <- "first Intron"
	
	m$annotation <- newAnn
	return(m)
}


getContagency <- function(test, ref, category)
{
	test.sum <- sum(test)
	test.pos <- test[category]
	if(is.na(test.pos)) test.pos <- 0
	test.neg <- test.sum - test.pos
	
	ref.sum <- sum(ref)
	ref.pos <- ref[category]
	ref.neg <- ref.sum - ref.pos

	return(matrix(c(test.pos, ref.pos, test.neg, ref.neg), nrow = 2, byrow = TRUE))
}

getOddMatrix <- function(test, ref)
{
	groups <- names(ref)
	OddMatrix <- matrix(NA, nrow = length(groups), ncol = 5)

	for(i in 1:length(groups))
		{
		ctg <- getContagency(test, ref, groups[i])
		fh <- fisher.test(ctg)
		pv <- fh$p.value
		odd <- fh$estimate["odds ratio"]
		confInt.down <- fh$conf.int[1]
		confInt.up <- fh$conf.int[2]
		
		OddMatrix[i, -1] <- c(pv, odd, confInt.down, confInt.up)
		}
	OddMatrix[,1] <- groups
	colnames(OddMatrix) <- c("Region", "P.value", "Odd", "Conf.int.low", "Conf.int.high")
	
	OddMatrix <- as.data.frame(OddMatrix)
	OddMatrix[,"P.value"] <- toNum(OddMatrix[,"P.value"])
	OddMatrix[,"Odd"] <- toNum(OddMatrix[,"Odd"])
	OddMatrix[,"Conf.int.low"] <- toNum(OddMatrix[,"Conf.int.low"])
	OddMatrix[,"Conf.int.high"] <- toNum(OddMatrix[,"Conf.int.high"])
	
	return(OddMatrix)
}

#toNum <- function(x) as.numeric(levels(x))[x]


getNbOverlap <- function(query, subject)
{
	gr.ovl <- findOverlaps(query = query, subject = subject, type = "any")
	hits <- queryHits(gr.ovl)
	return(length(unique(hits)))
}

getOverlap <- function(query, subject)
{
	gr.ovl <- findOverlaps(query = query, subject = subject, type = "any")
	hits <- queryHits(gr.ovl)
	return(hits)
}

getOddMatrixEnhancer <- function(contagList)
{
	groups <- names(contagList)
	OddMatrix <- matrix(NA, nrow = length(groups), ncol = 5)
	
	for(i in 1:length(groups))
		{
		ctg <- contagList[[i]]
		fh <- fisher.test(ctg)
		pv <- fh$p.value
		odd <- fh$estimate["odds ratio"]
		confInt.down <- fh$conf.int[1]
		confInt.up <- fh$conf.int[2]
		
		OddMatrix[i, -1] <- c(pv, odd, confInt.down, confInt.up)		
		}
	OddMatrix[,1] <- groups
	colnames(OddMatrix) <- c("Region", "P.value", "Odd", "Conf.int.low", "Conf.int.high")	
	
	OddMatrix <- as.data.frame(OddMatrix)
	OddMatrix[,"P.value"] <- toNum(OddMatrix[,"P.value"])
	OddMatrix[,"Odd"] <- toNum(OddMatrix[,"Odd"])
	OddMatrix[,"Conf.int.low"] <- toNum(OddMatrix[,"Conf.int.low"])
	OddMatrix[,"Conf.int.high"] <- toNum(OddMatrix[,"Conf.int.high"])
	
	return(OddMatrix)
}

getContagencyHistone <- function(test.pos, test.tot, ref.pos, ref.tot)
{
	test.neg <- test.tot - test.pos
	ref.neg <- ref.tot - ref.pos

	return(matrix(c(test.pos, ref.pos, test.neg, ref.neg), nrow = 2, byrow = TRUE))
}


histoneForestPlot <- function(inputF, randomF, histFiles)
{
	readMat <- read.xlsx(inputF, sheet = 1)

	groups <- unique(readMat$group)
	readMatList <- lapply(groups, function(i)
		readMat[readMat$group == i, ]
		)
	names(readMatList) <- groups

	readMatList.gr <- lapply(1:length(readMatList), function(i)
		makeGRangesFromDataFrame(readMatList[[i]][, c("chromosome", "start", "end", "alignment.id")],
		keep.extra.columns = TRUE, ignore.strand =TRUE)
		)
	names(readMatList.gr) <- groups	
	
	# LOAD RANDOM SEQUENCES (CONTROL)	
	# Random bed
	rdMat <- read.delim(randomF, header = FALSE)
	rdMat.gr <- makeGRangesFromDataFrame(rdMat,seqnames.field="V1", start.field = "V2",
		end.field = "V3", strand.field = "V6",
		keep.extra.columns = FALSE)

	# LOAD HISTONES
	#names(histFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histFiles))	
	histList.raw <- lapply(histFiles, read.delim, header = FALSE)

	lengthList <- c(500, 1000, 10000)

	ggmatList <- lapply(lengthList, function(l){
		myLength <- l
		histList <- list()
		for(i in 1:length(histList.raw))
			{
			hst <- histList.raw[[i]]
			hst.middle <- round((hst$V3 + hst$V2) / 2)
			histList[[i]] <- makeGRangesFromDataFrame(data.frame(chr = hst$V1, Start = hst.middle - myLength, End = hst.middle + myLength))
			}
		#names(histList) <- paste0(names(histFiles), "_+/-", myLength)
		names(histList) <- names(histFiles)

		# OVERLAP

		# Real peaks
		ovlList <- lapply(readMatList.gr, function(i)
			lapply(histList, function(j) getNbOverlap(query = i, subject = j))
			)

		# Random peaks	
		ovl.rd <- lapply(histList, getNbOverlap, query = rdMat.gr)	

		# CONTAGENCY
		contagList <- lapply(1:length(ovlList), function(i){
			ovl.name <- names(ovlList)[i]
			nbTest <- nrow(readMatList[[ovl.name]])

			lapply(1:length(ovlList[[i]]), function(j){
				getContagencyHistone(ovlList[[i]][[j]], nbTest, ovl.rd[[j]], nrow(rdMat))
				})
			})
		names(contagList) <- names(ovlList)

		for(i in 1:length(contagList))
			{
			names(contagList[[i]]) <- names(histList) 
			}

		# ODD MATRIX
		oddList <- lapply(contagList, getOddMatrixEnhancer)

		for(i in 1:length(oddList))
			{
			oddList[[i]]$Group <- names(oddList)[i]
			}
		oddMat <- do.call(rbind, oddList)
		oddMat$size <- myLength
		return(oddMat)
		})



	# PLOT
	ggmat <- do.call(rbind, ggmatList)

	# reverses the factor level ordering for labels after coord_flip()
	ggmat$Region <- factor(ggmat$Region, levels=rev(sort(unique(ggmat$Region))))

	fp <- ggplot(data=ggmat, aes(x=Region, y=Odd, ymin=Conf.int.low, ymax=Conf.int.high, color = Group)) +
			geom_pointrange(position = position_dodge(width = 0.5)) + 
			geom_hline(yintercept = 1, lty = 2) +  # add a dotted line at x=1 after flip
			coord_flip() +  # flip coordinates (puts labels on y axis)
			xlab("Label") + ylab("Odds ratio (95% CI)") +
			theme_bw()  # use a white background
	fp <- fp + facet_wrap(~size, ncol = length(lengthList))


	pdf(gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_histones_odd_ratio.pdf", inputF), width = 12)
	print(fp)
	dev.off()

	write.xlsx(ggmat,
		gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_histones_odd_ratio.xlsx", inputF), row.names = FALSE, overwrite = TRUE)
}


if(FALSE)
{

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################



# LOAD TRANSLOCATION SITES
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment"))
sampleName <- "G3_W250"
setwd(sampleName)
readMat <- read.xlsx(paste0(sampleName, "_aln_stat_FLANK_GROUP_GENES.xlsx"), sheet = 1)

groups <- unique(readMat$group)
readMatList <- lapply(groups, function(i)
	readMat[readMat$group == i, ]
	)
names(readMatList) <- groups

readMatList.gr <- lapply(1:length(readMatList), function(i)
	makeGRangesFromDataFrame(readMatList[[i]][, c("chromosome", "start", "end", "strand", "alignment.id")],
	keep.extra.columns = TRUE)
	)
names(readMatList.gr) <- groups	
	
# LOAD RANDOM SEQUENCES (CONTROL)	
# Random bed
setwd(file.path("/Volumes/home/Geoffroy/offTargets/Giando/randomBed"))
rdMat <- read.delim("random_hg38_w250.bed", header = FALSE)
rdMat.gr <- makeGRangesFromDataFrame(rdMat,seqnames.field="V1", start.field = "V2",
	end.field = "V3", strand.field = "V6",
	keep.extra.columns = FALSE)






##########################################################################################
# HISTONES

###############
# LOAD HISTONES
setwd("/Volumes/Home/Geoffroy/database/epigenomics_roadmap/E035")
histFiles <- c("E035-H3K36me3.broadPeak_hg38_homer.bed", "E035-H3K4me3.broadPeak_hg38_homer.bed",
	"E035-H3K27me3.broadPeak_hg38_homer.bed")
names(histFiles) <- gsub(".broadPeak_hg38_homer.bed", "", histFiles)	

histList.raw <- lapply(histFiles, read.delim, header = FALSE)

lengthList <- c(500, 1000, 10000)

ggmatList <- lapply(lengthList, function(l){
	myLength <- l
	histList <- list()
	for(i in 1:length(histList.raw))
		{
		hst <- histList.raw[[i]]
		hst.middle <- round((hst$V3 + hst$V2) / 2)
		histList[[i]] <- makeGRangesFromDataFrame(data.frame(chr = hst$V1, Start = hst.middle - myLength, End = hst.middle + myLength))
		}
	#names(histList) <- paste0(names(histFiles), "_+/-", myLength)
	names(histList) <- names(histFiles)

	#########
	# OVERLAP

	# Real peaks
	ovlList <- lapply(readMatList.gr, function(i)
		lapply(histList, function(j) getNbOverlap(query = i, subject = j))
		)

	# Random peaks	
	ovl.rd <- lapply(histList, getNbOverlap, query = rdMat.gr)	


	############
	# CONTAGENCY
	contagList <- lapply(1:length(ovlList), function(i){
		ovl.name <- names(ovlList)[i]
		nbTest <- nrow(readMatList[[ovl.name]])

		lapply(1:length(ovlList[[i]]), function(j){
			getContagencyHistone(ovlList[[i]][[j]], nbTest, ovl.rd[[j]], nrow(rdMat))
			})
		})
	names(contagList) <- names(ovlList)

	for(i in 1:length(contagList))
		{
		names(contagList[[i]]) <- names(histList) 
		}


	############
	# ODD MATRIX
	oddList <- lapply(contagList, getOddMatrixEnhancer)

	for(i in 1:length(oddList))
		{
		oddList[[i]]$Group <- names(oddList)[i]
		}
	oddMat <- do.call(rbind, oddList)
	oddMat$size <- myLength
	return(oddMat)
	})


######
# PLOT
ggmat <- do.call(rbind, ggmatList)

# reverses the factor level ordering for labels after coord_flip()
ggmat$Region <- factor(ggmat$Region, levels=rev(sort(unique(ggmat$Region))))

fp <- ggplot(data=ggmat, aes(x=Region, y=Odd, ymin=Conf.int.low, ymax=Conf.int.high, color = Group)) +
        geom_pointrange(position = position_dodge(width = 0.5)) + 
        geom_hline(yintercept = 1, lty = 2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Label") + ylab("Odds ratio (95% CI)") +
        theme_bw()  # use a white background
fp <- fp + facet_wrap(~size, ncol = length(lengthList))


setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment"))
setwd(sampleName)

pdf(paste0(sampleName, "_histones_odd_ratio.pdf"), width = 12)
print(fp)
dev.off()

write.xlsx(ggmat,
	paste0(sampleName, "_histones_odd_ratio.xlsx"), row.names = FALSE)


}




