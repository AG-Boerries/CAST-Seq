#library(openxlsx)
#library(ChIPseeker)
#library(clusterProfiler)
#library(GenomicFeatures)
#library(rtracklayer)
#library(org.Hs.eg.db)
#library(ggplot2)

#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

toNum <- function(x) as.numeric(levels(x))[x]

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


annotateShort <- function(m)
{
	rawAnn <- m$annotation
	newAnn <- rawAnn
	
	body <- grep("Exon|Intron", rawAnn)    
	promoter <- grep("Promoter", rawAnn)
	downstream <- grep("Downstream", rawAnn)
	
	newAnn[body] <- "Body"
	newAnn[promoter] <- "Promoter"
	newAnn[downstream] <- "Downstream"

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


annotateGene <- function(inputFile)
{
	# LOAD READMAT
	readMat <- read.xlsx(inputFile, sheet = 1)

	# CHANGE START AND END INTO gRNA ALIGNMENT CORRDINATES
	readMat.middle <-readMat
	#readMat.middle$start <- readMat.middle$middleCoord
	#readMat.middle$end <- readMat.middle$middleCoord

	readMat.gr <- makeGRangesFromDataFrame(readMat.middle[, c("chromosome", "start", "end", "alignment.id")],
		keep.extra.columns = TRUE, ignore.strand = TRUE)

	# ANNOTATE PEAKS (chipseeker)
	readMat.anno.gr <- annotatePeak(readMat.gr, tssRegion=c(-3000, 3000), 
							 TxDb=TXDB, annoDb=ORG.STR, overlap = "all", level = "gene")

	# Save annotation
	readMat.anno <-  as.data.frame(readMat.anno.gr)

	# MERGE READMAT AND ANNOTATION
	rownames(readMat) <- readMat$alignment.id
	rownames(readMat.anno) <- readMat.anno$alignment.id

	readMat.merge <- merge(readMat, readMat.anno[,-c(1:6)], by=0, all=TRUE)
	readMat.merge <- readMat.merge[match(readMat$alignment.id, readMat.merge$alignment.id), ]
	readMat.merge <- readMat.merge[, -1]# remove rownames

	write.xlsx(readMat.merge, gsub(".xlsx", "_GENES.xlsx", inputFile), row.names = FALSE)
}

annotateGeneTALEN <- function(inputFile)
{
	# LOAD READMAT
	readMat <- read.xlsx(inputFile, sheet = 1)

	# CHANGE START AND END INTO gRNA ALIGNMENT CORRDINATES
	readMat.middle <- readMat
	#readMat.middle$start <- readMat.middle$start + (readMat.middle$end - readMat.middle$start)
	#readMat.middle$end <- readMat.middle$start + (readMat.middle$end - readMat.middle$start)

	readMat.gr <- makeGRangesFromDataFrame(readMat.middle[, c("chromosome", "start", "end", "alignment.id")],
		keep.extra.columns = TRUE, ignore.strand = TRUE)

	# ANNOTATE PEAKS (chipseeker)
	readMat.anno.gr <- annotatePeak(readMat.gr, tssRegion=c(-3000, 3000), 
							 TxDb=TXDB, annoDb=ORG.STR, overlap = "all", level = "gene")

	# Save annotation
	readMat.anno <-  as.data.frame(readMat.anno.gr)

	# MERGE READMAT AND ANNOTATION
	rownames(readMat) <- readMat$alignment.id
	rownames(readMat.anno) <- readMat.anno$alignment.id

	readMat.merge <- merge(readMat, readMat.anno[,-c(1:6)], by=0, all=TRUE)
	readMat.merge <- readMat.merge[match(readMat$alignment.id, readMat.merge$alignment.id), ]
	readMat.merge <- readMat.merge[, -1]# remove rownames

	write.xlsx(readMat.merge, gsub(".xlsx", "_GENES.xlsx", inputFile), row.names = FALSE)
}



geneBarplot <- function(inputFile)
{
	# Read inputFile
	readMat.merge <- read.xlsx(inputFile, sheet = 1)
	readMat.merge <- annotateExonIntron(readMat.merge)

	groups <- unique(readMat.merge$group)
	peakAnnoList <- lapply(groups, function(i)
		readMat.merge[readMat.merge$group == i, ]
		)
	names(peakAnnoList) <- groups	
		
	# Plot
	ggmat <- lapply(1:length(peakAnnoList), function(i)
		data.frame(Regions = names(table(peakAnnoList[[i]]$annotation)),
				   Percentage = table(peakAnnoList[[i]]$annotation) / nrow(peakAnnoList[[i]]) * 100,
				   Group = names(peakAnnoList)[i])
		)
	ggmat <- do.call(rbind, ggmat)	
	ggmat$Regions <- factor(ggmat$Regions, levels = rev(c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)",
		"5' UTR", "first Exon", "other Exon", "first Intron", "other Intron", "3' UTR",
		"Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)",
		"Distal Intergenic"))
		)
	ggmat$Group <- factor(ggmat$Group, levels = rev(c("off.target", "hom.recomb", "CBS"))
		)		
			
	p <- ggplot(data=ggmat, aes(x=Regions, y=Percentage.Freq, fill=Group))
	p <- p + geom_bar(stat="identity", position=position_dodge())
	p <- p + theme_bw()	
	p <- p + coord_flip()

	pdf(gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_region_barplot.pdf", inputFile))
	#png(gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_region_barplot.png", inputFile), units="px", width=1600, height=1600, res=300)
	plot(p)
	dev.off()
}




geneForestPlot <- function(inputFile, randomFile)
{
	# Load input file
	readMat.merge <- read.xlsx(inputFile, sheet = 1)
	readMat.merge <- annotateExonIntron(readMat.merge)

	groups <- unique(readMat.merge$group)
	peakAnnoList <- lapply(groups, function(i)
		readMat.merge[readMat.merge$group == i, ]
		)
	names(peakAnnoList) <- groups

	# Load random file
	rdMat <- read.xlsx(randomFile, sheet = 1)
	rdMat <- annotateExonIntron(rdMat)

	# Forest plot
	tableList <- lapply(peakAnnoList, function(i) table(i$annotation))
	table.rd <- table(rdMat$annotation)

	oddList <- lapply(tableList, getOddMatrix, ref = table.rd)

	write.xlsx(oddList, gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_region_odd_ratio.xlsx", inputFile))

	for(i in 1:length(oddList))
		{
		oddList[[i]]$Group <- names(oddList)[i]
		}

	df <- do.call(rbind, oddList)

	# reverses the factor level ordering for labels after coord_flip()
	df$Region <- factor(df$Region, levels=rev(c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", "5' UTR",
		"first Exon", "first Intron", "other Exon", "other Intron", "3' UTR", "Downstream (<1kb)", "Downstream (1-2kb)",
		"Downstream (2-3kb)", "Distal Intergenic")))

	fp <- ggplot(data=df, aes(x=Region, y=Odd, ymin=Conf.int.low, ymax=Conf.int.high, color = Group)) +
			geom_pointrange(position = position_dodge(width = 0.5)) + 
			geom_hline(yintercept = 1, lty = 2) +  # add a dotted line at x=1 after flip
			coord_flip() +  # flip coordinates (puts labels on y axis)
			xlab("Label") + ylab("Odds ratio (95% CI)") +
			theme_bw()  # use a white background

	pdf(gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_region_odd_ratio.pdf", inputFile))
	#png(gsub("_aln_stat_FLANK_GROUP_GENES.xlsx", "_region_odd_ratio.png", inputFile), units="px", width=1600, height=1600, res=300)
	plot(fp)
	dev.off()
}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE)
{
############
# FANCF TEST

setwd(file.path("/home/gandri/offTargets/Giando/pipeline/FANCF_decoy_rep1/results/overlap_aln"))
readMat <- read.xlsx("FANCF_decoy_ovl_w250.xlsx")
readMat.gr <- makeGRangesFromDataFrame(readMat[, c("chromosome", "start", "end")],
	keep.extra.columns = TRUE)

readMat.anno.gr <- annotatePeak(readMat.gr, tssRegion=c(-3000, 3000), 
							 TxDb=txdb, annoDb="org.Hs.eg.db", overlap = "all", level = "gene")

#####################
# LOAD READMAT

#setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment"))
#sampleName <- "G3_W250"
#setwd(sampleName)
#readMat <- read.xlsx(paste0(sampleName, "_aln_stat_FLANK_GROUP.xlsx"), sheet = 1)


setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final/Giando"))
readMat <- read.xlsx("cas9_union.xlsx")

readMat.gr <- makeGRangesFromDataFrame(readMat[, c("chromosome", "start", "end")],
	keep.extra.columns = TRUE)


#############################
# ANNOTATE PEAKS (chipseeker)

# GET PROMOTER INFO
#promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# Peak Annotation
readMat.anno.gr <- annotatePeak(readMat.gr, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db", overlap = "all", level = "gene")

# Save annotation
readMat.anno <-  as.data.frame(readMat.anno.gr)

##############################
# MERGE READMAT AND ANNOTATION

rownames(readMat) <- readMat$alignment.id
rownames(readMat.anno) <- readMat.anno$alignment.id

readMat.merge <- merge(readMat, readMat.anno[,-c(1:6)], by=0, all=TRUE)
readMat.merge <- readMat.merge[match(readMat$alignment.id, readMat.merge$alignment.id), ]


#########################
# ADD ONCOGENE ANNOTATION
#setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/doc"))
#oncoMat <- read.delim("allOnco_May2018.tsv")
#onco <- as.character(oncoMat$symbol)

#readMat.merge$isOnco <- NA
#readMat.merge$isOnco[readMat.merge$SYMBOL %in% onco] <- "yes"
#readMat.merge <- readMat.merge[,-1]

# SAVE
#setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment"))
#setwd(sampleName)
#write.xlsx(readMat.merge, paste0(sampleName, "_aln_stat_FLANK_GROUP_GENES.xlsx"))


##########################################################################################
# GENE REGIONS PROPORTION

readMat.merge <- read.xlsx(paste0(sampleName, "_aln_stat_FLANK_GROUP_GENES.xlsx"), sheet = 1)
readMat.merge <- annotateExonIntron(readMat.merge)

groups <- unique(readMat.merge$group)
peakAnnoList <- lapply(groups, function(i)
	readMat.merge[readMat.merge$group == i, ]
	)
names(peakAnnoList) <- groups	
        


#########
# BARPLOT

ggmat <- lapply(1:length(peakAnnoList), function(i)
	data.frame(Regions = names(table(peakAnnoList[[i]]$annotation)),
			   Percentage = table(peakAnnoList[[i]]$annotation) / nrow(peakAnnoList[[i]]) * 100,
			   Group = names(peakAnnoList)[i])
	)
ggmat <- do.call(rbind, ggmat)	
ggmat$Regions <- factor(ggmat$Regions, levels = rev(c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)",
	"5' UTR", "first Exon", "other Exon", "first Intron", "other Intron", "3' UTR",
	"Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)",
	"Distal Intergenic"))
	)
ggmat$Group <- factor(ggmat$Group, levels = rev(c("off.target", "hom.recomb", "CBS"))
	)		
			
p <- ggplot(data=ggmat, aes(x=Regions, y=Percentage.Freq, fill=Group))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + theme_bw()	
p <- p + coord_flip()

pdf(paste0(sampleName, "_region_barplot.pdf"))
plot(p)
dev.off()

############################################################
# FOREST PLOT (COMPARE OUR REGIONS AGAINST RANDOM SEQUENCES)

# Real 
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
readMat <- read.xlsx("cas9_UNION.xlsx")

readMat.gr <- makeGRangesFromDataFrame(readMat[, c("chromosome", "start", "end")],
	keep.extra.columns = TRUE)
	
readMat.anno <- annotatePeak(readMat.gr, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db", overlap = "all")
readMat.anno <- annotateExonIntron(as.data.frame(readMat.anno))

# Random bed
setwd(file.path("/Volumes/home/Geoffroy/offTargets/Giando/doc"))
rdMat <- read.delim("random_w250.bed", header = FALSE)
rdMat.gr <- makeGRangesFromDataFrame(rdMat,seqnames.field="V1", start.field = "V2",
	end.field = "V3", strand.field = "V6",
	keep.extra.columns = FALSE)

# Peak Annotation
rdMat.anno.gr <- annotatePeak(rdMat.gr, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db", overlap = "all")
rdMat.anno <- as.data.frame(rdMat.anno.gr)
rdMat.anno <- annotateExonIntron(rdMat.anno)

# Forest plot
peakAnnoList <- list(cas9 = readMat.anno)
tableList <- lapply(peakAnnoList, function(i) table(i$annotation))
table.rd <- table(rdMat.anno$annotation)


oddList <- lapply(tableList, getOddMatrix, ref = table.rd)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/comparison/final_060519"))
write.xlsx(oddList, "cas9_region_odd_ratio.xlsx")

for(i in 1:length(oddList))
	{
	oddList[[i]]$Group <- names(oddList)[i]
	}

df <- do.call(rbind, oddList)

# reverses the factor level ordering for labels after coord_flip()
df$Region <- factor(df$Region, levels=rev(c("Promoter (2-3kb)", "Promoter (1-2kb)", "Promoter (<=1kb)", "5' UTR",
	"first Exon", "first Intron", "other Exon", "other Intron", "3' UTR", "Downstream (<1kb)", "Downstream (1-2kb)",
	"Downstream (2-3kb)", "Distal Intergenic")))
	
#df$Region <- factor(df$Region, levels=rev(c("Promoter", "5' UTR",
#	"Body", "3' UTR", "Downstream", "Distal Intergenic")))	

fp <- ggplot(data=df, aes(x=Region, y=Odd, ymin=Conf.int.low, ymax=Conf.int.high, color = "black")) +
        geom_pointrange(position = position_dodge(width = 0.5), color = "black") + 
        geom_hline(yintercept = 0, lty = 2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Label") + ylab("Odds ratio (95% CI)") +
        theme_bw()  # use a white background



pdf("cas9_region_odd_ratio.pdf")
plot(fp)
dev.off()


}