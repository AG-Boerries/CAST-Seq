#library(Biostrings)
#library(openxlsx)
#library(data.table)
#library(ggplot2)
#library(ggseqlogo)
#library(textreadr)
#library(parallel)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getEmpiricalPV <- function(x, y, type)
{
    x.ecdf <- stats::ecdf(x)
    if(type == "greater") return(1-x.ecdf(y))
    else if(type == "less") return(x.ecdf(y))
    return(NA)    
}

assignGroups <- function(realF, rdF, otsF, cutoff)
{
	ots <- read.delim(otsF, header = FALSE)

	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	realM.bf <- getBestFlanking(realM)
	rdM.bf <- getBestFlanking(rdM)
	
	score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
	#flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
	flanking.cutoff <- 25
	
	score.pv <- unlist(lapply(realM[, "score"], getEmpiricalPV, x = rdM[, "score"], type = "greater"))
	score.adj.pv <- p.adjust(score.pv, method = "BH")
	flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
	flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")
	#print(score.cutoff)
	#print(flanking.cutoff)	
	
	# plot score.cutoff (REAL)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (REAL)
	#xmax <- min(100, max(realM.bf))
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM.bf[realM.bf <= 100]), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot score.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM.bf), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	mygroups <- rep("CBS", nrow(realM))
	mygroups[realM.bf >= flanking.cutoff] <- "hom.recomb"
	mygroups[score.adj.pv < cutoff] <- "off.target"

	mygroups2 <- rep(NA, nrow(realM))
	#mygroups2[mygroups == "off.target" & (realM.bf > flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "off.target" & (realM.bf >= flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "hom.recomb"] <- "yes"
	
	mygroups3 <- rep(NA, nrow(realM))
	ld.idx <- isLargeDel(realM, ots)
	#print(length(ld.idx))
	if(length(ld.idx)!=0){
		if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
		}
	
	
	realM$group <- mygroups
	realM$is.hom.recomb. <- mygroups2
	realM$is.large.del. <- mygroups3
	
	realM$off.target.pvalue <- score.pv 
	realM$off.target.adj.pvalue <- score.adj.pv 
	realM$hom.recomb.pvalue <- flanking.pv
	realM$hom.recomb.adj.pvalue <- flanking.adj.pv
	
	
	write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE)	
}

assignGroups_2OT <- function(realF, rdF, otsF, cutoff)
{
	ots <- read.delim(otsF, header = FALSE)

	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	realM.bf <- getBestFlanking_2OT(realM)
	rdM.bf <- getBestFlanking_2OT(rdM)
	
	score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
	#flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
	flanking.cutoff <- 25
	
	score.pv <- unlist(lapply(realM[, "score"], getEmpiricalPV, x = rdM[, "score"], type = "greater"))
	score.adj.pv <- p.adjust(score.pv, method = "BH")
	flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
	flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")
	#print(score.cutoff)
	#print(flanking.cutoff)	
	
	# plot score.cutoff (REAL)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (REAL)
	#xmax <- min(100, max(realM.bf))
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM.bf[realM.bf <= 100]), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot score.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM.bf), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	mygroups <- rep("CBS", nrow(realM))
	mygroups[realM.bf >= flanking.cutoff] <- "hom.recomb"
	mygroups[score.adj.pv < cutoff] <- "off.target"

	mygroups2 <- rep(NA, nrow(realM))
	#mygroups2[mygroups == "off.target" & (realM.bf > flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "off.target" & (realM.bf >= flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "hom.recomb"] <- "yes"
	
	mygroups3 <- rep(NA, nrow(realM))
	ld.idx <- isLargeDel(realM, ots)
	#print(length(ld.idx))
	if(length(ld.idx)!=0){
		if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
		}
	
	
	realM$group <- mygroups
	realM$is.hom.recomb. <- mygroups2
	realM$is.large.del. <- mygroups3
	
	realM$off.target.pvalue <- score.pv 
	realM$off.target.adj.pvalue <- score.adj.pv 
	realM$hom.recomb.pvalue <- flanking.pv
	realM$hom.recomb.adj.pvalue <- flanking.adj.pv
	
	
	write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE)	
}

assignGroupsTALEN <- function(realF, rdF, otsF, cutoff)
{
	ots <- read.delim(otsF, header = FALSE)

	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	realM.bf <- getBestFlanking(realM)
	rdM.bf <- getBestFlanking(rdM)
	
	#score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
	#flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
	flanking.cutoff <- 25
	
	#print(score.cutoff)
	#print(flanking.cutoff)	
	
	bestCb <- realM$BestCB
	qv <- lapply(1:length(bestCb), function(i){
		if(is.na(bestCb[i])) return(1)
		qv.current <- realM[i, grep("_adj.pv$", colnames(realM))]
		qv.current <- as.numeric(qv.current[grep(bestCb[i], names(qv.current))])
		return(qv.current)
		})
	qv <- unlist(qv)
	
	
	flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
	flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")

	#nbSignif <- realM[, grep("_adj.pv", colnames(realM))]
	#nbSignif <- (colSums(nbSignif) / ncol(nbSignif)) * 100# percentage of significant pvalue
	
	mygroups <- rep("CBS", nrow(realM))
	mygroups[realM.bf >= flanking.cutoff] <- "hom.recomb"
	mygroups[qv < cutoff] <- "off.target"

	mygroups2 <- rep(NA, nrow(realM))
	mygroups2[mygroups == "off.target" & (realM.bf >= flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "hom.recomb"] <- "yes"
	
	mygroups3 <- rep(NA, nrow(realM))
	ld.idx <- isLargeDel(realM, ots)
	if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
	
	realM$group <- mygroups
	realM$is.hom.recomb. <- mygroups2
	realM$is.large.del. <- mygroups3
	
	realM$hom.recomb.pvalue <- flanking.pv
	realM$hom.recomb.adj.pvalue <- flanking.adj.pv
	
	write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE)	
	
	# flanking plot (REAL)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM.bf[realM.bf <= 100]), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM.bf), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()

}

getBestFlanking <- function(m)
{
	m.max <- apply(m[, c("flank.length", "flank.rev.length")],1,max)
	return(m.max)
}

getBestFlanking_2OT <- function(m)
{
	m.max <- apply(m[, c("flank.length", "flank.rev.length", "flank2.length", "flank2.rev.length")],1,max)
	return(m.max)
}

isLargeDel <- function(m, ots)
{
	m.sub <- m[, c("chromosome", "start", "end")]
	m.sub$strand <- "*"
	
	ots.sub <- data.frame(chr = as.character(ots[,1]),
						  start = as.numeric(ots[,2]),
						  end = as.numeric(ots[,3]),
						  strand = "*")
	
	
	# CONVERT TO GRANGES
	m.gr <- makeGRangesFromDataFrame(m.sub, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = FALSE)
		
	ots.gr <- makeGRangesFromDataFrame(ots.sub, seqnames.field = "chr",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = FALSE)
	
	# INTERSECT
	gr.ovl <- findOverlaps(query = m.gr, subject = ots.gr, type = "any", ignore.strand= TRUE)	
	ovl.idx <- unique(queryHits(gr.ovl))
	
	return(ovl.idx)
}


groupSummary <- function(inputF, outputF, clusters = NULL, score = NULL, pv = NULL)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	if(!is.null(clusters)) readMat <- readMat[readMat$collapseCluster > clusters, ]
	if(!is.null(score)) readMat <- readMat[readMat$score > score, ]
	if(!is.null(pv)) readMat <- readMat[readMat$adj.pvalue < pv, ]
	
	if(nrow(readMat)>0){
		write.xlsx(data.frame(table(readMat$group)), outputF, row.names = FALSE)
		}
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE)
{
################
# LOAD REAL DATA

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment/G3_W250"))
readMat.real <- read.xlsx("G3_W250_aln_stat_FLANK.xlsx", sheet = 1)

##################
# LOAD RANDOM DATA

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment/random_hg38_w250"))
readMat.rd <- read.xlsx("random_hg38_w250_aln_stat_FLANK.xlsx", sheet = 1)


#############################################################################
# ASSIGN GROUPS: off targets, homologous recombination, common breaking sites

cutoff <- 0.05
readMat.real <- assignGroups(readMat.real, readMat.rd, cutoff)


# SAVE
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment/G3_W250"))
write.xlsx(readMat.real, "G3_W250_aln_stat_FLANK_GROUP.xlsx")
}




