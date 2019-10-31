

getPAM <- function(x)
{
	length <- nchar(x)
	return(substring(x, length - 1, length))
}

isOnTargetOLD <- function(inputF, otsF, distance)
{
	ots <- read.delim(otsF, header = FALSE)
	chr <- as.character(ots[,1])
	start <- as.numeric(ots[,2])
	end <- as.numeric(ots[,3])

	readMat <- read.xlsx(inputF, sheet = 1)
	idx <- (readMat[,1] == chr) &
		( (readMat[,2] >= (start - distance) & readMat[,2] <= (end + distance)) | (readMat[,3] >= (start - distance) & readMat[,3] <= (end + distance)) )
	return(idx)
}

isOnTarget <- function(inputF, otsF, distance)
{
	readMat <- read.xlsx(inputF, sheet = 1)

	ots <- read.delim(otsF, header = FALSE)
	ots.sub <- data.frame(chr = as.character(ots[,1]),
						  start = as.numeric(ots[,2]),
						  end = as.numeric(ots[,3])
						  )
	
	
	# CONVERT TO GRANGES
	m.gr <- makeGRangesFromDataFrame(readMat, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = FALSE, ignore.strand = TRUE)
		
	ots.gr <- makeGRangesFromDataFrame(ots.sub, seqnames.field = "chr",
		start.field = "start", end.field = "end",
		keep.extra.columns = FALSE, ignore.strand = TRUE)
	
	# INTERSECT
	gr.ovl <- findOverlaps(query = m.gr, subject = ots.gr, type = "any", ignore.strand= TRUE)	
	ovl.idx <- unique(queryHits(gr.ovl))
	
	idx <- rep(FALSE, nrow(readMat))
	idx[ovl.idx] <- TRUE
	
	return(idx)
}

getOTscore <- function(inputF, randomF, cutoff, otsF, distance)
{
	# load files
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat.rd <- read.xlsx(randomF, sheet = 1)
	
	# Get thresholds
	score.thr <- sort(readMat.rd[, "score"], decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# score of the top XX%
	hom.thr <- sort(getBestFlanking(readMat.rd), decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# flanking length of the top XX%
	
	# on targets
	iot <- isOnTarget(inputF, otsF, distance)
	
	scoreList <- mclapply(1:nrow(readMat), function(i){
		score <- 0
		# adjusted pvalue
		if(readMat$adj.pvalue[i] <= 0.05) score <- score + 1 
		else if(readMat$adj.pvalue[i] <= 0.45) score <- score + 0.5
	
		# deduplicated reads
		if(readMat$collapseCluster[i] >= 3) score <- score + 1 
		else if(readMat$collapseCluster[i] == 2) score <- score + 0.5
	
		# cluster width
		if(readMat$width.raw[i] > 10) score <- score + 1
		else if(readMat$width.raw[i] > 0) score <- score + 0.5
	
		# PAM
		if(getPAM(readMat$pattern[i]) == "GG") score <- score + 1
		else if(getPAM(readMat$pattern[i]) == "AG") score <- score + 0.5
	
		# Alignment score
		if(readMat$score[i] > score.thr) score <- score + 1
		else if(readMat$score[i] == score.thr) score <- score + 0.5
	
		# Aln. start rel.
		aln.rel <- readMat$width.raw - abs(readMat$aln.start.rel)
		if(aln.rel[i] >= -50) score <- score + 1
		else if(aln.rel[i] >= -150) score <- score + 0.5
	
		# Homology
		hom <- max(c(readMat$flank.length[i], readMat$flank.rev.length[i]))
		if(hom < hom.thr) score <- score + 1
		else if(hom == hom.thr) score <- score + 0.5
	
		# Dist. On target
		if(!(iot[i])) score <- score + 1
		
		# MAPQ
		if(readMat$avg.MAPQ[i] >= 30) score <- score + 1 
	
		return(score)
		}, mc.cores = NBCPU)
	
	return(unlist(scoreList))
}


getHRscore <- function(inputF, randomF, cutoff, otsF, distance)
{
	# load files
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat.rd <- read.xlsx(randomF, sheet = 1)
	
	# Get thresholds
	score.thr <- sort(readMat.rd[, "score"], decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# score of the top XX%
	hom.thr <- sort(getBestFlanking(readMat.rd), decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# flanking length of the top XX%
	
	# on targets
	iot <- isOnTarget(inputF, otsF, distance)
	
	scoreList <- mclapply(1:nrow(readMat), function(i){
		score <- 0
		# adjusted pvalue
		if(readMat$adj.pvalue[i] <= 0.05) score <- score + 1 
		else if(readMat$adj.pvalue[i] <= 0.45) score <- score + 0.5
	
		# deduplicated reads
		if(readMat$collapseCluster[i] >= 3) score <- score + 1 
		else if(readMat$collapseCluster[i] == 2) score <- score + 0.5
	
		# cluster width
		if(readMat$width.raw[i] > 10) score <- score + 1
		else if(readMat$width.raw[i] > 0) score <- score + 0.5
	
		# Alignment score
		if(readMat$score[i] < score.thr) score <- score + 1
		else if(readMat$score[i] == score.thr) score <- score + 0.5
	
		# Aln. start rel.
		aln.rel <- readMat$width.raw - abs(readMat$aln.start.rel)
		if(aln.rel[i] <= -150) score <- score + 1
		else if(aln.rel[i] <= -50) score <- score + 0.5
	
		# Homology
		hom <- max(c(readMat$flank.length[i], readMat$flank.rev.length[i]))
		if(hom > 2*hom.thr) score <- score + 4
		else if(hom > hom.thr) score <- score + 2
		else if(hom == hom.thr) score <- score + 1
		
		# MAPQ
		if(readMat$avg.MAPQ[i] >= 30) score <- score + 1 
	
		return(score)
		}, mc.cores = NBCPU)
	
	return(unlist(scoreList))
}


getCBSscore <- function(inputF, randomF, cutoff, otsF, distance)
{
	# load files
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat.rd <- read.xlsx(randomF, sheet = 1)
	
	# Get thresholds
	score.thr <- sort(readMat.rd[, "score"], decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# score of the top XX%
	hom.thr <- sort(getBestFlanking(readMat.rd), decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# flanking length of the top XX%
	
	# on targets
	iot <- isOnTarget(inputF, otsF, distance)
	
	scoreList <- mclapply(1:nrow(readMat), function(i){
		score <- 0
		# adjusted pvalue
		if(readMat$adj.pvalue[i] <= 0.05) score <- score + 1 
		else if(readMat$adj.pvalue[i] <= 0.45) score <- score + 0.5
	
		# deduplicated reads
		if(readMat$collapseCluster[i] >= 3) score <- score + 1 
		else if(readMat$collapseCluster[i] == 2) score <- score + 0.5
	
		# cluster width
		if(readMat$width.raw[i] > 10) score <- score + 1
		else if(readMat$width.raw[i] > 0) score <- score + 0.5
	
		# Alignment score
		if(readMat$score[i] < score.thr) score <- score + 1
		else if(readMat$score[i] == score.thr) score <- score + 0.5
	
		# Aln. start rel.
		aln.rel <- readMat$width.raw - abs(readMat$aln.start.rel)
		if(aln.rel[i] <= -150) score <- score + 1
		else if(aln.rel[i] <= -50) score <- score + 0.5
	
		# Homology
		hom <- max(c(readMat$flank.length[i], readMat$flank.rev.length[i]))
		if(hom < hom.thr) score <- score + 1
		else if(hom == hom.thr) score <- score + 0.5

		# MAPQ
		if(readMat$avg.MAPQ[i] >= 30) score <- score + 1 
	
		return(score)
		}, mc.cores = NBCPU)
	
	return(unlist(scoreList))
}


addScore <- function(inputF, randomF, cutoff, otsF, distance)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat$OT.score <- getOTscore(inputF, randomF, cutoff, otsF, distance)
	readMat$HR.score <- getHRscore(inputF, randomF, cutoff, otsF, distance)
	readMat$CBS.score <- getCBSscore(inputF, randomF, cutoff, otsF, distance)

	write.xlsx(readMat, gsub(".xlsx", "_SCORE.xlsx", inputF), row.names = FALSE)
}

scoreDensity <- function(inputF)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	ggmat <- readMat[, c("OT.score", "HR.score", "CBS.score")]
	ggmat <- melt(ggmat)
	
	p <- ggplot(ggmat, aes(value, group = variable, fill = variable))
	p <- p + geom_density(alpha = 0.25)
	p <- p + scale_fill_manual(values = c("OT.score" = "red", "HR.score" = "blue", "CBS.score" = "grey"))
	p <- p + theme_bw()
	
	pdf(gsub(".xlsx", "_density.pdf", inputF))
	plot(p)
	dev.off()
}

###########################
# TALEN DEDICATED FUNCTIONS

getOTscoreTALEN <- function(inputF, randomF, cutoff, otsF, distance)
{
	# load files
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat.rd <- read.xlsx(randomF, sheet = 1)
	
	# Get thresholds
	score.thr <- sort(readMat.rd[, "score"], decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# score of the top XX%
	hom.thr <- sort(getBestFlanking(readMat.rd), decreasing = TRUE)[floor(nrow(readMat.rd) * cutoff)]# flanking length of the top XX%
	
	# on targets
	iot <- isOnTarget(inputF, otsF, distance)
	
	scoreList <- mclapply(1:nrow(readMat), function(i){
		score <- 0
		# adjusted pvalue
		if(readMat$adj.pvalue[i] <= 0.05) score <- score + 1 
		else if(readMat$adj.pvalue[i] <= 0.45) score <- score + 0.5
	
		# deduplicated reads
		if(readMat$collapseCluster[i] >= 3) score <- score + 1 
		else if(readMat$collapseCluster[i] == 2) score <- score + 0.5
	
		# cluster width
		if(readMat$width.raw[i] > 10) score <- score + 1
		else if(readMat$width.raw[i] > 0) score <- score + 0.5
	
		# PAM
		if(getPAM(readMat$pattern[i]) == "GG") score <- score + 1
		else if(getPAM(readMat$pattern[i]) == "AG") score <- score + 0.5
	
		# Alignment score
		if(readMat$score[i] > score.thr) score <- score + 1
		else if(readMat$score[i] == score.thr) score <- score + 0.5
	
		# Aln. start rel.
		aln.rel <- readMat$width.raw - abs(readMat$aln.start.rel)
		if(aln.rel[i] >= -50) score <- score + 1
		else if(aln.rel[i] >= -150) score <- score + 0.5
	
		# Homology
		hom <- max(c(readMat$flank.length[i], readMat$flank.rev.length[i]))
		if(hom < hom.thr) score <- score + 1
		else if(hom == hom.thr) score <- score + 0.5
	
		# Dist. On target
		if(!(iot[i])) score <- score + 1
		
		# MAPQ
		if(readMat$avg.MAPQ[i] >= 30) score <- score + 1 
	
		return(score)
		}, mc.cores = NBCPU)
	
	return(unlist(scoreList))
}


addScoreTALEN <- function(inputF, randomF, cutoff, otsF, distance)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat$OT.score <- getOTscore(inputF, randomF, cutoff, otsF, distance)
	readMat$HR.score <- getHRscore(inputF, randomF, cutoff, otsF, distance)
	readMat$CBS.score <- getCBSscore(inputF, randomF, cutoff, otsF, distance)

	write.xlsx(readMat, gsub(".xlsx", "_SCORE.xlsx", inputF), row.names = FALSE)
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE)
{
# DEFINE INPUTS
inputF <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/guide_aln", "FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")

randomF <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/random", "random_w250_aln_stat_FLANK.xlsx")

otsF <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/data/", "ots.bed")


# GET SCORE
getOTscore(inputF, randomF, 0.05, otsF, 250)
getHRscore(inputF, randomF, 0.05, otsF, 250)
getCBSscore(inputF, randomF, 0.05, otsF, 250)


addScore(inputF, randomF, 0.05, otsF, 250)


# PLOT DENSITY
scoreDensity(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/guide_aln", "FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx"))

}