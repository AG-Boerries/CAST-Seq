#library(tools)
#library(Biostrings)
#library(openxlsx)
#library(data.table)
#library(ggplot2)
#library(ggseqlogo)
#library(parallel)

#source("/home/gandri/offTargets/Giando/pipeline/script/bed2sequence.R")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

toNum <- function(x) as.numeric(levels(x))[x]

getRevComp <- function(x) as.character(reverseComplement(DNAString(x)))

strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

getBestAln <- function(alnList)
{
	scoreList <- unlist(lapply(alnList, score))
	return(alnList[[which.max(scoreList)]])
}

getBestAlnIdx <- function(alnList)
{
	scoreList <- unlist(lapply(alnList, score))
	return(which.max(scoreList))
}

getStrAlignment <- function(aln.file)
{
	#aln.str <- read_document(aln.file)
	aln.str <- read.delim(aln.file, sep = "\t", header = FALSE)
	aln.str <- as.character(t(aln.str))
	
	# get S1 index
	s1.idx <- grep("^S1", aln.str)
	
	s1.str <- lapply(s1.idx, function(i){
		ref.str <- aln.str[i]
		ref.str <- strsplit(ref.str, split = " ")
		ref.str <- unlist(ref.str)
		ref.str <- ref.str[ref.str != ""]
		return(ref.str[3])	
		})
	s1.str <- paste(s1.str, collapse = "")

	# get sequence idx
	sq.idx <- s1.idx - 2
	sq.str <- lapply(sq.idx, function(i){
		ref.str <- aln.str[i]
		ref.str <- strsplit(ref.str, split = " ")
		ref.str <- unlist(ref.str)
		ref.str <- ref.str[ref.str != ""]
		return(ref.str[3])	
		})
	sq.str <- paste(sq.str, collapse = "")
	
	return(c(sq.str, s1.str))
}


getStartPosition <- function(aln.file)
{
	#aln.str <- read_document(aln.file)
	aln.str <- read.delim(aln.file, sep = "\t", header = FALSE)
	aln.str <- as.character(t(aln.str))
	
	test.str <- aln.str[22]
	test.str <- strsplit(test.str, split = " ")
	test.str <- unlist(test.str)
	test.str <- test.str[test.str != ""]
	return(test.str[2])
}



getAlnStat <- function(aln)
{
	#pattern <- as.character(pattern(aln))
	#subject <- as.character(subject(aln))
	score <- score(aln)
	sq.id <- pid(aln)
	nb.mism <- nmismatch(aln)
	nb.indel <- nindel(aln)
	nb.ins <- insertion(nb.indel)[1]
	length.ins <- insertion(nb.indel)[2]
	nb.del <- deletion(nb.indel)[1]
	length.del <- deletion(nb.indel)[2]
	aln.chr <- compareStrings(aln)
	return(c(score, sq.id, nb.mism, nb.ins, length.ins, nb.del, length.del, aln.chr))
}

getAlnChar <- function(aln.test, aln.ref, ref, idn = NULL)
{
	alnChar <- rep(NA, nchar(ref))
	insChar <- rep(0, nchar(ref))
	index <- 1
	
	for(i in 1:nchar(aln.test))
		{

		currentRef <- substring(ref, index, index)
		
		#print(i)
		#print(index)
		#print(currentRef)
		#print(substring(aln.test, i, i))
		#print(substring(aln.ref, i, i))
		#print("")
	
		if(substring(aln.test, i, i) == currentRef & substring(aln.ref, i, i) == currentRef){
			if(is.null(idn)) alnChar[index] <- substring(aln.ref, i, i)
			else(alnChar[index] <- ".")
			index <- index + 1
			#print("match")
			} else if(substring(aln.test, i, i) != "-" & substring(aln.test, i, i) != currentRef & substring(aln.ref, i, i) == currentRef){
				alnChar[index] <- substring(aln.test, i, i)
				index <- index + 1
				#print("mismatch")
				} else if(substring(aln.test, i, i) == "-"){
					alnChar[index] <- "-1"
					index <- index + 1
					#print("deletion")
					}else if(substring(aln.ref, i, i) == "-"){
						insChar[index - 1] <- insChar[index - 1] + 1
						#print("insertion")		
						}else(print("else"))	
		}	
	#print(alnChar)
	#print(insChar)
	finalChar <- unlist(lapply(1:length(alnChar), function(i)
		c(alnChar[i], insChar[i])
		))
	return(finalChar)
}

getAlnSumStr <- function(aln.test, aln.ref)
{
	alnChar <- rep(NA, nchar(aln.test))

	for(i in 1:nchar(aln.test))
		{
		if(substring(aln.test, i, i) == substring(aln.ref, i, i)) alnChar[i] <- tolower(substring(aln.test, i, i))
		else if(substring(aln.test, i, i) == "-") alnChar[i] <- "-"
		else if(substring(aln.ref, i, i) == "-") alnChar[i] <- "+"
		else alnChar[i] <- toupper(substring(aln.test, i, i))
		}
	return(paste(alnChar, collapse=""))
}

getRangeAlignment <- function(sq, p1, p2, nR, sm)
{
	rAlnList <- lapply(nR, function(n){
		refSeq.current <- paste0(p1, strrep("N", n), p2)
		return(pairwiseAlignment(sq, subject = refSeq.current,
			type = "local-global", substitutionMatrix = sm,
			gapOpening = 1, gapExtension = 1))	
		})
	bestIdx <- getBestAlnIdx(rAlnList)
	return(list(rAlnList[[bestIdx]], nR[bestIdx]))
}


getGuideAlignmentDual <- function(inputF, guideLeft, guideRight, alnFolder)
{
	# DEFINE N RANGE
	nRange <- seq(5, 50)

	# DEFINE REFERENCE SEQUENCE
	refSeq.left <- guideLeft
	refSeq.right <- guideRight
	
	refSeq.left.rev <- getRevComp(refSeq.left)
	refSeq.right.rev <- getRevComp(refSeq.right)
	
	refSeqMat <- matrix(c(refSeq.left, refSeq.left.rev,
						  refSeq.left, refSeq.right,
						  refSeq.left, refSeq.right.rev,
						  refSeq.left.rev, refSeq.right,
						  refSeq.left.rev, refSeq.right.rev,
						  refSeq.right, refSeq.right.rev,
						  
						  refSeq.left.rev, refSeq.left,
						  refSeq.right, refSeq.left,
						  refSeq.right.rev, refSeq.left,
						  refSeq.right, refSeq.left.rev,
						  refSeq.right.rev, refSeq.left.rev,
						  refSeq.right.rev, refSeq.right,
						  
						  refSeq.left, refSeq.left,
						  refSeq.left.rev, refSeq.left.rev,
						  refSeq.right, refSeq.right,
						  refSeq.right.rev, refSeq.right.rev
						  ),
		ncol = 2, byrow = TRUE)
	
	rownames(refSeqMat) <- c("LF.LR", "LF.RF", "LF.RR", "LR.RF", "LR.RR", "RF.RR",
		"LR.LF", "RF.LF", "RR.LF", "RF.LR", "RR.LR", "RR.RF",
		"LF.LF", "LR.LR", "RF.RF", "RR.RR")

	# LOAD BED FILE
	if(file_ext(inputF) == "bed"){
		readMat <- read.delim(inputF, header = FALSE)
		colnames(readMat) <- c("chromosome", "start", "end", "ID_read", "size", "strand")
		}else{
			readMat <- read.xlsx(inputF)
			}	
	inputName <- gsub(paste0(".", file_ext(inputF)), "", basename(inputF))	
	print(inputName)
	
	# REGIONS 2 SEQUENCE
	sequences <- bed2sequence(readMat)
	#sequences <- c(refSeq, sequences)# add refSeq

	###########
	# ALIGNMENT
	sm <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
	sm[, "N"] <- 0
	sm["N", ] <- 0
	sm["N", "N"] <- 0

	alnList.raw <- lapply(1:nrow(refSeqMat), function(j){
		x <- mclapply(sequences, function(i){
			getRangeAlignment(i, refSeqMat[j, 1], refSeqMat[j, 2], nRange, sm)
			}, mc.cores = NBCPU)
		return(x)
		})
	names(alnList.raw) <- rownames(refSeqMat)	
	
	for(i in 1:length(alnList.raw))
		{
		names(alnList.raw[[i]]) <- paste0(inputName, "_aln", 1:length(alnList.raw[[i]]))	
		}
		
	alnList <- lapply(alnList.raw, function(i){
		lapply(i, function(j) j[[1]])
		})
	
	alnList.N <- lapply(alnList.raw, function(i){
		lapply(i, function(j) j[[2]])
		})
	
	scoreList.N <- lapply(alnList.raw, function(i){
		lapply(i, function(j) score(j[[1]]))
		})	
				
	alnN.mat <- matrix(unlist(alnList.N), nrow = length(sequences), ncol = nrow(refSeqMat), byrow = FALSE)		
	colnames(alnN.mat) <- paste0(rownames(refSeqMat), "_N")
	rownames(alnN.mat) <- paste0(inputName, "_aln", 1:nrow(readMat))
	write.table(alnN.mat, file.path(alnFolder, paste0(inputName, "_N_matrix.txt")), sep = "\t", quote = FALSE, row.names = TRUE)
	
	scoreN.mat <- matrix(unlist(scoreList.N), nrow = length(sequences), ncol = nrow(refSeqMat), byrow = FALSE)		
	colnames(scoreN.mat) <- paste0(rownames(refSeqMat), "_Score")
	rownames(scoreN.mat) <- paste0(inputName, "_aln", 1:nrow(readMat))
	write.table(scoreN.mat, file.path(alnFolder, paste0(inputName, "_score_matrix.txt")), sep = "\t", quote = FALSE, row.names = TRUE)
	
	# Plot score boxplot
	ggmat <- melt(scoreN.mat)
	ggmat$Var2 <- factor(ggmat$Var2, levels = rev(names(sort(colMeans(scoreN.mat)))))
	
	p <- ggplot(ggmat, aes(x = Var2, y = value))
	p <- p + geom_boxplot()
	p <- p + theme_bw()
	p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

	pdf(file.path(alnFolder, paste0(inputName, "_score_boxplot.pdf")))
	plot(p)
	dev.off()
	
	scoreN.mat <- scoreN.mat[order(-rowMeans(scoreN.mat)), ]
	ggmat <- melt(scoreN.mat[1:10, ])
	ggmat$Var1 <- factor(ggmat$Var1, levels = rev(names(sort(rowMeans(scoreN.mat)))))
		
	p <- ggplot(ggmat, aes(x = Var1, y = value))
	p <- p + geom_boxplot()
	p <- p + theme_bw()
	p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

	pdf(file.path(alnFolder, paste0(inputName, "_score_per_aln_boxplot.pdf")), width = 15)
	plot(p)
	dev.off()
		
	# Plot N boxplot
	ggmat <- melt(alnN.mat)
	ggmat$Var2 <- factor(ggmat$Var2, levels = rev(names(sort(colMeans(alnN.mat)))))
	
	p <- ggplot(ggmat, aes(x = Var2, y = value))
	p <- p + geom_boxplot()
	p <- p + theme_bw()
	p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

	pdf(file.path(alnFolder, paste0(inputName, "_N_boxplot.pdf")))
	plot(p)
	dev.off()
	

	# SAVE
	dir.create(alnFolder, showWarnings = FALSE)
	#alnList <- alnList[-1]# remove refSeq alignment

	lapply(1:length(alnList), function(j){
		dir.create(file.path(alnFolder, names(alnList)[j]), showWarnings = FALSE)
		lapply(1:length(alnList[[j]]), function(i){
			writePairwiseAlignments(alnList[[j]][[i]],
				file = file.path(alnFolder, names(alnList)[j], paste0(names(alnList[[j]])[i], "_TMP.txt")))
			})		
		})

	statMatList <- lapply(names(alnList), function(i){
		# GET PROPER ALIGNMENT (FIX ISSUE WITH INDEL IN FIRST / LAST POSITION)
		refSeqFolder <- i
		alnFiles <- list.files(file.path(alnFolder, refSeqFolder), pattern = "_TMP.txt")
		names(alnFiles) <- gsub("_TMP.txt", "", alnFiles)

		alnFiles <- alnFiles[names(alnList[[refSeqFolder]])]# Re-order
	
		alnList.str <- lapply(alnFiles, function(i)
			getStrAlignment(file.path(alnFolder, refSeqFolder, i))
			)

		alnMat <- do.call(rbind, alnList.str)# 1st column: test, 2nd column: ref
	
		# GET ALIGNMENT N
		alnN <- unlist(alnList.N[[i]])
			
		# GET ALIGNMENT SUMMARY STR
		alnSum <- unlist(lapply(1:nrow(alnMat), function(i)
			getAlnSumStr(alnMat[i, 1], alnMat[i, 2])
			))

		# START POSITION IN TEST SEQUENCE
		alnStart <- as.numeric(unlist(lapply(alnFiles, function(i)
			getStartPosition(file.path(alnFolder, refSeqFolder, i)))))

		# RELATIVE COORDINATES
		middleCoord <- floor((readMat$end - readMat$start) / 2)
		middleCoord.abs <- readMat$start + middleCoord
		alnStart.rel <- alnStart - middleCoord

		# STATISTIC
		statList <- lapply(alnList[[refSeqFolder]], getAlnStat)
		statMat <- do.call(rbind, statList)
		statMat <- cbind(alnMat, alnSum, alnN, statMat, middleCoord.abs, alnStart.rel)
		colnames(statMat) <- c("pattern", "subject", "aln.sum", "N", "score", "sq.hom", "nb.mism",
			"nb.ins", "length.ins", "nb.del", "length.del", "aln.chr", "middleCoord", "aln.start.rel")
		statMat <- as.data.frame(statMat)
		statMat[, c("N", "score", "sq.hom", "nb.mism", "nb.ins", "length.ins", "nb.del", "length.del", "middleCoord", "aln.start.rel")] <- apply(
					statMat[, c("N", "score", "sq.hom", "nb.mism", "nb.ins",
					"length.ins", "nb.del", "length.del", "middleCoord", "aln.start.rel")], 2, as.numeric)
	
		colnames(statMat) <- paste(refSeqFolder, colnames(statMat), sep = "_")
		return(statMat)
		})
	
	statMat <- do.call(cbind, statMatList)

	
	# MERGE READMAT AND STATMAT
	readMat.stat <- cbind(readMat,
						  alignment.id = paste0(inputName, "_aln", 1:nrow(readMat)),
						  statMat)	
	write.xlsx(readMat.stat, file.path(alnFolder, paste0(inputName, "_aln_stat.xlsx")))	

	# Density plot
	#pdf(file.path(alnFolder, paste0(inputName, "_aln_stat_density.pdf")))
	#par(mfrow=c(3,2))
	#plot(density(readMat.stat$score), main = "Score")
	#plot(density(readMat.stat$sq.hom), main = "Sequence Homology")
	#plot(density(readMat.stat$nb.mism), main = "Nb. mismatches")
	#plot(density(readMat.stat$nb.ins), main = "Nb. insertion")
	#plot(density(readMat.stat$nb.del), main = "Nb. deletion")
	#plot(density(readMat.stat$aln.start.rel), main = "Aln. relative start")
	#dev.off()

	#pdf(file.path(alnFolder, paste0(inputName, "_aln_stat_histogram.pdf")))
	#par(mfrow=c(3,2))
	#hist(readMat.stat$score, main = "Score")
	#hist(readMat.stat$sq.hom, main = "Sequence Homology")
	#hist(readMat.stat$nb.mism, main = "Nb. mismatches")
	#hist(readMat.stat$nb.ins, main = "Nb. insertion")
	#hist(readMat.stat$nb.del, main = "Nb. deletion")
	#hist(readMat.stat$aln.start.rel, main = "Aln. relative start")
	#dev.off()
}


assignPV <- function(realF, rdF)
{
	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	# remove random NNNNNN
	rdM <- rdM[-grep("N", rdM$LF.LR_pattern), ]
	
	# get score
	realM.score <- realM[, grep("score", colnames(realM))]
	rdM.score <- rdM[, grep("score", colnames(rdM))]
	
	# ecdf
	pvList <- lapply(1:ncol(realM.score), function(i){
		rd.ecdf <- stats::ecdf(rdM.score[, i])
		pv <- lapply(as.numeric(realM.score[, i]), function(j)
			1-rd.ecdf(j)
			)
		return(unlist(pv))
		})

	pvMat <- do.call(cbind, pvList)
	rownames(pvMat) <- realM$alignment.id
	colnames(pvMat) <- gsub("score", "pv", colnames(realM.score))
	
	qvMat <- t(apply(pvMat, 1, p.adjust))
	colnames(qvMat) <- gsub("pv", "adj.pv", colnames(qvMat))
	
	realM.pv <- cbind(realM, pvMat, qvMat)
	write.xlsx(realM.pv, realF, row.names = FALSE)
}


assignBestCombination <- function(inputFile)
{
	cmb <- c("LF.LR", "LF.RF", "LF.RR", "LR.RF", "LR.RR", "RF.RR",
		"LR.LF", "RF.LF", "RR.LF", "RF.LR", "RR.LR", "RR.RF",
		"LF.LF", "LR.LR", "RF.RF", "RR.RR")

	readMat <- read.xlsx(inputFile, sheet = 1)
	#lapply(cmb, function(i) grep(i, colnames(readMat)))

	dMat <- readMat[, grep("_aln.start.rel$", colnames(readMat))]
	nMat <- readMat[, grep("_N$", colnames(readMat))]
	cbName <- unlist(lapply(strsplit(colnames(nMat), split = "_"), function(i) i[1]))
	
	isMono <- paste0(substring(cbName, 2, 2), substring(cbName, 5, 5)) %in% c("FF", "RR")

	bestCb <- lapply(1:nrow(readMat), function(i){
		qv <- as.numeric(readMat[i, grep("_adj.pv$", colnames(readMat))])
		names(qv) <- cbName
		if(min(qv) > 0.1) return(NA)
		
		signifIdx <- which(qv <= 0.1)
		
		if(length(signifIdx) == 1) return(cbName[signifIdx])
		
		qv.sub <- qv[signifIdx]
		dMat.sub <- dMat[i,signifIdx]
		nMat.sub <- nMat[i,signifIdx]
		cbName.sub <- cbName[signifIdx]
		isMono.sub <- isMono[signifIdx]
		
		score <- unlist(lapply(1:length(signifIdx), function(j){
			score.current <- 0
			if(nMat.sub[, j] < 10){score.current <- score.current + 1
				}else if(nMat.sub[, j] < 20) score.current <- score.current + 2
		
			if(abs(dMat.sub[, j]) < 50) score.current <- score.current + 1
			
			if(!isMono.sub[j]) score.current <- score.current + 1
			#if(isMono.sub[j]) score.current <- score.current - 1
			
			return(score.current)		
			}))
		names(score) <- names(qv.sub)
		
		qv.sub.sub <- qv.sub[names(score)[score == max(score)]]
			
		return(names(qv.sub.sub[order(qv.sub.sub)])[1])
		})
	bestCb <- unlist(bestCb)
	
	readMat$BestCB <- bestCb
	write.xlsx(readMat, inputFile)
}


plotSites <- function(inputF, siteIdx, outputF)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	readMat <- readMat[siteIdx, ]# select sites
	rownames(readMat) <- apply(readMat[1:2], 1, paste, collapse = ";")

	scoreM <- readMat[,grep("_score", colnames(readMat))]
	scoreM$loc <- rownames(scoreM)

	distM <- readMat[,grep("_aln.start.rel", colnames(readMat))]
	distM$loc <- rownames(distM)

	nM <- readMat[,grep("_N", colnames(readMat))]
	nM$loc <- rownames(distM)

	scoreM.gg <- melt(scoreM)
	distM.gg <- melt(distM)
	nM.gg <- melt(nM)

	ggmat <- cbind(scoreM.gg, distance = distM.gg$value, N = nM.gg$value)
	ggmat$variable <- gsub("_score", "", ggmat$variable)

	p <- ggplot(ggmat, aes(x=distance, y=value, color = N))
	p <- p + geom_point(size = 4, alpha = 0.5) + scale_color_viridis_c(option = "B")
	p <- p + facet_wrap(~ loc)
	p <- p + theme_bw()#+ theme(legend.position="none")
	p <- p + geom_vline(xintercept=0, linetype = "dashed")
	p <- p + geom_text_repel(aes(label = variable), size = 3)
	
	pdf(outputF, width = 10, height = 8)
	plot(p)
	dev.off()	
}





############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################



if(FALSE)
{
library(BSgenome.Hsapiens.UCSC.hg38)# hg38!!!
NBCPU <- 12

inputF <- file.path("/home/gandri/offTargets/Giando/pipeline/HBBTal_D1_LEFT/results/guide_aln/HBBTalD1_S7_L001_w250.xlsx")
inputF <- file.path("/home/gandri/offTargets/Giando/pipeline/HBBTal_D1/results/overlap_aln/HBBTal_ovl_w250.xlsx")

guideLeft <- file.path("/home/gandri/offTargets/Giando/pipeline/HBBTal_D1_LEFT/data/gRNA_Left.fa")
guideLeft <- toupper(as.character(readDNAStringSet(guideLeft)))
guideRight <- file.path("/home/gandri/offTargets/Giando/pipeline/HBBTal_D1_LEFT/data/gRNA_Right.fa")
guideRight <- toupper(as.character(readDNAStringSet(guideRight)))
alnFolder <- file.path("/home/gandri/offTargets/Giando/test/")

#getGuideAlignmentDual(inputF, guideLeft, guideRight, alnFolder)


inputF <- file.path("/home/gandri/offTargets/Giando/pipeline/HBBTal_D1_LEFT/results/random/random_w250.bed")
alnFolder <- file.path("/home/gandri/offTargets/Giando/test_random/")

getGuideAlignmentDual(inputF, guideLeft, guideRight, alnFolder)


assignPV(file.path("~/cluster/master/offTargets/Giando/test/HBBTalD1_S7_L001_w250_aln_stat.xlsx"),
	file.path("~/cluster/master/offTargets/Giando/test_random/random_w250_aln_stat.xlsx"))# generate HBBTalD1_S7_L001_w250_aln_stat_PV.xlsx
assignBestCombination(file.path("~/cluster/master/offTargets/Giando/test/HBBTalD1_S7_L001_w250_aln_stat_PV.xlsx"))# generate HBBTalD1_S7_L001_w250_aln_stat_BC.xlsx


#########################
# PLOT SCORE FCT DISTANCE

readMat <- read.xlsx(file.path("~/cluster/master/offTargets/Giando/test/", "HBBTalD1_S7_L001_w250_aln_stat_PV.xlsx"))
readMat <- readMat[c(2, 10, 48, 51, 171, 190), ]# select sites
rownames(readMat) <- apply(readMat[1:2], 1, paste, collapse = ";")

scoreM <- readMat[,grep("_score", colnames(readMat))]
scoreM$loc <- rownames(scoreM)

distM <- readMat[,grep("_aln.start.rel", colnames(readMat))]
distM$loc <- rownames(distM)

nM <- readMat[,grep("_N", colnames(readMat))]
nM$loc <- rownames(distM)

scoreM.gg <- melt(scoreM)
distM.gg <- melt(distM)
nM.gg <- melt(nM)


ggmat <- cbind(scoreM.gg, distance = distM.gg$value, N = nM.gg$value)
ggmat$variable <- gsub("_score", "", ggmat$variable)

p <- ggplot(ggmat, aes(x=distance, y=value, color = N))
p <- p + geom_point(size = 4, alpha = 0.5) + scale_color_viridis_c(option = "B")
p <- p + facet_wrap(~ loc)
p <- p + theme_bw()#+ theme(legend.position="none")
p <- p + geom_vline(xintercept=0, linetype = "dashed")
p <- p + geom_text_repel(aes(label = variable), size = 3)
	
pdf("score_vs_distance.pdf", width = 10, height = 8)
plot(p)
dev.off()	
	
}







