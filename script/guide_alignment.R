#library(tools)
#library(Biostrings)
#library(openxlsx)
#library(data.table)
#library(ggplot2)
#library(ggseqlogo)
#library(textreadr)
#library(parallel)

#source("~/bitbucket/work/offTargets/Giando/bed2sequence.R")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

toNum <- function(x) as.numeric(levels(x))[x]

getRevComp <- function(x) as.character(reverseComplement(DNAString(x)))

strReverse <- function(x)
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

getBestAln <- function(aln1, aln2)
{
	if(score(aln2) > score(aln1)) return(aln2)
	return(aln1)
}

getStrAlignment <- function(aln.file)
{
	#aln.str <- read_document(aln.file)
	aln.str <- read.delim(aln.file, sep = "\t", header = FALSE)
	aln.str <- as.character(t(aln.str))

	test.str <- aln.str[22]
	test.str <- strsplit(test.str, split = " ")
	test.str <- unlist(test.str)
	test.str <- test.str[test.str != ""]
	test.str <- test.str[3]
	
	ref.str <- aln.str[24]
	ref.str <- strsplit(ref.str, split = " ")
	ref.str <- unlist(ref.str)
	ref.str <- ref.str[ref.str != ""]
	ref.str <- ref.str[3]
	
	return(c(test.str, ref.str))
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

getGuideAlignment <- function(inputF, guide, alnFolder, gnm = BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
{
	# DEFINE REFERENCE SEQUENCE
	refSeq <- guide

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
	sequences <- bed2sequence(readMat, g = gnm)
	sequences <- c(refSeq, sequences)# add refSeq


	###########
	# ALIGNMENT
	sm <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)

	alnList.norm <- mclapply(sequences, pairwiseAlignment, subject = refSeq,
		type = "local-global", substitutionMatrix = sm,
		gapOpening = 1, gapExtension = 1,
		mc.cores = NBCPU
		)

	alnList.rev <- mclapply(sequences, function(i)
		pairwiseAlignment(getRevComp(i), subject = refSeq,
			type = "local-global", substitutionMatrix = sm,
			gapOpening = 1, gapExtension = 1), mc.cores = NBCPU
		)

	print("TEST")
	print(score(alnList.norm[[1]]))
	print(score(alnList.rev[[1]]))

	# GET BEST ALIGNMENT
	alnList <- lapply(1:length(alnList.norm), function(i)
		getBestAln(alnList.norm[[i]], alnList.rev[[i]])
		)

	refSeqList <- rep("norm.", length(alnList))
	isRev <- unlist(lapply(1:length(alnList), function(i) identical(alnList[[i]], alnList.rev[[i]])))
	refSeqList[isRev] <- "rev. comp."
	refSeqList <- refSeqList[-1]# remove refSeq alignment	
	
	# SAVE
	dir.create(alnFolder, showWarnings = FALSE)

	dir1 <- file.path(tempdir(), "tmpDIR")
	dir.create(dir1)# TMP DIR
	
	guideTMP <- tempfile(tmpdir = dir1, fileext = ".txt")
	
	#writePairwiseAlignments(alnList[[1]],
  # 	file = file.path(alnFolder, "guide_aln_TMP.txt"))
	
	writePairwiseAlignments(alnList[[1]],
	                        file = guideTMP)

	alnList <- alnList[-1]# remove refSeq alignment
	names(alnList) <- paste0(inputName, "_aln", 1:length(alnList))	

	tmpFiles <- sapply(1:length(alnList), function(i) tempfile(tmpdir = dir1, fileext = ".txt"))
	
	lapply(1:length(alnList), function(i){
	  
	  #writePairwiseAlignments(alnList[[i]],
	  #                       file = file.path(alnFolder, paste0(names(alnList)[i], "_TMP.txt")))
	  writePairwiseAlignments(alnList[[i]], file = tmpFiles[i])
	  
	})
	

	# GET PROPER ALIGNMENT (FIX ISSUE WITH INDEL IN FIRST / LAST POSITION)
	#alnFiles <- list.files(alnFolder, pattern = "_TMP.txt")
	#names(alnFiles) <- gsub("_TMP.txt", "", alnFiles)

	#alnFiles <- alnFiles[names(alnList)]# Re-order
	
	#alnList.str <- lapply(alnFiles, function(i)
  #	getStrAlignment(file.path(alnFolder, i))
	#	)
	
	alnList.str <- lapply(tmpFiles, getStrAlignment)

	alnMat <- do.call(rbind, alnList.str)# 1st column: test, 2nd column: ref

	# GET ALIGNMENT SUMMARY STR
	alnSum <- unlist(lapply(1:nrow(alnMat), function(i)
		getAlnSumStr(alnMat[i, 1], alnMat[i, 2])
		))

	# START POSITION IN TEST SEQUENCE
	#alnStart <- as.numeric(unlist(lapply(alnFiles, function(i)
	#	getStartPosition(file.path(alnFolder, i)))))
	
	alnStart <- as.numeric(unlist(lapply(tmpFiles, function(i)
	  getStartPosition(i))))

	# RELATIVE COORDINATES
	middleCoord <- floor((readMat$end - readMat$start) / 2)
	middleCoord.abs <- readMat$start + middleCoord
	alnStart.rel <- alnStart - middleCoord
	alnStart.abs <- alnStart + readMat$start
	
	# STATISTIC
	statList <- lapply(alnList, getAlnStat)
	statMat <- do.call(rbind, statList)
	statMat <- cbind(names(alnList), alnMat, alnSum, statMat, refSeqList, middleCoord.abs, alnStart.abs, alnStart.rel)
	colnames(statMat) <- c("alignment.id", "pattern", "subject", "aln.sum", "score", "sq.hom", "nb.mism",
		"nb.ins", "length.ins", "nb.del", "length.del", "aln.chr", "aln.guide", "middleCoord", "aln.start.abs", "aln.start.rel")
	
	# MERGE READMAT AND STATMAT
	readMat.stat <- cbind(readMat, statMat)	
	readMat.stat[, c("score", "sq.hom", "nb.mism",
		"nb.ins", "length.ins", "nb.del", "length.del", "middleCoord", "aln.start.abs", "aln.start.rel")] <- apply(readMat.stat[, c("score", "sq.hom", "nb.mism", "nb.ins",
			"length.ins", "nb.del", "length.del", "middleCoord", "aln.start.abs", "aln.start.rel")], 2, as.numeric)	

	write.xlsx(readMat.stat, file.path(alnFolder, paste0(inputName, "_aln_stat.xlsx")), overwrite = TRUE)	

	# Density plot
	pdf(file.path(alnFolder, paste0(inputName, "_aln_stat_density.pdf")))
	par(mfrow=c(3,2))
	plot(density(readMat.stat$score), main = "Score")
	plot(density(readMat.stat$sq.hom), main = "Sequence Homology")
	plot(density(readMat.stat$nb.mism), main = "Nb. mismatches")
	plot(density(readMat.stat$nb.ins), main = "Nb. insertion")
	plot(density(readMat.stat$nb.del), main = "Nb. deletion")
	plot(density(readMat.stat$aln.start.rel), main = "Aln. relative start")
	dev.off()

	pdf(file.path(alnFolder, paste0(inputName, "_aln_stat_histogram.pdf")))
	par(mfrow=c(3,2))
	hist(readMat.stat$score, main = "Score")
	hist(readMat.stat$sq.hom, main = "Sequence Homology")
	hist(readMat.stat$nb.mism, main = "Nb. mismatches")
	hist(readMat.stat$nb.ins, main = "Nb. insertion")
	hist(readMat.stat$nb.del, main = "Nb. deletion")
	hist(readMat.stat$aln.start.rel, main = "Aln. relative start")
	dev.off()
	
	# DELETE TMP DIR
	unlink(dir1, recursive = T)
}

guidePlot <- function(inputFile, outputFile, hits = NULL, score = NULL, pv = NULL, ref = NULL, OMTonly = FALSE)
{
	###################
	# DISPLAY ALIGNMENT
	
	# READ INPUT FILE
	readMat.stat <- read.xlsx(inputFile)
	
	if(OMTonly) readMat.stat <- readMat.stat[readMat.stat$group == "OMT", ]

	# FILTER TOP SCORE
	readMat.stat <- readMat.stat[order(-readMat.stat$score), ]
	if(!is.null(hits)) readMat.stat <- readMat.stat[readMat.stat$hits > hits, ]
	if(!is.null(score)) readMat.stat <- readMat.stat[readMat.stat$score > score, ]# change score!!!
	if(!is.null(pv)) readMat.stat <- readMat.stat[readMat.stat$adj.pvalue < pv, ]

	if(nrow(readMat.stat)>1){
	
		alnDisplay <- lapply(1:nrow(readMat.stat), function(i)
			getAlnChar(as.character(readMat.stat[i, "pattern"]),
				as.character(readMat.stat[i, "subject"]),
				ref,
				".")
			)
		alnDisplay <- do.call(rbind, alnDisplay)
		#rownames(alnDisplay) <- paste(as.character(readMat.stat$chromosome), as.character(readMat.stat$start), as.character(readMat.stat$read), sep = ":")
		rownames(alnDisplay) <- paste(as.character(readMat.stat$chromosome), as.character(readMat.stat$start), as.character(readMat.stat$hits), sep = ":")
	
	
		# Add refseq on top
		ref.display <- rep(NA, ncol(alnDisplay))
		ref.display[seq(1, 46, by = 2)] <- unlist(strsplit(ref, split=""))
		ref.display[seq(2, 46, by = 2)] <- "0"
	
		alnDisplay <- rbind(GUIDE = ref.display, alnDisplay)
		alnDisplay <- data.frame(alnDisplay)
		
		write.xlsx(alnDisplay, gsub(".pdf", ".xlsx", outputFile), row.names = TRUE, overwrite = TRUE)
		#write.table(alnDisplay, gsub(".pdf", ".txt", outputFile), row.names = TRUE, sep = "\t", quote = FALSE)
		
		
		# HEATMAP
		ggmat <- data.table::melt(data.table(cbind(sq = rownames(alnDisplay), alnDisplay)), id.vars = "sq")

		ggmat$text <- ggmat$value
		ggmat$text[ggmat$text == 0] <- ""
		ggmat$value[ggmat$value == "."] <- "SAME"
		ggmat$value[ggmat$value == "0"] <- "zero"
		ggmat$value[grep("-", ggmat$value)] <- "neg"
		ggmat$value[grep("[0-9]", ggmat$value)] <- "pos"

		ggmat$sq <- factor(ggmat$sq, levels = rev(unique(ggmat$sq)))

		myColors <- c(SAME = "white",
					  A = "yellow",
					  C = "orange",
					  G = "green",
					  T = "pink",
					  N = "grey",
					  R = "snow2",
					  neg = "lightblue",
					  pos = "red",
					  zero = "white"
					  )

		p <- ggplot(ggmat, 
			aes(x = variable, y = sq, fill = factor(value))) + 
			geom_tile() + scale_fill_manual(values=myColors)
		p <- p + geom_text(aes(label=text), size=4)
		p <- p + geom_hline(yintercept=seq(1.5, length(ggmat$sq)+0.5))
		p <- p + theme(axis.text.x = element_blank(),
			axis.ticks = element_blank())
		p <- p + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())
		#p <- p + scale_fill_grey(na.value = "snow2")
	
		pdf.width = 8
		pdf.height = max(8*nrow(alnDisplay) / 25, 8)
		pdf(outputFile, width = pdf.width, height = pdf.height)
		#png(outputFile, units="in", width=pdf.width, height=pdf.height, res=300)
		plot(p)
		dev.off()
	
	}
}

logoPlot <- function(inputFile, outputFile, hits = NULL, score = NULL, pv = NULL, ref)
{
	# READ INPUT FILE
	readMat.stat <- read.xlsx(inputFile)

	# FILTER TOP SCORE
	readMat.stat <- readMat.stat[order(-readMat.stat$score), ]
	if(!is.null(hits)) readMat.stat <- readMat.stat[readMat.stat$hits > hits, ]
	if(!is.null(score)) readMat.stat <- readMat.stat[readMat.stat$score > score, ]# change score!!!
	if(!is.null(pv)) readMat.stat <- readMat.stat[readMat.stat$adj.pvalue < pv, ]
	
	if(nrow(readMat.stat)>1){
		alnDisplay <- lapply(1:nrow(readMat.stat), function(i)
		getAlnChar(as.character(readMat.stat[i, "pattern"]),
			as.character(readMat.stat[i, "subject"]),
			ref,
			NULL)
		)
	alnDisplay <- do.call(rbind, alnDisplay)
	rownames(alnDisplay) <- paste(as.character(readMat.stat$chromosome), as.character(readMat.stat$start), as.character(readMat.stat$read), sep = ":")
		
	alnDisplay.sub <- alnDisplay[, seq(1, ncol(alnDisplay), by = 2)]
	alnDisplay.ins <- alnDisplay[, seq(2, ncol(alnDisplay), by = 2)]


	alnDisplay.sub[grep("-", alnDisplay.sub)] <- "d"# deletion

	logoIn <- apply(alnDisplay.sub, 1, paste, collapse = "")
	#logoIn <- gsub("-", "", logoIn)

	p <- ggplot() +
		geom_logo(logoIn, method='p', 
		seq_type='other', namespace = unique(as.character(unlist(logoIn)))) + 
		theme_logo()
	p <- p + scale_fill_grey(na.value = "snow2")	

	pdf(outputFile, width = 12, height = 3)
	#png(outputFile, units="in", width=12, height=3, res=300)
	plot(p)
	dev.off()
	}


}



checkSignif <- function(readFile, hits = NULL, score = NULL, pv = NULL){
	readMat <- read.xslx(readFile, sheet = 1)
	hits.bool <- rep(TRUE, nrow(readMat))
	score.bool <- rep(TRUE, nrow(readMat))
	pv.bool <- rep(TRUE, nrow(readMat))
	
	if(!is.null(hits)) hits.bool[readMat$hits <= hits] <- FALSE
	if(!is.null(score)) score.bool[readMat$score <= score] <- FALSE
	if(!is.null(pv)) pv.bool[readMat$adj.pvalue <= pv] <- FALSE
	
	return(hits.bool & score.bool & pv.bool)
}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


if(FALSE)
{

# HEATMAP
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/test/results/Rpipeline/Translocome-G3_S12_L001_w250"))
inputFile <- "Translocome-G3_S12_L001_w250_aln_stat.xlsx"
score = NULL
pv = NULL
ref = "GTGAGTAGAGCGGAGGCAGGNRG"

# READ INPUT FILE
	readMat.stat <- read.xlsx(inputFile)

	# FILTER TOP SCORE
	readMat.stat <- readMat.stat[order(-readMat.stat$score), ]
	if(!is.null(score)) readMat.stat <- readMat.stat[readMat.stat$score > score, ]# change score!!!
	if(!is.null(pv)) readMat.stat <- readMat.stat[readMat.stat$adj.pvalue < pv, ]

	alnDisplay <- lapply(1:nrow(readMat.stat), function(i)
		getAlnChar(as.character(readMat.stat[i, "pattern"]),
			as.character(readMat.stat[i, "subject"]),
			ref,
			".")
		)
	alnDisplay <- do.call(rbind, alnDisplay)
	rownames(alnDisplay) <- paste(as.character(readMat.stat$chromosome), as.character(readMat.stat$start), as.character(readMat.stat$read), sep = ":")
	
	# Add refseq on top
	ref.display <- rep(NA, ncol(alnDisplay))
	ref.display[seq(1, 46, by = 2)] <- unlist(strsplit(ref, split=""))
	ref.display[seq(2, 46, by = 2)] <- "0"
	
	alnDisplay <- rbind(Ref = ref.display, alnDisplay)


# LOGO
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/test/results/Rpipeline/Translocome-G3_S12_L001_w250"))
inputFile <- "Translocome-G3_S12_L001_w250_aln_stat.xlsx"
score = NULL
pv = NULL
ref = "GTGAGTAGAGCGGAGGCAGGNRG"

# READ INPUT FILE
	readMat.stat <- read.xlsx(inputFile)

	# FILTER TOP SCORE
	readMat.stat <- readMat.stat[order(-readMat.stat$score), ]
	if(!is.null(score)) readMat.stat <- readMat.stat[readMat.stat$score > score, ]# change score!!!
	if(!is.null(pv)) readMat.stat <- readMat.stat[readMat.stat$adj.pvalue < pv, ]

	alnDisplay <- lapply(1:nrow(readMat.stat), function(i)
		getAlnCharLogo(as.character(readMat.stat[i, "pattern"]),
			as.character(readMat.stat[i, "subject"]),
			ref,
			NULL)
		)
	alnDisplay <- do.call(rbind, alnDisplay)
	rownames(alnDisplay) <- paste(as.character(readMat.stat$chromosome), as.character(readMat.stat$start), as.character(readMat.stat$read), sep = ":")
	
alnDisplay <- alnDisplay[1:25,]	
	
alnDisplay.sub <- alnDisplay[, seq(1, ncol(alnDisplay), by = 2)]
alnDisplay.ins <- alnDisplay[, seq(2, ncol(alnDisplay), by = 2)]


alnDisplay.sub[grep("-", alnDisplay.sub)] <- "d"

logoIn <- apply(alnDisplay.sub, 1, paste, collapse = "")
#logoIn <- gsub("-", "", logoIn)

p <- ggplot() +
	geom_logo(logoIn, method='p', 
    seq_type='other', namespace = unique(as.character(unlist(logoIn)))) + 
	theme_logo()
p <- p + scale_fill_grey(na.value = "snow2")

pdf(paste0(sampleName, "aln_logo.pdf"), width = 12, height = 3)
plot(p)
dev.off()




}







