#library(Biostrings)
#library(openxlsx)
#library(data.table)
#library(ggplot2)
#library(ggseqlogo)
#library(textreadr)
#library(parallel)

# set python and load longest_common_substring function
#require(reticulate)
#use_python("/usr/bin/python")
#source_python("~/bitbucket/work/offTargets/Giando/lcs.py")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getRevComp <- function(x) as.character(reverseComplement(DNAString(x)))

addFlanking <- function(inputFile, otBed, size)
{
	# inputFile: alignment output (.xlsx)
	# flankFile: bed file containing up-/downstream flank locations
	# outputFile: output file of the analysis (.xlsx)

	#########################
	# DEFINE FLANKING REGIONS

	flankMat <- read.delim(otBed, header = FALSE)
	flankMat[,2] <- flankMat[,2] - size
	flankMat[,3] <- flankMat[,3] + size
	flankSeq <- bed2sequence(flankMat, g = GNM)
	flankSeq <- flankSeq[[1]]
	
	#flankSeq.up <- as.character(readDNAStringSet(flankLeftF))
	#flankSeq.down <- as.character(readDNAStringSet(flankRightF))

	##########################
	# LOAD REGIONS OF INTEREST
	readMat <- read.xlsx(inputFile, sheet = 1)

	sequences <- bed2sequence(readMat, g = GNM)
	sequences.rev <- mclapply(sequences, getRevComp, mc.cores = NBCPU)

	###########
	# ALIGNMENT

	sq.match <- mclapply(sequences, longest_common_substring, flankSeq, mc.cores = NBCPU)
	sq.rev.match <- mclapply(sequences.rev, longest_common_substring, flankSeq, mc.cores = NBCPU)

	##############################
	# GET LENGTH OF BEST ALIGNMENT

	length <- unlist(mclapply(sq.match, nchar, mc.cores = NBCPU))
	length.rev <- unlist(mclapply(sq.rev.match, nchar, mc.cores = NBCPU))


	#########
	# SUMMARY
	flankSm <- data.frame(flank = unlist(sq.match),
						  flank.length = length,
						  flank.rev = unlist(sq.rev.match),
						  flank.rev.length = length.rev		  
						  ) 
	write.xlsx(cbind(readMat, flankSm), gsub(".xlsx", "_FLANK.xlsx", inputFile))
}


addFlankingFromSq <- function(inputFile, hom.sq, hom.sq2 = NULL)
{
	# inputFile: alignment output (.xlsx)
	# flankFile: bed file containing up-/downstream flank locations
	# outputFile: output file of the analysis (.xlsx)

	##########################
	# LOAD REGIONS OF INTEREST
	readMat <- read.xlsx(inputFile, sheet = 1)

	sequences <- bed2sequence(readMat, g = GNM)
	sequences.rev <- mclapply(sequences, getRevComp, mc.cores = NBCPU)

	###########
	# ALIGNMENT

	sq.match <- mclapply(sequences, longest_common_substring, hom.sq, mc.cores = NBCPU)
	sq.rev.match <- mclapply(sequences.rev, longest_common_substring, hom.sq, mc.cores = NBCPU)

	##############################
	# GET LENGTH OF BEST ALIGNMENT

	length <- unlist(mclapply(sq.match, nchar, mc.cores = NBCPU))
	length.rev <- unlist(mclapply(sq.rev.match, nchar, mc.cores = NBCPU))


	#########
	# SUMMARY
	flankSm <- data.frame(flank = unlist(sq.match),
						  flank.length = length,
						  flank.rev = unlist(sq.rev.match),
						  flank.rev.length = length.rev		  
						  ) 
						  
	if(!is.null(hom.sq2)){
		###########
		# ALIGNMENT

		sq.match2 <- mclapply(sequences, longest_common_substring, hom.sq2, mc.cores = NBCPU)
		sq.rev.match2 <- mclapply(sequences.rev, longest_common_substring, hom.sq2, mc.cores = NBCPU)

		##############################
		# GET LENGTH OF BEST ALIGNMENT

		length2 <- unlist(mclapply(sq.match2, nchar, mc.cores = NBCPU))
		length2.rev <- unlist(mclapply(sq.rev.match2, nchar, mc.cores = NBCPU))


		#########
		# SUMMARY
		flankSm2 <- data.frame(flank2 = unlist(sq.match2),
							  flank2.length = length2,
							  flank2.rev = unlist(sq.rev.match2),
							  flank2.rev.length = length2.rev		  
							  ) 
							  
		flankSm <- cbind(flankSm, flankSm2)					  
	}					  
						  
	write.xlsx(cbind(readMat, flankSm), gsub(".xlsx", "_FLANK.xlsx", inputFile))
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################



