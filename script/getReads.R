


getSubBam <- function(bam.in, bam.out, region.str)
{
	options(scipen =99)
	command=paste0("samtools view -h ", bam.in, " ", region.str, " > ", bam.out)
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

getFwdBam <- function(bam.in, bam.out)
{
	options(scipen =99)
	command=paste0("samtools view -F 0x10 -h ", bam.in, " > ", bam.out)
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

getRevBam <- function(bam.in, bam.out)
{
	options(scipen =99)
	command=paste0("samtools view -f 0x10 -h ", bam.in, " > ", bam.out)
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}



#samtools view -F 0x10# Fwd
#samtools view -f 0x10# Rev

bam2fa <- function(bam.in, fa.out)
{
	options(scipen =99)
	command=paste0("samtools bam2fq -F 260  ", bam.in, " | seqtk seq -A > ", fa.out)
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}


getReads <- function(inputFile, bamFile, outputDir){
	# LOAD READMAT
	readMat <- read.xlsx(inputFile, sheet = 1)
	#readMat <- head(readMat)# FOR TEST ONLY
	
	# CREATE TMP DIR
	dir.create(dir1 <- file.path(tempdir(), "bam_tmp"))
	
	# LOOP OVER SITES	
	lapply(1:nrow(readMat), function(i){
		pos.str <- paste0(as.character(readMat[i,1]), ":", as.character(readMat[i,2]), "-", as.character(readMat[i,3]))
		
		# both strands
		bam.tmp  <- tempfile(tmpdir = dir1, fileext = ".bam")
		getSubBam(bamFile, bam.tmp, pos.str)
	
		#faName <- file.path(outputDir, paste0(gsub("\\:", "_", pos.str), "_BOTH.fa"))
		#bam2fa(bam.tmp, faName)
		
		# Fwd strand
		bam.tmp.fwd  <- tempfile(tmpdir = dir1, fileext = ".bam")
		getFwdBam(bam.tmp, bam.tmp.fwd)
	
		faName <- file.path(outputDir, paste0(gsub("\\:", "_", pos.str), "_Fwd.fa"))
		bam2fa(bam.tmp.fwd, faName)
		
		# Rev strand
		bam.tmp.rev  <- tempfile(tmpdir = dir1, fileext = ".bam")
		getRevBam(bam.tmp, bam.tmp.rev)
	
		faName <- file.path(outputDir, paste0(gsub("\\:", "_", pos.str), "_Rev.fa"))
		bam2fa(bam.tmp.rev, faName)
		
	})

	# DELETE TMP DIR
	unlink(dir1, recursive = T)
	dir.exists(dir1)
}


addNbReads <- function(inputFile, bamFile){
	# LOAD READMAT
	readMat <- read.xlsx(inputFile, sheet = 1)
	#readMat <- head(readMat)# FOR TEST ONLY
	
	# CREATE TMP DIR
	dir.create(dir1 <- file.path(tempdir(), "bam_tmp"))
	
	# LOOP OVER SITES	
	nb <- lapply(1:nrow(readMat), function(i){
		pos.str <- paste0(as.character(readMat[i,1]), ":", as.character(readMat[i,2]), "-", as.character(readMat[i,3]))
	
		# both strands
		bam.tmp  <- tempfile(tmpdir = dir1, fileext = ".bam")
		getSubBam(bamFile, bam.tmp, pos.str)
	
		# Fwd
		nbFwd <- as.numeric(nbPlusReadBam(bam.tmp))
		
		# Rev
		nbRev <- as.numeric(nbMinusReadBam(bam.tmp))
	
		return(c(nbFwd, nbRev))
		})
	nb <- do.call(rbind, nb)
	colnames(nb) <- c("plus.strand", "minus.strand")

	# SAVE
	readMat <- cbind(readMat, nb)
	write.xlsx(readMat, gsub(".xlsx", "_READS.xlsx", inputFile), row.names = FALSE, overwrite = TRUE)

	# DELETE TMP DIR
	unlink(dir1, recursive = T)
	dir.exists(dir1)
}



if(FALSE)
{
# DO NOT RUN
library(openxlsx)

# BAM
bamFile <- "/home/gandrieux/offTargets/Giando/pipelineGit/samples/HBG1Rev_196/results/fastq_aln/Rev-196-1_S4_L001_AlignmentSort.bam"

# SITES
inputFile <- "/home/gandrieux/offTargets/Giando/pipelineGit/samples/HBG1Rev_196/results/guide_aln/Rev-196-1_S4_L001_w250_aln_stat.xlsx"

# OUTPUT DIRECTORY
outputDir <- "/home/gandrieux/offTargets/Giando/test/getReads"

getReads(inputFile, bamFile, outputDir)

addNbReads(inputFile, bamFile)

}