
getOnTargetReads <- function(inputF, otsF)
{
	m <- read.xlsx(inputF, sheet = 1)
	m.sub <- m[, c("chromosome", "start", "end", "read")]

	ots.sub <- read.delim(otsF, header = FALSE)
	
	# CONVERT TO GRANGES
	m.gr <- makeGRangesFromDataFrame(m.sub, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", ignore.strand = TRUE,
		keep.extra.columns = FALSE)
		
	ots.gr <- makeGRangesFromDataFrame(ots.sub, seqnames.field = "V1",
		start.field = "V2", end.field = "V3", ignore.strand = TRUE,
		keep.extra.columns = FALSE)
	
	# INTERSECT
	gr.ovl <- findOverlaps(query = m.gr, subject = ots.gr, type = "any", ignore.strand= TRUE)	
	
	# GET NUMBER OF READS ON ON-TARGET SITE
	nbReads <- sum(m.sub$read[unique(queryHits(gr.ovl))])
	
	return(nbReads)
}


getOGI <- function(rawFastq.de, filtFastq.de, groupF.de, rawFastq.u, filtFastq.u, groupF.u, otsF, outF)
{
	groupM.de <- read.xlsx(groupF.de, sheet = 1)
	groupM.u <- read.xlsx(groupF.de, sheet = 1)
	
	TOT.de <- as.numeric(nbReadFastqgz(rawFastq.de))
	RPF.de <- as.numeric(nbReadFastqgz(filtFastq.de))
	RON.de <- getOnTargetReads(groupF.de, otsF)

	TOT.u <- as.numeric(nbReadFastqgz(rawFastq.u))
	RPF.u <- as.numeric(nbReadFastqgz(filtFastq.u))
	RON.u <- getOnTargetReads(groupF.u, otsF)

	print(paste0("TOT.de: ", TOT.de))
	print(paste0("RPF.de: ", RPF.de))
	print(paste0("RON.de: ", RON.de))

	print(paste0("TOT.u: ", TOT.u))
	print(paste0("RPF.u: ", RPF.u))
	print(paste0("RON.u: ", RON.u))

	print(paste0("(RPF.de - RON.de) / TOT.de): ", (RPF.de - RON.de) / TOT.de))
	print(paste0("(RPF.u - RON.u) / TOT.u): ", (RPF.u - RON.u) / TOT.u))

	ogi <- ((RPF.de - RON.de) / TOT.de) / ((RPF.u - RON.u) / TOT.u)

	write(ogi, outF)
}






############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE)
{
# DEFINE INPUTS

rawFastq.de <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/data/fastq", "FANCF-withDecoyGel_S1_L001_R2_001.fastq")
filtFastq.de <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/fastq_aln", "FANCF-withDecoyGel_S1_L001_pos.fastq")
groupF.de <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/guide_aln", "FANCF-withDecoyGel_S1_L001_w250_aln_stat_FLANK_GROUP.xlsx")

rawFastq.u <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/data/fastq", "UT-FANCF-withDecoyGel_S2_L001_R2_001.fastq")
filtFastq.u <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/fastq_aln", "UT-FANCF-withDecoyGel_S2_L001_pos.fastq")
groupF.u <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/results/guide_aln", "UT-FANCF-withDecoyGel_S2_L001_w250_aln_stat_FLANK_GROUP.xlsx")

otsF <- file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/FANCF/data/", "ots.bed")

getOnTargetReads(groupF.u, otsF)

# GET OGI
getOGI(rawFastq.de, filtFastq.de, groupF.de, rawFastq.u, filtFastq.u, groupF.u, otsF)
}