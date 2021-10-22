#library(GenomicRanges)
#library(openxlsx)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################


addNormalizeCount <- function(m, libSize){
  #m <- read.delim(inputFile)
  reads <- m$read
  hits <- m$hits
  width <- m$width
  
  # PER MILLION (CPM like)
  reads.cpm <- (reads*10^6) / libSize
  hits.cpm <- (hits*10^6) / libSize

  # PER MILLION PER BP (TPM like)
  reads.tpm <- ((reads / width)*10^6) / libSize
  hits.tpm <- ((hits / width)*10^6) / libSize
  
  # ADD FEATURES
  m$read.perM <- log2(reads.cpm)
  m$read.perM.perBP <- reads.tpm
  m$hits.perM <- log2(hits.cpm)
  m$hits.perM.perBP <- hits.tpm
  
  #write.table(m, inputFile, sep = "\t", quote = FALSE, row.names = FALSE)
  return(m)
}


getLimits <- function(m, chrFile)
{
	chr.size <- read.delim(chrFile, stringsAsFactors = FALSE, header = FALSE)
	return(chr.size[match(m[,1], chr.size[,1]), 2])
}

expand <- function(readMat, size, chrFile)
{
	limits <- getLimits(readMat, chrFile)
	
	# remove sites on un-defined chr
	readMat <- readMat[!is.na(limits), ]
	limits <- limits[!is.na(limits)]
	
	readMat$start <- readMat$start - size
	readMat$start[readMat$start < 1] <- 1
	readMat$end <- readMat$end + size
	readMat$end[readMat$end > limits] <- limits[readMat$end > limits]
	return(readMat)
}


#results <- matrix(data=NA,ncol=9,nrow = 1)
#    colnames(results) <- c('Term','Count','Size.univ', 'Size.tot','p-value','adj.P.Val','odds ratio','Entrez','Symbol')
#    geneSet <- intersect(universe, geneSets[[i]])
#    a <- length(intersect(DEgenes,geneSet))
#    b <- length(setdiff(DEgenes,intersect(DEgenes,geneSet)))
#    c <- length(setdiff(geneSet,intersect(DEgenes,geneSet)))
#    d <- length(setdiff(universe,DEgenes)) - c
#    contigency.matrix <- cbind(c(a,b),c(c,d))
#    res <- fisher.test(contigency.matrix,alternative = 'greater')


doFisher <- function(tumor.reads, tumor.tot, ctl.reads, ctl.tot)
{
	a <- tumor.reads
	b <- ctl.reads
	c <- tumor.tot - tumor.reads
	d <- ctl.tot - ctl.reads
	contigency.matrix <- cbind(c(a,b),c(c,d))
    res <- fisher.test(contigency.matrix, alternative = 'greater')
	
	resMat <- data.frame(OddRatio = res$estimate[[1]], pvalue = res$p.value)
	
	return(resMat)

}


overlapFiltering <- function(m)
{
	m.unique <- m[, c("chromosome", "start", "end", "strand", "read", "hits")]
	m.unique <- m.unique[!duplicated(as.matrix(m.unique)), ]
	
	reads.ctl <- c()
	hits.ctl <- c()
	
	for(i in 1:nrow(m.unique))
		{
		idx <- m$chromosome == m.unique[i, "chromosome"] &
			   m$start == m.unique[i, "start"] &
			   m$end == m.unique[i, "end"]
			   
		reads.ctl <- c(reads.ctl, sum(m$"read.1"[idx]))
		hits.ctl <- c(hits.ctl, sum(m$"hits.1"[idx]))
		}
	
	m.unique <- cbind(m.unique, read.ctl = reads.ctl, hits.ctl = hits.ctl)
	return(m.unique)
}


doEnrichment <- function(testFile, refFile, nbTest, nbRef, size, chrFile)
{
	# LOAD HITS
	cl.test <- read.delim(testFile)
	cl.ref <- read.delim(refFile)
	
	print(head(cl.test))
	
	print(head(cl.ref))

	# EXPAND WINDOW
	cl.test <- expand(cl.test, size, chrFile)
	cl.ref <- expand(cl.ref, size, chrFile)

	# CONVERT TO GRANGES
	cl.test.gr <- makeGRangesFromDataFrame(cl.test, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = TRUE)
		
	cl.ref.gr <- makeGRangesFromDataFrame(cl.ref, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = TRUE)
	
	# INTERSECT
	gr.ovl <- findOverlaps(query = cl.test.gr, subject = cl.ref.gr, type = "any", ignore.strand= TRUE)
	
	if(length(queryHits(gr.ovl)) == 0){
		df.spe <- cl.test
		df.spe <- cbind(df.spe, read.ctl = 1, hits.ctl = 1)
		df.final <- df.spe	
	}else if(length(queryHits(gr.ovl)) == nrow(cl.test)){
	  df.ovl <- data.frame(cl.test[queryHits(gr.ovl),], cl.ref[subjectHits(gr.ovl),])
	  # MERGE MULTIPLE OVERLAP
	  df.ovl <- overlapFiltering(df.ovl)
	  df.spe <- data.frame()
	  df.final <- df.ovl
	}else{
		df.ovl <- data.frame(cl.test[queryHits(gr.ovl),], cl.ref[subjectHits(gr.ovl),])
		df.spe <- cl.test[-queryHits(gr.ovl),]

		# MERGE MULTIPLE OVERLAP
		df.ovl <- overlapFiltering(df.ovl)

		df.spe <- cbind(df.spe, read.ctl = 1, hits.ctl = 1)

		df.final <- rbind(df.spe, df.ovl)
	
	}
	
	# GET WIDTH
	df.final <- cbind(df.final, width.raw = (df.final$end - df.final$start) - (size * 2),
		width = df.final$end - df.final$start)

	# HYPERG TEST
	fhList <- lapply(1:nrow(df.final), function(i){
		doFisher(df.final[i, "read"], nbTest, df.final[i, "read.ctl"], nbRef)
		})
	fhList <- do.call(rbind, fhList)

	df.final <- cbind(df.final, fhList)

	# ADJUSTED PVALUE
	df.final$adj.pvalue <- p.adjust(df.final$pvalue, method = "BH")

	# MARK ARTIFICIAL READS
	df.final$artificial <- NA
	if(nrow(df.spe) != 0) df.final$artificial[1:nrow(df.spe)] <- "yes"

	# SORT ACCORDING TO PVALUE
	df.final <- df.final[order(df.final$pvalue, -df.final$OddRatio), ]

	# NORMALIZE COUNT
	df.final <- addNormalizeCount(m = df.final, libSize = nbTest)
	
	# SAVE
	outName <- paste0(gsub("_Alignment_hits.bed$", "", testFile), "_w", size, ".xlsx")
	write.xlsx(df.final, outName, overwrite = TRUE)
}


doEnrichmentDefault <- function(testFile, nbTest, size)
{
	nbRef <- nbTest

	# LOAD HITS
	cl.test <- read.delim(testFile)
	
	print(head(cl.test))

	# EXPAND WINDOW
	cl.test <- expand(cl.test, size)

	# CONVERT TO GRANGES
	cl.test.gr <- makeGRangesFromDataFrame(cl.test, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = TRUE)
	
	df.final <- cbind(cl.test, read.ctl = 1, hits.ctl = 1)
	
	# GET WIDTH
	df.final <- cbind(df.final, width = df.final$end - df.final$start)

	# HYPERG TEST
	fhList <- lapply(1:nrow(df.final), function(i)
		doFisher(df.final[i, "read"], nbTest, df.final[i, "read.ctl"], nbRef)
		)
	fhList <- do.call(rbind, fhList)

	df.final <- cbind(df.final, fhList)

	# ADJUSTED PVALUE
	df.final$adj.pvalue <- p.adjust(df.final$pvalue, method = "BH")

	# MARK ARTIFICIAL READS
	df.final$artificial <- "yes"

	# SORT ACCORDING TO PVALUE
	df.final <- df.final[order(df.final$pvalue, -df.final$OddRatio), ]

	# SAVE
	outName <- paste0(gsub("_Alignment_hits.bed$", "", testFile), "_w", size, ".xlsx")
	write.xlsx(df.final, outName, overwrite = TRUE)
}




############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


if(FALSE)
{
# TEST
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/results"))
testFile <- "G3_custom_Alignment_cluster.bed"
refFile <- "UT_custom_Alignment_cluster.bed"
nbTest <- nrow(read.delim("G3_custom_Alignment.bed", header = FALSE))
nbRef <- nrow(read.delim("UT_custom_Alignment.bed", header = FALSE))
size <- 250


testFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/GENEWIZ_90-556214738/EMD/EMD1_4/EMD-sample13-1/results/guide_aln/EMD-sample13-1_Alignment_hits.bed")
refFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/GENEWIZ_90-556214738/EMD/EMD1_4/EMD-sample13-1/results/guide_aln/EMD-sample23-1_Alignment_hits.bed")
chrFile <- file.path("~/Research/CASTSeq/pipelineGit/annotations/human/chrom.sizes")
nbTest <- 10^6
nbRef <- 10^6
size <- 250

###############################################
# OLD



# NUMBER OF READS
G3 <- 4095770
UT <- 3951966

# LOAD CLUSTERS
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/data"))
cl.G3 <- read.delim("G3_clusters.txt")
cl.UT <- read.delim("UT_clusters.txt")

# EXPAND WINDOW
my.size <- 250
cl.G3 <- expand(cl.G3, my.size)
cl.UT <- expand(cl.UT, my.size)


# CONVERT TO GRANGES
cl.G3.gr <- makeGRangesFromDataFrame(cl.G3, seqnames.field = "chromosome",
	start.field = "start", end.field = "end", strand.field = "strand",
	keep.extra.columns = TRUE)
		
cl.UT.gr <- makeGRangesFromDataFrame(cl.UT, seqnames.field = "chromosome",
	start.field = "start", end.field = "end", strand.field = "strand",
	keep.extra.columns = TRUE)
	
	
# INTERSECT
gr.ovl <- findOverlaps(query = cl.G3.gr, subject = cl.UT.gr, type = "any")	
df.ovl <- data.frame(cl.G3[queryHits(gr.ovl),], cl.UT[subjectHits(gr.ovl),])
df.spe <- cl.G3[-queryHits(gr.ovl),]

# MERGE MULTIPLE OVERLAP
df.ovl <- overlapFiltering(df.ovl)

df.spe <- cbind(df.spe, collapse.ctl = 1, collapseCluster.ctl = 1)

df.final <- rbind(df.spe, df.ovl)

# HYPERG TEST
fhList <- lapply(1:nrow(df.final), function(i)
	doFisher(df.final[i, "collapse"], G3, df.final[i, "collapse.ctl"], UT)
	)
fhList <- do.call(rbind, fhList)

df.final <- cbind(df.final, fhList)

# ADJUSTED PVALUE
df.final$adj.pvalue <- p.adjust(df.final$pvalue, method = "BH")

# MARK ARTIFICIAL READS
df.final$artificial <- NA
df.final$artificial[1:nrow(df.spe)] <- "yes"

# SORT ACCORDING TO PVALUE
df.final <- df.final[order(df.final$pvalue, -df.final$OddRatio), ]


# SAVE
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/data"))
write.xlsx(df.final, paste0("G3_clusters_W", my.size, "_Fisher.xlsx"), overwrite = TRUE)


}




