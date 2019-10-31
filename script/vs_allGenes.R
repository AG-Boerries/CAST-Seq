library(biomaRt)
library(org.Hs.eg.db)
library(openxlsx)
library(GenomicRanges)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")# 38


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

symbol2entrez <- function(symbol)
{
	entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
	entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
	entrez <- entrez[!is.na(entrez)]
	return(entrez)
}

entrez2symbol <- function(entrez)
{
	symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
	symbol <- unlist(lapply(symbol, function(i) return(i[1])))
	return(symbol)
}

getHumanEntrez <- function()
{
	x <- org.Hs.egSYMBOL
	mapped_genes <- mappedkeys(x)
	xx <- as.list(x[mapped_genes])
	return(names(xx))
}

gene2bed <- function(entrez, upstream = 0, downstream = 0)
{
	filterlist <- list(entrez)
	results <- getBM(attributes = c("entrezgene_id", "chromosome_name", "start_position", "end_position", "strand"),
		filters = c("entrezgene_id"),
	 	values = filterlist, mart = ensembl)
	results$chromosome_name <- paste("chr", results$chromosome_name, sep = "") 	
	chr <- c(paste("chr", 1:22, sep  =""), "chrX", "chrY")
	results <- results[results$chromosome_name %in% chr,]
	results <- results[!duplicated(results$entrezgene_id),]	
	ids <- paste("ID", 1:nrow(results), sep = "")
	
	geneStrand <- ifelse(results$strand == 1, "+", "-")
	bedMat <- data.frame(results$chromosome_name,
		results$start_position - upstream,
		results$end_position + downstream,
		results$entrezgene_id,
		geneStrand) 	
	colnames(bedMat) <- NULL	
	bedMat <- bedMat[!duplicated(bedMat[, 4]), ]
	#write.table(bedMat, outFile, sep = "\t", row.names = FALSE, quote = FALSE)
	return(bedMat)
}

addGenes <- function(inputFile, oncoFile, geneMat, genes.width = 0, site.width = 100000)
{
	# LOAD ONCOGENES
	oncoEntrez <- read.delim(oncoFile, header = FALSE, stringsAsFactors = FALSE)
	oncoEntrez <- as.character(t(oncoEntrez))

	# LOAD INPUT FILE
	readMat <- read.xlsx(inputFile, sheet = 1)
	readMat.sub <- readMat[, c("chromosome", "start", "end", "alignment.id")]
	
	# EXTEND SITES +/- site.width
	#middle <- readMat.sub$start + (readMat.sub$end - readMat.sub$start)
	middle <- readMat$middleCoord
	readMat.sub$start <- middle - site.width
	readMat.sub$end <- middle + site.width
	
	# LOAD ONCOGENES
	allentrez <- getHumanEntrez()
	#geneMat <- gene2bed(allentrez, upstream = genes.width, downstream = genes.width)
	colnames(geneMat) <- c("chromosome", "start", "end", "gene.entrez", "strand")
	geneMat$start <- geneMat$start - genes.width
	geneMat$end <- geneMat$end + genes.width
	geneMat <- geneMat[, -5]# remove strand
	geneMat$gene.symbol <- entrez2symbol(geneMat$gene.entrez)
	
	#write.xlsx(geneMat, "/home/gandri/org.Hs.eg.geneCoord.xlsx")

	# INTERSECT INPUT VS. ONCO
	readMat.gr <- makeGRangesFromDataFrame(readMat.sub, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = FALSE, ignore.strand = TRUE)
		
	geneMat.gr <- makeGRangesFromDataFrame(geneMat, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	
	gr.ovl <- findOverlaps(query = readMat.gr, subject = geneMat.gr, type = "any")	
	df.ovl <- data.frame(readMat.sub[queryHits(gr.ovl),],geneMat[subjectHits(gr.ovl),])
	
	print(head(df.ovl))

	geneChr <- lapply(readMat.sub$alignment.id, function(i){
		isAligned <- df.ovl$alignment.id == i
		if(sum(isAligned)==0) return(NA)
				
		symbol.current <- df.ovl$gene.symbol[isAligned]
		start.current <- df.ovl$start.1[isAligned]
		end.current <- df.ovl$end.1[isAligned]
		
		chrList <- lapply(1:length(symbol.current), function(j)
			paste0(symbol.current[j], " (", start.current[j], ":", end.current[j], ")")
			)
		return(paste(unlist(chrList), collapse = "; "))
		})
		
	geneList <- lapply(readMat.sub$alignment.id, function(i)
		df.ovl$gene.entrez[df.ovl$alignment.id == i]
		)
	oncoChr <- lapply(geneList, function(i){
		if(length(i)==0) return(NA)
		if(is.na(i)) return(NA)
		onco.entrez <- intersect(i, oncoEntrez)
		if(length(onco.entrez)==0) return(NA)
		return(paste(entrez2symbol(onco.entrez), collapse = "; "))
		})
	
	#write.xlsx(cbind(readMat, "Genes.within.100kb" = unlist(geneChr), "Oncogenes.within.100kb" = unlist(oncoChr)),
	#	gsub(".xlsx", "_final.xlsx", inputFile), row.names = FALSE)
	
	write.xlsx(cbind(readMat, "Genes.within.100kb" = unlist(geneChr), "Oncogenes.within.100kb" = unlist(oncoChr)),
		inputFile, row.names = FALSE)

}


addGenesTALEN <- function(inputFile, oncoFile, geneMat, genes.width = 0, site.width = 100000)
{
	# LOAD ONCOGENES
	oncoEntrez <- read.delim(oncoFile, header = FALSE, stringsAsFactors = FALSE)
	oncoEntrez <- as.character(t(oncoEntrez))

	# LOAD INPUT FILE
	readMat <- read.xlsx(inputFile, sheet = 1)
	readMat.sub <- readMat[, c("chromosome", "start", "end", "alignment.id")]
	
	# EXTEND SITES +/- site.width
	middle <- readMat.sub$start + (readMat.sub$end - readMat.sub$start)
	#middle <- readMat$middleCoord
	readMat.sub$start <- middle - site.width
	readMat.sub$end <- middle + site.width
	
	# LOAD ONCOGENES
	allentrez <- getHumanEntrez()
	#geneMat <- gene2bed(allentrez, upstream = genes.width, downstream = genes.width)
	colnames(geneMat) <- c("chromosome", "start", "end", "gene.entrez", "strand")
	geneMat$start <- geneMat$start - genes.width
	geneMat$end <- geneMat$end + genes.width
	geneMat <- geneMat[, -5]# remove strand
	geneMat$gene.symbol <- entrez2symbol(geneMat$gene.entrez)

	# INTERSECT INPUT VS. ONCO
	readMat.gr <- makeGRangesFromDataFrame(readMat.sub, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = FALSE, ignore.strand = TRUE)
		
	geneMat.gr <- makeGRangesFromDataFrame(geneMat, seqnames.field = "chromosome",
		start.field = "start", end.field = "end",
		keep.extra.columns = TRUE, ignore.strand = TRUE)
	
	gr.ovl <- findOverlaps(query = readMat.gr, subject = geneMat.gr, type = "any")	
	df.ovl <- data.frame(readMat.sub[queryHits(gr.ovl),],geneMat[subjectHits(gr.ovl),])
	
	print(head(df.ovl))

	geneChr <- lapply(readMat.sub$alignment.id, function(i){
		isAligned <- df.ovl$alignment.id == i
		if(sum(isAligned)==0) return(NA)
				
		symbol.current <- df.ovl$gene.symbol[isAligned]
		start.current <- df.ovl$start.1[isAligned]
		end.current <- df.ovl$end.1[isAligned]
		
		chrList <- lapply(1:length(symbol.current), function(j)
			paste0(symbol.current[j], " (", start.current[j], ":", end.current[j], ")")
			)
		return(paste(unlist(chrList), collapse = "; "))
		})
		
	geneList <- lapply(readMat.sub$alignment.id, function(i)
		df.ovl$gene.entrez[df.ovl$alignment.id == i]
		)
	oncoChr <- lapply(geneList, function(i){
		if(length(i)==0) return(NA)
		if(is.na(i)) return(NA)
		onco.entrez <- intersect(i, oncoEntrez)
		if(length(onco.entrez)==0) return(NA)
		return(paste(entrez2symbol(onco.entrez), collapse = "; "))
		})
	
	#write.xlsx(cbind(readMat, "Genes.within.100kb" = unlist(geneChr), "Oncogenes.within.100kb" = unlist(oncoChr)),
	#	gsub(".xlsx", "_final.xlsx", inputFile), row.names = FALSE)
	
	write.xlsx(cbind(readMat, "Genes.within.100kb" = unlist(geneChr), "Oncogenes.within.100kb" = unlist(oncoChr)),
		inputFile, row.names = FALSE)

}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


if(FALSE)
{

inputFile <- file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D1/results/overlap_aln/", "G3_WT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx")
genes.width <- 0
site.width <- 100000
oncoFile <- file.path("~/cluster/master/offTargets/Giando/pipeline/annotations/CancerGenesList_ENTREZ.txt")

# G3 WT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D1/results/overlap_aln/"))
addGenes("G3_WT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)


# G3 HiFi
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_HiFi_D1/results/overlap_aln/"))
addGenes("G3_HiFi_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)

# 399 WT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/WT399_D1/results/overlap_aln/"))
addGenes("WT399_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)


# 399 HiFi
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/HiFi399_D1/results/overlap_aln/"))
addGenes("HiFi399_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)


# FANCF
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/FANCF_decoy_rep1/results/overlap_aln/"))
addGenes("FANCF_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)

# VEGFA
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/VEGFA_decoy_rep1/results/overlap_aln/"))
addGenes("VEGFA_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)


# OT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/OT6/results/overlap_aln/"))
addGenes("OT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)


# HBB
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/HBBTal_D1/results/overlap_aln/"))
addGenes("HBBTal_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx", oncoFile = oncoFile,
	genes.width = 0, site.width = 100000)

















}









