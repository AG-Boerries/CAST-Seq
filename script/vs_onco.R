#library(biomaRt)
#library(org.Hs.eg.db)
#library(openxlsx)
#library(GenomicRanges)

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

#gene2bed <- function(entrez, upstream = 0, downstream = 0)
#{
#	filterlist <- list(entrez)
#	results <- getBM(attributes = c("entrezgene", "chromosome_name", "start_position", "end_position", "strand"),
#		filters = c("entrezgene"),
#	 	values = filterlist, mart = ensembl)
#	results$chromosome_name <- paste("chr", results$chromosome_name, sep = "") 	
#	chr <- c(paste("chr", 1:22, sep  =""), "chrX", "chrY")
#	results <- results[results$chromosome_name %in% chr,]
#	results <- results[!duplicated(results$entrezgene),]	
#	ids <- paste("ID", 1:nrow(results), sep = "")
#	
#	geneStrand <- ifelse(results$strand == 1, "+", "-")
#	bedMat <- data.frame(results$chromosome_name,
#		results$start_position - upstream,
#		results$end_position + downstream,
#		results$entrezgene,
#		geneStrand) 	
#	colnames(bedMat) <- NULL	
#	bedMat <- bedMat[!duplicated(bedMat[, 4]), ]
#	#write.table(bedMat, outFile, sep = "\t", row.names = FALSE, quote = FALSE)
#	return(bedMat)
#}

addOnco <- function(inputFile, oncoFile, width)
{
	# LOAD INPUT FILE
	readMat <- read.xlsx(inputFile, sheet = 1)
	readMat.sub <- readMat[, c("chromosome", "start", "end", "alignment.id")]
	readMat.sub$strand <- "*"
	
	# LOAD ONCOGENES
	oncoEntrez <- read.delim(oncoFile, header = FALSE, stringsAsFactors = FALSE)
	oncoEntrez <- as.character(t(oncoEntrez))
	oncoMat <- gene2bed(oncoEntrez, upstream = width, downstream = width)
	colnames(oncoMat) <- c("chromosome", "start", "end", "onco.entrez", "strand")
	

	# INTERSECT INPUT VS. ONCO
	readMat.gr <- makeGRangesFromDataFrame(readMat.sub, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = FALSE)
		
	oncoMat.gr <- makeGRangesFromDataFrame(oncoMat, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = TRUE)
	
	gr.ovl <- findOverlaps(query = readMat.gr, subject = oncoMat.gr, type = "any")	
	df.ovl <- data.frame(readMat.sub[queryHits(gr.ovl),],oncoMat[subjectHits(gr.ovl),])
	
	print(head(df.ovl))

	oncoList <- lapply(readMat.sub$alignment.id, function(i)
		df.ovl$onco.entrez[match(i, df.ovl$alignment.id)]
		)
	oncoList.symbol <- lapply(oncoList, function(i){
		if(is.na(i)) return(NA)
		return(entrez2symbol(i))
		})
	
	# RELATIVE DISTANCE
	startList <- lapply(readMat.sub$alignment.id, function(i)
		df.ovl$start.1[match(i, df.ovl$alignment.id)]
		)	
	startList <- unlist(lapply(startList, function(i) i + width))	
	
	rel.dist <- startList - (readMat.sub$start + (floor(readMat.sub$end - readMat.sub$start)))
		

	# SAVE
	oncoList.str <- lapply(oncoList, function(i){
		if(is.na(i)) return(NA)
		return(paste(i, collapse = ";"))
		})
	oncoList.str <- unlist(oncoList.str)	

	oncoList.symbol.str <- lapply(oncoList.symbol, function(i){
		if(is.na(i)) return(NA)
		return(paste(i, collapse = ";"))
		})
	oncoList.symbol.str <- unlist(oncoList.symbol.str)
	
	write.xlsx(cbind(readMat, onco.start.rel = rel.dist, onco.symbol = oncoList.symbol.str),
		gsub(".xlsx", "_ONCO.xlsx", inputFile), row.names = FALSE)

}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

