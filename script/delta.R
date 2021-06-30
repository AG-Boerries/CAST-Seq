


toNUM <- function(x) as.numeric(levels(x))[x]



# TO DO: get DELTA, use grange to intersect ots and sites, remove everything +/- otsD of ots


getDelta <- function(inputFile, otsF = NULL, otsD = 50, distance = 2500)
{
	# inputFile: bed file (output from fastq alignment pipeline)
	# output: generate a _delta.bed file and a _delta_density.pdf plot
	
	inputName <- gsub(".bed$", "", inputFile)
	
	# read the bed file
	bed.raw <-  as.data.frame(read.table(inputFile, header=FALSE,sep="\t", stringsAsFactors=FALSE, quote="", fill = TRUE))
	
	bed <- bed.raw	
	
	# merge identical coordinates
	print("Deduplicate")
	
	bed.plus <- bed[bed$V6 == "+", ]
	bed.minus <- bed[bed$V6 == "-", ]

	start.cutoff <- 3
	end.cutoff <- 3

	chr <- unique(bed.minus$V1)
	bed.minus.chr <- mclapply(chr, function(i){
		tempbed <- bed.minus[bed.minus$V1 == i, c(1:3, 6)]
		tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

		bedCD <- tempbed[1,]
		bedCD.reads <- 1
		
		if(nrow(tempbed) > 1){
			start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
			end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

			for(n in 2:nrow(tempbed)) # start for loop for each row
				{ 
				if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
					bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
					}else{
						bedCD <- rbind(bedCD, tempbed[n,])
						bedCD.reads <- c(bedCD.reads, 1)
						}	
				}
			}
		
		bedCD <- cbind(bedCD, reads = bedCD.reads)
		rownames(bedCD) <- NULL	
		colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
		return(bedCD)
		}, mc.cores = NBCPU
		)

	chr <- unique(bed.plus$V1)
	bed.plus.chr <- mclapply(chr, function(i){
		tempbed <- bed.plus[bed.plus$V1 == i, c(1:3, 6)]
		tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

		bedCD <- tempbed[1,]
		bedCD.reads <- 1
		if(nrow(tempbed) > 1){
			start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
			end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

			for(n in 2:nrow(tempbed)) # start for loop for each row
				{ 
				if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
					bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
					}else{
						bedCD <- rbind(bedCD, tempbed[n,])
						bedCD.reads <- c(bedCD.reads, 1)
						}	
				}
			}
		bedCD <- cbind(bedCD, reads = bedCD.reads)
		rownames(bedCD) <- NULL	
		colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
		return(bedCD)
		}, mc.cores = NBCPU
		)

	bed2 <- do.call(rbind, c(bed.minus.chr, bed.plus.chr))

	# change start for minus strand
	bed3 <- bed2   # create a copy of the bed file in bed2
	bed3[bed3[,4]== "-",2] <- bed3[bed3[,4]== "-",3] # if the strand (col 6) is negative copy the position in col 3 to col 2 
	bed3[,3] <- bed3[,2] # change col 3 to initial nucleotide 

	bed4 <- bed3[order(bed3[,1],bed3[,2],decreasing=F),]  # create bed 4 and reorder the start site

	# Calculate delta
	print("Calculate delta")
	delta <- c()   # empty object to fill during the for loop
	for(i in unique(bed4[,1])){ # start for loop for each chomosome
		tempbed <- subset(bed4, chromosome==i) # take the data for one choromose 
		chrDelta <- as.numeric(tempbed[1,2])
		if(nrow(tempbed) > 1) chrDelta <- c(chrDelta, as.numeric(tempbed[2:nrow(tempbed),2]) - as.numeric(tempbed[1:(nrow(tempbed)-1),2]))
		delta <- c(delta, chrDelta)
		}

	# Plot delta density
	#pdf(paste0(inputName, "_delta_density_RAW.pdf"))
	#tiff(paste0(inputName, "_delta_density.tiff"), units="in", width=5, height=5, res = 800)
	#png(paste0(inputName, "_delta_density.png"), units="px", width=1600, height=1600, res=300)
	#plot(density(log10(delta+1)), main = "", xlab = "log10(read base distance)")
	#dev.off()
	
	ggmat <- data.frame(dist = log10(delta+1))
	p <- ggplot(ggmat) + 
	  geom_density(aes(x=dist))
	p <- p + annotate("rect", xmin = -Inf,
	                  xmax = log10(distance),
	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	p <- p + geom_vline(xintercept=log10(distance), 
	                    color = "black", size=2)
	p <- p + theme_bw(base_size = 16)
	p <- p + xlab("log10(read base distance)")
	p <- p + ggtitle(paste0("Cutoff: ", distance))
	ggsave(plot = p, filename = paste0(inputName, "_delta_density_RAW.pdf"),
	       width = 6, height = 6)
	
	
	# Percentage of reads below / above the distance threshold
	bellow <- sum(delta <= distance)
	above <- sum(delta > distance)
	bellow.pc <- round(bellow / (bellow + above) * 100, digits = 2)
	above.pc <- round(above / (bellow + above) * 100, digits = 2)
	toxlsx <- data.frame(distance = distance,
	                     bellow.percentage = bellow.pc,
	                     above.percentage = above.pc)
	write.xlsx(toxlsx, paste0(inputName, "_delta_density_RAW.xlsx"), overwrite = TRUE)
	
	bed5 <- cbind(bed4, delta)
	colnames(bed5) <- c("chromosome", "start", "end", "strand", "read", "delta")
	write.table(bed5, paste0(inputName, "_delta_RAW.bed"), quote = FALSE,  sep = "\t", row.names = FALSE)
	
	
	#################################
	# REMOVE READS CLOSE TO ON TARGET
	
	# Remove ON-target reads OLD
	#if(!is.null(otsF)){
    #	ots <- read.delim(otsF, header = FALSE)
	#	chr <- as.character(ots[,1])
	#	start <- as.numeric(ots[,2])
	#	end <- as.numeric(ots[,3])
	#	strand <- as.character(ots[,6])
#
#		idx <- (bed[,1] == chr) & (bed[,6] == strand) &
#			( (bed[,2] >= (start - otsD) & bed[,2] <= (end + otsD)) | (bed[,3] >= (start - otsD) & bed[,3] <= (end + otsD)) )
#		print(sum(idx))
#		bed <- bed[!idx, ]
#	}
	
	
	# Remove ON-target reads NEW
	if(!is.null(otsF)){
		bed <- bed.raw
		
		print(paste0("Remove ON-target reads +/- ", otsD, "bp"))
		ots <- read.delim(otsF, header = FALSE)
		
		gr1 <- makeGRangesFromDataFrame(bed, seqnames.field = "V1",
			start.field = "V2", end.field = "V3", strand = "V6",
			keep.extra.columns = FALSE, ignore.strand = FALSE)
			
		gr2 <- makeGRangesFromDataFrame(ots, seqnames.field = "V1",
			start.field = "V2", end.field = "V3", strand = "V6",
			keep.extra.columns = FALSE, ignore.strand = FALSE)
	
		gr.ovl <- findOverlaps(query = gr1, subject = gr2, type = "any", maxgap = otsD, ignore.strand=TRUE) # REMOVE BOTH STRANDS !!!
		if(length(queryHits(gr.ovl)) != 0){
			bed <- bed[-queryHits(gr.ovl), ]	
			}
	print(paste0(" Remove +/- ", otsD))
	print(paste0("bed.raw: ", nrow(bed.raw)))
	print(paste0("bed (wo +/- ", otsD, "): ", nrow(bed)))
	
	
	# merge identical coordinates
	print("Deduplicate")
	
	bed.plus <- bed[bed$V6 == "+", ]
	bed.minus <- bed[bed$V6 == "-", ]

	start.cutoff <- 3
	end.cutoff <- 3

	chr <- unique(bed.minus$V1)
	bed.minus.chr <- mclapply(chr, function(i){
		tempbed <- bed.minus[bed.minus$V1 == i, c(1:3, 6)]
		tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

		bedCD <- tempbed[1,]
		bedCD.reads <- 1
		
		if(nrow(tempbed) > 1){
			start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
			end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

			for(n in 2:nrow(tempbed)) # start for loop for each row
				{ 
				if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
					bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
					}else{
						bedCD <- rbind(bedCD, tempbed[n,])
						bedCD.reads <- c(bedCD.reads, 1)
						}	
				}
			}
		
		bedCD <- cbind(bedCD, reads = bedCD.reads)
		rownames(bedCD) <- NULL	
		colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
		return(bedCD)
		}, mc.cores = NBCPU
		)

	chr <- unique(bed.plus$V1)
	bed.plus.chr <- mclapply(chr, function(i){
		tempbed <- bed.plus[bed.plus$V1 == i, c(1:3, 6)]
		tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

		bedCD <- tempbed[1,]
		bedCD.reads <- 1
		if(nrow(tempbed) > 1){
			start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
			end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

			for(n in 2:nrow(tempbed)) # start for loop for each row
				{ 
				if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
					bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
					}else{
						bedCD <- rbind(bedCD, tempbed[n,])
						bedCD.reads <- c(bedCD.reads, 1)
						}	
				}
			}
		bedCD <- cbind(bedCD, reads = bedCD.reads)
		rownames(bedCD) <- NULL	
		colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
		return(bedCD)
		}, mc.cores = NBCPU
		)

	bed2 <- do.call(rbind, c(bed.minus.chr, bed.plus.chr))

	# change start for minus strand
	bed3 <- bed2   # create a copy of the bed file in bed2
	bed3[bed3[,4]== "-",2] <- bed3[bed3[,4]== "-",3] # if the strand (col 6) is negative copy the position in col 3 to col 2 
	bed3[,3] <- bed3[,2] # change col 3 to initial nucleotide 

	bed4 <- bed3[order(bed3[,1],bed3[,2],decreasing=F),]  # create bed 4 and reorder the start site

	# Calculate delta
	print("Calculate delta")
	delta <- c()   # empty object to fill during the for loop
	for(i in unique(bed4[,1])){ # start for loop for each chomosome
		tempbed <- subset(bed4, chromosome==i) # take the data for one choromose 
		chrDelta <- as.numeric(tempbed[1,2])
		if(nrow(tempbed) > 1) chrDelta <- c(chrDelta, as.numeric(tempbed[2:nrow(tempbed),2]) - as.numeric(tempbed[1:(nrow(tempbed)-1),2]))
		delta <- c(delta, chrDelta)
		}
	
	
	}
	# Plot delta density
	#pdf(paste0(inputName, "_delta_density.pdf"))
	#plot(density(log10(delta+1)), main = "", xlab = "log10(read base distance)")
	#dev.off()
	
	ggmat <- data.frame(dist = log10(delta+1))
	p <- ggplot(ggmat) + 
	  geom_density(aes(x=dist))
	p <- p + annotate("rect", xmin = -Inf,
	                  xmax = log10(distance),
	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	p <- p + geom_vline(xintercept=log10(distance), 
	                    color = "black", size=2)
	p <- p + theme_bw(base_size = 16)
	p <- p + xlab("log10(read base distance)")
	p <- p + ggtitle(paste0("Cutoff: ", distance))
	ggsave(plot = p, filename = paste0(inputName, "_delta_density.pdf"),
	       width = 6, height = 6)
	
	# Percentage of reads below / above the distance threshold
	bellow <- sum(delta <= distance)
	above <- sum(delta > distance)
	bellow.pc <- round(bellow / (bellow + above) * 100, digits = 2)
	above.pc <- round(above / (bellow + above) * 100, digits = 2)
	toxlsx <- data.frame(distance = distance,
	                     bellow.percentage = bellow.pc,
	                     above.percentage = above.pc)
	write.xlsx(toxlsx, paste0(inputName, "_delta_density.xlsx"), overwrite = TRUE)
	
	
	bed5 <- cbind(bed4, delta)
	colnames(bed5) <- c("chromosome", "start", "end", "strand", "read", "delta")
	write.table(bed5, paste0(inputName, "_delta.bed"), quote = FALSE,  sep = "\t", row.names = FALSE)
	
}




getDeltaShuffle <- function(inputFile, nb, distance, genome.size)
{
	# create tmp dir
	dir.create(dir1 <- file.path(tempdir(), "shuffle"))
	
	# create "nb" shuffle bed files
	bedtmp <- tempfile(tmpdir = dir1, fileext = ".bed")
	write.table(read.delim(inputFile), bedtmp, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	
	mclapply(1:nb, function(i){
		shuffleBed(inFile = bedtmp,
		   myGenome = genome.size,
		   outFile = tempfile(tmpdir = dir1, fileext = ".bed"),
		   opt.string="")
		}, mc.cores = NBCPU)
	file.remove(bedtmp)
		
	# get delta
	shuffleFiles <- list.files(dir1, pattern = ".bed$", full.names = TRUE)
	
	deltaList <- mclapply(shuffleFiles, function(i){
		# read the bed file
		bed <-  as.data.frame(read.table(i, header=FALSE,sep="\t", stringsAsFactors=FALSE, quote="", fill = TRUE))

		# merge identical coordinates
		bed.plus <- bed[bed$V4 == "+", ]
		bed.minus <- bed[bed$V4 == "-", ]

		start.cutoff <- 3
		end.cutoff <- 2

		chr <- unique(bed.minus$V1)
		bed.minus.chr <- mclapply(chr, function(i){
			tempbed <- bed.minus[bed.minus$V1 == i, c(1:4)]
			tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

			bedCD <- tempbed[1,]
			bedCD.reads <- 1
			if(nrow(tempbed) > 1){
				start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
				end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

				for(n in 2:nrow(tempbed)) # start for loop for each row
					{ 
					if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
						bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
						}else{
							bedCD <- rbind(bedCD, tempbed[n,])
							bedCD.reads <- c(bedCD.reads, 1)
							}	
					}
				}
			bedCD <- cbind(bedCD, reads = bedCD.reads)
			rownames(bedCD) <- NULL	
			colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
			return(bedCD)
			}, mc.cores = 1
			)

		chr <- unique(bed.plus$V1)
		bed.plus.chr <- mclapply(chr, function(i){
			tempbed <- bed.plus[bed.plus$V1 == i, c(1:4)]
			tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

			bedCD <- tempbed[1,]
			bedCD.reads <- 1
			if(nrow(tempbed) > 1){
				start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
				end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

				for(n in 2:nrow(tempbed)) # start for loop for each row
					{ 
					if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
						bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
						}else{
							bedCD <- rbind(bedCD, tempbed[n,])
							bedCD.reads <- c(bedCD.reads, 1)
							}	
					}
				}
			bedCD <- cbind(bedCD, reads = bedCD.reads)
			rownames(bedCD) <- NULL	
			colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
			return(bedCD)
			}, mc.cores = 1
			)

		bed2 <- do.call(rbind, c(bed.minus.chr, bed.plus.chr))
	
		# change start for minus strand
		bed3 <- na.omit(bed2)   # create a copy of the bed file in bed2

		bed3[bed3[,4]== "-",2] <- bed3[bed3[,4]== "-",3] # if the strand (col 6) is negative copy the position in col 3 to col 2 
		bed3[,3] <- bed3[,2] # change col 3 to initial nucleotide 

		bed4 <- bed3[order(bed3[,1],bed3[,2],decreasing=F),]  # create bed 4 and reorder the start site

		# Calculate delta
		delta <- c()   # empty object to fill during the for loop
		for(i in unique(bed4[,1])){ # start for loop for each chomosome
			tempbed <- subset(bed4, chromosome==i) # take the data for one choromose 
			chrDelta <- as.numeric(tempbed[1,2])
			if(nrow(tempbed) > 1) chrDelta <- c(chrDelta, as.numeric(tempbed[2:nrow(tempbed),2]) - as.numeric(tempbed[1:(nrow(tempbed)-1),2]))
			delta <- c(delta, chrDelta)
			}
		return(delta)	
		}, mc.cores = NBCPU)
	
	# density	
	#pdf(gsub(".bed", ".shuffle_delta_density.pdf", inputFile))
	#tiff(gsub(".bed", ".shuffle_delta_density.pdf", inputFile), units="in", width=5, height=5, res = 800)
	#png(gsub(".bed", ".shuffle_delta_density.png", inputFile), units="px", width=1600, height=1600, res=300)
	#plot(density(log10(unlist(deltaList)+1)), main = paste0("Cutoff: ", distance), xlab = "log10(read base distance)")
	#abline(v=log10(distance), col="black", lwd=3)
	#dev.off()	
	
	
	# density ggplot
  ggmat <- data.frame(dist = log10(unlist(deltaList)+1))
	p <- ggplot(ggmat) + 
	  geom_density(aes(x=dist))
	p <- p + annotate("rect", xmin = -Inf,
	              xmax = log10(distance),
	              ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	p <- p + geom_vline(xintercept=log10(distance), 
	              color = "black", size=2)
	p <- p + theme_bw(base_size = 16)
	p <- p + xlab("log10(read base distance)")
	p <- p + ggtitle(paste0("Cutoff: ", distance))
	ggsave(plot = p, filename = gsub(".bed", ".shuffle_delta_density.pdf", inputFile),
	       width = 6, height = 6)
	
	# Percentage of reads below / above the distance threshold
	bellow <- sum(unlist(deltaList) <= distance)
	above <- sum(unlist(deltaList) > distance)
	bellow.pc <- round(bellow / (bellow + above) * 100, digits = 2)
	above.pc <- round(above / (bellow + above) * 100, digits = 2)
	toxlsx <- data.frame(distance = distance,
	                     bellow.percentage = bellow.pc,
	                     above.percentage = above.pc)
	write.xlsx(toxlsx, gsub(".bed", ".shuffle_delta_density.xlsx", inputFile), overwrite = TRUE)
	
	
	# Percentage of pre-defined delta
	deltaPC <- lapply(deltaList, function(d){
		return(sum(d < distance) / length(d))
		})
	deltaPC.m <- mean(unlist(deltaPC))	
	write(deltaPC.m, gsub(".bed", "_PERCENTAGE.txt", inputFile))
	
	# delete tmp dir
	unlink(dir1, recursive = T)
	dir.exists(dir1)
}






deltaDensity <- function(deltaF, otsF, distance)
{
	readMat <- read.delim(deltaF)

	ots <- read.delim(otsF, header = FALSE)
	#chr <- as.character(ots[,1])
	#start <- as.numeric(ots[,2])
	#end <- as.numeric(ots[,3])

	#idx <- (readMat[,1] == chr) & (readMat[,2] >= (start - distance)) & (readMat[,3] <= (end + distance))
	#idx <- (readMat[,1] == chr) &
	#	( (readMat[,2] >= (start - distance) & readMat[,2] <= (end + distance)) | (readMat[,3] >= (start - distance) & readMat[,3] <= (end + distance)) )
	#readMat.sub <- readMat[!idx, ]

	
	# Remove ON-target reads NEW
	if(!is.null(otsF)){
		#print(paste0("Remove ON-target reads +/- ", otsD, "bp"))
		ots <- read.delim(otsF, header = FALSE)
		
		gr1 <- makeGRangesFromDataFrame(readMat, seqnames.field = "chromosome",
			start.field = "start", end.field = "end",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
			
		gr2 <- makeGRangesFromDataFrame(ots, seqnames.field = "V1",
			start.field = "V2", end.field = "V3",
			keep.extra.columns = FALSE, ignore.strand = TRUE)
	
		gr.ovl <- findOverlaps(query = gr1, subject = gr2, type = "any", maxgap = otsD)
		if(length(queryHits(gr.ovl)) != 0){
			readMat.sub <- readMat[-queryHits(gr.ovl), ]
		}else readMat.sub <- readMat
	
	}else readMat.sub <- readMat


	pdf(gsub(".bed", "_density_minus_target.pdf", deltaF))
	plot(density(log10(readMat.sub[, "delta"]+1)), main = "", xlab = "log10(read base distance)")
	dev.off()
}



if(FALSE)
{
NBCPU <- 1

# TEST getDelta
inputFile <- file.path("/home/gandrieux/offTargets/Giando/test/Rev-196-1_S4_L001_Alignment.bed")
otsF <- file.path("/home/gandrieux/offTargets/Giando/pipelineGit/samples/HBG1Rev_196/data/ots.bed")
otsD <- 50
getDelta(inputFile, otsF, otsD)

deltaF <- file.path("/home/gandrieux/offTargets/Giando/test/Rev-196-1_S4_L001_Alignment_delta.bed")
deltaDensity(deltaF, otsF, distance)


# TEST getDeltaShuffle <- function(inputFile, nb, distance, genome.size)
inputFile <- "/home/gandrieux/offTargets/Giando/test/deltaShuffle/Rev-UT-2_S6_L001_Alignment_delta.bed" 
distance <- 1500
genome.size <- file.path("/home/gandrieux/offTargets/Giando/pipelineGit/annotations/human", "chrom.sizes")
nb <- 2
getDeltaShuffle(inputFile, nb, distance, genome.size)



############################
# new deduplication
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/results/")
bed <-  as.data.frame(read.table("G3_cat_Alignment.bed",header=FALSE,sep="\t", stringsAsFactors=FALSE, quote="", fill = TRUE))

bed.plus <- bed[bed$V6 == "+", ]
bed.minus <- bed[bed$V6 == "-", ]

start.cutoff <- 3
end.cutoff <- 2

chr <- unique(bed.minus$V1)
bed.minus.chr <- mclapply(chr, function(i){
	tempbed <- bed.minus[bed.minus$V1 == i, c(1:3, 6)]
	tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

	start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
	end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

	bedCD <- tempbed[1,]
	bedCD.reads <- 1
	for(n in 2:nrow(tempbed)) # start for loop for each row
		{ 
		if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
			bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
			}else{
				bedCD <- rbind(bedCD, tempbed[n,])
				bedCD.reads <- c(bedCD.reads, 1)
				}	
		}
	bedCD <- cbind(bedCD, reads = bedCD.reads)
	rownames(bedCD) <- NULL	
	colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
	return(bedCD)
	}, mc.cores = 1
	)

bed.plus.chr <- mclapply(chr, function(i){
	tempbed <- bed.plus[bed.plus$V1 == i, c(1:3, 6)]
	tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

	start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
	end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

	bedCD <- tempbed[1,]
	bedCD.reads <- 1
	for(n in 2:nrow(tempbed)) # start for loop for each row
		{ 
		if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
			bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
			}else{
				bedCD <- rbind(bedCD, tempbed[n,])
				bedCD.reads <- c(bedCD.reads, 1)
				}	
		}
	bedCD <- cbind(bedCD, reads = bedCD.reads)
	rownames(bedCD) <- NULL	
	colnames(bedCD) <- c("chromosome", "start", "end", "strand", "read")	
	return(bedCD)
	}, mc.cores = 1
	)


bed2 <- do.call(rbind, c(bed.minus.chr, bed.plus.chr))



bed.minus.chr22 <- bed.minus[bed.minus$V1 == "chr22", ]
bed.plus.chr22 <- bed.plus[bed.plus$V1 == "chr22", ]


tempbed <- bed.minus[bed.minus$V1 == chr, c(1:3, 6)]
	tempbed <- tempbed[order(tempbed$V2, tempbed$V3), ]

start.diff <- c(0, tempbed$V2[2:length(tempbed$V2)] - tempbed$V2[1:(length(tempbed$V2)-1)])
end.diff <- c(0, tempbed$V3[2:length(tempbed$V3)] - tempbed$V3[1:(length(tempbed$V3)-1)])

bedCD <- tempbed[1,]
bedCD.reads <- 1
for(n in 2:nrow(tempbed)) # start for loop for each row
	{ 
	if(start.diff[n] < start.cutoff & abs(end.diff[n]) < end.cutoff){
		bedCD.reads[length(bedCD.reads)] <- bedCD.reads[length(bedCD.reads)] + 1
		}else{
			bedCD <- rbind(bedCD, tempbed[n,])
			bedCD.reads <- c(bedCD.reads, 1)
			}	
	}
bedCD <- cbind(bedCD, reads = bedCD.reads)
rownames(bedCD) <- NULL	
colnames(bedCD) <- NULL	
	

##########################################################################################
# set the working directory
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/results/")

readMat <- read.delim("G3_cat_Alignment_delta.bed")

# On-target site
chr.ots <- "chr3"
start.ots <- 46372985 - 20000
end.ots <- 46373015 + 20000

idx <- (readMat[,1] == chr.ots) & (readMat[,2] >= start.ots) & (readMat[,3] <= end.ots)
readMat.sub <- readMat[!idx, ]

pdf("test_delta_minus_off_target.pdf")
plot(density(log10(readMat.sub[, "delta"])))
dev.off()



# UNTR

setwd("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/test/results/")

readMat <- read.delim("Transloc-UTFs_S25_L001_Alignment_delta.bed")

# On-target site
chr.ots <- "chr3"
start.ots <- 46372985 - 20000
end.ots <- 46373015 + 20000

idx <- (readMat[,1] == chr.ots) & (readMat[,2] >= start.ots) & (readMat[,3] <= end.ots)
readMat.sub <- readMat[!idx, ]

pdf("Transloc-UTFs_S25_L001_delta_minus_off_target.pdf")
plot(density(log10(readMat.sub[, "delta"])))
dev.off()



}