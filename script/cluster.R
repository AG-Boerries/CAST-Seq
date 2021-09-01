
getThreshold <- function(deltaFile, cutoff = 0.05) 
{
	bed <- read.delim(deltaFile)
	d <- bed$delta
	d <- sort(d)
	idx <- floor(length(d) * cutoff)
	return(d[idx])
}




getHits <- function(testFile, threshold)# get hits
{
	print("getHits")
	CUTOFF <- threshold
	#CUTOFF  <- 1500

	bed <- read.delim(testFile, stringsAsFactors = FALSE)
	clusterTab <- c() # for each chomosome calculate the difference between one start site and the previous one
	colTOkeep <- c(1:5)
	colClust <- length(colTOkeep) +1
	for(i in unique(bed[,1])){ # start for loop for each chomosome
		tempbed <- subset(bed,  chromosome==i) # take the data for one choromose 
		for(n in 1:nrow(tempbed)){   # start for loop for each row
			if( n == 1){ 
				clusterTab <- rbind(clusterTab, unlist(c(tempbed[1,colTOkeep], 1)))  # report read in the final table
				}
			if(n >= 2 ){
				if( tempbed$delta[n] <= CUTOFF ){
					clusterTab[nrow(clusterTab), 3] <- tempbed$end[n]
					clusterTab[nrow(clusterTab), colClust-1 ] <- as.numeric(clusterTab[nrow(clusterTab), colClust-1]) +  tempbed$read[n] 
					clusterTab[nrow(clusterTab), colClust] <- as.numeric(clusterTab[nrow(clusterTab), colClust]) +1  # add collapsed sequence
				}
				if( tempbed$delta[n] >= CUTOFF ){ 
					clusterTab <- rbind(clusterTab, unlist(c(tempbed[n,colTOkeep], 1)))
				}
			}
		}
	}
	colnames(clusterTab)[colClust] <- "hits"
	#clusterTab <- cbind(clusterTab, width = as.numeric(clusterTab[, "end"]) - as.numeric(clusterTab[, "start"]))
	write.table(clusterTab, gsub("delta", "hits", testFile), sep = "\t", row.names=F, quote = FALSE)
}






if(FALSE)
{

# set the working directory
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/results/")
bedCD3 <- read.delim("G3_all_Alignment_delta.bed")

#bedCD3 <- read.table("UT_delta.txt", header=T, stringsAsFactors=F)

CUTOFF <- 222
clusterTab <- c() # for each chomosome calculate the difference between one start site and the previous one
colTOkeep <- c(1,2,3,4,5,6)
colClust <- length(colTOkeep) +1
for(i in unique(bedCD3[,1])){ # start for loop for each chomosome
	tempbed <- subset(bedCD3,  chromosome==i) # take the data for one choromose 
	for(n in 1:nrow(tempbed)){   # start for loop for each row
		if( n == 1){ 
			clusterTab <- rbind(clusterTab, unlist(c(tempbed[1,colTOkeep], 1)))  # report read in the final table
			}
		if(n >= 2 ){
			if( tempbed$delta[n] <= CUTOFF ){
				clusterTab[nrow(clusterTab), 3] <- tempbed$end[n]
				clusterTab[nrow(clusterTab), colClust-1 ] <- as.numeric(clusterTab[nrow(clusterTab), colClust-1]) +  tempbed$collapse[n] 
				clusterTab[nrow(clusterTab), colClust] <- as.numeric(clusterTab[nrow(clusterTab), colClust]) +1  # add collapsed sequence
			}
			if( tempbed$delta[n] >= CUTOFF ){ 
				clusterTab <- rbind(clusterTab, unlist(c(tempbed[n,colTOkeep], 1)))
			}
		}
	}
}


colnames(clusterTab)[colClust] <-"collapseCluster"
write.table(clusterTab, "G3_clusters.txt",sep = "\t", row.names=F)
#write.table(clusterTab, "UT_clusters.txt",sep = "\t", row.names=F)

}


