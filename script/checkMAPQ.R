# bash script (bai file is needed)
#samtools view G3-WT-d1_S1_L001_AlignmentSort.bam chr22:20000000-21000000 | \
#    awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
    
    
getAvgMAPQ <- function(functionstring="samtools view", bamFile, chr, start, end)
{
	options(scipen =99)
	# check if bai file exist, if not create it
	if(!file.exists(paste0(bamFile, ".bai"))){
		command=paste("samtools index", bamFile, sep = " ")
		try(system(command))
		}
	
	# Get average MAPQ
	region <- paste0(chr, ":", start, "-", end)
	command=paste(functionstring, bamFile, region, "| awk '{sum+=$5} END { print sum/NR}' ", sep = " ")
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

# DO NOT RUN
#nbReadFastq(bamFile = "G3-WT-d1_S1_L001_AlignmentSort.bam", chr = "chr22", start = 20000000, end = 21000000)



addMAPQ <- function(inputF, bamFile)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	
	mpq <- lapply(1:nrow(readMat), function(i){
		getAvgMAPQ(bamFile = bamFile, chr = readMat$chromosome[i], start = readMat$start[i], end = readMat$end[i])
		})
	
	readMat$avg.MAPQ <- as.numeric(unlist(mpq))
	write.xlsx(readMat, inputF, row.names = FALSE)
}
