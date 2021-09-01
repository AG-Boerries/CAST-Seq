
############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################




readLength <- function(fastq, out = NULL)
{
	# zcat G3-WT-d1_S1_L001_trim2.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' > G3-WT-d1_S1_L001_trim2.readLength.txt
	options(scipen =99)
	if(is.null(out)){
		outFile <- gsub(".fastq.gz", ".readLength.txt", fastq)}else{outFile <- out}
	command=paste("gunzip -c",fastq,"| awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' >", outFile,sep=" ")
  	cat(command,"\n")
  	try(system(command))
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE){
  library(ggplot2)
  library(parallel)
  
  # Get all fastq.gz
  #setwd("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results/fastq_aln")
  #setwd("/home/gandri/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results/fastq_aln")
  setwd("/home/gandri/offTargets/Giando/pipelineGit/samples/")
  #fastqFiles <- list.files(pattern = ".fastq.gz$", recursive = TRUE)
  #fastqFiles <- fastqFiles[!grepl("Headtoheadtrimmed", fastqFiles)]
  #mclapply(fastqFiles, readLength, mc.cores = 12)
  
  # Get read distribution
  lengthFiles <- list.files(pattern = ".readLength.txt$", recursive = TRUE)
  names(lengthFiles) <- gsub(".txt", "", basename(lengthFiles))
  lengthList <- lapply(lengthFiles, read.delim, header = FALSE, colClasses = c("integer", "integer"), sep = " ")
  
  # Plot read length density distribution
  #setwd("~/cluster/master/offTargets/Giando/readLength")
  setwd("/home/gandri/offTargets/Giando/readLength")
  lapply(1:length(lengthList), function(i){
    p <- ggplot(data=lengthList[[i]], aes(x=V1, y=V2))
    p <- p + geom_bar(stat="identity")
    p <- p + theme_bw()
    p <- p + xlab("read length (bp)") + ylab("number of reads")
    ggsave(plot = p, filename = paste0(names(lengthList)[i], ".pdf"), width = 7, height = 7)
  })
  
  
  
}



