
getReadLength <- function(fastqGZ, outFile){
  
  #fastq <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".fastq")
  fastq <- gsub(".gz$", "", fastqGZ)
  if(!file.exists(fastq)){
    gzipCommand <- paste0("gzip -k -d ", fastqGZ)
    system(gzipCommand)
  }
  
  mycommand <- paste0("awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ", fastq,  " > ", outFile)
  system(mycommand)
  
  unlink(fastq)
}


plotReadLength <- function(fastqGZ, outName, bin = TRUE){

  getReadLength(fastqGZ, paste0(outName, ".txt"))
  
  # Read count matrix and create a scatterplot
  countMat <- read.delim(paste0(outName, ".txt"), header = FALSE, sep = " ")
  colnames(countMat) <- c("length", "count")
  countMat <- countMat[order(countMat$length), ]
  
  countMat$pass <- ifelse(countMat$length<30, "no", "yes")
  
  if(!bin){
    p <- ggplot(countMat, aes(x=length, y=count)) + 
      geom_point(aes(colour =pass), alpha = 0.5) + 
      scale_color_manual(values=c(yes = "green4", no = "red4"))
    p <- p + theme(legend.position = "none")
    p <- p + theme_bw(base_size = 14) +
      xlab("read length (bp)") + ylab("number of reads")
    
    ggsave(plot = p, filename = paste0(outName, "_read_length_distribution.pdf"),
           width = 6, height = 4)
  }else{
    # BIN
    bins <- c(0,5,10,15,20,25,30,50,100,500)
    
    countMat.bin <- lapply(seq(1,length(bins)-1), function(i){
      idx.min <- bins[i]
      idx.max <- bins[i+1]
      
      idx.str <- paste0(idx.min, "_", idx.max)
      idx.count <- sum(countMat$count[countMat$length >= idx.min & countMat$length < idx.max])
      return(data.frame(bin = idx.str, count = idx.count))
    })
    countMat.bin <- do.call(rbind, countMat.bin)
    countMat.bin$bin <- factor(countMat.bin$bin, levels = unique(countMat.bin$bin))
    
    countMat.bin$pass <- c(rep("no", 6), rep("yes", 3))
    
    p<-ggplot(data=countMat.bin, aes(x=bin, y=count)) +
      geom_bar(aes(fill = pass), stat="identity") +
      scale_fill_manual(values=c(yes = "green4", no = "red4"))
    p <- p + theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      xlab("read length (bp)") + ylab("number of reads")
    p <- p + theme(legend.position = "none")
    
    ggsave(plot = p, filename = paste0(outName, "_read_length_distribution_BIN.pdf"),
           width = 6, height = 4)
  }
}


if(FALSE){
  library(ggplot2)
  
  # MAIN
  # test
  setwd(file.path("~/Tmp"))
  fastqGZ <- "BCL11A-BE-1_trim3_soft.fastq.gz"
  
  outName <- "BCL11A-BE-1_trim3_soft"
  
  
  # CYBB
  fastqGZ <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-887623654/CYBB/CYBB_crispr/Cpf1-CYBB/T-Cpf1CYBB1/fastq_aln/T-Cpf1CYBB1_trim3_soft.fastq.gz")
  dir.create("~/Research/CASTSeq/Test/read_length")
  outName <- "~/Research/CASTSeq/Test/read_length/T-Cpf1CYBB1_trim3_soft"
  
  plotReadLength(fastqGZ, outName)
  
  # COL17A1
  fastqGZ <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-851042618/COL17A1_longer_pos/COL17A1_g1/COL17A1-D10Ag1/COL17A1-D10Ag1-1/fastq_aln/COL17A1-D10Ag1-1_trim3_soft.fastq.gz")
  outName <- "~/Research/CASTSeq/Test/read_length/COL17A1-D10Ag1-1_trim3_soft"
  
  plotReadLength(fastqGZ, outName)
  
  # COL17A1 g3
  fastqGZ <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-851042618/COL17A1/COL17A1_g3_crispr/COL17A1-D10Ag3/COL17A1-D10Ag2-1/fastq_aln/COL17A1-D10Ag2-1_trim3_soft.fastq.gz")
  outName <- "~/Research/CASTSeq/Test/read_length/COL17A1-D10Ag2-1_trim3_soft"
  
  plotReadLength(fastqGZ, outName)
  
}

