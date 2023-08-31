
############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

#nbReadFastq <- function(inFile)
#{
#  options(scipen =99)
#  command=paste0("echo $(cat ", inFile, "|wc -l)/4|bc")
#  cat(command,"\n")
#  try(system(command, intern = TRUE))
#}

#nbReadFastqgz <- function(inFile)
#{
#  options(scipen =99)
#  command=paste0("zcat < ", inFile, " | echo $((`wc -l`/4))")
#  cat(command,"\n")
#  try(system(command, intern = TRUE))
#}


nbReadBAM <- function(inFile){
  options(scipen =99)
  command=paste0("samtools view -c -F 260 ", inFile)
  cat(command,"\n")
  try(system(command, intern = TRUE))
}

nbReadBED <- function(inFile){
  options(scipen =99)
  command=paste0("echo $(wc -l ", inFile, ")")
  cat(command,"\n")
  try(system(command, intern = TRUE))
}


countReads <- function(fastqFolder, sampleFolder, testSample, testLabel, ext1, ext2, outputFile){
  
  countMatList <- lapply(1:length(testSample), function(i){
    
    sample.current <- testSample[i]
    label.current <- testLabel[i]
    fastq.current <- fastqFolder[i]
    
    # R1
    f1 <- file.path(fastq.current, paste0(sample.current, ext1))
    r1 <- as.numeric(nbReadFastqgz(f1))
    
    r1.pc <- 100
    
    # R2
    f2 <- file.path(fastq.current, paste0(sample.current, ext2))
    r2 <- as.numeric(nbReadFastqgz(f2))
    
    r2.pc <- 100
    
    # ASSEMBLED
    assembled <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current,  "fastq_aln", paste0(sample.current, "_assembled.fastq.gz"))))
    if(r1 == 0){assembled.pc <- 0}else assembled.pc <- round(assembled / r1 * 100, digits = 2)
    
    # MERGED
    merged <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_merged.fastq.gz"))))
    if(r1 == 0){merged.pc <- 0}else merged.pc <- round(merged / r1 * 100, digits = 2)
    
    # POS
    pos <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_pos.fastq.gz"))))
    if(merged == 0){pos.pc <- 0}else pos.pc <- round(pos / merged * 100, digits = 2)
    
    # NEG (mispriming)
    neg <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_Filt2.fastq.gz"))))
    if(pos == 0){neg.pc <- 0}else neg.pc <- round(neg / pos * 100, digits = 2)
    
    # TRIM1
    trim1 <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_trim1.fastq.gz"))))
    if(neg == 0){trim1.pc <- 0}else trim1.pc <- round(trim1 / neg * 100, digits = 2)
    
    # TRIM2
    trim2 <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_trim2.fastq.gz"))))
    if(trim1 == 0){trim2.pc <- 0}else trim2.pc <- round(trim2 / trim1 * 100, digits = 2)
    
    # TRIM3
    trim3 <- as.numeric(nbReadFastqgz(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_trim3.fastq.gz"))))
    if(trim2 == 0){trim3.pc <- 0}else trim3.pc <- round(trim3 / trim2 * 100, digits = 2)
    
    # BAM
    bam <- as.numeric(nbReadBAM(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_AlignmentSort.bam"))))
    if(trim3 == 0){bam.pc <- 0}else bam.pc <- round(bam / trim3 * 100, digits = 2)
    
    # BED
    bed <- nbReadBED(file.path(sampleFolder, label.current, "fastq_aln", paste0(sample.current, "_Alignment.bed")))
    bed <- as.numeric(sapply(strsplit(bed, split = " "), function(i) i[1]))
    if(bam == 0){bed.pc <- 0}else bed.pc <- round(bed / bam * 100, digits = 2)
    
    # SAVE
    toxlsx <- data.frame(count = c(r1, r2, assembled, merged, pos, neg, trim1, trim2, trim3, bam, bed),
                         percentage = c(r1.pc, r2.pc, assembled.pc, merged.pc, pos.pc, neg.pc, trim1.pc, trim2.pc, trim3.pc, bam.pc, bed.pc)
    )
    colnames(toxlsx) <- paste0(label.current, ".", colnames(toxlsx))
    return(toxlsx)
  })
  countMat <- do.call(cbind, countMatList)
  countMat <- data.frame(File = c("R1", "R2", "assembled", "merged", "pos", "mispriming", "trim1", "trim2", "trim3", "bam", "bed"),
                         countMat)
  write.xlsx(countMat, outputFile)

}

countReadsPlot <- function(inputFile, outputFile){
  qcMat <- read.xlsx(inputFile, sheet = 1)
  rownames(qcMat) <- qcMat$File
  qcMat <- qcMat[, grepl ("percentage", colnames(qcMat))]
  colnames(qcMat) <- gsub(".percentage", "", colnames(qcMat))
  
  pheatmap(qcMat, color = magma(10), border_color = "grey60", fontsize = 12,
           cellwidth = 14, cellheight = 14, cluster_rows = FALSE, cluster_cols = FALSE,
           breaks = seq(0,100, by=10),
           filename = outputFile)
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# DO NOT RUN
if(FALSE){
  sampleFolder <- file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/")
  testSample <- "EMD101-sample1-1"
  utSample <- "EMD101-sample3-1"
  
  countReads(sampleFolder = sampleFolder,
             testSample = testSample,
             utSample = utSample)
  
}