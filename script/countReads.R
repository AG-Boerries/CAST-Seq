
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


countReads <- function(fastqFolder, sampleFolder, testSample, utSample){
  # R1
  r1 <- as.numeric(nbReadFastqgz(file.path(fastqFolder, paste0(testSample, "_R1_001.fastq.gz"))))
  r1.UT <- as.numeric(nbReadFastqgz(file.path(fastqFolder, paste0(utSample, "_R1_001.fastq.gz"))))
  
  r1.pc <- 100
  r1.UT.pc <- 100
  
  # R2
  r2 <- as.numeric(nbReadFastqgz(file.path(fastqFolder, paste0(testSample, "_R2_001.fastq.gz"))))
  r2.UT <- as.numeric(nbReadFastqgz(file.path(fastqFolder, paste0(utSample, "_R2_001.fastq.gz"))))
  
  r2.pc <- 100
  r2.UT.pc <- 100
  
  # ASSEMBLED
  assembled <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_assembled.fastq.gz"))))
  assembled.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_assembled.fastq.gz"))))
  
  if(r1 == 0){assembled.pc <- 0}else assembled.pc <- round(assembled / r1 * 100, digits = 2)
  if(r1.UT == 0){assembled.UT.pc <- 0}else assembled.UT.pc <- round(assembled.UT / r1.UT * 100, digits = 2)
  
  # MERGED
  merged <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_merged.fastq.gz"))))
  merged.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_merged.fastq.gz"))))
  
  if(r1 == 0){merged.pc <- 0}else merged.pc <- round(merged / r1 * 100, digits = 2)
  if(r1.UT == 0){merged.UT.pc <- 0}else merged.UT.pc <- round(merged.UT / r1.UT * 100, digits = 2)
  
  # POS
  pos <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_pos.fastq.gz"))))
  pos.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_pos.fastq.gz"))))
  
  if(merged == 0){pos.pc <- 0}else pos.pc <- round(pos / merged * 100, digits = 2)
  if(merged.UT == 0){pos.UT.pc <- 0}else pos.UT.pc <- round(pos.UT / merged.UT * 100, digits = 2)
  
  # NEG (mispriming)
  neg <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_Filt2.fastq.gz"))))
  neg.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_Filt2.fastq.gz"))))
  
  if(pos == 0){neg.pc <- 0}else neg.pc <- round(neg / pos * 100, digits = 2)
  if(pos.UT == 0){neg.UT.pc <- 0}else neg.UT.pc <- round(neg.UT / pos.UT * 100, digits = 2)
  
  # TRIM1
  trim1 <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_trim1.fastq.gz"))))
  trim1.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_trim1.fastq.gz"))))
  
  if(neg == 0){trim1.pc <- 0}else trim1.pc <- round(trim1 / neg * 100, digits = 2)
  if(neg.UT == 0){trim1.UT.pc <- 0}else trim1.UT.pc <- round(trim1.UT / neg.UT * 100, digits = 2)
  
  # TRIM2
  trim2 <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_trim2.fastq.gz"))))
  trim2.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_trim2.fastq.gz"))))
  
  if(trim1 == 0){trim2.pc <- 0}else trim2.pc <- round(trim2 / trim1 * 100, digits = 2)
  if(trim1.UT == 0){trim2.UT.pc <- 0}else trim2.UT.pc <- round(trim2.UT / trim1.UT * 100, digits = 2)
  
  # TRIM3
  trim3 <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_trim3.fastq.gz"))))
  trim3.UT <- as.numeric(nbReadFastqgz(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_trim3.fastq.gz"))))
  
  if(trim2 == 0){trim3.pc <- 0}else trim3.pc <- round(trim3 / trim2 * 100, digits = 2)
  if(trim2.UT == 0){trim3.UT.pc <- 0}else trim3.UT.pc <- round(trim3.UT / trim2.UT * 100, digits = 2)
  
  # BAM
  bam <- as.numeric(nbReadBAM(file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_AlignmentSort.bam"))))
  bam.UT <- as.numeric(nbReadBAM(file.path(sampleFolder, "results/fastq_aln", paste0(utSample, "_AlignmentSort.bam"))))
  
  if(trim3 == 0){bam.pc <- 0}else bam.pc <- round(bam / trim3 * 100, digits = 2)
  if(trim3.UT == 0){bam.UT.pc <- 0}else bam.UT.pc <- round(bam.UT / trim3.UT * 100, digits = 2)
  
  # BED
  bed <- nbReadBED(file.path(sampleFolder, "results/guide_aln", paste0(testSample, "_Alignment.bed")))
  bed <- as.numeric(sapply(strsplit(bed, split = " "), function(i) i[1]))
  bed.UT <- nbReadBED(file.path(sampleFolder, "results/guide_aln", paste0(utSample, "_Alignment.bed")))
  bed.UT <- as.numeric(sapply(strsplit(bed.UT, split = " "), function(i) i[1]))
  
  if(bam == 0){bed.pc <- 0}else bed.pc <- round(bed / bam * 100, digits = 2)
  if(bam.UT == 0){bed.UT.pc <- 0}else bed.UT.pc <- round(bed.UT / bam.UT * 100, digits = 2)
  
  # SAVE
  toxlsx <- data.frame(File = c("R1", "R2", "assembled", "merged", "pos", "mispriming", "trim1", "trim2", "trim3", "bam", "bed"),
                       Treated.count = c(r1, r2, assembled, merged, pos, neg, trim1, trim2, trim3, bam, bed),
                       Treated.percentage = c(r1.pc, r2.pc, assembled.pc, merged.pc, pos.pc, neg.pc, trim1.pc, trim2.pc, trim3.pc, bam.pc, bed.pc),
                       UN.Treated.count = c(r1.UT, r2.UT, assembled.UT, merged.UT, pos.UT, neg.UT, trim1.UT, trim2.UT, trim3.UT, bam.UT, bed.UT),
                       UN.Treated.percentage = c(r1.UT.pc, r2.UT.pc, assembled.UT.pc, merged.UT.pc, pos.UT.pc, neg.UT.pc, trim1.UT.pc, trim2.UT.pc, trim3.UT.pc, bam.UT.pc, bed.UT.pc)
                       )
  
  write.xlsx(toxlsx, file.path(sampleFolder, "results/fastq_aln", paste0(testSample, "_QC.xlsx")), overwrite = TRUE)
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