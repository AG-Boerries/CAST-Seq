




bamCoverage <- function(functionstring="~/Programs/bedtools2/bin/bedtools coverage", bam, bed, outFile, opt.string="")
{
  options(scipen =99)
  
  command=paste(functionstring, "-a", bed, "-b", bam, opt.string, " | gzip >", outFile)
  cat(command,"\n")
  try(system(command))
}



ONreadout <- function(bamFile, otsFile, gRNA.orientation, window.size = 5000, sampleName, outDir){
  
  # READ INPUT FILES
  ots <- read.delim(otsFile, header = FALSE)
  ots.title <- paste0("ON: ", ots$V1, ": ", ots$V2, " - ", ots$V3, "; ", ots$V6, " strand")
  
  ots.site <- ots$V2 + floor((ots$V3 - ots$V2)/2)
  
  ots$V6 <- "+"# always on positive strand
  
  dir.create(dir1 <- file.path(tempdir(), "cleavageDir"))
  
  # extend ots
  ots$V2 <- ots$V2 - window.size
  ots$V3 <- ots$V3 + window.size
  ots.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  write.table(ots, ots.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # calculate coverage per base
  cov.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.pos.tmp, opt.string="-bed -d -s")
  cov.pos <- read.delim(cov.pos.tmp, header = FALSE)
  
  cov.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.neg.tmp, opt.string="-bed -d -S")
  cov.neg <- read.delim(cov.neg.tmp, header = FALSE)
  
  # calculate total coverage
  cov.tot.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.pos.tmp, opt.string="-bed -s")
  cov.tot.pos <- read.delim(cov.tot.pos.tmp, header = FALSE)
  cov.tot.pos <- cov.tot.pos$V7
  
  cov.tot.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.neg.tmp, opt.string="-bed -S")
  cov.tot.neg <- read.delim(cov.tot.neg.tmp, header = FALSE)
  cov.tot.neg <- cov.tot.neg$V7
  
  cov.tot <- cov.tot.pos + cov.tot.neg
  cov.df <- data.frame(TOTAL = cov.tot,
                       POS = cov.tot.pos,
                       NEG = cov.tot.neg,
                       POS.PC = round(cov.tot.pos / cov.tot * 100, digits = 2),
                       NEG.PC = round(cov.tot.neg / cov.tot * 100, digits = 2)
  )
  rownames(cov.df) <- "Read coverage"
  
  if(gRNA.orientation == "forward"){
    cov.df$DEL <- cov.df$POS
    cov.df$INV <- cov.df$NEG
    cov.df$DEL.PC <- cov.df$POS.PC
    cov.df$INV.PC <- cov.df$NEG.PC
  }else if(gRNA.orientation == "reverse"){
    cov.df$DEL <- cov.df$NEG
    cov.df$INV <- cov.df$POS
    cov.df$DEL.PC <- cov.df$NEG.PC
    cov.df$INV.PC <- cov.df$POS.PC
  }else{
    print("Wrong gRNA orientation. Must be forward or reverse")
  }
  write.xlsx(cov.df, file = file.path(outDir, paste0(sampleName, "_", gsub("000$", "kb", window.size), ".xlsx")),
             row.names = TRUE, overwrite = TRUE)
  
  
  # plot coverage
  ggmat <- data.frame(CHR = c(cov.pos$V1, cov.neg$V1),
                      POSITION = c(cov.pos$V2 + cov.pos$V7, cov.neg$V2 + cov.neg$V7),
                      STRAND = c(rep("POS", nrow(cov.pos)), rep("NEG", nrow(cov.neg))),
                      COV = c(cov.pos$V8, cov.neg$V8)
  )
  ggmat$DIST <- ggmat$POSITION - ots.site
  
  ggmat$EVENT <- NA
  if(gRNA.orientation == "forward"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "DEL" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "INV" 
  }else if(gRNA.orientation == "reverse"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "INV" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "DEL" 
  }else{
    print("Wrong gRNA orientation. Must be forward or reverse")
  }
  ggmat$EVENT <- factor(ggmat$EVENT, levels = c("DEL", "INV"))
  
  # COVERAGE
  p <- ggplot(ggmat, aes(x= POSITION))
  p <- p + geom_line(aes(y=COV, color = STRAND), size = 1, alpha = 0.75)
  p <- p + geom_ribbon(aes(ymin = 0, ymax = COV, fill = STRAND), alpha = 0.2)
  p <- p + scale_color_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                       NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + scale_fill_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                      NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) 
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + annotation_logticks() 
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + xlab("") + ylab("read coverage")
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  ggsave(plot = p, filename = file.path(outDir, paste0(sampleName, "_coverage_ON_", gsub("000$", "kb", window.size), ".pdf")), 
         width = 8, height = 5)
  
  # DEL / INV
  p <- ggplot(ggmat, aes(x= POSITION))
  p <- p + geom_line(aes(y=COV, color = EVENT), size = 1, alpha = 0.75)
  p <- p + geom_ribbon(aes(ymin = 0, ymax = COV, fill = EVENT), alpha = 0.2)
  p <- p + scale_color_manual(values=c(DEL = "orange2",
                                       INV = "orchid2"))
  p <- p + scale_fill_manual(values=c(DEL = "orange2",
                                      INV = "orchid2"))
  p <- p + scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) 
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + annotation_logticks() 
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + xlab("") + ylab("read coverage")
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  ggsave(plot = p, filename = file.path(outDir, paste0(sampleName, "_coverage_ON_", gsub("000$", "kb", window.size), "_event.pdf")),
         width = 8, height = 5)
  
  # RATIO
  cov.tot <- cov.pos$V8 + cov.neg$V8
  ggmat <- data.frame(CHR = c(cov.pos$V1, cov.neg$V1),
                      POSITION = c(cov.pos$V2 + cov.pos$V7, cov.neg$V2 + cov.neg$V7),
                      STRAND = c(rep("POS", nrow(cov.pos)), rep("NEG", nrow(cov.neg))),
                      PC = c(cov.pos$V8 / cov.tot * 100, cov.neg$V8 / cov.tot * 100))
  
  ggmat$EVENT <- NA
  if(gRNA.orientation == "forward"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "DEL" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "INV" 
  }else if(gRNA.orientation == "reverse"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "INV" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "DEL" 
  }else{
    print("Wrong gRNA orientation. Must be forward or reverse")
  }
  ggmat$EVENT <- factor(ggmat$EVENT, levels = c("DEL", "INV"))
  
  p <- ggplot(data=ggmat, aes(x=POSITION, y=PC, fill=STRAND)) +
    geom_bar(stat="identity")+
    #geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
    #          color="white", size=3.5)+
    scale_fill_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                               NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  p <- p + xlab("") + ylab("strand percentage")
  ggsave(plot = p, file.path(outDir, filename = paste0(sampleName, "_percentage_ON_", gsub("000$", "kb", window.size), ".pdf")),
         width = 8, height = 5)
  
  
  p <- ggplot(data=ggmat, aes(x=POSITION, y=PC, fill=EVENT)) +
    geom_bar(stat="identity")+
    #geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
    #          color="white", size=3.5)+
    scale_fill_manual(values=c(DEL = "orange2",
                               INV = "orchid2"))
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  p <- p + xlab("") + ylab("strand percentage")
  ggsave(plot = p, filename = file.path(outDir, paste0(sampleName, "_percentage_ON_", gsub("000$", "kb", window.size), "_event.pdf")),
         width = 8, height = 5)
  
  # Delete tmp dir
  unlink(dir1, recursive = T)
}





# DO NOT RUN
if(FALSE){
  library(ggplot2)
  library(openxlsx)
  
  # PARAMTERS
  window.size <- 5000
  outDir <- file.path("~/Research/CASTSeq/ONreadout")
  dir.create(outDir)
  
  bamFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/CASTseq_manuscript/G3_WT_D1/results/fastq_aln/G3-WT-d1_S1_L001_AlignmentSort.bam")
  otsFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/CASTseq_manuscript/G3_WT_D1/data/ots.bed")
  sampleName <- "G3_WT_D1"
  gRNA.orientation <- "forward"
  ONreadout(bamFile, otsFile, gRNA.orientation, window.size = 5000, sampleName, outDir)
  
  
  bamFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/CASTseq_manuscript/G3Rev_WT_D1/results/fastq_aln/G3Rev_UT_D1_L001_AlignmentSort.bam")
  otsFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/CASTseq_manuscript/G3Rev_WT_D1/data/ots.bed")
  sampleName <- "G3Rev_WT_D1"
  gRNA.orientation <- "reverse"
  ONreadout(bamFile, otsFile, gRNA.orientation, window.size = 5000, sampleName, outDir)
  
  
  bamFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-1-g1g3/results/fastq_aln/COL17A1-1-g1g3_AlignmentSort.bam")
  otsFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-1-g1g3/data/ots.bed")
  sampleName <- "COL17A1-1-g1g3"
  
  bamFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/results/fastq_aln/EMD101-sample1-1_AlignmentSort.bam")
  otsFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/data/ots.bed")
  sampleName <- "EMD101-sample1-1"
  

  

  
  
  # READ INPUT FILES
  ots <- read.delim(otsFile, header = FALSE)
  ots.title <- paste0("ON: ", ots$V1, ": ", ots$V2, " - ", ots$V3, "; ", ots$V6, " strand")
  
  ots.site <- ots$V2 + floor((ots$V3 - ots$V2)/2)
  
  ots$V6 <- "+"# always on positive strand
  
  dir.create(dir1 <- file.path(tempdir(), "cleavageDir"))
  
  # extend ots
  ots$V2 <- ots$V2 - window.size
  ots$V3 <- ots$V3 + window.size
  ots.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  write.table(ots, ots.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ots5 <- ots
  ots5$V3 <- ots.site
  ots5.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  write.table(ots5, ots5.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ots3 <- ots
  ots3$V2 <- ots.site
  ots3.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  write.table(ots3, ots3.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # calculate coverage per base
  cov.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.pos.tmp, opt.string="-bed -d -s")
  cov.pos <- read.delim(cov.pos.tmp, header = FALSE)

  cov.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.neg.tmp, opt.string="-bed -d -S")
  cov.neg <- read.delim(cov.neg.tmp, header = FALSE)
  
  # calculate total coverage
  cov.tot.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.pos.tmp, opt.string="-bed -s")
  cov.tot.pos <- read.delim(cov.tot.pos.tmp, header = FALSE)
  cov.tot.pos <- cov.tot.pos$V7
  
  cov.tot.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.neg.tmp, opt.string="-bed -S")
  cov.tot.neg <- read.delim(cov.tot.neg.tmp, header = FALSE)
  cov.tot.neg <- cov.tot.neg$V7
  
  cov.tot <- cov.tot.pos + cov.tot.neg
  cov.df <- data.frame(TOTAL = cov.tot,
                       POS = cov.tot.pos,
                       NEG = cov.tot.neg,
                       POS.PC = round(cov.tot.pos / cov.tot * 100, digits = 2),
                       NEG.PC = round(cov.tot.neg / cov.tot * 100, digits = 2)
                       )
  rownames(cov.df) <- "Read coverage"
  
  if(gRNA.orientation == "forward"){
    cov.df$DEL <- cov.df$POS
    cov.df$INV <- cov.df$NEG
    cov.df$DEL.PC <- cov.df$POS.PC
    cov.df$INV.PC <- cov.df$NEG.PC
  }else if(gRNA.orientation == "reverse"){
    cov.df$DEL <- cov.df$NEG
    cov.df$INV <- cov.df$POS
    cov.df$DEL.PC <- cov.df$NEG.PC
    cov.df$INV.PC <- cov.df$POS.PC
  }else{
    print("Wrong gRNA orientation. Must be forward or reverse")
  }
  write.xlsx(cov.df, file = paste0("~/tmp/", sampleName, "_", gsub("000$", "kb", window.size), ".xlsx"),
             row.names = TRUE, overwrite = TRUE)
  
  
  # plot coverage
  ggmat <- data.frame(CHR = c(cov.pos$V1, cov.neg$V1),
                      POSITION = c(cov.pos$V2 + cov.pos$V7, cov.neg$V2 + cov.neg$V7),
                      STRAND = c(rep("POS", nrow(cov.pos)), rep("NEG", nrow(cov.neg))),
                      COV = c(cov.pos$V8, cov.neg$V8)
                      )
  ggmat$DIST <- ggmat$POSITION - ots.site
  
  ggmat$EVENT <- NA
  if(gRNA.orientation == "forward"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "DEL" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "INV" 
  }else if(gRNA.orientation == "reverse"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "INV" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "DEL" 
  }else{
    print("Wrong gRNA orientation. Must be forward or reverse")
  }
  ggmat$EVENT <- factor(ggmat$EVENT, levels = c("DEL", "INV"))
  
  # COVERAGE
  p <- ggplot(ggmat, aes(x= POSITION))
  p <- p + geom_line(aes(y=COV, color = STRAND), size = 1, alpha = 0.75)
  p <- p + geom_ribbon(aes(ymin = 0, ymax = COV, fill = STRAND), alpha = 0.2)
  p <- p + scale_color_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                       NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + scale_fill_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                       NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) 
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + annotation_logticks() 
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + xlab("") + ylab("read coverage")
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  ggsave(plot = p, filename = paste0("~/tmp/", sampleName, "_coverage_ON_", gsub("000$", "kb", window.size), ".pdf"), width = 8, height = 5)

  # DEL / INV
  p <- ggplot(ggmat, aes(x= POSITION))
  p <- p + geom_line(aes(y=COV, color = EVENT), size = 1, alpha = 0.75)
  p <- p + geom_ribbon(aes(ymin = 0, ymax = COV, fill = EVENT), alpha = 0.2)
  p <- p + scale_color_manual(values=c(DEL = "orange2",
                                       INV = "orchid2"))
  p <- p + scale_fill_manual(values=c(DEL = "orange2",
                                       INV = "orchid2"))
  p <- p + scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) 
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + annotation_logticks() 
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + xlab("") + ylab("read coverage")
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  ggsave(plot = p, filename = paste0("~/tmp/", sampleName, "_coverage_ON_", gsub("000$", "kb", window.size), "_event.pdf"), width = 8, height = 5)
  
  # RATIO
  cov.tot <- cov.pos$V8 + cov.neg$V8
  ggmat <- data.frame(CHR = c(cov.pos$V1, cov.neg$V1),
                      POSITION = c(cov.pos$V2 + cov.pos$V7, cov.neg$V2 + cov.neg$V7),
                      STRAND = c(rep("POS", nrow(cov.pos)), rep("NEG", nrow(cov.neg))),
                      PC = c(cov.pos$V8 / cov.tot * 100, cov.neg$V8 / cov.tot * 100))
  
  ggmat$EVENT <- NA
  if(gRNA.orientation == "forward"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "DEL" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "INV" 
  }else if(gRNA.orientation == "reverse"){
    ggmat$EVENT[ggmat$STRAND == "POS"] <- "INV" 
    ggmat$EVENT[ggmat$STRAND == "NEG"] <- "DEL" 
  }else{
    print("Wrong gRNA orientation. Must be forward or reverse")
  }
  ggmat$EVENT <- factor(ggmat$EVENT, levels = c("DEL", "INV"))
  
  p <- ggplot(data=ggmat, aes(x=POSITION, y=PC, fill=STRAND)) +
    geom_bar(stat="identity")+
    #geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
    #          color="white", size=3.5)+
    scale_fill_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                               NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  p <- p + xlab("") + ylab("strand percentage")
  ggsave(plot = p, filename = paste0("~/tmp/", sampleName, "_percentage_ON_", gsub("000$", "kb", window.size), ".pdf"), width = 8, height = 5)
  
  
  p <- ggplot(data=ggmat, aes(x=POSITION, y=PC, fill=EVENT)) +
    geom_bar(stat="identity")+
    #geom_text(aes(y=label_ypos, label=len), vjust=1.6, 
    #          color="white", size=3.5)+
    scale_fill_manual(values=c(DEL = "orange2",
                               INV = "orchid2"))
  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  p <- p + xlab("") + ylab("strand percentage")
  ggsave(plot = p, filename = paste0("~/tmp/", sampleName, "_percentage_ON_", gsub("000$", "kb", window.size), "_event.pdf"), width = 8, height = 5)
  
  
  
  
  
  
  
  
  
  
  # IDENTIFY BLOCKS
  cov.pos.block <- cov.pos[cov.pos$V8 != 0, ]

  delta <- as.numeric(cov.pos.block[1,7])
  if(nrow(cov.pos.block) > 1) delta <- c(delta, as.numeric(cov.pos.block[2:nrow(cov.pos.block),7]) - as.numeric(cov.pos.block[1:(nrow(cov.pos.block)-1),7]))
  cov.pos.block$DELTA <- delta
  
  tempbed <- cov.pos.block
  tempbed$V2 <- tempbed$V2 + tempbed$V7-1
  tempbed$V3 <- tempbed$V2 + tempbed$V7
  clusterTab <- c()

  CUTOFF <- 1
  for(n in 1:nrow(tempbed)){   # start for loop for each row
    if( n == 1){ 
      clusterTab <- rbind(clusterTab, unlist(tempbed[1,]))  # report read in the final table
    }
    if(n >= 2 ){
      if( tempbed$DELTA[n] <= CUTOFF ){
        clusterTab[nrow(clusterTab), 3] <- tempbed[n, 3]
        clusterTab[nrow(clusterTab), 7] <- tempbed[n, 7]
        clusterTab[nrow(clusterTab), 8] <- as.numeric(clusterTab[nrow(clusterTab), 8]) +  tempbed$V8[n] 
      }
      if( tempbed$DELTA[n] > CUTOFF ){ 
        clusterTab <- rbind(clusterTab, unlist(tempbed[n,]))
      }
    }  
  }
  clusterTab <- data.frame(clusterTab)
  clusterTab$SIZE <- as.numeric(clusterTab$V3) - as.numeric(clusterTab$V2)
  
  
  unlink(dir1, recursive = T)
  
}









