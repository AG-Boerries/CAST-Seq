




bamCoverage <- function(functionstring="bedtools coverage", bam, bed, outFile, opt.string="")
{
  options(scipen =99)
  
  command=paste(functionstring, "-a", bed, "-b", bam, opt.string, " | gzip >", outFile)
  cat(command,"\n")
  try(system(command))
}



chunkBed <- function(min, max, size){
  limits <- ceiling(seq(min, max, length.out = size+1))
  bins <- lapply(2:length(limits), function(k){
    return(c(limits[k-1], limits[k]))
  })
  return(do.call(rbind, bins))
}

ONreadout.old <- function(bamFile, otsFile, gRNA.orientation, window.size = 5000, sampleName, outDir){
  
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
  
  # calculate coverage per bins
  bins <- chunkBed(ots$V2, ots$V3, 100)
  bins <- data.frame(ots$V1,
                     bins,
                     paste0("bin_", seq(1:nrow(bins))),
                     1000,
                     "+")
  colnames(bins) <- colnames(ots)
  bins.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  write.table(bins, bins.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cov.bin.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = bins.tmp, outFile = cov.bin.pos.tmp, opt.string="-bed -s")
  cov.bin.pos <- read.delim(cov.bin.pos.tmp, header = FALSE)
  
  cov.bin.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile, bed = bins.tmp, outFile = cov.bin.neg.tmp, opt.string="-bed -S")
  cov.bin.neg <- read.delim(cov.bin.neg.tmp, header = FALSE)
  
  cov.bin.df <- cbind(cov.bin.pos[, c(1:4, 7)], cov.bin.neg$V7, cov.bin.pos$V7 + cov.bin.neg$V7)
  colnames(cov.bin.df) <- c("chr", "start", "end", "bin", "POS", "NEG", "TOT")
  write.xlsx(cov.bin.df, file = file.path(outDir, paste0(sampleName, "_", gsub("000$", "kb", window.size), "_bins.xlsx")),
             row.names = FALSE, overwrite = TRUE)
  
  ###############
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

ONreadout <- function(bamFile, otsFile, gRNA.orientation, window.size = 5000, sampleName, outDir){
  # bamFile can be of length 1 or > 1
  bamFile.str <- paste(bamFile, collapse = " ")
  
  # TOTAL COVERAGE (USED FOR CPM)
  bamFile.size <- sapply(bamFile, function(x) as.numeric(nbReadBAM(x)))
  
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
  
  # Initialize data with first bam file
  # calculate coverage per base
  cov.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile.str, bed = ots.tmp, outFile = cov.pos.tmp, opt.string="-bed -d -s")
  cov.pos <- read.delim(cov.pos.tmp, header = FALSE)
  
  cov.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile.str, bed = ots.tmp, outFile = cov.neg.tmp, opt.string="-bed -d -S")
  cov.neg <- read.delim(cov.neg.tmp, header = FALSE)
  
  # calculate total coverage
  cov.tot.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile.str, bed = ots.tmp, outFile = cov.tot.pos.tmp, opt.string="-bed -s")
  cov.tot.pos <- read.delim(cov.tot.pos.tmp, header = FALSE)
  cov.tot.pos <- cov.tot.pos$V7
  
  cov.tot.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile.str, bed = ots.tmp, outFile = cov.tot.neg.tmp, opt.string="-bed -S")
  cov.tot.neg <- read.delim(cov.tot.neg.tmp, header = FALSE)
  cov.tot.neg <- cov.tot.neg$V7
  
  cov.tot <- cov.tot.pos + cov.tot.neg
  cov.df <- data.frame(TOTAL = cov.tot,
                       POS = cov.tot.pos,
                       NEG = cov.tot.neg)
  rownames(cov.df) <- "Read coverage"
  #rownames(cov.df) <- basename(bamFile)
  
  # Adjust percentages
  cov.df$POS.PC <- round(cov.df$POS / cov.df$TOTAL * 100, digits = 2)
  cov.df$NEG.PC <- round(cov.df$NEG / cov.df$TOTAL * 100, digits = 2)
  
  
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
             rowNames = TRUE, overwrite = TRUE)
  
  # calculate coverage per bins
  bins <- chunkBed(ots$V2, ots$V3, 100)
  bins <- data.frame(ots$V1,
                     bins,
                     paste0("bin_", seq(1:nrow(bins))),
                     1000,
                     "+")
  colnames(bins) <- colnames(ots)
  bins.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  write.table(bins, bins.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cov.bin.pos.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile.str, bed = bins.tmp, outFile = cov.bin.pos.tmp, opt.string="-bed -s")
  cov.bin.pos <- read.delim(cov.bin.pos.tmp, header = FALSE)
  
  cov.bin.neg.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
  bamCoverage(bam = bamFile.str, bed = bins.tmp, outFile = cov.bin.neg.tmp, opt.string="-bed -S")
  cov.bin.neg <- read.delim(cov.bin.neg.tmp, header = FALSE)
  
  cov.bin.df <- cbind(cov.bin.pos[, c(1:4, 7)], cov.bin.neg$V7, cov.bin.pos$V7 + cov.bin.neg$V7)
  colnames(cov.bin.df) <- c("chr", "start", "end", "bin", "POS", "NEG", "TOT")
  write.xlsx(cov.bin.df, file = file.path(outDir, paste0(sampleName, "_", gsub("000$", "kb", window.size), "_bins.xlsx")),
             rowNames = FALSE, overwrite = TRUE)
  
  ###############
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
  
  # COUNT 2 CPM
  ggmat$COV.RAW <- ggmat$COV
  ggmat$COV <- log2(((ggmat$COV / sum(bamFile.size))*10^6)+1)
  
  # COVERAGE
  p <- ggplot(ggmat, aes(x= POSITION))
  p <- p + geom_line(aes(y=COV, color = STRAND), size = 1, alpha = 0.75)
  p <- p + geom_ribbon(aes(ymin = 0, ymax = COV, fill = STRAND), alpha = 0.2)
  p <- p + scale_color_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                       NEG = rgb(175, 174, 239, maxColorValue = 255)))
  p <- p + scale_fill_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                      NEG = rgb(175, 174, 239, maxColorValue = 255)))
  if(max(ggmat$COV)>10000){
    p <- p + scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) 
    p <- p + annotation_logticks() 
  }
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())

  p <- p + geom_vline(xintercept=ots.site, linetype="dashed", size = 1)
  p <- p + xlab("") + ylab("read coverage")
  p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
  p <- p + ggtitle(ots.title)
  ggsave(plot = p, filename = file.path(outDir, paste0(sampleName, "_coverage_ON_", gsub("000$", "kb", window.size), ".pdf")), 
         width = 8, height = 5)
  
  write.table(ggmat, file.path(outDir, paste0(sampleName, "_coverage_ON_", gsub("000$", "kb", window.size), ".data.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # DEL / INV
  p <- ggplot(ggmat, aes(x= POSITION))
  p <- p + geom_line(aes(y=COV, color = EVENT), size = 1, alpha = 0.75)
  p <- p + geom_ribbon(aes(ymin = 0, ymax = COV, fill = EVENT), alpha = 0.2)
  p <- p + scale_color_manual(values=c(DEL = "orange2",
                                       INV = "orchid2"))
  p <- p + scale_fill_manual(values=c(DEL = "orange2",
                                      INV = "orchid2"))
  if(max(ggmat$COV)>10000){
    p <- p + scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) 
    p <- p + annotation_logticks() 
  }
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
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

coverage_single <- function(inputFile, bamFile, window.size = 100, outDir){
  # LOAD READMAT
  readMat <- read.xlsx(inputFile, sheet = 1)  
  readMat <- readMat[readMat$hits > 1,]
  
  if(nrow(readMat)==0){
    print(paste0("no site with hits > 1 in ", inputFile))
    return(NA)
  }
  
  mclapply(1:nrow(readMat), function(i){
    subDir <- readMat$group[i]
    dir.create(file.path(outDir, subDir), recursive = TRUE, showWarnings = FALSE)
    
    # CREATE SITE BED
    ots <- data.frame(V1 = readMat$chromosome[i],
                      V2 = readMat$start[i],
                      V3 = readMat$end[i],
                      V4 = readMat$SYMBOL[i],
                      V5 = 1000,
                      V6 = "*")
    ots.title <- paste0(ots$V1, ": ", ots$V2, " - ", ots$V3, "; ", ots$V4)
    ots.name <- paste0(ots$V1, "_", ots$V2, "_", ots$V3, "_", ots$V4)
    
    ots.site <- ots$V2 + floor((ots$V3 - ots$V2)/2)
    
    dir.create(dir1 <- file.path(tempdir(), paste0("cleavageDir", i)))
    
    # extend ots
    ots$V2 <- ots$V2 - window.size
    ots$V3 <- ots$V3 + window.size
    ots.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    write.table(ots, ots.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # calculate coverage per base
    cov.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tmp, opt.string="-bed -d")
    cov <- read.delim(cov.tmp, header = FALSE)
    
    # calculate total coverage
    cov.tot.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.tmp, opt.string="-bed")
    cov.tot <- read.delim(cov.tot.tmp, header = FALSE)
    cov.tot <- cov.tot$V7
    
    # calculate coverage per bins
    bins <- chunkBed(ots$V2, ots$V3, 100)
    bins <- data.frame(ots$V1,
                       bins,
                       paste0("bin_", seq(1:nrow(bins))),
                       1000,
                       "*")
    colnames(bins) <- colnames(ots)
    bins.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    write.table(bins, bins.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cov.bin.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    bamCoverage(bam = bamFile, bed = bins.tmp, outFile = cov.bin.tmp, opt.string="-bed")
    cov.bin <- read.delim(cov.bin.tmp, header = FALSE)
    
    cov.bin.df <- data.frame(cov.bin[, c(1:4, 7)],
                             SYMBOL = readMat$SYMBOL[i])
    colnames(cov.bin.df) <- c("chr", "start", "end", "bin", "TOT", "SYMBOL")
    write.xlsx(cov.bin.df, file = file.path(outDir, subDir, paste0(ots.name, "_", gsub("000$", "kb", window.size), "_bins.xlsx")),
               row.names = FALSE, overwrite = TRUE)
    
    ###############
    # plot coverage
    ggmat <- data.frame(CHR = cov$V1,
                        POSITION = cov$V2 + cov$V7,
                        COV = cov$V8)
    
    # COVERAGE
    p <- ggplot(ggmat, aes(x= POSITION))
    p <- p + geom_line(aes(y=COV), size = 1, alpha = 0.75)
    
    p <- p + geom_vline(xintercept = readMat$start[i], 
                        color = "orangered", size=1)
    p <- p + geom_vline(xintercept = readMat$end[i], 
                        color = "orangered", size=1)
    
    if(readMat$group[i] != "NBS"){
      p <- p + geom_vline(xintercept = readMat$aln.start.abs[i], linetype="dashed", 
                          color = "limegreen", size=1.5)
    }
    
    if(max(ggmat$COV >= 10000)){
      p <- p + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
      p <- p + annotation_logticks() 
    }
    
    p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
    p <- p + xlab("") + ylab("read coverage")
    p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
    p <- p + ggtitle(ots.title)
    ggsave(plot = p, filename = file.path(outDir, subDir, paste0(ots.name, "_", gsub("000$", "kb", window.size), "_bins.pdf")), 
           width = 6, height = 4)
    
    # Delete tmp dir
    unlink(dir1, recursive = T)
    
    
  }, mc.cores = NBCPU)
}

coverage_double <- function(inputFile, bamFile, window.size = 100, outDir){
    # LOAD READMAT
    readMat <- read.xlsx(inputFile, sheet = 1)  
    readMat <- readMat[readMat$hits > 1,]
    
    if(nrow(readMat)==0){
      print(paste0("no site with hits > 1 in ", inputFile))
      return(NA)
    }
    
    mclapply(1:nrow(readMat), function(i){
      subDir <- readMat$group[i]
      dir.create(file.path(outDir, subDir), recursive = TRUE, showWarnings = FALSE)
      
      # CREATE SITE BED
      ots <- data.frame(V1 = readMat$chromosome[i],
                        V2 = readMat$start[i],
                        V3 = readMat$end[i],
                        V4 = readMat$SYMBOL[i],
                        V5 = 1000,
                        V6 = "*")
      ots.title <- paste0(ots$V1, ": ", ots$V2, " - ", ots$V3, "; ", ots$V4)
      ots.name <- paste0(ots$V1, "_", ots$V2, "_", ots$V3, "_", ots$V4)
      
      ots.site <- ots$V2 + floor((ots$V3 - ots$V2)/2)
      
      dir.create(dir1 <- file.path(tempdir(), paste0("cleavageDir", i)))
      
      # extend ots
      ots$V2 <- ots$V2 - window.size
      ots$V3 <- ots$V3 + window.size
      ots.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
      write.table(ots, ots.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      # calculate coverage per base
      cov.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
      bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tmp, opt.string="-bed -d")
      cov <- read.delim(cov.tmp, header = FALSE)
      
      # calculate total coverage
      cov.tot.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
      bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.tmp, opt.string="-bed")
      cov.tot <- read.delim(cov.tot.tmp, header = FALSE)
      cov.tot <- cov.tot$V7
      
      # calculate coverage per bins
      bins <- chunkBed(ots$V2, ots$V3, 100)
      bins <- data.frame(ots$V1,
                         bins,
                         paste0("bin_", seq(1:nrow(bins))),
                         1000,
                         "*")
      colnames(bins) <- colnames(ots)
      bins.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
      write.table(bins, bins.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
      cov.bin.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
      bamCoverage(bam = bamFile, bed = bins.tmp, outFile = cov.bin.tmp, opt.string="-bed")
      cov.bin <- read.delim(cov.bin.tmp, header = FALSE)
      
      cov.bin.df <- data.frame(cov.bin[, c(1:4, 7)],
                               SYMBOL = readMat$SYMBOL[i])
      colnames(cov.bin.df) <- c("chr", "start", "end", "bin", "TOT", "SYMBOL")
      write.xlsx(cov.bin.df, file = file.path(outDir, subDir, paste0(ots.name, "_", gsub("000$", "kb", window.size), "_bins.xlsx")),
                 row.names = FALSE, overwrite = TRUE)
      
      ###############
      # plot coverage
      ggmat <- data.frame(CHR = cov$V1,
                          POSITION = cov$V2 + cov$V7,
                          COV = cov$V8)
      
      # COVERAGE
      p <- ggplot(ggmat, aes(x= POSITION))
      p <- p + geom_line(aes(y=COV), size = 1, alpha = 0.75)
      
      p <- p + geom_vline(xintercept = readMat$start[i], 
                          color = "orangered", size=1)
      p <- p + geom_vline(xintercept = readMat$end[i], 
                          color = "orangered", size=1)
      
      if(readMat$group[i] != "NBS"){
        p <- p + geom_vline(xintercept = readMat$aln.start.abs.gRNA1[i], linetype="dashed", 
                            color = "limegreen", alpha = 0.5, size=1.5)
        p <- p + geom_vline(xintercept = readMat$aln.start.abs.gRNA2[i], linetype="dashed", 
                            color = "mediumpurple", alpha = 0.5, size=1.5)
      }
      
      if(max(ggmat$COV >= 10000)){
        p <- p + scale_y_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
        p <- p + annotation_logticks() 
      }
      
      p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
      p <- p + xlab("") + ylab("read coverage")
      p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
      p <- p + ggtitle(ots.title)
      ggsave(plot = p, filename = file.path(outDir, subDir, paste0(ots.name, "_", gsub("000$", "kb", window.size), "_bins.pdf")), 
             width = 6, height = 4)
      
      # Delete tmp dir
      unlink(dir1, recursive = T)
      
      
    }, mc.cores = NBCPU)
}

coverage_talenOLD <- function(inputFile, bamFile, window.size = 100, outDir){
  # LOAD READMAT
  readMat <- read.xlsx(inputFile, sheet = 1)  
  readMat <- readMat[readMat$hits > 1,]
  
  if(nrow(readMat)==0){
    print(paste0("no site with hits > 1 in ", inputFile))
    return(NA)
  }
  
  mclapply(1:nrow(readMat), function(i){
    subDir <- readMat$group[i]
    dir.create(file.path(outDir, subDir), recursive = TRUE, showWarnings = FALSE)
    
    # CREATE SITE BED
    ots <- data.frame(V1 = readMat$chromosome[i],
                      V2 = readMat$start[i],
                      V3 = readMat$end[i],
                      V4 = readMat$SYMBOL[i],
                      V5 = 1000,
                      V6 = "*")
    ots.title <- paste0(ots$V1, ": ", ots$V2, " - ", ots$V3, "; ", ots$V4)
    ots.name <- paste0(ots$V1, "_", ots$V2, "_", ots$V3, "_", ots$V4)
    
    ots.site <- ots$V2 + floor((ots$V3 - ots$V2)/2)
    
    dir.create(dir1 <- file.path(tempdir(), paste0("cleavageDir", i)))
    
    # extend ots
    ots$V2 <- ots$V2 - window.size
    ots$V3 <- ots$V3 + window.size
    ots.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    write.table(ots, ots.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # calculate coverage per base
    cov.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tmp, opt.string="-bed -d")
    cov <- read.delim(cov.tmp, header = FALSE)
    
    # calculate total coverage
    cov.tot.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    bamCoverage(bam = bamFile, bed = ots.tmp, outFile = cov.tot.tmp, opt.string="-bed")
    cov.tot <- read.delim(cov.tot.tmp, header = FALSE)
    cov.tot <- cov.tot$V7
    
    # calculate coverage per bins
    bins <- chunkBed(ots$V2, ots$V3, 100)
    bins <- data.frame(ots$V1,
                       bins,
                       paste0("bin_", seq(1:nrow(bins))),
                       1000,
                       "*")
    colnames(bins) <- colnames(ots)
    bins.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    write.table(bins, bins.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cov.bin.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    bamCoverage(bam = bamFile, bed = bins.tmp, outFile = cov.bin.tmp, opt.string="-bed")
    cov.bin <- read.delim(cov.bin.tmp, header = FALSE)
    
    cov.bin.df <- data.frame(cov.bin[, c(1:4, 7)],
                             SYMBOL = readMat$SYMBOL[i])
    colnames(cov.bin.df) <- c("chr", "start", "end", "bin", "TOT", "SYMBOL")
    write.xlsx(cov.bin.df, file = file.path(outDir, subDir, paste0(ots.name, "_", gsub("000$", "kb", window.size), "_bins.xlsx")),
               row.names = FALSE, overwrite = TRUE)
    
    ###############
    # plot coverage
    ggmat <- data.frame(CHR = cov$V1,
                        POSITION = cov$V2 + cov$V7,
                        COV = cov$V8)
    
    # COVERAGE
    p <- ggplot(ggmat, aes(x= POSITION))
    p <- p + geom_line(aes(y=COV), size = 1, alpha = 0.75)
    
    p <- p + geom_vline(xintercept = readMat$start[i], 
                        color = "orangered", size=1)
    p <- p + geom_vline(xintercept = readMat$end[i], 
                        color = "orangered", size=1)
    

    
    if(max(ggmat$COV >= 10000)){
      p <- p + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
      p <- p + annotation_logticks() 
    }
    
    if(!is.na(readMat$BestCB[i]) & readMat$group[i] != "NBS"){
      cb <- readMat$BestCB[i]
      gRNA1.x <- readMat[i, paste0(cb, "_aln.start.abs")]
      aln.length <- nchar(readMat[i, paste0(cb, "_subject")])
      gRNA2.x <- gRNA1.x + aln.length
      
      if(substring(cb, 1, 1) == "L"){
        p <- p + geom_vline(xintercept = gRNA1.x, linetype="dashed", 
                            color = "limegreen", size=1.5)
        p <- p + geom_vline(xintercept = gRNA2.x, linetype="dashed", 
                            color = "mediumpurple", size=1.5)
      }else{
        p <- p + geom_vline(xintercept = gRNA2.x, linetype="dashed", 
                            color = "limegreen", size=1.5)
        p <- p + geom_vline(xintercept = gRNA1.x, linetype="dashed", 
                            color = "mediumpurple", size=1.5)
      }
      

    }
    
    p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
    p <- p + xlab("") + ylab("read coverage")
    p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
    p <- p + ggtitle(ots.title)
    ggsave(plot = p, filename = file.path(outDir, subDir, paste0(ots.name, "_", gsub("000$", "kb", window.size), "_bins.pdf")), 
           width = 6, height = 4)
    
    # Delete tmp dir
    unlink(dir1, recursive = T)
    
    
  }, mc.cores = NBCPU)
}

coverage_cleavage_crispr <- function(inputFile, covDir, bin.size = 100){
  # LOAD READMAT
  readMat <- read.xlsx(inputFile, sheet = 1)  
  readMat <- readMat[readMat$hits > 1,]
  
  if(nrow(readMat)==0){
    print(paste0("no site with hits > 1 in ", inputFile))
    return(NA)
  }
  
  mclapply(1:nrow(readMat), function(i){
    chr <- readMat$chromosome[i]
    start <- readMat$start.site[i]
    end <- readMat$end.site[i]
    
    covFile <- file.path(covDir, paste(chr, start, end, bin.size, "bins.xlsx", sep = "_"))
    cov.bin.df <- read.xlsx(covFile)
    
    # COVERAGE
    bamLabels <- colnames(cov.bin.df)[seq(5, ncol(cov.bin.df))]
    ggmatList <- lapply(bamLabels, function(j){
      return(data.frame(CHR = cov.bin.df$chr,
                        POSITION = cov.bin.df$start,
                        COV = cov.bin.df[, j],
                        SAMPLE = j))
    })
    ggmat <- do.call(rbind, ggmatList)
    
    # COVERAGE
    ots.title <- paste0(chr, ": ", start, " - ", end)
    ots.name <- paste0(chr, "_", start, "_", end)
    
    p <- ggplot(ggmat, aes(x= POSITION))
    p <- p + geom_line(aes(y=COV, color = SAMPLE), size = 1, alpha = 0.75)
    
    p <- p + geom_vline(xintercept = readMat$start.site[i], 
                        color = "orangered", size=1)
    p <- p + geom_vline(xintercept = readMat$end.site[i], 
                        color = "orangered", size=1)
    
    p <- p + geom_vline(xintercept = readMat$start.bin[i], 
                        color = "grey20", size=0.75, linetype = "dashed")
    p <- p + geom_vline(xintercept = readMat$end.bin[i], 
                        color = "grey20", size=0.75, linetype = "dashed")
    
    if(max(ggmat$COV >= 10000)){
      p <- p + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
      p <- p + annotation_logticks() 
    }
    
    #cb.color <- brewer.pal(n = 10, name = "Set3")
    cb.color <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
    #cb.color <- c("red", "blue", "green")
    
    gRNA1.x <- readMat[i, "aln.start.abs"]
    aln.length <- nchar(readMat[i,  "subject"])
    gRNA2.x <- gRNA1.x + aln.length
    p <- p + annotate("rect",
                      xmin = gRNA1.x, xmax = gRNA2.x,
                      ymin = 0, ymax = max(ggmat$COV),
                      alpha = 0.8, fill = cb.color[1])
    
    p <- p + theme_bw(base_size = 16) + theme(panel.grid.minor = element_blank())
    p <- p + xlab("") + ylab("read coverage")
    p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold"))
    p <- p + ggtitle(ots.title)
    ggsave(plot = p, filename = file.path(covDir, paste0(ots.name, "_", bin.size, "_bins_cleavage.pdf")), 
           width = 6, height = 4)
    
  }, mc.cores = NBCPU)
}

coverage_cleavage_double <- function(inputFile, covDir, bin.size = 100){
  # LOAD READMAT
  readMat <- read.xlsx(inputFile, sheet = 1)  
  readMat <- readMat[readMat$hits > 1,]
  
  if(nrow(readMat)==0){
    print(paste0("no site with hits > 1 in ", inputFile))
    return(NA)
  }
  
  mclapply(1:nrow(readMat), function(i){
    chr <- readMat$chromosome[i]
    start <- readMat$start.site[i]
    end <- readMat$end.site[i]
    
    covFile <- file.path(covDir, paste(chr, start, end, bin.size, "bins.xlsx", sep = "_"))
    cov.bin.df <- read.xlsx(covFile)
    
    # COVERAGE
    bamLabels <- colnames(cov.bin.df)[seq(5, ncol(cov.bin.df))]
    ggmatList <- lapply(bamLabels, function(j){
      return(data.frame(CHR = cov.bin.df$chr,
                        POSITION = cov.bin.df$start,
                        COV = cov.bin.df[, j],
                        SAMPLE = j))
    })
    ggmat <- do.call(rbind, ggmatList)
    
    # COVERAGE
    ots.title <- paste0(chr, ": ", start, " - ", end)
    ots.name <- paste0(chr, "_", start, "_", end)
    
    p <- ggplot(ggmat, aes(x= POSITION))
    p <- p + geom_line(aes(y=COV, color = SAMPLE), size = 1, alpha = 0.75)
    
    p <- p + geom_vline(xintercept = readMat$start.site[i], 
                        color = "orangered", size=1)
    p <- p + geom_vline(xintercept = readMat$end.site[i], 
                        color = "orangered", size=1)
    
    p <- p + geom_vline(xintercept = readMat$start.bin[i], 
                        color = "grey20", size=0.75, linetype = "dashed")
    p <- p + geom_vline(xintercept = readMat$end.bin[i], 
                        color = "grey20", size=0.75, linetype = "dashed")
    
    if(max(ggmat$COV >= 10000)){
      p <- p + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
      p <- p + annotation_logticks() 
    }

    p <- p + geom_vline(xintercept = readMat$aln.start.abs.gRNA1[i], linetype="dashed", 
                        color = "limegreen", alpha = 0.5, size=1.5)
    p <- p + geom_vline(xintercept = readMat$aln.start.abs.gRNA2[i], linetype="dashed", 
                        color = "mediumpurple", alpha = 0.5, size=1.5)
    
    p <- p + theme_bw(base_size = 16) + theme(panel.grid.minor = element_blank())
    p <- p + xlab("") + ylab("read coverage")
    p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold"))
    p <- p + ggtitle(ots.title)
    ggsave(plot = p, filename = file.path(covDir, paste0(ots.name, "_", bin.size, "_bins_cleavage.pdf")), 
           width = 6, height = 4)
    
  }, mc.cores = NBCPU)
}


coverage_cleavage_talen <- function(inputFile, covDir, bin.size = 100){
  # LOAD READMAT
  readMat <- read.xlsx(inputFile, sheet = 1)  
  readMat <- readMat[readMat$hits > 1,]
  
  if(nrow(readMat)==0){
    print(paste0("no site with hits > 1 in ", inputFile))
    return(NA)
  }
  
  mclapply(1:nrow(readMat), function(i){
    chr <- readMat$chromosome[i]
    start <- readMat$start.site[i]
    end <- readMat$end.site[i]
    
    covFile <- file.path(covDir, paste(chr, start, end, bin.size, "bins.xlsx", sep = "_"))
    cov.bin.df <- read.xlsx(covFile)
    
    # COVERAGE
    bamLabels <- colnames(cov.bin.df)[seq(5, ncol(cov.bin.df))]
    ggmatList <- lapply(bamLabels, function(j){
      return(data.frame(CHR = cov.bin.df$chr,
                        POSITION = cov.bin.df$start,
                        COV = cov.bin.df[, j],
                        SAMPLE = j))
    })
    ggmat <- do.call(rbind, ggmatList)
    
    # COVERAGE
    ots.title <- paste0(chr, ": ", start, " - ", end)
    ots.name <- paste0(chr, "_", start, "_", end)
    
    p <- ggplot(ggmat, aes(x= POSITION))
    p <- p + geom_line(aes(y=COV, color = SAMPLE), size = 1, alpha = 0.75)
    
    p <- p + geom_vline(xintercept = readMat$start.site[i], 
                        color = "orangered", size=1)
    p <- p + geom_vline(xintercept = readMat$end.site[i], 
                        color = "orangered", size=1)
    
    p <- p + geom_vline(xintercept = readMat$start.bin[i], 
                        color = "grey20", size=0.75, linetype = "dashed")
    p <- p + geom_vline(xintercept = readMat$end.bin[i], 
                        color = "grey20", size=0.75, linetype = "dashed")
    
    if(max(ggmat$COV >= 10000)){
      p <- p + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
      p <- p + annotation_logticks() 
    }
    
    #cb.color <- brewer.pal(n = 10, name = "Set3")
    cb.color <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
    #cb.color <- c("red", "blue", "green")
    
    for(x in 1:sum(grepl("_COMB$", colnames(readMat)))){
      cb <- readMat[i, paste0("top", x, "_COMB")]
      gRNA1.x <- readMat[i, paste0("top", x, "_aln.start.abs")]
      aln.length <- nchar(readMat[i,  paste0("top", x, "_subject")])
      gRNA2.x <- gRNA1.x + aln.length
      p <- p + annotate("rect",
                        xmin = gRNA1.x, xmax = gRNA2.x,
                        ymin = 0, ymax = max(ggmat$COV),
                        alpha = .5, fill = cb.color[x])
    }
    
    p <- p + theme_bw(base_size = 16) + theme(panel.grid.minor = element_blank())
    p <- p + xlab("") + ylab("read coverage")
    p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold"))
    p <- p + ggtitle(ots.title)
    ggsave(plot = p, filename = file.path(covDir, paste0(ots.name, "_", bin.size, "_bins_cleavage.pdf")), 
           width = 6, height = 4)
    
  }, mc.cores = NBCPU)
}

coverage_bin <- function(inputFile, outputFile, bamFiles, bamLabels, bin.size = 100, top = "sum", d = 100, plot = FALSE, outDir){
  # LOAD READMAT
  readMat <- read.xlsx(inputFile, sheet = 1)  
  readMat <- readMat[readMat$hits > 1,]
  
  if(nrow(readMat)==0){
    print(paste0("no site with hits > 1 in ", inputFile))
    return(NA)
  }
  
  bestCoordList <- lapply(1:nrow(readMat), function(i){
    # CREATE SITE BED
    ots <- data.frame(V1 = readMat$chromosome[i],
                      V2 = readMat$start[i],
                      V3 = readMat$end[i],
                      V4 = readMat$siteID[i],
                      V5 = 1000,
                      V6 = "*")
    ots.title <- paste0(ots$V1, ": ", ots$V2, " - ", ots$V3)
    ots.name <- paste0(ots$V1, "_", ots$V2, "_", ots$V3)
    
    ots.site <- ots$V2 + floor((ots$V3 - ots$V2)/2)
    
    dir.create(dir1 <- file.path(tempdir(), paste0("cleavageDir", i)))
    
    # calculate coverage per bins
    bins <- chunkBed(ots$V2, ots$V3, bin.size)
    bins <- data.frame(ots$V1,
                       bins,
                       paste0("bin_", seq(1:nrow(bins))),
                       1000,
                       "*")
    colnames(bins) <- colnames(ots)
    bins.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    write.table(bins, bins.tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cov.bin.tmp <- tempfile(tmpdir = dir1, fileext = ".bed")
    
    cov.bin.df.list <- lapply(bamFiles, function(bam){
      bamCoverage(bam = bam, bed = bins.tmp, outFile = cov.bin.tmp, opt.string="-bed")
      cov.bin <- read.delim(cov.bin.tmp, header = FALSE)
      
      cov.bin.df <- data.frame(cov.bin[, c(1:4, 7)])
      colnames(cov.bin.df) <- c("chr", "start", "end", "bin", "TOT")
      return(cov.bin.df)
    })
    
    cov.bin.df <- do.call(cbind, lapply(cov.bin.df.list, function(j) j$TOT))
    colnames(cov.bin.df) <- bamLabels
    cov.bin.df <- data.frame(cov.bin.df.list[[1]][, 1:4], cov.bin.df)
    
    write.xlsx(cov.bin.df, file = file.path(outDir, paste0(ots.name, "_", bin.size, "_bins.xlsx")),
               row.names = FALSE, overwrite = TRUE)
    
    # SELECT BEST COVERAGE
    rel.dist <- abs(cov.bin.df$start - median(cov.bin.df$start))
    
    idx <- NA
    if(top == "sum"){
      score <- apply(cov.bin.df[, -c(1:4), drop = FALSE], 1, sum)
    }else if(top == "mean"){  
      score <- apply(cov.bin.df[, -c(1:4), drop = FALSE], 1, mean)
    }else("Unkown selection metric, please select between sum and mean")
    
    score.max <- max(score)
    idx <- which(score == score.max)
    if(length(idx) != 1){
      idx <- which(score == score.max)
      idx <- idx[which.min(rel.dist[idx])]
    }
    bc <- floor(cov.bin.df$start[idx] + ((cov.bin.df$end[idx] - cov.bin.df$start[idx]) / 2))

    if(plot){
      ###############
      # plot coverage
      ggmatList <- lapply(bamLabels, function(j){
        return(data.frame(CHR = cov.bin.df$chr,
                          POSITION = cov.bin.df$start,
                          COV = cov.bin.df[, j],
                          SAMPLE = j))
      })
      ggmat <- do.call(rbind, ggmatList)
      
      # COVERAGE
      p <- ggplot(ggmat, aes(x= POSITION))
      p <- p + geom_line(aes(y=COV, color = SAMPLE), size = 1, alpha = 0.75)
      
      p <- p + geom_vline(xintercept = readMat$start[i], 
                          color = "orangered", size=1)
      p <- p + geom_vline(xintercept = readMat$end[i], 
                          color = "orangered", size=1)
      p <- p + geom_vline(xintercept = bc, 
                          color = "blue", size=1)
      
      if(max(ggmat$COV >= 10000)){
        p <- p + scale_y_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
        p <- p + annotation_logticks() 
      }
      
      p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())
      p <- p + xlab("") + ylab("read coverage")
      p <- p + theme(legend.position = "bottom", plot.title = element_text(size = 18, face = "bold"))
      p <- p + ggtitle(ots.title)
      ggsave(plot = p, filename = file.path(outDir, paste0(ots.name, "_", bin.size, "_bins.pdf")), 
             width = 6, height = 4)
    }

    # Delete tmp dir
    unlink(dir1, recursive = T)
    
    return(bc)
  })
  
  bestCoord <- unlist(bestCoordList)
  bestCoord.start <- bestCoord - d
  bestCoord.end <- bestCoord + d
  
  # save site coordinates
  readMat$start.site <- readMat$start
  readMat$end.site <- readMat$end
  
  # save bin coordinates
  readMat$start.bin <- bestCoord.start
  readMat$end.bin <- bestCoord.end
  
  # set bin coordinate as default 
  readMat$start <- readMat$start.bin
  readMat$end <- readMat$end.bin
  
  write.xlsx(readMat, outputFile, overwrite = TRUE)
}

extend_bin <- function(inputFile, outputFile, otsF, ext.size = 200){
  readMat <- read.xlsx(inputFile, sheet = 1) 
  
  isON <- getONidx(readMat, otsF)
  isON[is.na(isON)] <- "no"
  
  if(sum(isON == "yes") != 0){
    idx <- which(isON == "yes")
    readMat$start.bin[idx] <- readMat$start.bin[idx] - ext.size
    readMat$end.bin[idx] <- readMat$end.bin[idx] + ext.size
    
    readMat$start <- readMat$start.bin
    readMat$end <- readMat$end.bin
  }
  
  write.xlsx(readMat, outputFile, overwrite = TRUE)
}

# DO NOT RUN
if(FALSE){
  library(ggplot2)
  library(openxlsx)
  library(parallel)
  
  bamCoverage <- function(functionstring="/Users/geoffroyandrieux/Programs/bedtools2-2.30.0/bin/bedtools coverage", bam, bed, outFile, opt.string="")
  {
    options(scipen =99)
    
    command=paste(functionstring, "-a", bed, "-b", bam, opt.string, " | gzip >", outFile)
    cat(command,"\n")
    try(system(command))
  }
  
  nbReadBAM <- function(inFile){
         options(scipen =99)
         command=paste0("/Users/geoffroyandrieux/Programs/samtools-1.16.1/samtools view -c -F 260 ", inFile)
         cat(command,"\n")
         try(system(command, intern = TRUE))
    }
  
  ########################################################################
  # MULTI BAM TEST
  
  setwd("~/Research/CASTSeq/Test/ONreadout/")
  bamFiles <- c("GE-285-1_AlignmentSort.bam", "GE-285-2_AlignmentSort.bam")
  otsFile <- c("ots.bed")
  gRNA.orientation <- "forward"
  window.size <- 5000
  sampleName <- "GE-285"
  outDir <- "~/Research/CASTSeq/Test/ONreadout/"
  
  ONreadout(bamFile = bamFiles,
            otsFile = otsFile,
            gRNA.orientation = gRNA.orientation,
            window.size = window.size,
            sampleName = sampleName,
            outDir = outDir)
  
  ########################################################################
  # COVERAGE PLOT TEST
  
  inputFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9/KRT9-T1_g1/results/guide_aln/KRT9-T1_w250_FINAL.xlsx")
  bamFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9/KRT9-T1_g1/results/fastq_aln/KRT9-T1_AlignmentSort.bam")
  window.size = 1000
  outDir <- file.path("~/Research/CASTSeq/test/")
  
  
  start_time <- Sys.time()
  
  coverage_double(inputFile, bamFile, window.size = 1000, outDir)
    
  end_time <- Sys.time()
  end_time - start_time
  
  # TALEN
  inputFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/STAT3/STAT3-T3edited-1/results/guide_aln/STAT3-T3edited-1_w250_FINAL.xlsx")
  bamFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/STAT3/STAT3-T3edited-1/results/fastq_aln/STAT3-T3edited-1_AlignmentSort.bam")
  window.size = 100
  outDir <- file.path("~/Research/CASTSeq/test/")
  
  coverage_talen(inputFile, bamFile, window.size = 100, outDir)
  
  ########################################################################
  # ON READOUT TEST
  
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









