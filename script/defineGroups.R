#library(Biostrings)
#library(openxlsx)
#library(data.table)
#library(ggplot2)
#library(ggseqlogo)
#library(textreadr)
#library(parallel)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getEmpiricalPV <- function(x, y, type)
{
    x.ecdf <- stats::ecdf(x)
    if(type == "greater") return(1-x.ecdf(y))
    else if(type == "less") return(x.ecdf(y))
    return(NA)    
}


getONidx <- function(realM, otsF){
  
  ots.df <- read.delim(otsF, header = FALSE)
  ots.df <- ots.df[1,,drop = FALSE]# only take the first ON target into account
  ots.gr <- makeGRangesFromDataFrame(ots.df,
                                     seqnames.field = "V1", start.field = "V2", end.field = "V3",
                                     keep.extra.columns = FALSE, ignore.strand = TRUE)
  realM.gr <- makeGRangesFromDataFrame(realM,
                                       seqnames.field = "chromosome", start.field = "start", end.field = "end",
                                       keep.extra.columns = FALSE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = ots.gr, subject = realM.gr, type = "any", maxgap = 0)

  isON <- rep(NA, nrow(realM))
  isON[subjectHits(gr.ovl)] <- "yes"
  
  return(isON)
}

assignGroups <- function(realF, rdF, otsF, cutoff)
{
	ots <- read.delim(otsF, header = FALSE)
	print(realF)
	print(file.exists(realF))
	
	print(rdF)
	print(file.exists(rdF))
	
	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	realM.bf <- getBestFlanking(realM)
	rdM.bf <- getBestFlanking(rdM)
	
	score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
	#flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
	flanking.cutoff <- 25
	
	score.pv <- unlist(lapply(realM[, "score"], getEmpiricalPV, x = rdM[, "score"], type = "greater"))
	score.adj.pv <- p.adjust(score.pv, method = "BH")
	flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
	flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")
	#print(score.cutoff)
	#print(flanking.cutoff)	
	
	if(nrow(realM) > 10){
	  # plot score.cutoff (REAL)
	  pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", realF))
	  #png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	  plot(density(realM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	  abline(v=score.cutoff, col="black", lwd=3)
	  dev.off()
	  
	  
	  ggmat <- data.frame(dist = realM[, "score"])
	  p <- ggplot(ggmat) + 
	    geom_density(aes(x=dist))
	  p <- p + annotate("rect", xmin = score.cutoff,
	                    xmax = Inf,
	                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	  p <- p + geom_vline(xintercept=score.cutoff, 
	                      color = "black", size=2)
	  p <- p + theme_bw(base_size = 16)
	  p <- p + xlab("substring length (bp)")
	  p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff))
	  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", realF),
	         width = 6, height = 6)
	  
	}

	if(length(realM.bf) > 10){
  	# plot flanking.cutoff (REAL)
  	#xmax <- min(100, max(realM.bf))
  	#pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF))
  	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
  	#plot(density(realM.bf[realM.bf <= 100]), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
  	#abline(v=flanking.cutoff, col="black", lwd=3)
  	#dev.off()
  	
  
  	ggmat <- data.frame(dist = realM.bf[realM.bf <= 100])
  	p <- ggplot(ggmat) + 
  	  geom_density(aes(x=dist))
  	p <- p + annotate("rect", xmin = flanking.cutoff,
  	                  xmax = Inf,
  	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  	p <- p + geom_vline(xintercept=flanking.cutoff, 
  	                    color = "black", size=2)
  	p <- p + theme_bw(base_size = 16)
  	p <- p + xlab("substring length (bp)")
  	p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
  	ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF),
  	       width = 6, height = 6)
	}
	
	# plot score.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	ggmat <- data.frame(dist = rdM[, "score"])
	p <- ggplot(ggmat) + 
	  geom_density(aes(x=dist))
	p <- p + annotate("rect", xmin = score.cutoff,
	                  xmax = Inf,
	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	p <- p + geom_vline(xintercept=score.cutoff, 
	                    color = "black", size=2)
	p <- p + theme_bw(base_size = 16)
	p <- p + xlab("guide alignment score")
	p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff))
	ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", rdF),
	       width = 6, height = 6)
	
	
	# plot flanking.cutoff (RANDOM)
	#pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	#plot(density(rdM.bf), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	#abline(v=flanking.cutoff, col="black", lwd=3)
	#dev.off()
	
	ggmat <- data.frame(dist = rdM.bf)
	p <- ggplot(ggmat) + 
	  geom_density(aes(x=dist))
	p <- p + annotate("rect", xmin = flanking.cutoff,
	                  xmax = Inf,
	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	p <- p + geom_vline(xintercept=flanking.cutoff, 
	                    color = "black", size=2)
	p <- p + theme_bw(base_size = 16)
	p <- p + xlab("substring length (bp)")
	p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
	ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF),
	       width = 6, height = 6)
	
	
	#
	mygroups <- rep("NBS", nrow(realM))
	mygroups[realM.bf >= flanking.cutoff] <- "HMT"
	mygroups[score.pv < cutoff] <- "OMT"# QV or PV

	mygroups2 <- rep(NA, nrow(realM))
	#mygroups2[mygroups == "off.target" & (realM.bf > flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "OMT" & (realM.bf >= flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "HMT"] <- "yes"
	
	mygroups3 <- rep(NA, nrow(realM))
	ld.idx <- isLargeDel(realM, ots)
	#print(length(ld.idx))
	if(length(ld.idx)!=0){
		if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
		}
	
	realM$group <- mygroups
	
	realM$is.ON <- getONidx(realM, otsF)
	
	realM$is.HMT <- mygroups2
	realM$is.large.del. <- mygroups3
	
	realM$OMT.pvalue <- score.pv 
	realM$OMT.adj.pvalue <- score.adj.pv 
	realM$HMT.pvalue <- flanking.pv
	realM$HMT.adj.pvalue <- flanking.adj.pv
	
	
	write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE, overwrite = TRUE)	
}

assignGroups2 <- function(realF, rdF, otsF, cutoff)
{
  ots <- read.delim(otsF, header = FALSE)
  print(realF)
  print(file.exists(realF))
  
  print(rdF)
  print(file.exists(rdF))
  
  realM <- read.xlsx(realF, sheet = 1)
  rdM <- read.xlsx(rdF, sheet = 1)
  
  realM.bf <- getBestFlanking(realM)
  rdM.bf <- getBestFlanking(rdM)
  
  score.cutoff1 <- sort(rdM[, "score.gRNA1"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
  score.cutoff2 <- sort(rdM[, "score.gRNA2"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
  
  score <- sapply(1:nrow(realM), function(x) max(c(realM[x, "score.gRNA1"], realM[x, "score.gRNA2"])))
  
  #flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
  flanking.cutoff <- 25
  
  score.pv1 <- unlist(lapply(realM[, "score.gRNA1"], getEmpiricalPV, x = rdM[, "score.gRNA1"], type = "greater"))
  score.adj.pv1 <- p.adjust(score.pv1, method = "BH")
  score.pv2 <- unlist(lapply(realM[, "score.gRNA2"], getEmpiricalPV, x = rdM[, "score.gRNA2"], type = "greater"))
  score.adj.pv2 <- p.adjust(score.pv2, method = "BH")
  
  flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
  flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")
  #print(score.cutoff)
  #print(flanking.cutoff)	
  
  if(nrow(realM) > 10){
    # plot score.cutoff (REAL)
    ggmat <- data.frame(dist = realM[, "score.gRNA1"])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = score.cutoff1,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=score.cutoff1, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff1))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA1.pdf", realF),
           width = 6, height = 6)
    
    ggmat <- data.frame(dist = realM[, "score.gRNA2"])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = score.cutoff2,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=score.cutoff2, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff2))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA2.pdf", realF),
           width = 6, height = 6)
    
  }
  
  if(length(realM.bf) > 10){
    # plot flanking.cutoff (REAL)
    ggmat <- data.frame(dist = realM.bf[realM.bf <= 100])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = flanking.cutoff,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=flanking.cutoff, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF),
           width = 6, height = 6)
  }
  
  # plot score.cutoff (RANDOM)
  ggmat <- data.frame(dist = rdM[, "score.gRNA1"])
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = score.cutoff1,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=score.cutoff1, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("guide alignment score")
  p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff1))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA1.pdf", rdF),
         width = 6, height = 6)
  
  ggmat <- data.frame(dist = rdM[, "score.gRNA2"])
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = score.cutoff2,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=score.cutoff2, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("guide alignment score")
  p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff2))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA2.pdf", rdF),
         width = 6, height = 6)
  
  
  # plot flanking.cutoff (RANDOM)
  ggmat <- data.frame(dist = rdM.bf)
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = flanking.cutoff,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=flanking.cutoff, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("substring length (bp)")
  p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF),
         width = 6, height = 6)
  
  
  score.pv <- sapply(1:length(score.pv1), function(x) min(c(score.pv1[x], score.pv2[x])))
  score.adj.pv <- p.adjust(score.pv, method = "BH")
  
  gRNA <- sapply(1:length(score.pv), function(x){
    grnas <- c("gRNA1", "gRNA2")
    return(toString(grnas[which(c(score.pv1[x], score.pv2[x]) == score.pv[x])]))
  })
  
  #
  mygroups <- rep("NBS", nrow(realM))
  mygroups[realM.bf >= flanking.cutoff] <- "HMT"
  mygroups[score.pv < cutoff] <- "OMT"# QV or PV
  
  mygroups2 <- rep(NA, nrow(realM))
  #mygroups2[mygroups == "off.target" & (realM.bf > flanking.cutoff)] <- "yes"
  mygroups2[mygroups == "OMT" & (realM.bf >= flanking.cutoff)] <- "yes"
  mygroups2[mygroups == "HMT"] <- "yes"
  
  mygroups3 <- rep(NA, nrow(realM))
  ld.idx <- isLargeDel(realM, ots)
  #print(length(ld.idx))
  if(length(ld.idx)!=0){
    if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
  }
  
  realM$group <- mygroups
  realM$gRNA <- gRNA
  realM$score <- score
  
  realM$is.ON <- getONidx(realM, otsF)
  
  realM$is.HMT <- mygroups2
  realM$is.large.del. <- mygroups3
  
  realM$OMT.pvalue <- score.pv 
  realM$OMT.adj.pvalue <- score.adj.pv 
  realM$HMT.pvalue <- flanking.pv
  realM$HMT.adj.pvalue <- flanking.adj.pv
  
  write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE, overwrite = TRUE)	
}

assignGroupsDoubleNickase <- function(realF, rdF, otsF, cutoff)
{
  ots <- read.delim(otsF, header = FALSE)
  print(realF)
  print(file.exists(realF))
  
  print(rdF)
  print(file.exists(rdF))
  
  realM <- read.xlsx(realF, sheet = 1)
  rdM <- read.xlsx(rdF, sheet = 1)
  
  realM.bf <- getBestFlanking(realM)
  rdM.bf <- getBestFlanking(rdM)
  
  score.cutoff1 <- sort(rdM[, "score.gRNA1"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
  score.cutoff2 <- sort(rdM[, "score.gRNA2"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
  score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
  
  #flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
  flanking.cutoff <- 25
  
  score.pv1 <- unlist(lapply(realM[, "score.gRNA1"], getEmpiricalPV, x = rdM[, "score.gRNA1"], type = "greater"))
  score.adj.pv1 <- p.adjust(score.pv1, method = "BH")
  score.pv2 <- unlist(lapply(realM[, "score.gRNA2"], getEmpiricalPV, x = rdM[, "score.gRNA2"], type = "greater"))
  score.adj.pv2 <- p.adjust(score.pv2, method = "BH")
  
  score.pv <- unlist(lapply(realM[, "score"], getEmpiricalPV, x = rdM[, "score"], type = "greater"))
  score.adj.pv <- p.adjust(score.pv, method = "BH")
  
  flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
  flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")
  #print(score.cutoff)
  #print(flanking.cutoff)	
  
  if(nrow(realM) > 10){
    # plot score.cutoff (REAL)
    ggmat <- data.frame(dist = realM[, "score.gRNA1"])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = score.cutoff1,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=score.cutoff1, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff1))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA1.pdf", realF),
           width = 6, height = 6)
    
    ggmat <- data.frame(dist = realM[, "score.gRNA2"])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = score.cutoff2,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=score.cutoff2, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff2))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA2.pdf", realF),
           width = 6, height = 6)
    
    ggmat <- data.frame(dist = realM[, "score"])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = score.cutoff,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=score.cutoff, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", realF),
           width = 6, height = 6)
    
  }
  
  if(length(realM.bf) > 10){
    # plot flanking.cutoff (REAL)
    ggmat <- data.frame(dist = realM.bf[realM.bf <= 100])
    p <- ggplot(ggmat) + 
      geom_density(aes(x=dist))
    p <- p + annotate("rect", xmin = flanking.cutoff,
                      xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
    p <- p + geom_vline(xintercept=flanking.cutoff, 
                        color = "black", size=2)
    p <- p + theme_bw(base_size = 16)
    p <- p + xlab("substring length (bp)")
    p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
    ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF),
           width = 6, height = 6)
  }
  
  # plot score.cutoff (RANDOM)
  ggmat <- data.frame(dist = rdM[, "score.gRNA1"])
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = score.cutoff1,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=score.cutoff1, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("guide alignment score")
  p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff1))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA1.pdf", rdF),
         width = 6, height = 6)
  
  ggmat <- data.frame(dist = rdM[, "score.gRNA2"])
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = score.cutoff2,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=score.cutoff2, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("guide alignment score")
  p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff2))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.gRNA2.pdf", rdF),
         width = 6, height = 6)
  
  ggmat <- data.frame(dist = rdM[, "score"])
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = score.cutoff,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=score.cutoff, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("guide alignment score")
  p <- p + ggtitle(paste0("Cutoff score: ", score.cutoff))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", rdF),
         width = 6, height = 6)
  
  
  # plot flanking.cutoff (RANDOM)
  ggmat <- data.frame(dist = rdM.bf)
  p <- ggplot(ggmat) + 
    geom_density(aes(x=dist))
  p <- p + annotate("rect", xmin = flanking.cutoff,
                    xmax = Inf,
                    ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  p <- p + geom_vline(xintercept=flanking.cutoff, 
                      color = "black", size=2)
  p <- p + theme_bw(base_size = 16)
  p <- p + xlab("substring length (bp)")
  p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
  ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF),
         width = 6, height = 6)
  
  which_signif <- sapply(1:length(score.pv), function(x){
    mystr <- c()
    if(score.pv[x] < cutoff) mystr <- c(mystr, "CUMULATIVE") 
    if(score.pv1[x] < cutoff) mystr <- c(mystr, "gRNA1") 
    if(score.pv2[x] < cutoff) mystr <- c(mystr, "gRNA2")
    return(toString(mystr))  
  })
  
  #
  mygroups <- rep("NBS", nrow(realM))
  mygroups[realM.bf >= flanking.cutoff] <- "HMT"
  mygroups[score.pv < cutoff] <- "OMT"# QV or PV
  
  mygroups2 <- rep(NA, nrow(realM))
  #mygroups2[mygroups == "off.target" & (realM.bf > flanking.cutoff)] <- "yes"
  mygroups2[mygroups == "OMT" & (realM.bf >= flanking.cutoff)] <- "yes"
  mygroups2[mygroups == "HMT"] <- "yes"
  
  mygroups3 <- rep(NA, nrow(realM))
  ld.idx <- isLargeDel(realM, ots)
  #print(length(ld.idx))
  if(length(ld.idx)!=0){
    if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
  }
  
  realM$group <- mygroups

  realM$is.ON <- getONidx(realM, otsF)
  realM$which.signif <- which_signif
  
  realM$is.HMT <- mygroups2
  realM$is.large.del. <- mygroups3
  
  realM$OMT.pvalue <- score.pv 
  realM$OMT.adj.pvalue <- score.adj.pv 
  
  realM$OMT.pvalue.gRNA1 <- score.pv1 
  realM$OMT.adj.pvalue.gRNA1 <- score.adj.pv1
  
  realM$OMT.pvalue.gRNA2 <- score.pv2 
  realM$OMT.adj.pvalue.gRNA2 <- score.adj.pv2
  
  realM$HMT.pvalue <- flanking.pv
  realM$HMT.adj.pvalue <- flanking.adj.pv
  
  write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE, overwrite = TRUE)	
}

assignGroups_2OT <- function(realF, rdF, otsF, cutoff)
{
	ots <- read.delim(otsF, header = FALSE)

	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	realM.bf <- getBestFlanking_2OT(realM)
	rdM.bf <- getBestFlanking_2OT(rdM)
	
	score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
	#flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
	flanking.cutoff <- 25
	
	score.pv <- unlist(lapply(realM[, "score"], getEmpiricalPV, x = rdM[, "score"], type = "greater"))
	score.adj.pv <- p.adjust(score.pv, method = "BH")
	flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
	flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")
	#print(score.cutoff)
	#print(flanking.cutoff)	
	
	# plot score.cutoff (REAL)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (REAL)
	#xmax <- min(100, max(realM.bf))
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
	plot(density(realM.bf[realM.bf <= 100]), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot score.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_score_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM[, "score"]), main = paste0("Cutoff score: ", score.cutoff), xlab = "guide alignment score")
	abline(v=score.cutoff, col="black", lwd=3)
	dev.off()
	
	# plot flanking.cutoff (RANDOM)
	pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	plot(density(rdM.bf), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	abline(v=flanking.cutoff, col="black", lwd=3)
	dev.off()
	
	mygroups <- rep("NBS", nrow(realM))
	mygroups[realM.bf >= flanking.cutoff] <- "HMT"
	mygroups[score.adj.pv < cutoff] <- "OMT"

	mygroups2 <- rep(NA, nrow(realM))
	#mygroups2[mygroups == "off.target" & (realM.bf > flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "OMT" & (realM.bf >= flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "HMT"] <- "yes"
	
	mygroups3 <- rep(NA, nrow(realM))
	ld.idx <- isLargeDel(realM, ots)
	#print(length(ld.idx))
	if(length(ld.idx)!=0){
		if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
		}
	
	
	realM$group <- mygroups
	
	realM$is.ON <- getONidx(realM, otsF)
	
	realM$is.HMT <- mygroups2
	realM$is.large.del. <- mygroups3
	
	realM$OMT.pvalue <- score.pv 
	realM$OMT.adj.pvalue <- score.adj.pv 
	realM$HMT.pvalue <- flanking.pv
	realM$HMT.adj.pvalue <- flanking.adj.pv
	
	
	write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE, overwrite = TRUE)	
}

assignGroupsTALEN <- function(realF, rdF, otsF, cutoff)
{
	ots <- read.delim(otsF, header = FALSE)

	realM <- read.xlsx(realF, sheet = 1)
	rdM <- read.xlsx(rdF, sheet = 1)
	
	realM.bf <- getBestFlanking(realM)
	rdM.bf <- getBestFlanking(rdM)
	
	#score.cutoff <- sort(rdM[, "score"], decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# score of the top XX%
	#flanking.cutoff <- sort(rdM.bf, decreasing = TRUE)[floor(nrow(rdM) * cutoff)]# flanking length of the top XX%
	flanking.cutoff <- 25
	
	#print(score.cutoff)
	#print(flanking.cutoff)	
	
	bestCb <- realM$BestCB
	qv <- lapply(1:length(bestCb), function(i){
		if(is.na(bestCb[i])) return(1)
		qv.current <- realM[i, grep("_adj.pv$", colnames(realM))]
		qv.current <- as.numeric(qv.current[grep(bestCb[i], names(qv.current))])
		return(qv.current)
		})
	qv <- unlist(qv)
	
	
	flanking.pv <- unlist(lapply(realM.bf, getEmpiricalPV, x = rdM.bf, type = "greater"))
	flanking.adj.pv <- p.adjust(flanking.pv, method = "BH")

	#nbSignif <- realM[, grep("_adj.pv", colnames(realM))]
	#nbSignif <- (colSums(nbSignif) / ncol(nbSignif)) * 100# percentage of significant pvalue
	
	mygroups <- rep("NBS", nrow(realM))
	mygroups[realM.bf >= flanking.cutoff] <- "HMT"
	mygroups[qv < cutoff] <- "OMT"

	mygroups2 <- rep(NA, nrow(realM))
	mygroups2[mygroups == "OMT" & (realM.bf >= flanking.cutoff)] <- "yes"
	mygroups2[mygroups == "HMT"] <- "yes"
	
	mygroups3 <- rep(NA, nrow(realM))
	ld.idx <- isLargeDel(realM, ots)
	if(!is.na(ld.idx)) mygroups3[ld.idx] <- "yes"
	
	realM$group <- mygroups
	
	realM$is.ON <- getONidx(realM, otsF)
	
	realM$is.HMT <- mygroups2
	realM$is.large.del. <- mygroups3
	
	realM$HMT.pvalue <- flanking.pv
	realM$HMT.adj.pvalue <- flanking.adj.pv
	
	write.xlsx(realM, gsub(".xlsx", "_GROUP.xlsx", realF), row.names = FALSE, overwrite = TRUE)	
	
	if(length(realM.bf) > 10){
  	# flanking plot (REAL)
  	#pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF))
  	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", realF), units="px", width=1600, height=1600, res=300)
  	#plot(density(realM.bf[realM.bf <= 100]), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
  	#abline(v=flanking.cutoff, col="black", lwd=3)
  	#dev.off()
  	
  	ggmat <- data.frame(dist = realM.bf[realM.bf <= 100])
  	p <- ggplot(ggmat) + 
  	  geom_density(aes(x=dist))
  	p <- p + annotate("rect", xmin = flanking.cutoff,
  	                  xmax = Inf,
  	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
  	p <- p + geom_vline(xintercept=flanking.cutoff, 
  	                    color = "black", size=2)
  	p <- p + theme_bw(base_size = 16)
  	p <- p + xlab("substring length (bp)")
  	p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
  	ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", realF),
  	       width = 6, height = 6)
	}
	
	# plot flanking.cutoff (RANDOM)
	#pdf(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF))
	#png(gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.png", rdF), units="px", width=1600, height=1600, res=300)
	#plot(density(rdM.bf), main = paste0("Cutoff length: ", flanking.cutoff), xlab = "substring length (bp)")
	#abline(v=flanking.cutoff, col="black", lwd=3)
	#dev.off()
	
	ggmat <- data.frame(dist = rdM.bf)
	p <- ggplot(ggmat) + 
	  geom_density(aes(x=dist))
	p <- p + annotate("rect", xmin = flanking.cutoff,
	                  xmax = Inf,
	                  ymin = -Inf, ymax = Inf, fill = "blue3", alpha = .2)
	p <- p + geom_vline(xintercept=flanking.cutoff, 
	                    color = "black", size=2)
	p <- p + theme_bw(base_size = 16)
	p <- p + xlab("substring length (bp)")
	p <- p + ggtitle(paste0("Cutoff length: ", flanking.cutoff))
	ggsave(plot = p, filename = gsub("_aln_stat_FLANK.xlsx", "_flanking_cutoff.pdf", rdF),
	       width = 6, height = 6)

}

getBestFlanking <- function(m)
{
	m.max <- apply(m[, c("flank.length", "flank.rev.length")],1,max)
	return(m.max)
}

getBestFlanking_2OT <- function(m)
{
	m.max <- apply(m[, c("flank.length", "flank.rev.length", "flank2.length", "flank2.rev.length")],1,max)
	return(m.max)
}

isLargeDel <- function(m, ots)
{
	m.sub <- m[, c("chromosome", "start", "end")]
	m.sub$strand <- "*"
	
	ots.sub <- data.frame(chr = as.character(ots[,1]),
						  start = as.numeric(ots[,2]),
						  end = as.numeric(ots[,3]),
						  strand = "*")
	
	
	# CONVERT TO GRANGES
	m.gr <- makeGRangesFromDataFrame(m.sub, seqnames.field = "chromosome",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = FALSE)
		
	ots.gr <- makeGRangesFromDataFrame(ots.sub, seqnames.field = "chr",
		start.field = "start", end.field = "end", strand.field = "strand",
		keep.extra.columns = FALSE)
	
	# INTERSECT
	gr.ovl <- findOverlaps(query = m.gr, subject = ots.gr, type = "any", ignore.strand= TRUE)	
	ovl.idx <- unique(queryHits(gr.ovl))
	
	return(ovl.idx)
}


groupSummary <- function(inputF, outputF, hits = NULL, score = NULL, pv = NULL)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	if(!is.null(hits)) readMat <- readMat[readMat$hits > hits, ]
	if(!is.null(score)) readMat <- readMat[readMat$score > score, ]
	if(!is.null(pv)) readMat <- readMat[readMat$adj.pvalue < pv, ]
	
	if(nrow(readMat)>0){
		write.xlsx(data.frame(table(readMat$group)), outputF, row.names = FALSE, overwrite = TRUE)
		}
}

hitsBarplot <- function(inputF, pv = NULL, top = NULL, showNBS = TRUE, log = TRUE, outName = NULL){
  realM <- read.xlsx(inputF, sheet = 1)
  if(!is.null(pv)) realM <- realM[realM$adj.pvalue < pv, ]
  if(!showNBS) realM <- realM[realM$group != "NBS", ]
  if(!is.null(top)){
    realM <- realM[order(-realM$hits), ]
    if(nrow(realM) > top) realM <- realM[1:top, ]
  }
  
  
  ggmat <- realM[, c("hits", "group")]
  ggmat$group[realM$group == "OMT" & realM$is.HMT == "yes"] <- "OMT/HMT"
  ggmat$group[1] <- "ON"
  ggmat$group <- factor(ggmat$group, levels = c("ON", "OMT", "HMT", "OMT/HMT", "NBS"))
  
  ggmat$Name <- paste(realM$chromosome, realM$start, realM$hits, sep = ":")
  ggmat$Name <- factor(ggmat$Name, levels = ggmat$Name[order(-ggmat$hits)])
  ggmat <- ggmat[order(ggmat$Name), ]
  
  p <- ggplot(data=ggmat, aes(x=Name, y=hits, fill = group)) +
    geom_bar(stat="identity", width = 0.75)
  p <- p + geom_text(aes(label=hits, color = group),  vjust=0.5, hjust = -0.25, size=3.5, angle = 90, fontface = "bold")
  #p <- p + geom_text(aes(label=hits, color = group),  vjust=-0.5, hjust = -0.25, size=3, angle = 45, fontface = "bold")
  #p <- p + ylim(c(0, max(ggmat$hits) + (max(ggmat$hits)*0.5)))
  p <- p + theme_bw(base_size = 18) + theme(panel.grid.minor = element_blank())

  p <- p + scale_fill_manual(values=c(ON = "green3",
                                      OMT = "red3",
                                      HMT = "blue3",
                                      "OMT/HMT" = "goldenrod1",
                                      NBS = "grey"))
  p <- p + scale_color_manual(values=c(ON = "green3",
                                      OMT = "red3",
                                      HMT = "blue3",
                                      "OMT/HMT" = "goldenrod1",
                                      NBS = "grey"))
  p <- p + xlab("") + ylab("Hits")
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(log){
    p <- p + scale_y_log10(limits = c(1, max(ggmat$hits) + (max(ggmat$hits)*0.5)))
    p <- p + annotation_logticks(sides = "l", size = 0.25) 
  }else{
    p <- p + ylim(c(0, (max(ggmat$hits) + max(ggmat$hits*0.1))))
  }

  
  mywidth <- nrow(ggmat) * 7.5 /10
  mywidth <- max(c(mywidth, 5))
  
  if(is.null(outName)){
    pdfName <- gsub("_FINAL.xlsx", "_hits_barplot.pdf", inputF)
  }else{
    pdfName <- outName
  }
  ggsave(plot = p, filename = pdfName,
         width = mywidth, height = 7)
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE)
{
  library(openxlsx)
  library(ggplot2)
################
# LOAD REAL DATA

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment/G3_W250"))
readMat.real <- read.xlsx("G3_W250_aln_stat_FLANK.xlsx", sheet = 1)

##################
# LOAD RANDOM DATA

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment/random_hg38_w250"))
readMat.rd <- read.xlsx("random_hg38_w250_aln_stat_FLANK.xlsx", sheet = 1)


#############################################################################
# ASSIGN GROUPS: off targets, homologous recombination, common breaking sites

cutoff <- 0.05
readMat.real <- assignGroups(readMat.real, readMat.rd, cutoff)


# SAVE
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/alignment/G3_W250"))
write.xlsx(readMat.real, "G3_W250_aln_stat_FLANK_GROUP.xlsx")


###########################
# TEST hitsBarplot FUNCTION


inputF <- file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/results/guide_aln/EMD101-sample3-1_w250_FINAL.xlsx")
pv <- 0.05
top <- 10
log <- TRUE

setwd(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221"))

inputList <- list.files(pattern = "_FINAL.xlsx", recursive = TRUE, full.names = TRUE)
lapply(inputList, function(i){
  print(i)
  hitsBarplot(i, pv = 0.05, top = 50, showNBS = TRUE, log = TRUE,
              outName = gsub("_FINAL.xlsx", "_hits_barplot.pdf", i))
  hitsBarplot(i, pv = 0.05, top = 50, showNBS = FALSE, log = TRUE,
              outName = gsub("_FINAL.xlsx", "_hits_barplot_woNBS.pdf", i))
})

# EMENDO OVERLAP
setwd(file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/"))

inputList <- list.files(pattern = "_FINAL.xlsx", recursive = TRUE, full.names = TRUE)
lapply(inputList, function(i){
  print(i)
  hitsBarplot(i, pv = NULL, top = NULL, showNBS = TRUE, log = TRUE,
              outName = gsub("_FINAL.xlsx", "_hits_barplot.pdf", i))
  hitsBarplot(i, pv = NULL, top = NULL, showNBS = FALSE, log = TRUE,
              outName = gsub("_FINAL.xlsx", "_hits_barplot_woNBS.pdf", i))
})

}




