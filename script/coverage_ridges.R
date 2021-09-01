library(ggplot2)
library(openxlsx)
library(ggridges)



setwd(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-556214738/EMD/EMD4/"))
binFiles <- list.files(pattern = "_5kb_bins.xlsx", recursive = TRUE, full.names = TRUE)
names(binFiles) <- gsub("_5kb_bins.xlsx", "", basename(binFiles))

binList <- lapply(binFiles, read.xlsx, sheet = 1)

binList <- binList[1:4]

ggList <- lapply(1:length(binList), function(j){
  binMat <- binList[[j]]
  
  #binMat$POS <- ceiling(binMat$POS / 100)
  #binMat$NEG <- ceiling(binMat$NEG / 100)
  #binMat$TOT <- ceiling(binMat$TOT / 100)
  
  count.pos <- lapply(1:nrow(binMat), function(i){
    if(binMat$POS[i] == 0) return(NA)
    return(data.frame(start = binMat$start[i], strand = rep("POS", binMat$POS[i])))
  })
  count.pos <- count.pos[!is.na(count.pos)]
  count.pos <- do.call(rbind, count.pos)
  
  count.neg <- lapply(1:nrow(binMat), function(i){
    if(binMat$NEG[i] == 0) return(NA)
    return(data.frame(start = binMat$start[i], strand = rep("NEG", binMat$NEG[i])))
  })
  count.neg <- count.neg[!is.na(count.neg)]
  count.neg <- do.call(rbind, count.neg)
  
  count.tot <- lapply(1:nrow(binMat), function(i){
    if(binMat$TOT[i] == 0) return(NA)
    return(data.frame(start = binMat$start[i], strand = rep("TOT", binMat$TOT[i])))
  })
  count.tot <- count.tot[!is.na(count.tot)]
  count.tot <- do.call(rbind, count.tot)
  
  ggmat <- do.call(rbind, list(count.pos, count.neg, count.tot))
  ggmat$SAMPLE = names(binList)[j]
  
  #
  ggmat <- do.call(rbind, list(ggmat[sample(which(ggmat$strand == "POS"), 1000), ],
                               ggmat[sample(which(ggmat$strand == "NEG"), 1000), ],
                               ggmat[sample(which(ggmat$strand == "TOT"), 1000), ]))
  
  
  return(ggmat)
})


ggmat <- do.call(rbind, ggList)

p <- ggplot(ggmat[ggmat$strand != "TOT", ], aes(x = start, y = SAMPLE, fill = strand))
p <- p + geom_density_ridges(scale = 1, bandwidth = 100, alpha = 0.5)
p <- p + theme_bw(base_size = 14)
p <- p + xlab("") + ylab("")

ggsave(plot = p, filename = "~/tmp/test.pdf", height = 15, width = 10)

ggplot(count.tot, aes(start)) +
  geom_density(adjust = 5)