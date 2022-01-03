library(ggplot2)
library(openxlsx)
library(ggridges)


sample.current <- "Cas9"

#setwd(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-576086333/EMD/", sample.current))
#setwd(file.path("~/tmp/EMD_cluster/EMD/", sample.current))
setwd(file.path("~/Research/CASTSeq/test/coverage_ridge/", sample.current))


binFiles <- list.files(pattern = "_10kb_bins.xlsx", recursive = TRUE, full.names = TRUE)
names(binFiles) <- gsub("_10kb_bins.xlsx", "", basename(binFiles))
#names(binFiles) <- gsub("EMD-sample", "", names(binFiles))

binList <- lapply(binFiles, read.xlsx, sheet = 1)

ggList <- lapply(1:length(binList), function(j){
  binMat <- binList[[j]]
  binMat.tot <- sum(binMat$TOT)
  
  count.pos <- binMat$POS / binMat.tot 
  count.neg <- binMat$NEG / binMat.tot 

  ggmat <- data.frame(start = rep(binMat$start, 2),
                      count = c(count.pos, count.neg),
                      strand = c(rep("POS", length(count.pos)), rep("NEG", length(count.neg))),
                      event = c(rep("DEL", length(count.pos)), rep("INV", length(count.neg))),
                      SAMPLE = names(binList)[j]
                      )
  
  return(ggmat)
})

ggmat <- do.call(rbind, ggList)

#mylevels <- unlist(lapply(c(9:10), function(i) paste0(i , c("-1", "-2"))))
mylevels <- c("Cas9-KO-G3", "DF-KO-G3")
ggmat$SAMPLE <- factor(ggmat$SAMPLE, levels = mylevels)
ggmat <- ggmat[!is.na(ggmat$SAMPLE), ]

# STRAND COLOR
p <- ggplot(ggmat, aes(x = start, y = count, fill = strand, color = strand))
p <- p + geom_area(alpha = 0.5, position = "identity")
p <- p + theme_bw(base_size = 14)
p <- p + xlab("") + ylab("")
p <- p + facet_grid(rows = vars(ggmat$SAMPLE))
p <- p + theme(strip.text.y = element_text(size=10, face = "bold", color="black"))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + scale_fill_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                    NEG = rgb(175, 174, 239, maxColorValue = 255)))
p <- p + scale_color_manual(values=c(POS = rgb(245, 171, 173, maxColorValue = 255),
                                     NEG = rgb(175, 174, 239, maxColorValue = 255)))

myheight <- 2
myheight <-myheight +  0.75*length(binList)
ggsave(plot = p, filename = paste0("~/tmp/",sample.current ,"_strand.pdf"), width = 5, height = myheight)

# EVENT COLOR
p <- ggplot(ggmat, aes(x = start, y = count, fill = event, color = event))
p <- p + geom_area(alpha = 0.5, position = "identity")
p <- p + theme_bw(base_size = 14)
p <- p + xlab("") + ylab("")
p <- p + facet_grid(rows = vars(ggmat$SAMPLE))
p <- p + theme(strip.text.y = element_text(size=10, face = "bold", color="black"))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + scale_color_manual(values=c(DEL = "orange2",
                                     INV = "orchid2"))
p <- p + scale_fill_manual(values=c(DEL = "orange2",
                                    INV = "orchid2"))

ggsave(plot = p, filename = paste0("~/tmp/",sample.current ,"_event.pdf"), width = 5, height = myheight)



################################################################
# V2

ggList <- lapply(1:length(binList), function(j){
  binMat <- binList[[j]]
  ggmat <- data.frame(start = binMat$LOC,
                      count = log2(binMat$NORM),
                      SAMPLE = names(binList)[j]
  )
  
  return(ggmat)
})

ggmat <- do.call(rbind, ggList)

#mylevels <- unlist(lapply(c(9:10), function(i) paste0(i , c("-1", "-2"))))
mylevels <- c("Cas9-KO-G3", "DF-KO-G3")
ggmat$SAMPLE <- factor(ggmat$SAMPLE, levels = mylevels)
ggmat <- ggmat[!is.na(ggmat$SAMPLE), ]


p<-ggplot(ggmat, aes(x=start, y=count, group=SAMPLE)) +
  geom_line(aes(color=SAMPLE))+
  geom_point(aes(color=SAMPLE))
p <- p + theme_bw(base_size = 16)
p <- p + xlab("Distance from ON-target (bp)") + ylab("Normalized read count (log2 CPM)")
p <- p + scale_color_manual(values=c("Cas9-KO-G3" = "black",
                                     "DF-KO-G3" = rgb(191, 181, 123, maxColorValue = 255)))
p <- p + xlim(-1200, 1200)
ggsave(plot = p, filename = paste0(sample.current, "_log2_CPM.pdf"),
       width = 8, height = 4)


################################################################
# V3

ggList <- lapply(1:length(binList), function(j){
  binMat <- binList[[j]]
  binMat.tot <- sum(binMat$TOT)
  
  count.norm <- binMat$TOT / binMat.tot 
  
  ggmat <- data.frame(start = binMat$LOC,
                      count = count.norm,
                      SAMPLE = names(binList)[j]
  )
  
  return(ggmat)
})

ggmat <- do.call(rbind, ggList)

mylevels <- c("Cas9-KO-G3", "DF-KO-G3")
ggmat$SAMPLE <- factor(ggmat$SAMPLE, levels = mylevels)
ggmat <- ggmat[!is.na(ggmat$SAMPLE), ]


p<-ggplot(ggmat, aes(x=start, y=count, group=SAMPLE)) +
  geom_line(aes(color=SAMPLE), size = 1.25, alpha = 0.75)+
  geom_point(aes(color=SAMPLE), size = 2.5, alpha = 0.75)
p <- p + theme_bw(base_size = 16)
p <- p + xlab("Distance from ON-target (bp)") + ylab("Normalized read count")
p <- p + scale_color_manual(values=c("Cas9-KO-G3" = rgb(212, 212, 212, maxColorValue = 255),
                                     "DF-KO-G3" = rgb(191, 181, 123, maxColorValue = 255)))
p <- p + xlim(-1200, 1200)
ggsave(plot = p, filename = paste0(sample.current, "_Norm.pdf"),
       width = 8, height = 4)

write.xlsx(ggmat, file = paste0(sample.current, "_Norm.xlsx"))

################################################
# V3
ggList <- lapply(1:length(binList), function(j){
  binMat <- binList[[j]]
  ggmat <- data.frame(start = binMat$LOC,
                      count = log2(binMat$NORM),
                      SAMPLE = names(binList)[j]
  )
  
  return(ggmat)
})

ggmat <- do.call(rbind, ggList)

mylevels <- c("Cas9-KO-G3", "DF-KO-G3")
ggmat$SAMPLE <- factor(ggmat$SAMPLE, levels = mylevels)
ggmat <- ggmat[!is.na(ggmat$SAMPLE), ]


p<-ggplot(ggmat, aes(x=start, y=count, group=SAMPLE)) +
  geom_line(aes(color=SAMPLE), size = 1.25, alpha = 0.75)+
  geom_point(aes(color=SAMPLE), size = 2.5, alpha = 0.75)
p <- p + theme_bw(base_size = 16)
p <- p + xlab("Distance from ON-target (bp)") + ylab("Normalized read count (log2 CPM)")
p <- p + scale_color_manual(values=c("Cas9-KO-G3" = rgb(212, 212, 212, maxColorValue = 255),
                                     "DF-KO-G3" = rgb(191, 181, 123, maxColorValue = 255)))
p <- p + xlim(-1200, 1200)
ggsave(plot = p, filename = paste0(sample.current, "_log2_CPM.pdf"),
       width = 8, height = 4)



######################################
# BARPLOT



