library(ggplot2)
library(openxlsx)
library(ggridges)


sample.current <- "EMD3"

#setwd(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-576086333/EMD/", sample.current))
setwd(file.path("~/tmp/EMD_cluster/EMD/", sample.current))

binFiles <- list.files(pattern = "_5kb_bins.xlsx", recursive = TRUE, full.names = TRUE)
names(binFiles) <- gsub("_5kb_bins.xlsx", "", basename(binFiles))
names(binFiles) <- gsub("EMD-sample", "", names(binFiles))

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

mylevels <- unlist(lapply(c(9:10), function(i) paste0(i , c("-1", "-2"))))
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
myheight <-myheight +  0.5*length(binList)
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





