




barcodeHoppingFilter <- function(inputFile, outputFile, hardCoef = 2.5, softCoef = 2){
  if(softCoef > hardCoef) stop("barcodeHoppingFilter: softCoef cannot be higher than hardCoef")
  
  cast <- read.xlsx(inputFile, sheet = 1)
  
  if(nrow(cast) < 10){
    print("less than 10 sites, skip barcode hoping filter")
    cast$BarcodeHoping <- NA
    write.xlsx(cast, gsub("_RAW.xlsx$", ".xlsx", inputFile))
  }
  
  log10RATIO <- log10(cast$read / cast$hits)
  
  Q <- quantile(log10RATIO, probs=c(.25, .75), na.rm = FALSE)
  iqr <- IQR(log10RATIO)
  
  # hard Tukey coef
  hard.up <-  Q[2]+hardCoef*iqr # Upper Range  
  hard.low<- Q[1]-hardCoef*iqr # Lower Range
  
  soft.up <-  Q[2]+softCoef*iqr # Upper Range  
  soft.low<- Q[1]-softCoef*iqr # Lower Range
  
  cast$log10read_hits <- log10RATIO
  cast$BarcodeHopping <- "no"
  #cast$BarcodeHopping[log10RATIO < soft.low | log10RATIO > soft.up] <- "likely"
  #cast$BarcodeHopping[log10RATIO < hard.low | log10RATIO > hard.up] <- "yes"
  
  cast$BarcodeHopping[log10RATIO < soft.low] <- "likely"
  cast$BarcodeHopping[log10RATIO < hard.low] <- "yes"
  
  #cast <- cast[cast$BarcodeHopping != "yes", ]
  
  write.xlsx(cast, outputFile, overwrite = TRUE)
}











if(FALSE){
  library(openxlsx)
  library(ggplot2)
  library(ggstatsplot)
  
  
  # LOAD FINAL FILES
  setwd("~/Research/CASTSeq/report/EMD_Jan22/wetransfer_emd_ovl_new-zip_2022-01-17_0800/EMD_OVL_new/EMD9_dual/")
  castFiles <- list.files(pattern = "OVL1_FINAL.xlsx", recursive = TRUE)
  names(castFiles) <- dirname(castFiles)
  names(castFiles) <- gsub("MEG01_EMD-", "", names(castFiles))
  names(castFiles) <- gsub("_OVL1", "", names(castFiles))
  
  castList <- lapply(castFiles, read.xlsx)
  
  # REMOVE ON-TARGET
  #castList <- lapply(castList, function(i) i[i$is.ON != "yes" | is.na(i$is.ON), ])
  
  # DEFINE OUTLIERS
  tukeyCoef <- 2.5
  
  castList <- lapply(castList, function(cast){
    
    log10RATIO <- log10(cast$read / cast$hits)
    
    Q <- quantile(log10RATIO, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(log10RATIO)
    
    up <-  Q[2]+tukeyCoef*iqr # Upper Range  
    low<- Q[1]-tukeyCoef*iqr # Lower Range
    cast$BarcodeHoping <- "yes"
    cast$BarcodeHoping[log10RATIO > low & log10RATIO < up] <- "no"
    
    return(cast)
  })
  
  
  # PLOT
  ggList <- lapply(1:length(castList), function(i){
    cast <- castList[[i]]
    siteName <- paste(cast$chromosome, cast$start, cast$SYMBOL, sep = "_")
    rh <- cast$read / cast$hits
    rh.UT <- cast$read.ctl / cast$hits.ctl
    
    ggmat <- data.frame(SITE = siteName,
                        READS = cast$read,
                        HITS = cast$hits,
                        RATIO = rh,
                        READS.UT = cast$read.ctl,
                        HITS.UT = cast$hits.ctl,
                        RATIO.UT = rh.UT,
                        SAMPLE = names(castList)[i])
    ggmat$LOG10READS <- log10(ggmat$READS)
    ggmat$lOG10HITS <- log10(ggmat$HITS)
    ggmat$LOG10RATIO <- log10(ggmat$RATIO)
    
    ggmat$OUTLIER <- cast$BarcodeHoping
    return(ggmat)
  })
  ggmat <- do.call(rbind, ggList)
  
  # REMOVE LOW HITS
  #ggmat <- ggmat[ggmat$HITS >= 3, ]
  
  set.seed(123456)
  p <- ggbetweenstats(ggmat,
                      SAMPLE, LOG10RATIO,
                      outlier.tagging = TRUE, outlier.label = SITE,
                      outlier.coef = tukeyCoef,
                      outlier.label.args = c(max.overlaps = 20, cex = 3),
                      pairwise.comparisons = FALSE, results.subtitle = FALSE)
  
  setwd("~/Research/CASTSeq/report/EMD_Jan22")
  ggsave(plot = p, filename = paste0("EMD9_dual_coef", tukeyCoef,"_barcode_hoping_boxplot.pdf"),
         width = 20, height = 5)
  
  
  
  
  p <- ggplot(ggmat, aes(x=SAMPLE, y=RATIO))
  p <- p + geom_boxplot()
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave( plot = p, filename = "~/tmp/test.pdf",
          width = 8, height = 5)
  
  ggmat$LOGRATIO <- log10(ggmat$RATIO)
  
  
  
  p<-ggplot(ggmat, aes(x=log10(RATIO), fill=SAMPLE)) +
    geom_density(alpha=0.4)
  p <- p + theme_bw()
  
  p<-ggplot(ggmat, aes(x=READS, fill=SAMPLE)) +
    geom_density(alpha=0.4)
  p <- p + theme_bw()
  
  p<-ggplot(ggmat, aes(x=HITS, fill=SAMPLE)) +
    geom_density(alpha=0.4)
  p <- p + theme_bw()
  
  p<-ggplot(ggmat, aes(x=LOGREADS, fill=SAMPLE)) +
    geom_density(alpha=0.4)
  p <- p + theme_bw()
  
  
  #######################################################
  # KNN
  ggmat$LOGREADS <- log10(ggmat$READS)
  ggmat$LOGHITS <- log10(ggmat$HITS)
  
  p <- ggplot(ggmat, aes(x=LOGREADS, y =HITS, color=SAMPLE))
  p <- p + geom_point(size = 4)
  
  library(dbscan)
  
  WKNN_Outlier <- apply(kNNdist(x=ggmat[, c("LOGREADS", "HITS")], k = 2, all = T), 1, mean)  # Weighted KNN outlier score (mean)
  
  
  rank_WKNN_Outlier <- order(x=WKNN_Outlier, decreasing = TRUE) 
  WKNN_Result <- data.frame(ID = rank_WKNN_Outlier, score = WKNN_Outlier[rank_WKNN_Outlier])
  
  head(WKNN_Result, top_n)
  
  p <- ggplot() + geom_point(data=ggmat, mapping=aes(x=LOGREADS, y=HITS), shape = 19)
  
  p <- p +
    geom_point(data=ggmat[rank_WKNN_Outlier[1:top_n],], mapping=aes(x=LOGREADS,y=HITS), shape=19, 
               color="red", size=2) +
    geom_text(data=ggmat[rank_WKNN_Outlier[1:top_n],],
              mapping=aes(x=(LOGREADS-0.5), y=HITS, label=rank_WKNN_Outlier[1:top_n]), size=2.5)
  
  
  
}










