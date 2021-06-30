
library(GenomicRanges)
library(openxlsx)
library(rtracklayer)
library(ChIPpeakAnno)
library(UpSetR)

#source("~/Scripts/work/tools/venn_script.r")

hg19to38 <- function(m)
{
  chain <- import.chain("~/Research/database/ucsc/hg19ToHg38.over.chain")
  
  m$strand <- "*"
  m.gr <- makeGRangesFromDataFrame(m, seqnames.field = "V1",
                                   start.field = "V2", end.field = "V3", strand.field = "strand",
                                   keep.extra.columns = TRUE)
  
  m38.gr = liftOver(m.gr, chain)
  
  m38.gr = unlist(m38.gr)
  genome(m38.gr) = "hg38"
  return(data.frame(m38.gr))
}

getReads <- function(dfA, dfB, width)
{
  reads <- rep(0, nrow(dfA))
  
  grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
  
  idxA <- queryHits(gr.ovl)
  if(length(idxA) == 0) return(NA)
  
  idxA.unique <- unique(idxA)
  idxB <- subjectHits(gr.ovl)
  
  reads[idxA.unique] <- sapply(idxA.unique, function(i){
    currentIndex <- idxB[idxA == i]
    return(sum(dfB$read[currentIndex]))
  })
  
  return(reads)
}


getHits <- function(dfA, dfB, width)
{
  hits <- rep(0, nrow(dfA))
  
  grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
  
  idxA <- queryHits(gr.ovl)
  if(length(idxA) == 0) return(NA)
  
  idxA.unique <- unique(idxA)
  idxB <- subjectHits(gr.ovl)
  
  hits[idxA.unique] <- sapply(idxA.unique, function(i){
    currentIndex <- idxB[idxA == i]
    return(sum(dfB$hits[currentIndex]))
  })
  
  return(hits)
}

getPVs <- function(dfA, dfB, width)
{
  pvs <- rep(0, nrow(dfA))
  
  grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
  
  idxA <- queryHits(gr.ovl)
  if(length(idxA) == 0) return(NA)
  
  idxA.unique <- unique(idxA)
  idxB <- subjectHits(gr.ovl)
  
  pvs[idxA.unique] <- sapply(idxA.unique, function(i){
    currentIndex <- idxB[idxA == i]
    return(min(dfB$pvalue[currentIndex]))
  })
  
  return(pvs)
}

getQVs <- function(dfA, dfB, width)
{
  qvs <- rep(0, nrow(dfA))
  
  grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
  
  idxA <- queryHits(gr.ovl)
  if(length(idxA) == 0) return(NA)
  
  idxA.unique <- unique(idxA)
  idxB <- subjectHits(gr.ovl)
  
  qvs[idxA.unique] <- sapply(idxA.unique, function(i){
    currentIndex <- idxB[idxA == i]
    return(min(dfB$adj.pvalue[currentIndex]))
  })
  
  return(qvs)
}

getCommon <- function(dfA, dfB, width, nb.signif = 1)
{
  grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns = TRUE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = width)
  if(length(queryHits(gr.ovl)) == 0) return(NA)
  
  df.ovl <- data.frame(dfA[queryHits(gr.ovl),], dfB[subjectHits(gr.ovl),])
  
  pv.nameA <- grep("^adj.pvalue", colnames(df.ovl))[1]
  pv.nameB <- grep("^adj.pvalue", colnames(df.ovl))[2]
  if(nb.signif == 1){
    toKeep <- df.ovl[, pv.nameA] < 0.05 | df.ovl[, pv.nameB] < 0.05# AND= 2 significnat; OR=1significant
  }else if(nb.signif == 2){
    toKeep <- df.ovl[, pv.nameA] < 0.05 & df.ovl[, pv.nameB] < 0.05# AND= 2 significnat; OR=1significant
  }else{
    print("WRONG nb.signif argument. Can only be 1 or 2")
    toKeep <- rep(FALSE, nrow(df.ovl))
  }
  df.ovl.filt <- df.ovl[toKeep, ]
  
  df.ovl.filt1 <- df.ovl.filt[, 1:ncol(dfA)]
  colnames(df.ovl.filt1) <- colnames(dfA)
  df.ovl.filt2 <- df.ovl.filt[, (ncol(dfA)+1):ncol(df.ovl.filt)]
  colnames(df.ovl.filt2) <- colnames(dfB)
  
  df.ovl.filt <- rbind(df.ovl.filt1, df.ovl.filt2)
  df.ovl.filt <- df.ovl.filt[order(df.ovl.filt$chromosome, df.ovl.filt$start, df.ovl.filt$end), ]
  
  df.ovl.filt <- df.ovl.filt[!duplicated(df.ovl.filt), ]# remove duplicates
  
  return(df.ovl.filt)
}


getCommon_fromList <- function(dfList, width, nb.signif)
{
  combMat <- combn(1:length(dfList), 2)
  df.ovl.filt <- lapply(1:ncol(combMat), function(i)
    getCommon(dfList[[combMat[1,i]]], dfList[[combMat[2,i]]], width, nb.signif)
  )
  df.ovl.filt <- do.call(rbind, df.ovl.filt)	
  df.ovl.filt <- df.ovl.filt[order(df.ovl.filt$chromosome, df.ovl.filt$start, df.ovl.filt$end), ]
  df.ovl.filt <- df.ovl.filt[!duplicated(df.ovl.filt), ]# remove duplicates	
  
  return(df.ovl.filt)
}



dfComparisonList <- function(fList, nList, width = 2000, nb.signif = 1, NBS = TRUE)
{
  dfList <- lapply(fList, read.xlsx, sheet = 1)
  dfList <- lapply(dfList, function(i) i[,-4])
  #grList <- lapply(dfList, function(i) makeGRangesFromDataFrame(i, seqnames.field = "chromosome",
  #	start.field = "start", end.field = "end",
  #	keep.extra.columns = TRUE, ignore.strand = TRUE))
  
  # Remove NBS
  if(!NBS) dfList <- lapply(dfList, function(i) i[i$group != "NBS", ])
  
  # Overlaping sites	
  df.ovl.filt.raw <- getCommon_fromList(dfList, width, nb.signif)	
  
  # merge overlap (+/- 2000bp)
  df.ovl.filt <- df.ovl.filt.raw[, 1:7]
  
  idx <- 0
  df.final <- c()
  for(chr in unique(df.ovl.filt[,1]))
  {
    tempdf <- subset(df.ovl.filt,  chromosome==chr)
    df.final <- rbind(df.final, tempdf[1,])
    idx <- idx + 1
    if(nrow(tempdf) > 1){
      for(i in 2:nrow(tempdf))
      {
        if(tempdf[i, "start"] <= (df.final[idx, "end"] + width))
        {
          df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] <-
            df.final[idx, c("read", "hits", "read.ctl", "hits.ctl")] + tempdf[i, c("read", "hits", "read.ctl", "hits.ctl")]# add read, hits...
          if(tempdf[i, "end"] > df.final[idx, "end"]) df.final[idx, "end"] <- tempdf[i, "end"]# change end if new end is higher
          
        }else{
          df.final <- rbind(df.final, tempdf[i, ])
          idx <- idx + 1
        }
      }
    }
  }
  gr.final <- makeGRangesFromDataFrame(df.final, seqnames.field = "chromosome",
                                       start.field = "start", end.field = "end",
                                       keep.extra.columns = TRUE, ignore.strand = TRUE)	
  
  # Specific sites
  speList <- lapply(dfList, function(i){
    gr <- makeGRangesFromDataFrame(i, seqnames.field = "chromosome",
                                   start.field = "start", end.field = "end",
                                   keep.extra.columns = TRUE, ignore.strand = TRUE)
    gr.ovl <- findOverlaps(query = gr, subject = gr.final, type = "any", maxgap = width)
    df.spe <- i[-queryHits(gr.ovl),]
    return(df.spe)
  })
  names(speList) <- nList	
  
  # Intersect common sites vs. individual sites
  readList <- lapply(dfList, getReads, dfA = df.final, width = 0)
  readM <- do.call(cbind, readList)
  colnames(readM) <- paste0(nList, "_read")
  
  hitList <- lapply(dfList, getHits, dfA = df.final, width = 0)
  hitM <- do.call(cbind, hitList)
  colnames(hitM) <- paste0(nList, "_hits")
  
  pvList <- lapply(dfList, getPVs, dfA = df.final, width = 0)
  pvM <- do.call(cbind, pvList)
  colnames(pvM) <- paste0(nList, "_pvalue")
  
  qvList <- lapply(dfList, getQVs, dfA = df.final, width = 0)
  qvM <- do.call(cbind, qvList)
  colnames(qvM) <- paste0(nList, "_adj.pvalue")
  
  df.final <- cbind(df.final, readM, hitM, pvM, qvM)	
  
  # Set the final number of "read" and "hits"	
  df.final$read <- rowSums(df.final[, grep("_read", colnames(df.final))])	
  df.final$hits <- rowSums(df.final[, grep("_hits", colnames(df.final))])	
  
  # Re-order rows
  df.final <- df.final[order(-df.final$hits, -df.final$read), ]
  
  # return
  compL <- c(list(common.filtered = df.final), list(common = df.ovl.filt.raw), speList)
  return(compL)
  
}


makeUpset <- function(ovlFile){
  mysheets <- getSheetNames(ovlFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = ovlFile)
  names(mList) <- mysheets
  
  ovlMat <- mList[["common.filtered"]]
  
  # get nb of reads
  nbMat <- ovlMat[, grepl("*_read$", colnames(ovlMat))]
  
  siteNames <- paste(ovlMat$chromosome, ovlMat$start, ovlMat$end, sep = "_")
  sampleNames <- gsub("_read$", "", colnames(nbMat))
  rownames(nbMat) <- siteNames
  colnames(nbMat) <- sampleNames
  
  # GET siteList from nbMat
  siteList <- lapply(1:ncol(nbMat), function(i){
    return(rownames(nbMat)[as.numeric(nbMat[, i]) != 0])	
  })
  names(siteList) <- colnames(nbMat)
  
  # ADD specific sites
  for(i in names(siteList)){
    speMat <- mList[[i]]
    speMat <- speMat[speMat$adj.pvalue < 0.05, ]
    siteList[[i]] <- c(siteList[[i]], paste(speMat$chromosome, speMat$start, speMat$end, sep = "_"))
  }
  
  # UpSet plot
  
  pdf(gsub(".xlsx", "_UpSetR.pdf", ovlFile), width = 10, onefile = FALSE)
  print(upset(fromList(siteList), order.by = "freq", nsets = length(siteList),
              sets = names(siteList), keep.order = TRUE,
              text.scale = 2,
              mainbar.y.label = "Intersection", sets.x.label = "sites per sample"))
  dev.off()	
  
  
  # Save Venn
  #vennOut <- doVenn(siteList)
}


if(FALSE){
##################################################################
# DO NOT RUN


# MASAKO MK FILES 18.02.21
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/MK_Kapa/results/guide_aln/MK-T-Kapa_w250_aln_stat_FLANK_GROUP_GENES.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/MK_Rep2/results/guide_aln/MK-T2-Rep2_w250_aln_stat_FLANK_GROUP_GENES.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])

ovlDir <- file.path("~/Research/CASTSeq/revision/overlap/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "MK_Kapa_Rep2.xlsx")

makeUpset("MK_Kapa_Rep2.xlsx")

}
