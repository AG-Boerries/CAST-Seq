#library(tools)
#library(Biostrings)
#library(openxlsx)
#library(data.table)
#library(ggplot2)
#library(ggseqlogo)
#library(parallel)

#source("/home/gandri/offTargets/Giando/pipeline/script/bed2sequence.R")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

#toNum <- function(x) as.numeric(levels(x))[x]

getRevComp <- function(x) as.character(reverseComplement(DNAString(x)))

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

getBestAln <- function(alnList)
{
  scoreList <- unlist(lapply(alnList, score))
  return(alnList[[which.max(scoreList)]])
}

getBestAlnIdx <- function(alnList)
{
  scoreList <- unlist(lapply(alnList, score))
  return(which.max(scoreList))
}

getStrAlignmentOLD <- function(aln.file)
{
  #aln.str <- read_document(aln.file)
  aln.str <- read.delim(aln.file, sep = "\t", header = FALSE)
  aln.str <- as.character(t(aln.str))
  
  # get S1 index
  s1.idx <- grep("^S1", aln.str)
  
  s1.str <- lapply(s1.idx, function(i){
    ref.str <- aln.str[i]
    ref.str <- strsplit(ref.str, split = " ")
    ref.str <- unlist(ref.str)
    ref.str <- ref.str[ref.str != ""]
    return(ref.str[3])	
  })
  s1.str <- paste(s1.str, collapse = "")
  
  # get sequence idx
  sq.idx <- s1.idx - 2
  sq.str <- lapply(sq.idx, function(i){
    ref.str <- aln.str[i]
    ref.str <- strsplit(ref.str, split = " ")
    ref.str <- unlist(ref.str)
    ref.str <- ref.str[ref.str != ""]
    return(ref.str[3])	
  })
  sq.str <- paste(sq.str, collapse = "")
  
  return(c(sq.str, s1.str))
}

getStrAlignment <- function(aln){
  return(c(as.character(pattern(aln)), as.character(subject(aln))))
}

getStartPosition <- function(aln.file)
{
  #aln.str <- read_document(aln.file)
  aln.str <- read.delim(aln.file, sep = "\t", header = FALSE)
  aln.str <- as.character(t(aln.str))
  
  test.str <- aln.str[22]
  test.str <- strsplit(test.str, split = " ")
  test.str <- unlist(test.str)
  test.str <- test.str[test.str != ""]
  return(test.str[2])
}

getAlnStat <- function(aln)
{
  #pattern <- as.character(pattern(aln))
  #subject <- as.character(subject(aln))
  score <- score(aln)
  sq.id <- pid(aln)
  nb.mism <- nmismatch(aln)
  nb.indel <- nindel(aln)
  nb.ins <- insertion(nb.indel)[1]
  length.ins <- insertion(nb.indel)[2]
  nb.del <- deletion(nb.indel)[1]
  length.del <- deletion(nb.indel)[2]
  aln.chr <- compareStrings(aln)
  aln.stat <- c(score, sq.id, nb.mism, nb.ins, length.ins, nb.del, length.del, aln.chr)
  names(aln.stat) <- c("score", "sq.hom", "nb.mism","nb.ins", "length.ins", "nb.del", "length.del", "aln.chr")
  return(aln.stat)
}

getAlnChar <- function(aln.test, aln.ref, ref, idn = NULL)
{
  alnChar <- rep(NA, nchar(ref))
  insChar <- rep(0, nchar(ref))
  index <- 1
  
  for(i in 1:nchar(aln.test))
  {
    
    currentRef <- substring(ref, index, index)
    
    #print(i)
    #print(index)
    #print(currentRef)
    #print(substring(aln.test, i, i))
    #print(substring(aln.ref, i, i))
    #print("")
    
    if(substring(aln.test, i, i) == currentRef & substring(aln.ref, i, i) == currentRef){
      if(is.null(idn)) alnChar[index] <- substring(aln.ref, i, i)
      else(alnChar[index] <- ".")
      index <- index + 1
      #print("match")
    } else if(substring(aln.test, i, i) != "-" & substring(aln.test, i, i) != currentRef & substring(aln.ref, i, i) == currentRef){
      alnChar[index] <- substring(aln.test, i, i)
      index <- index + 1
      #print("mismatch")
    } else if(substring(aln.test, i, i) == "-"){
      alnChar[index] <- "-1"
      index <- index + 1
      #print("deletion")
    }else if(substring(aln.ref, i, i) == "-"){
      insChar[index - 1] <- insChar[index - 1] + 1
      #print("insertion")		
    }else(print("else"))	
  }	
  #print(alnChar)
  #print(insChar)
  finalChar <- unlist(lapply(1:length(alnChar), function(i)
    c(alnChar[i], insChar[i])
  ))
  return(finalChar)
}

getAlnSumStr <- function(aln.test, aln.ref)
{
  alnChar <- rep(NA, nchar(aln.test))
  
  for(i in 1:nchar(aln.test))
  {
    if(substring(aln.test, i, i) == substring(aln.ref, i, i)) alnChar[i] <- tolower(substring(aln.test, i, i))
    else if(substring(aln.test, i, i) == "-") alnChar[i] <- "-"
    else if(substring(aln.ref, i, i) == "-") alnChar[i] <- "+"
    else alnChar[i] <- toupper(substring(aln.test, i, i))
  }
  return(paste(alnChar, collapse=""))
}

getRangeAlignment <- function(sq, p1, p2, nR, sm)
{
  rAlnList.raw <- mclapply(nR, function(n){
    refSeq.current <- paste0(p1, strrep("N", n), p2)
    return(getTAlignment(sq, refSeq.current, sm))	
  }, mc.cores = NBCPU)
  
  rAlnList <- lapply(rAlnList.raw, function(x) x[[1]])#ALN
  names(rAlnList) <- paste0("N", nR)
  
  posList <-  lapply(rAlnList.raw, function(x) x[[2]])#POS
  names(posList) <- paste0("N", nR)
  
  return(list(ALN = rAlnList,
              POS = posList))
}

getTAlignment <- function(sq, ref, sm){
  sq.split <- splitSequence(sq, size = nchar(ref)+5)
  sq.split.aln <- lapply(sq.split, function(sqi){
    return(pairwiseAlignment(sqi, subject = ref,
                             type = "local-global", substitutionMatrix = sm,
                             gapOpening = 1, gapExtension = 1))	
  })
  
  firstBase <- sapply(sq.split.aln, function(j) substring(pattern(j), 1, 1))
  firstBase[firstBase != "T"] <- "notT"
  firstBase <- factor(firstBase, levels = c("T", "notT"))
  
  lastBase <- sapply(sq.split.aln, function(j) substring(pattern(j), nchar(pattern(j)), nchar(pattern(j))))
  lastBase[lastBase != "A"] <- "notA"
  lastBase <- factor(lastBase, levels = c("A", "notA"))

  sq.split.aln.score <- sapply(sq.split.aln, score)
  
  sq.split.aln <- sq.split.aln[order(firstBase, lastBase, -sq.split.aln.score)]
  sq.aln.selected <- sq.split.aln[[1]]
  pos <- as.numeric(strsplit(names(sq.split.aln)[1], split = "_")[[1]][2]) + start(pattern(sq.aln.selected))
  
  return(list(ALN = sq.aln.selected,
              POS = pos))
}


splitSequence <- function(sq, size = 30){
  sq.str <- strsplit(sq, split = "")[[1]]
  tindex <- which(sq.str == "T")
  tindex <- tindex[tindex < length(sq.str) - size]
  sq.split <- lapply(tindex, function(x){
    paste0(sq.str[seq(x, x+size)], collapse = "")
  })
  names(sq.split) <- paste0("index_", tindex)
  return(sq.split)
}


getGuideAlignmentDual <- function(inputF, outputF, guideLeft, guideRight, gnm = BSgenome.Hsapiens.UCSC.hg38::Hsapiens, plot = TRUE, binCoord = FALSE)
{
  # DEFINE N RANGE
  nRange <- seq(8, 28)
  
  # DEFINE REFERENCE SEQUENCE
  refSeq.left <- guideLeft
  refSeq.right <- guideRight
  
  refSeq.left.rev <- getRevComp(refSeq.left)
  refSeq.right.rev <- getRevComp(refSeq.right)
  
  refSeqMat <- matrix(c(refSeq.left, refSeq.left.rev,
                        refSeq.left, refSeq.right.rev,
                        refSeq.right, refSeq.right.rev,
                        refSeq.right, refSeq.left.rev
  ),
  ncol = 2, byrow = TRUE)
  
  rownames(refSeqMat) <- c("LF.LR", "LF.RR", "RF.RR","RF.LR")
  
  # LOAD FILE
  if(file_ext(inputF) == "bed"){
    readMat <- read.delim(inputF, header = FALSE)
    colnames(readMat) <- c("chromosome", "start", "end", "ID_read", "size", "strand")
  }else{
    readMat <- read.xlsx(inputF)
  }	
  
  # USE BIN COORDINATES (DEFAULT); OTHERWISE:
  if(!binCoord & "start.site" %in% colnames(readMat)){
    readMat$start <- readMat$start.site
    readMat$end <- readMat$end.site
  }
  
  inputName <- gsub(paste0(".", file_ext(inputF)), "", basename(inputF))	
  print(inputName)
  alnFolder <- dirname(inputF)
  
  # REGIONS 2 SEQUENCE
  sequences <- bed2sequence(readMat, g = gnm)
  if(file_ext(inputF) == "bed") sequences <- sequences[!grepl("N", sequences)]
  print(length(sequences))
  #sequences <- c(refSeq, sequences)# add refSeq
  
  #print(sequences[[1]])
  
  ###########
  # ALIGNMENT
  #sm <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
  #sm[, "N"] <- 0
  #sm["N", ] <- 0
  #sm["N", "N"] <- 0
  
  # Toni's matrix
  #sm <- matrix(c(1,0.6,0.7,0.3,0,0.9,1,0.8,0.8,0,0.3,0.2,1,0.4,0,0.5,0.5,0.7,1,0,0,0,0,0,0),
  #             nrow = 5, ncol = 5, byrow = TRUE)
  sm <- matrix(c(2,0.6,0.7,0.3,0,0.9,2,0.8,0.8,0,0.3,0.2,2,0.4,0,0.5,0.5,0.7,2,0,0,0,0,0,0),
               nrow = 5, ncol = 5, byrow = TRUE)
  rownames(sm) <- c("G", "A", "T", "C", "N")
  colnames(sm) <- c("G", "A", "T", "C", "N")
  
  print(sm)
  
  # ALIGNMENT
  alnList.raw <- lapply(sequences, function(i){
    x.raw <- lapply(1:nrow(refSeqMat), function(j){
      getRangeAlignment(i, refSeqMat[j, 1], refSeqMat[j, 2], nRange, sm)
    })
    
    x <- lapply(x.raw, function(k) k[[1]])
    names(x) <- rownames(refSeqMat)
    
    pos <- lapply(x.raw, function(k) k[[2]])
    names(pos) <- rownames(refSeqMat)
    
    return(list(ALN = x,
                POS = pos))
  })
  
  alnList <- lapply(alnList.raw, function(k) k[[1]])
  names(alnList) <- paste0(inputName, "_aln", 1:length(alnList.raw))
  
  posList <- lapply(alnList.raw, function(k) k[[2]])
  names(posList) <- paste0(inputName, "_aln", 1:length(alnList.raw))
  
  # SAVE RDS
  saveRDS(alnList, file.path(alnFolder, "alnList.raw.RDS"))
  saveRDS(posList, file.path(alnFolder, "posList.raw.RDS"))
  
  # score matrix per site
  scoreList.raw <- mclapply(1:length(alnList), function(i){
    alnList.current <- alnList[[i]]
    scores <- lapply(alnList.current, function(j) sapply(j, score))
    scores <- do.call(rbind, scores)
    return(scores)
  }, mc.cores = NBCPU)
  names(scoreList.raw) <- names(alnList)
  
  # SAVE HEATMAPS AND TXT FILES
  nDir <- file.path(alnFolder, "N_alignment")
  dir.create(nDir, showWarnings = FALSE)
  if(plot) dir.create(file.path(nDir, "score_heatmaps"))
  
  lapply(1:length(scoreList.raw), function(i){
    write.table(scoreList.raw[[i]], sep = "\t", quote = FALSE,
                file.path(nDir, paste0(names(scoreList.raw)[i], "_score_heatmap.txt")))
    if(plot){
      pheatmap(scoreList.raw[[i]], color = magma(10),
               filename = file.path(nDir, "score_heatmaps", paste0(names(scoreList.raw)[i], "_score_heatmap.pdf")),
               show_rownames = TRUE, show_colnames = TRUE,
               cluster_cols = FALSE,
               cellwidth = 12, cellheight = 12, fontsize = 10,
               width = length(nRange)*6/16, height = 5)
    }
    
  })
  
  # 
  write.xlsx(readMat, outputF, overwrite = TRUE)	
  
}


assignPV <- function(inputF, rdF, outputF, nbTop = 3, plot = TRUE){
  alnFolder <- dirname(inputF)
  inputName <- gsub(paste0(".", file_ext(inputF)), "", basename(inputF))
  readMat <- read.xlsx(inputF, sheet = 1)
  
  alnFolder.rd <- dirname(rdF)
  
  # LOAD CASTseq AND RANDOM
  alnMatFiles <- list.files(file.path(alnFolder, "N_alignment"), pattern = "_score_heatmap.txt", full.names = TRUE)
  names(alnMatFiles) <- gsub("_score_heatmap.txt", "", basename(alnMatFiles))
  alnMatList <- lapply(alnMatFiles, read.delim)
  
  alnMatFiles.rd <- list.files(file.path(alnFolder.rd, "N_alignment"), pattern = "_score_heatmap.txt", full.names = TRUE)
  names(alnMatFiles.rd) <- gsub("_score_heatmap.txt", "", basename(alnMatFiles.rd))
  alnMatList.rd <- lapply(alnMatFiles.rd, read.delim)
  
  nbRow <- nrow(alnMatList[[1]])
  nbCol <- ncol(alnMatList[[1]])
  myRowname <- rownames(alnMatList[[1]])
  myColname <- colnames(alnMatList[[1]])
  pvMatList <- lapply(1:length(alnMatList), function(i){
    pvMat <- matrix(0, nrow = nbRow, ncol = nbCol)
    rownames(pvMat) <- myRowname
    colnames(pvMat) <- myColname
    return(pvMat)
  })
  names(pvMatList) <- names(alnMatList)
  
  # FILL pvMatList
  for(i in 1:nbRow){
    for(j in 1:nbCol){
      score.rd <- sapply(alnMatList.rd, function(alnMat.rd) alnMat.rd[i, j])
      rd.ecdf <- stats::ecdf(score.rd)
      
      for(idx in 1:length(alnMatList)){
        pv <- 1 - rd.ecdf(alnMatList[[idx]][i, j])
        pvMatList[[idx]][i, j] <- pv
      }
    }
  }
  
  # plot pvalue heatmap
  if(plot){
    dir.create(file.path(alnFolder, "N_alignment", "pv_heatmaps"))
    lapply(1:length(pvMatList), function(i){
      mymat <- -log10(pvMatList[[i]])
      mymat[mymat>5] <- 5
      pheatmap(mymat, color = magma(10),
               filename = file.path(alnFolder, "N_alignment", "pv_heatmaps", gsub("_score_heatmap.txt", "_pv_heatmap.pdf", basename(alnMatFiles[i]))),
               show_rownames = TRUE, show_colnames = TRUE,
               cluster_cols = FALSE,
               cellwidth = 12, cellheight = 12, fontsize = 10,
               width = ncol(mymat)*6/16, height = 5)
    })
  }

  
  
  # SELECT THE TOP 10 BEST SITES (BASED ON PVALUE AND SCORE)
  alnMatList.long <- lapply(1:length(alnMatList), function(i){
    wide <- data.table(COMB = rownames(alnMatList[[i]]),
                       alnMatList[[i]])
    long <- melt(setDT(wide), id.vars = c("COMB"), variable.name = "N")
    
    wide.pv <- data.table(COMB = rownames(pvMatList[[i]]),
                          pvMatList[[i]])
    long.pv <- melt(setDT(wide.pv), id.vars = c("COMB"), variable.name = "N")
    
    long$PV <- long.pv$value
    long <- long[order(-long$value, long$PV), ]# score then pvalue
    long$TOP <- paste0("top", 1:nrow(long))
    return(long)
  })
  names(alnMatList.long) <- names(alnMatList)
  
  alnMatList.long <- lapply(alnMatList.long, function(i) i[1:nbTop])
  
  alnList.raw <- readRDS(file.path(alnFolder, "alnList.raw.RDS"))
  posList.raw <- readRDS(file.path(alnFolder, "posList.raw.RDS"))
  alnMatList.long <- alnMatList.long[match(names(alnList.raw), names(alnMatList.long))]
  
  alnList <- lapply(1:length(alnMatList.long), function(i){
    alnMat.long <- alnMatList.long[[i]]
    topList <- lapply(1:nrow(alnMat.long), function(j){
      comb.current <- alnMat.long$COMB[j]
      n.current <- alnMat.long$N[j]
      return(alnList.raw[[i]][[comb.current]][[n.current]])
    })
    names(topList) <- alnMat.long$TOP
    return(topList)
  })
  names(alnList) <- names(alnMatList.long)
  
  posList <- lapply(1:length(alnMatList.long), function(i){
    alnMat.long <- alnMatList.long[[i]]
    topList <- lapply(1:nrow(alnMat.long), function(j){
      comb.current <- alnMat.long$COMB[j]
      n.current <- alnMat.long$N[j]
      return(posList.raw[[i]][[comb.current]][[n.current]])
    })
    names(topList) <- alnMat.long$TOP
    return(topList)
  })
  names(posList) <- names(alnMatList.long)
  
  statMatList <- lapply(1:length(alnList), function(j){
    # GET PROPER ALIGNMENT (FIX ISSUE WITH INDEL IN FIRST / LAST POSITION)

    alnList.sub <- alnList[[j]]
    
    statMat.current <- lapply(1:length(alnList.sub), function(i){
      alnList.sub.sub <- alnList.sub[[i]]
      
      alnStr <- getStrAlignment(alnList.sub.sub)
      names(alnStr) <- c("pattern", "subject")
      
      alnSum <- getAlnSumStr(alnStr[1], alnStr[2])
      names(alnSum) <- "aln.sum"
      
      alnStart <- as.numeric(posList[[j]][[i]])
      middleCoord <- floor((readMat$end[j] - readMat$start[j]) / 2)
      middleCoord.abs <- readMat$start[j] + middleCoord
      names(middleCoord.abs) <- "middleCoord"
      alnStart.rel <- alnStart - middleCoord
      names(alnStart.rel) <- "aln.start.rel"
      alnStart.abs <- alnStart + readMat$start[j]
      names(alnStart.abs) <- "aln.start.abs"
      
      alnStat <- getAlnStat(alnList[[j]][[i]])
      
      statMat <- c(alnStr, alnSum, alnStat, middleCoord.abs, alnStart.abs, alnStart.rel)
      statMat <- as.data.frame(t(statMat))
      statMat[, c("score", "sq.hom", "nb.mism", "nb.ins", "length.ins", "nb.del", "length.del", "middleCoord", "aln.start.abs", "aln.start.rel")] <- apply(
        statMat[, c("score", "sq.hom", "nb.mism", "nb.ins",
                    "length.ins", "nb.del", "length.del", "middleCoord", "aln.start.abs", "aln.start.rel")], 2, as.numeric)
      colnames(statMat) <- paste0("top", i, "_", colnames(statMat))
      statMat <- statMat[, !grepl("score", colnames(statMat)), drop = FALSE]# remove scroe column
      return(statMat)
    })
    statMat.current <- do.call(cbind, statMat.current)
    return(statMat.current)
  })
  
  statMat <- do.call(rbind, statMatList)
  
  topMat <- lapply(alnMatList.long, function(i){
    topComb <- i$COMB
    names(topComb) <- paste0(i$TOP, "_COMB")
    return(topComb)
  })
  topMat <- do.call(rbind, topMat)
  
  nMat <- lapply(alnMatList.long, function(i){
    topN <- as.numeric(gsub("N", "", i$N))
    names(topN) <- paste0(i$TOP, "_N")
    return(topN)
  })
  nMat <- do.call(rbind, nMat)
  
  scoreMat <- lapply(alnMatList.long, function(i){
    topScore <- i$value
    names(topScore) <- paste0(i$TOP, "_score")
    return(topScore)
  })
  scoreMat <- do.call(rbind, scoreMat)
  
  pvMat <- lapply(alnMatList.long, function(i){
    topN <- i$PV
    names(topN) <- paste0(i$TOP, "_pv")
    return(topN)
  })
  pvMat <- do.call(rbind, pvMat)
  
  topMat.sum <- lapply(alnMatList.long, function(i){
    topstr <- apply(i, 1, function(j){
      paste(j[1:4], collapse = ";")
    })
    names(topstr) <- paste0(i$TOP, "_STR")
    return(data.frame(t(topstr)))
  })
  topMat.sum <- do.call(rbind, topMat.sum)
  
  # MERGE READMAT, TOPMAT AND STATMAT
  readMat.stat <- cbind(readMat,
                        alignment.id = names(alnList),
                        topMat.sum,
                        topMat,
                        nMat,
                        scoreMat,
                        pvMat,
                        statMat)	
  write.xlsx(readMat.stat, outputF, overwrite = TRUE)	
}





assignBestCombinationDoubleNickase <- function(inputFile)
{
  cmb <- c("LF.LR", "LF.RF", "LF.RR", "LR.RF", "LR.RR", "RF.RR",
           "LR.LF", "RF.LF", "RR.LF", "RF.LR", "RR.LR", "RR.RF",
           "LF.LF", "LR.LR", "RF.RF", "RR.RR")
  
  readMat <- read.xlsx(inputFile, sheet = 1)
  #lapply(cmb, function(i) grep(i, colnames(readMat)))
  
  dMat <- readMat[, grep("_aln.start.rel$", colnames(readMat))]
  nMat <- readMat[, grep("_N$", colnames(readMat))]
  cbName <- unlist(lapply(strsplit(colnames(nMat), split = "_"), function(i) i[1]))
  
  isMono <- paste0(substring(cbName, 2, 2), substring(cbName, 5, 5)) %in% c("FF", "RR")
  
  bestCb <- lapply(1:nrow(readMat), function(i){
    
    qv <- as.numeric(readMat[i, grep("_adj.pv$", colnames(readMat))])
    names(qv) <- cbName
    
    gRNAscore <- as.numeric(readMat[i, grep("_score$", colnames(readMat))])
    names(gRNAscore) <- cbName
    
    if(min(qv) > 0.1) return(NA)
    
    signifIdx <- which(qv <= 0.1)
    
    if(length(signifIdx) == 1) return(cbName[signifIdx])
    
    qv.sub <- qv[signifIdx]
    gRNAscore.sub <- gRNAscore[signifIdx]
    dMat.sub <- dMat[i,signifIdx]
    nMat.sub <- nMat[i,signifIdx]
    cbName.sub <- cbName[signifIdx]
    isMono.sub <- isMono[signifIdx]
    
    maxScore <- max(gRNAscore.sub)
    maxScoreIdx <- which(gRNAscore.sub == maxScore)
    
    if(length(maxScoreIdx) == 1) return(cbName.sub[maxScoreIdx])
    
    qv.sub <- qv.sub[maxScoreIdx]
    gRNAscore.sub <- gRNAscore.sub[maxScoreIdx]
    dMat.sub <- dMat.sub[maxScoreIdx]
    nMat.sub <- nMat.sub[maxScoreIdx]
    cbName.sub <- cbName.sub[maxScoreIdx]
    isMono.sub <- isMono.sub[maxScoreIdx]
    
    score <- unlist(lapply(1:length(maxScoreIdx), function(j){
      score.current <- 0
      if(nMat.sub[, j] < 10){score.current <- score.current + 1
      }else if(nMat.sub[, j] < 20) score.current <- score.current + 2
      
      if(abs(dMat.sub[, j]) < 50) score.current <- score.current + 1
      
      if(!isMono.sub[j]) score.current <- score.current + 1
      #if(isMono.sub[j]) score.current <- score.current - 1
      
      return(score.current)		
    }))
    names(score) <- names(qv.sub)
    
    qv.sub.sub <- qv.sub[names(score)[score == max(score)]]
    
    return(names(qv.sub.sub[order(qv.sub.sub)])[1])
  })
  bestCb <- unlist(bestCb)
  
  readMat$BestCB <- bestCb
  
  bestCbN <- sapply(1:length(bestCb), function(i){
    if(is.na(bestCb[i])) return(NA)
    
    return(as.numeric(readMat[i, paste0(bestCb[i], "_N")]))
  })
  readMat$N.BestCB <- bestCbN
  
  write.xlsx(readMat, inputFile, overwrite = TRUE)
}



plotSites <- function(inputF, siteIdx, outputF)
{
  readMat <- read.xlsx(inputF, sheet = 1)
  readMat <- readMat[siteIdx, ]# select sites
  rownames(readMat) <- apply(readMat[1:2], 1, paste, collapse = ";")
  
  scoreM <- readMat[,grep("_score", colnames(readMat))]
  scoreM$loc <- rownames(scoreM)
  
  distM <- readMat[,grep("_aln.start.rel", colnames(readMat))]
  distM$loc <- rownames(distM)
  
  nM <- readMat[,grep("_N", colnames(readMat))]
  nM$loc <- rownames(distM)
  
  scoreM.gg <- melt(scoreM)
  distM.gg <- melt(distM)
  nM.gg <- melt(nM)
  
  ggmat <- cbind(scoreM.gg, distance = distM.gg$value, N = nM.gg$value)
  ggmat$variable <- gsub("_score", "", ggmat$variable)
  
  p <- ggplot(ggmat, aes(x=distance, y=value, color = N))
  p <- p + geom_point(size = 4, alpha = 0.5) + scale_color_viridis_c(option = "B")
  p <- p + facet_wrap(~ loc)
  p <- p + theme_bw()#+ theme(legend.position="none")
  p <- p + geom_vline(xintercept=0, linetype = "dashed")
  p <- p + geom_text_repel(aes(label = variable), size = 3)
  
  pdf(outputF, width = 10, height = 8)
  plot(p)
  dev.off()	
}





############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE){
  # TEST
  require(BSgenome.Hsapiens.UCSC.hg38)
  require(openxlsx)
  NBCPU <- 4
  inputF <- file.path("/Users/geoffroyandrieux/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9/KRT9-T3_g1/results/guide_aln/KRT9-T3_w250_FINAL.xlsx")
  guideLeft <- toupper("CCGAGAATTGAGTTCCTGCANRG")
  guideRight <- toupper("CTCTTACTTGGATAAGGTGCNRG")
  alnFolder <- file.path("~/Research/CASTSeq/test")
  getGuideAlignmentDual(inputF, guideLeft, guideRight, alnFolder, BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
  
  # TEST MASAKO STAT3 24.08.22
  ovlD <- "~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/STAT3/TALEN/STAT3-CD34-T3-GSE56/OVL2_SIGNIF1/"
  projectName <- "STAT3-CD34-T3-GSE56"
  randomD <- "~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/STAT3/TALEN/STAT3-CD34-T3-GSE56/random/"
  randomName <- "random_w200"
  
  inputF <- file.path(ovlD, paste0(projectName, "_FINAL.xlsx"))
  rdF <- file.path(randomD, paste0(randomName, "_ALN.xlsx"))
  outputF <- file.path(ovlD, paste0(projectName, "_ALN_PV.xlsx"))
  nbTop = 3
  plot = FALSE
  
  assignPV(inputF,
          rdF,
          outputF,
          nbTop = 3,
          plot = FALSE)
}

