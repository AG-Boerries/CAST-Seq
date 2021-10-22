
# load model

ddPCR_linear <- function(inputFile, ciCutoff, modelRDS){
  lmObject <- readRDS(modelRDS)
  
  readMat <- read.xlsx(inputFile, sheet = 1)
  
  newData <- data.frame(CAST = readMat$hits)
  newData.log <- newData
  newData.log[newData.log < 0] <- 0
  newData.log <- log10(newData.log+1)
  
  pred <- predict(lmObject, newData.log, se.fit = TRUE, interval = "confidence", level = ciCutoff)
  #pred <- predict(lmObject, newData, se.fit = TRUE, interval = "prediction", level = 0.95)
  
  pred$fit[pred$fit < 0] <- 0
  pred.raw <- round((10^pred$fit[,1])-1)
  pred.lwr <- (10^pred$fit[,2])-1
  pred.upr <- (10^pred$fit[,3])-1
  pred.str <- sapply(1:length(pred.lwr), function(i){
    paste0("(", round(pred.lwr[i], 0), " : ", round(pred.upr[i], 0), ")")
  })
  
  newData <- cbind(newData, pred.raw, pred.str)
  colnames(newData) <- c("CASTSeq.hits", "ddPCR.estimate", paste0("CI(",ciCutoff*100, "%)"))
  
  newData <- cbind(readMat[, c("chromosome", "start", "end", "group")], newData)
  
  write.xlsx(newData, gsub("_FINAL.xlsx", "_ddPCR_estimate.xlsx", inputFile),
             row.names = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), overwrite = TRUE)
}


if(FALSE){
  library(openxlsx)
  
  inputFile <- file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/results/guide_aln/EMD101-sample1-1_w250_FINAL.xlsx")
  ciCutoff = 0.9
  modelRDS <- file.path("~/Research/CASTSeq/pipelineGit/script/lmObject.RDS")
  
  ddPCR_linear(inputFile, ciCutoff, modelRDS)
  
}