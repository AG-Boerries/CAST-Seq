
############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################


posTransform.fun = function(region) {
  return(region)
}

extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  
  print(zoom_bed)
  
  rbind(bed, zoom_bed)
}

extend_cytoband = function(bed, chromosome, start, end, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  
  idx <- which(zoom_bed[[2]] <= start & zoom_bed[[3]] >= start
               | zoom_bed[[2]] <= end & zoom_bed[[3]] >= end)
  if(length(idx) > 1){
    zoom_bed = zoom_bed[idx[1]:idx[2], ]
  }else zoom_bed = zoom_bed[idx, ]
  
  zoom_bed[[2]][1] <- start
  zoom_bed[[3]][nrow(zoom_bed)] <- end
  #print(zoom_bed)
  
  rbind(bed, zoom_bed)
}

extend_bed = function(bed, chromosome, start, end, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  
  idx <- which(zoom_bed[[2]] >= start & zoom_bed[[3]] <= end)
  zoom_bed = zoom_bed[idx, ]
  
  #zoom_bed[[2]][1] <- start
  #zoom_bed[[3]][nrow(zoom_bed)] <- end
  #print(zoom_bed)
  
  rbind(bed, zoom_bed)
}



circlizePipeline <- function(siteFile, zoom.size = 50000, label = FALSE,
                             bestScore.cutoff = 9, bestFlank.cutoff = 25,
                             gene.bed = NULL, ots.bed = NULL, 
                             realigned = FALSE,
                             outFile,
                             cytobandFile = NULL){
  # READ INPUT FILE
  siteM <- read.xlsx(siteFile)
  
  siteM.chr <- factor(gsub("chr", "", siteM$chromosome), levels = c(1:22, "X", "Y"))
  siteM <- siteM[order(siteM.chr, siteM$start), ]
  
  # CHANGE MIDDLE COORDINATE OF THE ON TARGET
  if(!is.null(ots.bed) & realigned){
    ots.df <- read.delim(ots.bed, header = FALSE)
    ots.df <- ots.df[1,,drop = FALSE]# only take the first ON target into account
    
    on.middle <- ots.df[,2] + round((ots.df[,3] - ots.df[,2]) / 2)
    
    # intersect
    ots.gr <- makeGRangesFromDataFrame(ots.df,
                                       seqnames.field = "V1", start.field = "V2", end.field = "V3",
                                       keep.extra.columns = FALSE, ignore.strand = TRUE)
    siteM.gr <- makeGRangesFromDataFrame(siteM,
                                         seqnames.field = "chromosome", start.field = "start", end.field = "end",
                                         keep.extra.columns = FALSE, ignore.strand = TRUE)
    
    gr.ovl <- findOverlaps(query = ots.gr, subject = siteM.gr, type = "any", maxgap = 0)
    
    
    siteM$middleCoord[subjectHits(gr.ovl)] <- on.middle
  }
  
  # DEFINE ZOOM
  if(is.null(ots.bed)){
    idx <- which.max(siteM$read)
    zoom.chr <- siteM$chromosome[idx]
    zoom.start <- siteM$middleCoord[idx] - zoom.size
    zoom.end <- siteM$middleCoord[idx] + zoom.size
  }else{
    ots.df <- read.delim(ots.bed, header = FALSE)
    zoom.chr <- ots.df[1,1]
    zoom.start <- ots.df[1,2] - zoom.size
    zoom.end <- ots.df[1,3] + zoom.size
  }
  
  # ADD ZOOM
  siteM <- extend_bed(siteM, zoom.chr, start = zoom.start, end = zoom.end)
  
  # INIT (MUST BE ONLINE !!!)
  cytoband = read.cytoband(cytobandFile)
  cytoband_df = cytoband$df
  chromosome = cytoband$chromosome
  
  new_cytoband_df <- extend_cytoband(cytoband_df, zoom.chr, start = zoom.start, end = zoom.end)
  
  xrange = c(cytoband$chr.len, cytoband$chr.len[zoom.chr])
  normal_chr_index = 1:length(cytoband$chr.len)
  zoomed_chr_index = length(cytoband$chr.len)+1
  
  sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                   xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index]))
  sector.width[normal_chr_index] <- sector.width[normal_chr_index] 
  sector.width[zoomed_chr_index] <- sector.width[zoomed_chr_index] / 1.25
  
  pdf(outFile, width = 3, height = 3)
  circos.clear()
  circos.par(start.degree = 90, "gap.degree" = c(rep(3, (length(normal_chr_index) - 1)), 10, 20))
  
  
  circos.initializeWithIdeogram(new_cytoband_df, 
                                sector.width = sector.width,
                                plotType = NULL)
  
  ##############
  # GROUP LABELS
  zoom.idx <- grep("^zoom_", siteM$chromosome)
  zoom.coord <- sapply(zoom.idx, function(i) paste(siteM[i, c("start", "end")], collapse = "-"))
  all.coord <- sapply(1:nrow(siteM), function(i) paste(siteM[i, c("start", "end")], collapse = "-"))
  toRemove <- setdiff(which(all.coord %in% zoom.coord), zoom.idx)
  
  siteM.sub <- siteM[-toRemove, ]# remove duplicated sites, keep only the ones in zoom area
  group <- siteM.sub[, "group"]
  #group[group == "off.target"] <- "OMT"
  #group[group == "HR"] <- "HMT"
  group[group == "OMT" & siteM.sub$is.hom.recomb. == "yes"] <- "OMT/HMT"
  #group[group == "CBS"] <- "NBS"
  
  # ON target
  if(!is.null(ots.bed)){
    ots.df <- read.delim(ots.bed, header = FALSE)
    ots.df <- ots.df[1,,drop = FALSE]# only take the first ON target into account
    ots.df[,1] <- paste0("zoom_", ots.df[,1])
    ots.gr <- makeGRangesFromDataFrame(ots.df,
                                       seqnames.field = "V1", start.field = "V2", end.field = "V3",
                                       keep.extra.columns = FALSE, ignore.strand = TRUE)
    siteM.sub.gr <- makeGRangesFromDataFrame(siteM.sub,
                                             seqnames.field = "chromosome", start.field = "start", end.field = "end",
                                             keep.extra.columns = FALSE, ignore.strand = TRUE)
    
    gr.ovl <- findOverlaps(query = ots.gr, subject = siteM.sub.gr, type = "any", maxgap = 0)
    
    group[subjectHits(gr.ovl)] <- "ON"
  }else{# Use the on with higher reads as ON target
    group[which.max(siteM.sub$read)] <- "ON"
  }
  
  group.bed <- data.frame(chr = siteM.sub$chromosome,
                          start =  siteM.sub$middleCoord,
                          end = siteM.sub$middleCoord,
                          value = group)
  group.color <- rep("grey", length(group))
  group.color[group == "ON"] <- "green"
  group.color[group == "OMT"] <- "red"
  group.color[group == "HMT"] <- "blue"
  group.color[group == "OMT/HMT"] <- "goldenrod1"
  
  
  # REMOVE NBS LABELS
  isNBS <- group.bed$value == "NBS"
  group.bed <- group.bed[!isNBS, ]
  group.color <- group.color[!isNBS]
  
  if(label){
    circos.genomicLabels(group.bed, labels.column = 4, side = "outside",
                         cex = 0.5, font = 1,
                         col = group.color, line_col = group.color,
                         line_lwd = par("lwd")*1.5, 
                         connection_height = convert_height(0.5, "mm"),
                         labels_height = convert_height(0.5, "cm"),
                         track.margin = c(0.025,0))
  }
  
  
  # ADD CHR NAMES
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,1, gsub("chr", "", CELL_META$sector.index), cex = 0.8, col = "black",
                facing = "clockwise", niceFacing = TRUE)
  }, track.height = 0.05, bg.border = NA, track.margin = c(0,0.015))
  
  # ADD IDEOGRAM
  circos.genomicIdeogram(new_cytoband_df,  track.margin = c(0,0.015))
  
  ##################
  # BEST SCORE TRACK
  bestScore <- siteM[, "score"]
  bestScore.chr <- siteM[, "chromosome"]
  bestScore.x <- siteM$middleCoord
  
  # fill missing chr
  missing.chr <- setdiff(unique(new_cytoband_df[[1]]), bestScore.chr)
  bestScore <- c(bestScore, rep(0, length(missing.chr)))
  bestScore.chr <- c(bestScore.chr, missing.chr)
  bestScore.x  <- c(bestScore.x, rep(0, length(missing.chr)))
  
  circos.track(factors = bestScore.chr, x = bestScore.x, y = bestScore, bg.col = "#EEEEEE",
               bg.border = NA, track.height = 0.1, panel.fun = function(x, y) {
                 
                 cell.xlim = get.cell.meta.data("cell.xlim")
                 cell.ylim = get.cell.meta.data("cell.ylim")
                 
                 print(cell.xlim)
                 print(cell.ylim)
                 
                 # reference lines
                 for(yi in floor(seq(0, max(bestScore), length.out = 3))) {
                   circos.lines(cell.xlim, c(yi, yi), lty = 2, col = "white") 
                 }
                 
                 xlim = get.cell.meta.data("xlim")
                 ylim = get.cell.meta.data("ylim")
                 
                 #print(xlim)
                 #print(ylim)
                 circos.rect(xlim[1], bestScore.cutoff, xlim[2], ceiling(ylim[2])+3, col = "#FF000020", border = NA)
                 
                 circos.points(x[y < bestScore.cutoff], y[y < bestScore.cutoff & y !=0], pch = 16, cex = 0.75, col = "grey")
                 circos.points(x[y >= bestScore.cutoff], y[y >= bestScore.cutoff], pch = 16, cex = 0.75, col = "red")
               })
  
  circos.yaxis(side = "left", at = c(0, bestScore.cutoff, max(bestScore)), #floor(seq(0, max(bestScore), length.out = 3))
               sector.index = get.all.sector.index()[1], labels.cex = 0.5,
               tick.length = convert_x(0.75, "mm", get.all.sector.index()[1],get.cell.meta.data("track.index")),
               lwd = par("lwd") * 1.25)
  
  # BEST SCORE (ZOOM ONLY)
  idx <- siteM$chromosome == paste0("zoom_", zoom.chr)
  bestScore.zoom <- siteM[idx, "score"]
  bestScore.zoom.chr <- siteM[idx, "chromosome"]
  bestScore.zoom.x1 <- siteM[idx, "start"]
  bestScore.zoom.x2 <- siteM[idx, "end"]
  
  bestScoreM.zoom <- data.frame(chr = bestScore.zoom.chr,
                                start = bestScore.zoom.x1,
                                end= bestScore.zoom.x2,
                                value1 = bestScore.zoom)
  
  circos.genomicTrackPlotRegion(bestScoreM.zoom, stack = TRUE, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = "black", border = NA, ytop = 0, ybottom = -1)
    #circos.genomicText(region, value, labels = c(2,5))
  }, bg.border = NA, bg.col = NA, track.index = get.cell.meta.data("track.index") - 2)
  
  ##################
  # BEST FLANK TRACK
  bestFlank <- apply(as.matrix(siteM[, c("flank.length", "flank.rev.length")]), 1, max)
  bestFlank[bestFlank > 30] <- 30
  bestFlank.chr <- siteM[, "chromosome"]
  #bestFlank.cutoff <- 25
  bestFlank.x <- siteM$middleCoord
  #bestFlank.x <-siteM$start + round((siteM$end - siteM$start) / 2)
  
  # fill missing chr
  missing.chr <- setdiff(unique(new_cytoband_df[[1]]), bestFlank.chr)
  bestFlank <- c(bestFlank, rep(0, length(missing.chr)))
  bestFlank.chr <- c(bestFlank.chr, missing.chr)
  bestFlank.x  <- c(bestFlank.x, rep(0, length(missing.chr)))
  
  
  circos.track(factors = bestFlank.chr, x = bestFlank.x, y = bestFlank, bg.col = "#EEEEEE",
               bg.border = NA, track.height = 0.1, panel.fun = function(x, y) {
                 
                 cell.xlim = get.cell.meta.data("cell.xlim")
                 cell.ylim = get.cell.meta.data("cell.ylim")
                 
                 print(cell.xlim)
                 print(cell.ylim)
                 
                 # reference lines
                 #for(xi in seq(cell.xlim[1], cell.xlim[2], length.out = 3)) {
                 # circos.lines(c(xi, xi), cell.ylim, lty = 2, col = "white") 
                 #}
                 for(yi in floor(seq(0, max(bestFlank), length.out = 3))) {
                   circos.lines(cell.xlim, c(yi, yi), lty = 2, col = "white") 
                 }
                 
                 xlim = get.cell.meta.data("xlim")
                 ylim = get.cell.meta.data("ylim")
                 
                 print(xlim)
                 print(ylim)
                 circos.rect(xlim[1], bestFlank.cutoff, xlim[2], ceiling(ylim[2])+5, col = "cornflowerblue", border = NA)
                 
                 circos.points(x[y < bestFlank.cutoff], y[y < bestFlank.cutoff & y !=0], pch = 16, cex = 0.75, col = "grey")
                 circos.points(x[y >= bestFlank.cutoff], y[y >= bestFlank.cutoff], pch = 16, cex = 0.75, col = "blue")
                 
               },  track.margin = c(0,0))
  
  circos.yaxis(side = "left", at = c(0, bestFlank.cutoff , max(bestFlank)), #floor(seq(0, max(bestFlank), length.out = 3))
               sector.index = get.all.sector.index()[1], labels.cex = 0.5,
               tick.length = convert_x(0.75, "mm", get.all.sector.index()[1], get.cell.meta.data("track.index")),
               lwd = par("lwd") * 1.25)
  
  #################
  # GENE ANNOTATION
  if(!is.null(gene.bed)){
    geneM <- read.delim(gene.bed, header = FALSE)# must be 4 columns: chr, start, end, symbol
    geneM[,1] <- paste0("zoom_", geneM[,1])
    
    if(class(geneM[,2]) == "factor") geneM[,2] <- as.numeric(levels(geneM[,2]))[geneM[,2]]
    if(class(geneM[,3]) == "factor") geneM[,3] <- as.numeric(levels(geneM[,3]))[geneM[,3]]
    if(class(geneM[,4]) == "factor") geneM[,4] <- as.character(geneM[,4])
    color <- c("grey40")
    
    circos.genomicTrackPlotRegion(geneM, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color, border = NA)
      #circos.genomicText(region, value, labels = c(2,5))
    }, bg.border = NA, bg.col = NA, track.height = 0.05, track.index = get.cell.meta.data("track.index") - 2, track.margin = c(0,0))
    
  }
  

  
  ###########
  # LINK ZOOM
  
  # intersect zoom and original chr indices
  zoom.df <- new_cytoband_df[grepl("^zoom_", new_cytoband_df[,1]), ]
  zoom.df[,1] <- gsub("zoom_", "", zoom.df[,1])
  zoom.gr <- makeGRangesFromDataFrame(zoom.df,
                                      seqnames.field = "V1", start.field = "V2", end.field = "V3",
                                      keep.extra.columns = FALSE, ignore.strand = TRUE)
  ori.gr <- makeGRangesFromDataFrame(new_cytoband_df[!grepl("^zoom_", new_cytoband_df[,1]), ],
                                     seqnames.field = "V1", start.field = "V2", end.field = "V3",
                                     keep.extra.columns = FALSE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = zoom.gr, subject = ori.gr, type = "any", maxgap = 0)
  
  circos.genomicLink(new_cytoband_df[grepl("^zoom_", new_cytoband_df[,1]),],
                     new_cytoband_df[subjectHits(gr.ovl),],
                     col = rgb(169, 169, 169, max = 255, alpha = 25),
                     border = NA)
  
  ################################
  # LINK ON TARGETS TO ALL TARGETS
  bed2 <- siteM.sub[, 1:5]
  #bed2$start <- bed2$start + round((bed2$end - bed2$start) / 2)
  bed2$start <- siteM.sub$middleCoord
  
  # Re-Arange layout
  bed2$end <- bed2$start + 1000
  bed2$start <- bed2$start - 1000
  
  bed2$start[bed2$chromosome != "zoom_chr3"] <- bed2$start[bed2$chromosome != "zoom_chr3"] - 20000000
  bed2$end[bed2$chromosome != "zoom_chr3"] <- bed2$end[bed2$chromosome != "zoom_chr3"] + 20000000
  
  
  # INCREASE LINE THINKNESS OF HIGH COVERED SITES
  #isHigh <- bed2$hits >= 5
  #bed2$start[isHigh] <- bed2$start[isHigh] - 500
  #bed2$end[isHigh] <- bed2$end[isHigh] + 500
  #bed2$start[isHigh & bed2$chromosome != "zoom_chr3"] <- bed2$start[isHigh& bed2$chromosome != "zoom_chr3"] - 20000000
  #bed2$end[isHigh& bed2$chromosome != "zoom_chr3"] <- bed2$end[isHigh & bed2$chromosome != "zoom_chr3"] + 20000000
  
  # Manually define the ON-TARGET
  if(!is.null(ots.bed)){
    ots.df <- read.delim(ots.bed, header = FALSE)
    ots.df <- ots.df[1, , drop = FALSE]
    ots.df[,1] <- gsub("zoom_", "", ots.df[,1])
    ots.df[,2] <- ots.df[,2] + round((ots.df[,3] - ots.df[,2]) /2)
    ots.df[,3] <- ots.df[,2]
    
    ots.chr <- paste0("zoom_", ots.df[,1])
    ots.start <- ots.df[,2] - 1000
    ots.end <- ots.df[,3] + 1000
    
    bed1 <- data.frame(chr = rep(ots.chr, nrow(bed2)),
                       start = rep(ots.start, nrow(bed2)),# 46359961
                       end = rep(ots.end, nrow(bed2)),# 46380781
                       value1 = 0)
  }else{# use middle coordinates
    ots.idx <- which.max(siteM.sub$read)
    bed1 <- data.frame(chr = rep(siteM.sub$chromosome[ots.idx], nrow(bed2)),
                       start = rep(siteM.sub$middleCoord[ots.idx] - 1000, nrow(bed2)),# 46359961
                       end = rep(siteM.sub$middleCoord[ots.idx] + 1000, nrow(bed2)),# 46380781
                       value1 = 0)
  }
  
  
  link.color <- rep(rgb(190, 190, 190, max = 255, alpha = 175), length(group))# grey
  link.color[group == "ON"] <- rgb(0, 255, 0, max = 255, alpha = 225)# green
  link.color[group == "OMT"] <- rgb(255, 0, 0, max = 255, alpha = 175)# red
  link.color[group == "HMT"] <- rgb(0, 0, 255, max = 255, alpha = 175)# blue
  link.color[group == "OMT/HMT"] <- rgb(255, 193, 37, max = 255, alpha = 175)# goldenrod1
  
  # ON to ON
  idx <- which.max(siteM.sub$read)# ON TARGET IDX
  bed2.sub <- bed2[idx, ]# SELECT ON
  bed1.sub <- bed1[idx, ]
  
  bed1.sub$end <- bed1.sub$start + 2000
  bed2.sub$start <- bed2.sub$end - round(((bed2.sub$end - bed2.sub$start) / 2.5))
  link.color.sub = link.color[idx]
  
  # ON to OFF
  bed2 <- bed2[-idx, ]# REMOVE ON
  bed1 <- bed1[-idx, ]
  link.color = link.color[-idx]
  group <- group[-idx]
  
  circos.genomicLink(bed1[group != "NBS",], bed2[group != "NBS",], col = link.color[group != "NBS"], border = NA, h.ratio = 0.5)
  circos.genomicLink(bed1.sub, bed1.sub, col = link.color.sub, border = NA, h.ratio = 0.1)# ON -ON
  #circos.genomicLink(bed1.sub, bed2.sub, col = link.color.sub, border = NA, h.ratio = 0.1)# ON -ON
  
  dev.off()
  
}





############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE){
  
  library(circlize)
  library(openxlsx)
  library(GenomicRanges)
  
  # TEST FUNCTION
  siteFile <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/MMUL_OVL/MMUL_OVL_aln_stat_FLANK_GROUP_GENES.xlsx")
  #zoom.chr <- "chr3"
  #zoom.start <- 46350000
  #zoom.end <- 46390000
  
  zoom.size <- 25000
  label = FALSE
  
  bestScore.cutoff <- 9  
  bestFlank.cutoff <- 25
  gene.bed <- file.path("~/Research/CASTSeq/circlize/genes/CD33_rheMac8.bed")
  ots.bed <- NULL
  
  cytobandFile <- file.path("~/Research/CASTSeq/pipelineGit/annotations/Macaca_mulatta/cytoBandIdeo_rheMac8.txt")
  
  outFile <- "~/tmp/test_circlize.pdf"
  circlizePipeline(siteFile, zoom.size, label, 
                   bestScore.cutoff = 9, bestFlank.cutoff = 25,
                   gene.bed = gene.bed, ots.bed = NULL, 
                   realigned = FALSE,
                   outFile,
                   cytobandFile = cytobandFile)
  
  
  
  
  # WITH OTS
  ots.bed <- file.path("~/cluster/cluster/offTargets/Giando/pipelineGit/samples/G3_WT_D1/data/ots.bed")
  outFile <- "home/gandrieux/tmp/test_circlize_OTS.pdf"
  circlizePipeline(siteFile, zoom.size, label, 
                   bestScore.cutoff = 9, bestFlank.cutoff = 25,
                   gene.bed = NULL, ots.bed, 
                   realigned = FALSE,
                   outFile,
                   species = "hg38")
  
  ##########
  # G3_WT_D1  
  
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples",
                        "G3_WT_D1/results/guide_aln/G3-WT-d1_S1_L001_w250_FINAL.xlsx")
  ots.bed <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples",
                       "G3_WT_D1/data/ots.bed")
  gene.bed <- file.path("~/Research/CASTSeq/circlize/genes/CCR5_CRR5.bed")
  outFile <- "~/tmp/test_circlize.pdf"
  
  zoom.size <- 25000
  label = FALSE 
  
  circlizePipeline(siteFile, zoom.size, label, 
                   bestScore.cutoff = 9, bestFlank.cutoff = 25,
                   gene.bed, ots.bed, 
                   realigned = FALSE,
                   outFile,
                   species = "hg38")
  
  
  
}
