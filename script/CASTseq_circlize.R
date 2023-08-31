
set.seed(1234)

start_time <- Sys.time()


#####################################################
###                                               ###
###         CAST-SEQ CIRCLIZE STAND ALONE         ###
###                                               ###
#####################################################




##########################################################################################
################                     PARSE                              ################## 
##########################################################################################
require(optparse)

option_list = list(
  make_option(c("--siteFile"), type="character", default=NULL, 
              help="Path of the CASTseq XXX_FINAL.xlsx file", metavar="character"),
  make_option(c("--zoom"), type="integer", default=25000, 
              help="Size of the blow-up area (if 0, then no blow-up) [default= %default]", metavar="integer"),   
  make_option(c("--label"), type="logical", default=FALSE, 
              help="show (or not) the labels for OMT and HMT sites [default= %default]", metavar="logical"),   
  make_option(c("--pv"), type="double", default=NULL, 
              help="(optional) pvalue threhold [default= %default]", metavar="double"),
  make_option(c("--bestScore"), type="double", default=NULL, 
              help="(optional) threhold characterizing the OMT. Otherwise it is set as the minimum OMTs score. [default= %default]", metavar="double"),
  make_option(c("--bestFlank"), type="double", default=25, 
              help="Threhold characterizing the HMT. [default= %default]", metavar="double"),
  make_option(c("--NBS"), type="logical", default=TRUE, 
              help="show (or not) the NBS. [default= %default]", metavar="logical"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="(optional) path of the bed file containing the gene coordinate for the blow-up. If NULL, use the coordinates from siteFile [default= %default]", metavar="character"),
  make_option(c("--ots"), type="character", default=NULL, 
              help="path of the ON-target bed file containing the ON-target coordinate (tab separated file with chromosome start end, no header)[default= %default]", metavar="character"),
  make_option(c("--realigned"), type="logical", default=TRUE, 
              help="realigned (or not) the ON-target coordinate (for aesthetic purpose) [default= %default]", metavar="logical"),
  make_option(c("--outFile"), type="character", default=NULL, 
              help="Path of the output circos plot (.pdf)", metavar="character"),
  make_option(c("--species"), type="character", default="hg38", 
              help="hg38 or mm10 are the only supported species [default= %default]", metavar="character"),
  make_option(c("--nbHits"), type="integer", default=NULL, 
              help="(optional) display sites with at least this number of hits [default= %default]", metavar="integer"),
  make_option(c("--top"), type="double", default=NULL, 
              help="(optional) display only the top XXX sites [default= %default]", metavar="integer")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$siteFile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (siteFile).n", call.=FALSE)
}

if(is.null(opt$ots)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (ots).n", call.=FALSE)
}

if(is.null(opt$outFile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (outFile).n", call.=FALSE)
}


##########################################################################################
################                   REQUIRED LIBRARIES                 ####################
##########################################################################################

require(circlize)
require(openxlsx)
require(GenomicRanges)

##########################################################################################
################                       FUNCTIONS                      ####################
##########################################################################################


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
  
  # merge zoom cytobands
  if(nrow(zoom_bed)>1){
    zoom_bed[[3]][1] <- end
    zoom_bed <- zoom_bed[1, ,drop = FALSE]
  }
  
  rbind(bed, zoom_bed)
}

extend_bed = function(bed, chromosome, start, end, prefix = "zoom_") {
  zoom.gr <- makeGRangesFromDataFrame(data.frame(chromosome = chromosome, start = start, end = end),
                                      seqnames.field = "chromosome", start.field = "start", end.field = "end",
                                      keep.extra.columns = FALSE, ignore.strand = TRUE)
  bed.gr <- makeGRangesFromDataFrame(bed,
                                     seqnames.field = "chromosome", start.field = "start", end.field = "end",
                                     keep.extra.columns = FALSE, ignore.strand = TRUE)
  
  gr.ovl <- findOverlaps(query = zoom.gr, subject = bed.gr, type = "any", maxgap = 0)
  
  zoom_bed = bed
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  
  idx <- subjectHits(gr.ovl)
  zoom_bed = zoom_bed[idx, ]
  
  rbind(bed, zoom_bed)
}


circlizePipeline <- function(siteFile, zoom.size = 25000, label = FALSE,
                             PV.cutoff = NULL,
                             bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                             showNBS = TRUE,
                             gene.bed = NULL, ots.bed = NULL, 
                             realigned = FALSE,
                             outFile,
                             species = "hg38",
                             nbHits = NULL,
                             top = NULL){
  # READ INPUT FILE
  siteM <- read.xlsx(siteFile)
  siteM <- siteM[siteM$chromosome %in% paste0("chr", c(1:22, "X", "Y")), ]
  
  print(paste0("siteM: ", nrow(siteM), " sites"))
  
  # SHOW NBS (TRUE / FALSE)
  if(!(showNBS)){
    siteM <- siteM[siteM[, "group"] != "NBS", ]
    print(paste0("siteM: ", nrow(siteM), " NO NBS sites"))
  }
  
  if(!is.null(nbHits)){
    # ONLY SHOW THE TOP XXX SITES
    keep <- siteM$hits >= nbHits
    if(sum(keep) == 0){
      
      keep <- rep(TRUE, nrow(siteM))
    }
    siteM <- siteM[keep, , drop = FALSE]
    print(paste0("siteM: ", nrow(siteM), " sites with at least ", nbHits, " hits."))
  }
  
  if(!is.null(top)){
    # ONLY SHOW THE TOP XXX SITES
    top <- min(c(top, nrow(siteM)))
    siteM <- siteM[1:top, , drop = FALSE]
    print(paste0("siteM: ", nrow(siteM), " top sites"))
  }
  
  siteM.chr <- factor(gsub("chr", "", siteM$chromosome), levels = c(1:22, "X", "Y"))
  siteM <- siteM[order(siteM.chr, siteM$start), ]
  
  # SELECT SIGNIFICANT (DEPRECATED)
  siteM$is.ON[is.na(siteM$is.ON)] <- "no"
  if(!is.null(PV.cutoff)) siteM <- siteM[siteM$adj.pvalue < PV.cutoff | siteM$is.ON == "yes", ]# keep signif and/or ON target
  
  # CHECK IF ON TARGET IS PRESENT
  onOverlap <- sum(siteM$is.ON == "yes")!=0
  
  # IF NOT ADD A FAKE ON TARGET SITE FOR AESTHETIC
  if(!onOverlap){
    print("no ON target site detected, use best OMT instead (MAX READS)")
    siteM$is.ON[which.max(siteM$hits)] <- "yes"
  }
  
  # FIX TALEN-SPECIFIC COLUMNS
  # score
  if(!("score" %in% colnames(siteM))){
    siteM$score <- siteM$top1_score
  }
  if(!("middleCoord" %in% colnames(siteM))){
    siteM$middleCoord <- siteM$top1_middleCoord
  }
  
  # FIX OLD TALEN-SPECIFIC COLUMNS
  # score
  if(is.null(siteM$score)){
    siteM$score <- sapply(1:nrow(siteM), function(i){
      if(siteM$group[i] != "OMT") return(5)
      bestCB <- siteM$BestCB[i]
      return(siteM[i, paste0(bestCB, "_score")])
    })
  }
  if(is.null(siteM$middleCoord)){
    siteM$middleCoord <- siteM$LF.LR_middleCoord
  }
  
  
  # CHANGE MIDDLE COORDINATE OF THE ON TARGET
  if(!is.null(ots.bed) & realigned & onOverlap){
    ots.df <- read.delim(ots.bed, header = FALSE)
    ots.df <- ots.df[1,,drop = FALSE]# only take the first ON target into account
    
    on.middle <- ots.df[,2] + round((ots.df[,3] - ots.df[,2]) / 2)
    
    siteM$middleCoord[siteM$is.ON == "yes"] <- on.middle
  }
  
  # DEFINE ZOOM
  idx <- siteM$is.ON == "yes"
  zoom.chr <- siteM$chromosome[idx]
  zoom.start <- siteM$middleCoord[idx] - zoom.size
  zoom.end <- siteM$middleCoord[idx] + zoom.size
  
  # ADD ZOOM
  siteM <- extend_bed(siteM, zoom.chr, start = zoom.start, end = zoom.end)
  bestScore.raw <- siteM$score
  
  # INIT (MUST BE ONLINE !!!)
  if(species == "hg38") cytoband = read.cytoband(species = species, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  if(species == "mm10") cytoband = read.cytoband(species = species, chromosome.index = paste0("chr", c(1:19, "X", "Y")))
  cytoband_df = cytoband$df
  chromosome = cytoband$chromosome
  
  new_cytoband_df <- extend_cytoband(cytoband_df, zoom.chr, start = zoom.start, end = zoom.end)
  
  # DEFINE ZOOM LIMITS
  zoom_lower_limit <- min(new_cytoband_df$V2[grepl("zoom_", new_cytoband_df$V1)])
  zoom_upper_limit <- max(new_cytoband_df$V3[grepl("zoom_", new_cytoband_df$V1)])
  
  xrange = c(cytoband$chr.len, cytoband$chr.len[zoom.chr])
  
  if(species == "hg38"){
    normal_chr_index = 1:24
    zoomed_chr_index = 25
    
    sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                     xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index]))
    sector.width[1:24] <- sector.width[1:24] 
    sector.width[25] <- sector.width[25] / 1.25
  }
  if(species == "mm10"){
    normal_chr_index = 1:21
    zoomed_chr_index = 22
    
    sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                     xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index]))
    sector.width[1:21] <- sector.width[1:21] 
    sector.width[22] <- sector.width[22] / 1.25
  }
  
  pdf(outFile, width = 3.5, height = 3.5)
  circos.clear()
  if(species == "hg38") circos.par(start.degree = 90, "gap.degree" = c(rep(3, 23), 10, 20))
  if(species == "mm10") circos.par(start.degree = 90, "gap.degree" = c(rep(3, 20), 10, 20))
  
  
  circos.initializeWithIdeogram(new_cytoband_df, 
                                sector.width = sector.width,
                                plotType = NULL)
  
  ##############
  # GROUP LABELS
  zoom.idx <- grep("^zoom_", siteM$chromosome)
  zoom.coord <- sapply(zoom.idx, function(i) paste(siteM[i, c("start", "end")], collapse = "-"))
  all.coord <- sapply(1:nrow(siteM), function(i) paste(siteM[i, c("start", "end")], collapse = "-"))
  toRemove <- setdiff(which(all.coord %in% zoom.coord), zoom.idx)
  
  # ADJUST ZOOM LIMITS
  siteM$start[grepl("zoom_", siteM$chromosome) & siteM$start < zoom_lower_limit] <- zoom_lower_limit
  siteM$end[grepl("zoom_", siteM$chromosome) & siteM$end > zoom_upper_limit] <- zoom_upper_limit
  siteM$middleCoord[grepl("zoom_", siteM$chromosome) & siteM$middleCoord < zoom_lower_limit] <- zoom_lower_limit
  siteM$middleCoord[grepl("zoom_", siteM$chromosome) & siteM$middleCoord > zoom_upper_limit] <- zoom_upper_limit
  
  
  
  # REMOVE DUPLICATED SITES
  siteM.sub <- siteM[-toRemove, ]# remove duplicated sites, keep only the ones in zoom area
  group <- siteM.sub[, "group"]
  group[group == "OMT" & siteM.sub$is.HMT == "yes"] <- "OMT/HMT"
  
  # ON target
  #group[siteM.sub$is.ON == "yes"] <- "ON"
  group[siteM.sub$is.ON == "yes" & (group == "OMT" | group == "OMT/HMT")] <- "ON"
  
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
  bestScore <- siteM$score
  bestScore.chr <- siteM[, "chromosome"]
  bestScore.x <- siteM$middleCoord
  bestScore.OMT <- ifelse(siteM$group == "OMT", TRUE, FALSE)
  
  bestScore.cutoff <- floor(min(bestScore[siteM$group == "OMT"]))
  if(is.null(bestScore.cutoff) | bestScore.cutoff == Inf) bestScore.cutoff <- max(bestScore)
  if(nrow(siteM)==2) bestScore.cutoff <- bestScore.cutoff * 0.75# fix a visual bug when only one site
  
  # fill missing chr
  missing.chr <- setdiff(unique(new_cytoband_df[[1]]), bestScore.chr)
  bestScore <- c(bestScore, rep(0, length(missing.chr)))
  bestScore.chr <- c(bestScore.chr, missing.chr)
  bestScore.x  <- c(bestScore.x, rep(0, length(missing.chr)))
  bestScore.OMT <- c(bestScore.OMT, rep(FALSE, length(missing.chr)))
  
  #
  bestScore[bestScore >= bestScore.cutoff & !bestScore.OMT] <- bestScore.cutoff -1 
  bestScore[bestScore < bestScore.cutoff & bestScore.OMT] <- bestScore.cutoff
  
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
                 circos.rect(xlim[1], bestScore.cutoff, xlim[2], (ceiling(ylim[2])+ceiling(ylim[2])*0.2), col = "#FF000020", border = NA)
                 
                 circos.points(x[y < bestScore.cutoff & y !=0], y[y < bestScore.cutoff & y !=0], pch = 16, cex = 0.75, col = "grey")
                 circos.points(x[y >= bestScore.cutoff & y !=0], y[y >= bestScore.cutoff & y !=0], pch = 16, cex = 0.75, col = "red")
               })
  
  circos.yaxis(side = "left", at = c(bestScore.cutoff, floor(max(bestScore))), #floor(seq(0, max(bestScore), length.out = 3))
               sector.index = get.all.sector.index()[1], labels.cex = 0.5,
               tick.length = convert_x(0.75, "mm", get.all.sector.index()[1],get.cell.meta.data("track.index")),
               lwd = par("lwd") * 1.25)
  
  # BEST SCORE (ZOOM ONLY)
  idx <- siteM$chromosome == paste0("zoom_", zoom.chr)
  bestScore.zoom <- bestScore.raw[idx]
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
  
  bestFlank.x <- siteM$middleCoord
  
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
                 
                 for(yi in floor(seq(0, max(bestFlank), length.out = 3))) {
                   circos.lines(cell.xlim, c(yi, yi), lty = 2, col = "white") 
                 }
                 
                 xlim = get.cell.meta.data("xlim")
                 ylim = get.cell.meta.data("ylim")
                 
                 print(xlim)
                 print(ylim)
                 circos.rect(xlim[1], bestFlank.cutoff, xlim[2], ceiling(ylim[2])+5, col = "cornflowerblue", border = NA)
                 
                 circos.points(x[y < bestFlank.cutoff & y !=0], y[y < bestFlank.cutoff & y !=0], pch = 16, cex = 0.75, col = "grey")
                 circos.points(x[y >= bestFlank.cutoff], y[y >= bestFlank.cutoff], pch = 16, cex = 0.75, col = "blue")
                 
               },  track.margin = c(0,0))
  
  circos.yaxis(side = "left", at = c(0, bestFlank.cutoff), #floor(seq(0, max(bestFlank), length.out = 3))
               sector.index = get.all.sector.index()[1], labels.cex = 0.5,
               tick.length = convert_x(0.75, "mm", get.all.sector.index()[1], get.cell.meta.data("track.index")),
               lwd = par("lwd") * 1.25)
  
  #################
  # GENE ANNOTATION
  plotGene <- FALSE
  if(!is.null(gene.bed)){
    geneM <- read.delim(gene.bed, header = FALSE)# must be 4 columns: chr, start, end, symbol
    geneM[,1] <- paste0("zoom_", geneM[,1])
    plotGene <- TRUE
  }else if("yes" %in% siteM.sub$is.ON){# define gene based on ON-target coordinates
    geneM <- siteM.sub[siteM.sub$is.ON == "yes", c("chromosome", "geneStart", "geneEnd", "SYMBOL")]
    plotGene <- TRUE
  }
  
  if(plotGene){
    if(class(geneM[,2]) == "factor") geneM[,2] <- as.numeric(levels(geneM[,2]))[geneM[,2]]
    if(class(geneM[,3]) == "factor") geneM[,3] <- as.numeric(levels(geneM[,3]))[geneM[,3]]
    if(class(geneM[,4]) == "factor") geneM[,4] <- as.character(geneM[,4])
    
    # adjust gene coordinates
    zoom.min <- new_cytoband_df[grepl("^zoom_", new_cytoband_df[,1]), "V2"]
    zoom.max <- new_cytoband_df[grepl("^zoom_", new_cytoband_df[,1]), "V3"]
    geneM$geneStart <- sapply(geneM$geneStart, function(x) max(c(x, zoom.min)))
    geneM$geneEnd <- sapply(geneM$geneEnd, function(x) min(c(x, zoom.max)))
    
    isMinusGene <- sapply(1:nrow(geneM), function(x) as.numeric(geneM$geneEnd[x]) < as.numeric(geneM$geneStart[x]))
    geneM$geneStartRAW <- geneM$geneStart
    geneM$geneEndRAW <- geneM$geneEnd
    geneM$geneStart[isMinusGene] <- geneM$geneEndRAW[isMinusGene]
    geneM$geneEnd[isMinusGene] <- geneM$geneStartRAW[isMinusGene]
    
    #if(geneM$geneEnd < geneM$geneStart) geneM$geneEnd <- geneM$geneStart
    
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
                     new_cytoband_df[subjectHits(gr.ovl)[1],],
                     col = "grey50",
                     border = NA)
  
  ################################
  # LINK ON TARGETS TO ALL TARGETS
  bed2 <- siteM.sub[, 1:4]
  #bed2$start <- bed2$start + round((bed2$end - bed2$start) / 2)
  bed2$start <- siteM.sub$middleCoord
  
  # Re-Arange layout
  mysize <- round((zoom.size * 2000) / 50000)
  bed2$end <- bed2$start + mysize
  bed2$start <- bed2$start - mysize
  
  bed2$start[bed2$chromosome != paste0("zoom_", zoom.chr)] <- bed2$start[bed2$chromosome != paste0("zoom_", zoom.chr)] - 30000000
  bed2$end[bed2$chromosome != paste0("zoom_", zoom.chr)] <- bed2$end[bed2$chromosome != paste0("zoom_", zoom.chr)] + 30000000
  
  # Manually define the ON-TARGET
  ots.idx <- which(group == "ON")
  if(length(ots.idx) == 0) ots.idx <- which(siteM.sub$is.ON == "yes")
  bed1 <- data.frame(chr = rep(siteM.sub$chromosome[ots.idx], nrow(bed2)),
                     start = rep(siteM.sub$middleCoord[ots.idx] - mysize, nrow(bed2)),# 46359961
                     end = rep(siteM.sub$middleCoord[ots.idx] + mysize, nrow(bed2)),# 46380781
                     value1 = 0)
  
  link.color <- rep(rgb(190, 190, 190, max = 255, alpha = 175), length(group))# grey
  link.color[group == "ON"] <- rgb(0, 255, 0, max = 255, alpha = 225)# green
  link.color[group == "OMT"] <- rgb(255, 0, 0, max = 255, alpha = 175)# red
  link.color[group == "HMT"] <- rgb(0, 0, 255, max = 255, alpha = 175)# blue
  link.color[group == "OMT/HMT"] <- rgb(255, 193, 37, max = 255, alpha = 175)# goldenrod1
  
  # ON to ON
  idx <- ots.idx
  bed1.sub <- bed1[idx, ]
  
  bed1.sub$start <- bed1.sub$start + round(bed1.sub$end - bed1.sub$start) 
  bed1.sub$start <- bed1.sub$start - (mysize*2)
  
  link.color.sub = link.color[idx]
  
  # ON to OFF
  bed2 <- bed2[-idx, ]# REMOVE ON
  bed1 <- bed1[-idx, ]
  link.color = link.color[-idx]
  
  circos.genomicLink(bed1, bed2, col = link.color, border = NA, h.ratio = 0.5)
  circos.genomicLink(bed1.sub, bed1.sub, col = link.color.sub, border = NA, h.ratio = 0.1)# ON -ON
  
  dev.off()
  
}

circlizePipelineNOZOOM <- function(siteFile, label = FALSE,
                                   PV.cutoff = NULL,
                                   bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                                   showNBS = TRUE,
                                   ots.bed = NULL, 
                                   realigned = FALSE,
                                   outFile,
                                   species = "hg38",
                                   nbHits = NULL,
                                   top = NULL){
  # READ INPUT FILE
  siteM <- read.xlsx(siteFile)
  siteM <- siteM[siteM$chromosome %in% paste0("chr", c(1:22, "X", "Y")), ]
  
  print(paste0("siteM: ", nrow(siteM), " sites"))
  
  # SHOW NBS (TRUE / FALSE)
  if(!(showNBS)){
    siteM <- siteM[siteM[, "group"] != "NBS", ]
    print(paste0("siteM: ", nrow(siteM), " NO NBS sites"))
  }
  
  if(!is.null(nbHits)){
    # ONLY SHOW THE TOP XXX SITES
    keep <- siteM$hits >= nbHits
    if(sum(keep) == 0){
      
      keep <- rep(TRUE, nrow(siteM))
    }
    siteM <- siteM[keep, , drop = FALSE]
    print(paste0("siteM: ", nrow(siteM), " sites with at least ", nbHits, " hits."))
  }
  
  if(!is.null(top)){
    # ONLY SHOW THE TOP XXX SITES
    top <- min(c(top, nrow(siteM)))
    siteM <- siteM[1:top, , drop = FALSE]
    print(paste0("siteM: ", nrow(siteM), " top sites"))
  }
  
  siteM.chr <- factor(gsub("chr", "", siteM$chromosome), levels = c(1:22, "X", "Y"))
  siteM <- siteM[order(siteM.chr, siteM$start), ]
  
  # SELECT SIGNIFICANT (DEPRECATED)
  siteM$is.ON[is.na(siteM$is.ON)] <- "no"
  if(!is.null(PV.cutoff)) siteM <- siteM[siteM$adj.pvalue < PV.cutoff | siteM$is.ON == "yes", ]# keep signif and/or ON target
  
  # CHECK IF ON TARGET IS PRESENT
  onOverlap <- sum(siteM$is.ON == "yes")!=0
  
  # IF NOT ADD A FAKE ON TARGET SITE FOR AESTHETIC
  if(!onOverlap){
    print("no ON target site detected, use best OMT instead (MAX READS)")
    siteM$is.ON[which.max(siteM$hits)] <- "yes"
  }
  
  # FIX TALEN-SPECIFIC COLUMNS
  # score
  if(!("score" %in% colnames(siteM))){
    siteM$score <- siteM$top1_score
  }
  if(!("middleCoord" %in% colnames(siteM))){
    siteM$middleCoord <- siteM$top1_middleCoord
  }
  
  # FIX OLD TALEN-SPECIFIC COLUMNS
  # score
  if(is.null(siteM$score)){
    siteM$score <- sapply(1:nrow(siteM), function(i){
      if(siteM$group[i] != "OMT") return(5)
      bestCB <- siteM$BestCB[i]
      return(siteM[i, paste0(bestCB, "_score")])
    })
  }
  if(is.null(siteM$middleCoord)){
    siteM$middleCoord <- siteM$LF.LR_middleCoord
  }
  
  # CHANGE MIDDLE COORDINATE OF THE ON TARGET
  if(!is.null(ots.bed) & realigned & onOverlap){
    ots.df <- read.delim(ots.bed, header = FALSE)
    ots.df <- ots.df[1,,drop = FALSE]# only take the first ON target into account
    
    on.middle <- ots.df[,2] + round((ots.df[,3] - ots.df[,2]) / 2)
    
    siteM$middleCoord[siteM$is.ON == "yes"] <- on.middle
  }
  
  
  # INIT (MUST BE ONLINE !!!)
  if(species == "hg38") cytoband = read.cytoband(species = species, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  if(species == "mm10") cytoband = read.cytoband(species = species, chromosome.index = paste0("chr", c(1:19, "X", "Y")))
  cytoband_df = cytoband$df
  chromosome = cytoband$chromosome
  
  xrange = cytoband$chr.len
  if(species == "hg38"){
    normal_chr_index = 1:24
    
    sector.width = xrange[normal_chr_index] / sum(xrange[normal_chr_index])
  }
  if(species == "mm10"){
    normal_chr_index = 1:21
    
    sector.width = xrange[normal_chr_index] / sum(xrange[normal_chr_index])
  }
  
  pdf(outFile, width = 3.5, height = 3.5)
  circos.clear()
  if(species == "hg38") circos.par(start.degree = 90, "gap.degree" = c(rep(3, 23), 12))
  if(species == "mm10") circos.par(start.degree = 90, "gap.degree" = c(rep(3, 20), 12))
  
  circos.initializeWithIdeogram(cytoband_df, 
                                sector.width = sector.width,
                                plotType = NULL)
  
  ##############
  # GROUP LABELS
  siteM.sub <- siteM
  group <- siteM.sub[, "group"]
  group[group == "OMT" & siteM.sub$is.HMT == "yes"] <- "OMT/HMT"
  
  # ON target
  group[siteM.sub$is.ON == "yes" & (group == "OMT" | group == "OMT/HMT")] <- "ON"
  
  
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
  circos.genomicIdeogram(cytoband_df,  track.margin = c(0,0.015))
  
  ##################
  # BEST SCORE TRACK
  bestScore <- siteM[, "score"]
  bestScore.chr <- siteM[, "chromosome"]
  bestScore.x <- siteM$middleCoord
  bestScore.OMT <- ifelse(siteM$group == "OMT", TRUE, FALSE)
  
  bestScore.cutoff <- floor(min(bestScore[siteM$group == "OMT"]))
  if(is.null(bestScore.cutoff) | bestScore.cutoff == Inf) bestScore.cutoff <- max(bestScore)
  if(nrow(siteM)==2) bestScore.cutoff <- bestScore.cutoff * 0.75# fix a visual bug when only one site
  
  # fill missing chr
  missing.chr <- setdiff(unique(cytoband_df[[1]]), bestScore.chr)
  bestScore <- c(bestScore, rep(0, length(missing.chr)))
  bestScore.chr <- c(bestScore.chr, missing.chr)
  bestScore.x  <- c(bestScore.x, rep(0, length(missing.chr)))
  bestScore.OMT <- c(bestScore.OMT, rep(FALSE, length(missing.chr)))
  
  #
  bestScore[bestScore >= bestScore.cutoff & !bestScore.OMT] <- bestScore.cutoff -1 
  bestScore[bestScore < bestScore.cutoff & bestScore.OMT] <- bestScore.cutoff
  
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
                 
                 circos.points(x[y < bestScore.cutoff & y !=0], y[y < bestScore.cutoff & y !=0], pch = 16, cex = 0.75, col = "grey")
                 circos.points(x[y >= bestScore.cutoff], y[y >= bestScore.cutoff], pch = 16, cex = 0.75, col = "red")
               })
  
  circos.yaxis(side = "left", at = c(bestScore.cutoff, floor(max(bestScore))), #floor(seq(0, max(bestScore), length.out = 3))
               sector.index = get.all.sector.index()[1], labels.cex = 0.5,
               tick.length = convert_x(0.75, "mm", get.all.sector.index()[1],get.cell.meta.data("track.index")),
               lwd = par("lwd") * 1.25)
  
  ##################
  # BEST FLANK TRACK
  bestFlank <- apply(as.matrix(siteM[, c("flank.length", "flank.rev.length")]), 1, max)
  bestFlank[bestFlank > 30] <- 30
  bestFlank.chr <- siteM[, "chromosome"]
  #bestFlank.cutoff <- 25
  bestFlank.x <- siteM$middleCoord
  #bestFlank.x <-siteM$start + round((siteM$end - siteM$start) / 2)
  
  # fill missing chr
  missing.chr <- setdiff(unique(cytoband_df[[1]]), bestFlank.chr)
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
                 
                 circos.points(x[y < bestFlank.cutoff & y !=0], y[y < bestFlank.cutoff & y !=0], pch = 16, cex = 0.75, col = "grey")
                 circos.points(x[y >= bestFlank.cutoff], y[y >= bestFlank.cutoff], pch = 16, cex = 0.75, col = "blue")
                 
               },  track.margin = c(0,0))
  
  circos.yaxis(side = "left", at = c(0, bestFlank.cutoff), #floor(seq(0, max(bestFlank), length.out = 3))
               sector.index = get.all.sector.index()[1], labels.cex = 0.5,
               tick.length = convert_x(0.75, "mm", get.all.sector.index()[1], get.cell.meta.data("track.index")),
               lwd = par("lwd") * 1.25)
  
  ################################
  # LINK ON TARGETS TO ALL TARGETS
  bed2 <- siteM.sub[, 1:4]
  #bed2$start <- bed2$start + round((bed2$end - bed2$start) / 2)
  bed2$start <- siteM.sub$middleCoord
  
  # Re-Arange layout
  mysize <- 30000000
  bed2$end <- bed2$start + mysize
  bed2$start <- bed2$start - mysize
  
  # Manually define the ON-TARGET
  if(!is.null(ots.bed)){
    ots.idx <- which(group == "ON")
    if(length(ots.idx) == 0) ots.idx <- which(siteM.sub$is.ON == "yes")
    bed1 <- data.frame(chr = rep(siteM.sub$chromosome[ots.idx], nrow(bed2)),
                       start = rep(siteM.sub$middleCoord[ots.idx] - mysize, nrow(bed2)),# 46359961
                       end = rep(siteM.sub$middleCoord[ots.idx] + mysize, nrow(bed2)),# 46380781
                       value1 = 0)
  }else{# use middle coordinates
    ots.idx <- which.max(siteM.sub$score)
    bed1 <- data.frame(chr = rep(siteM.sub$chromosome[ots.idx], nrow(bed2)),
                       start = rep(siteM.sub$middleCoord[ots.idx] - mysize, nrow(bed2)),# 46359961
                       end = rep(siteM.sub$middleCoord[ots.idx] + mysize, nrow(bed2)),# 46380781
                       value1 = 0)
  }
  
  
  link.color <- rep(rgb(190, 190, 190, max = 255, alpha = 175), length(group))# grey
  link.color[group == "ON"] <- rgb(0, 255, 0, max = 255, alpha = 225)# green
  link.color[group == "OMT"] <- rgb(255, 0, 0, max = 255, alpha = 175)# red
  link.color[group == "HMT"] <- rgb(0, 0, 255, max = 255, alpha = 175)# blue
  link.color[group == "OMT/HMT"] <- rgb(255, 193, 37, max = 255, alpha = 175)# goldenrod1
  
  # ON to ON
  #idx <- which.max(siteM.sub$score)# ON TARGET IDX
  idx <- ots.idx
  #bed2.sub <- bed2[idx, ]# SELECT ON
  bed1.sub <- bed1[idx, ]
  
  bed1.sub$start <- bed1.sub$start + round(bed1.sub$end - bed1.sub$start) 
  bed1.sub$start <- bed1.sub$start - (mysize*2)
  #bed1.sub$end <- bed1.sub$end + mysize
  #bed1.sub$end <- bed1.sub$start + 4000
  
  #bed1.sub$end <- bed1.sub$start + round(((bed1.sub$end - bed1.sub$start) / 2.5))
  #bed2.sub$start <- bed2.sub$end - round(((bed2.sub$end - bed2.sub$start) / 2.5))
  link.color.sub = link.color[idx]
  
  # ON to OFF
  bed2 <- bed2[-idx, ]# REMOVE ON
  bed1 <- bed1[-idx, ]
  link.color = link.color[-idx]
  
  circos.genomicLink(bed1, bed2, col = link.color, border = NA, h.ratio = 0.5)
  circos.genomicLink(bed1.sub, bed1.sub, col = link.color.sub, border = NA, h.ratio = 0.1)# ON -ON
  #circos.genomicLink(bed1.sub, bed2.sub, col = link.color.sub, border = NA, h.ratio = 0.1)# ON -ON
  
  dev.off()
  
}

##########################################################################################
################                         MAIN                         ####################
##########################################################################################

# RUN
siteFile = opt$siteFile
zoom.size = opt$zoom
label = opt$label
PV.cutoff = opt$pv
bestScore.cutoff = opt$bestScore
bestFlank.cutoff = opt$bestFlank
showNBS = opt$NBS
gene.bed = opt$gene
ots.bed = opt$ots
realigned = opt$realigned
outFile = opt$outFile
species = opt$species
nbHits = opt$nbHits
top = opt$top



if(zoom.size != 0){
  circlizePipeline(siteFile = siteFile,
                   zoom.size = zoom.size,
                   label = label,
                   PV.cutoff = PV.cutoff,
                   bestScore.cutoff = bestScore.cutoff,
                   bestFlank.cutoff = bestFlank.cutoff,
                   showNBS = showNBS,
                   gene.bed = gene.bed,
                   ots.bed = ots.bed,
                   realigned = realigned,
                   outFile = outFile,
                   species = species,
                   nbHits = nbHits,
                   top = top
  )
}else{
  circlizePipelineNOZOOM(siteFile,
                         label = label,
                         PV.cutoff = PV.cutoff,
                         bestScore.cutoff = bestScore.cutoff,
                         bestFlank.cutoff = bestFlank.cutoff,
                         showNBS = showNBS,
                         ots.bed = ots.bed, 
                         realigned = realigned,
                         outFile = outFile,
                         species = species,
                         nbHits = nbHits,
                         top = top
  )
}


end_time <- Sys.time()
end_time - start_time




if(FALSE){
  # HO TO RUN IT
  #Rscript CASTseq_circlize.R --siteFile "~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/GENEWIZ_90-737856818/Masako/ELANDE-KO_T3T4_OVL1/ELANDE-KO_T3T4_OVL1_FINAL.xlsx" \
  # --ots "~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/Masako/ELANE_KO/ELANE-KO-Treat3/data/ots.bed" \
  # --outFile "~/Research/CASTSeq/circlize/Masako_300822/ELANDE-KO_T3T4_OVL1/ELANDE-KO_T3T4_OVL1_circlize_25k_TEST.pdf"
  
  
  # DO NOT RUN
  library(circlize)
  library(openxlsx)
  library(GenomicRanges)
  
  
  # TEST CYBB (June 23)
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-887623654/CYBB/CYBB_crispr/Cpf1-CYBB/OVL2_SIGNIF2/Cpf1-CYBB_FINAL.xlsx")
  zoom.size <- 25000
  label <- FALSE
  PV.cutoff <- NULL
  bestScore.cutoff <- NULL
  bestFlank.cutoff <- 25
  showNBS <- TRUE
  gene.bed <- NULL
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-887623654/CYBB/CYBB_crispr/Cpf1-CYBB/data/ots.bed"
  realigned <- TRUE
  outFile <- file.path("~/Research/CASTSeq/Test/circlize/CYBB_FINAL_circlize_25k.pdf")
  species <- "hg38"
  nbHits <- NULL
  top <- NULL
  
  
  # TEST CD3E KI (June 23)
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-873723755/CD3E/CD3E_crispr/CD3E_KI/OVL2_SIGNIF2/CD3E_KI_FINAL_fake.xlsx")
  zoom.size <- 25000
  label <- FALSE
  PV.cutoff <- NULL
  bestScore.cutoff <- NULL
  bestFlank.cutoff <- 25
  showNBS <- TRUE
  gene.bed <- NULL
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-873723755/CD3E/CD3E_crispr/CD3E_KI/data/ots.bed"
  realigned <- TRUE
  outFile <- file.path("~/Research/CASTSeq/Test/circlize/CD3E_KI_FINAL_circlize_25k.pdf")
  species <- "hg38"
  nbHits <- NULL
  top <- NULL
  
  
  # TEST FANCF (NO ON OMT ON TARGET SITE)
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-825768046/FANCF/FANCF_crispr/FANCF-ABE/OVL2_SIGNIF2/FANCF-ABE_FINAL.xlsx")
  zoom.size <- 25000
  label <- FALSE
  PV.cutoff <- NULL
  bestScore.cutoff <- NULL
  bestFlank.cutoff <- 25
  showNBS <- TRUE
  gene.bed <- NULL
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-825768046/FANCF/FANCF_crispr/FANCF-ABE/data/ots.bed"
  realigned <- TRUE
  outFile <- file.path("~/Research/CASTSeq/Test/circlize/FANCF-ABE_FINAL_circlize_25k.pdf")
  species <- "hg38"
  nbHits <- NULL
  top <- NULL
  
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-830103735/EMX1/EMX1_crispr/EMX1-ABE/OVL2_SIGNIF2/EMX1-ABE_FINAL.xlsx")
  zoom.size <- 25000
  label <- FALSE
  PV.cutoff <- NULL
  bestScore.cutoff <- NULL
  bestFlank.cutoff <- 25
  showNBS <- TRUE
  gene.bed <- NULL
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples_V2/GENEWIZ_90-830103735/EMX1/EMX1_crispr/EMX1-ABE/data/ots.bed"
  realigned <- TRUE
  outFile <- file.path("~/Research/CASTSeq/Test/circlize/EMX1-ABE_FINAL_circlize_25k.pdf")
  species <- "hg38"
  nbHits <- NULL
  top <- NULL
  
  circlizePipeline(siteFile = siteFile,
                   zoom.size = zoom.size,
                   label = label,
                   PV.cutoff = PV.cutoff,
                   bestScore.cutoff = bestScore.cutoff,
                   bestFlank.cutoff = bestFlank.cutoff,
                   showNBS = showNBS,
                   gene.bed = gene.bed,
                   ots.bed = ots.bed,
                   realigned = realigned,
                   outFile = outFile,
                   species = species,
                   nbHits = nbHits,
                   top = top
  )
  
  
  # TEST OLD TALEN 04.11.22
  siteFile <- file.path("~/Research/CASTSeq/circlize/test_TALEN_old/CCR5-WT_OVL2_FINAL.xlsx")
  zoom.size <- 25000
  label <- FALSE
  PV.cutoff <- NULL
  bestScore.cutoff <- NULL
  bestFlank.cutoff <- 25
  showNBS <- TRUE
  gene.bed <- NULL
  ots.bed <- "~/Research/CASTSeq/circlize/test_TALEN_old/ots.bed"
  realigned <- TRUE
  outFile <- file.path("~/Research/CASTSeq/circlize/test_TALEN_old/CCR5-WT_OVL2_FINAL_circlize_25k.pdf")
  species <- "hg38"
  nbHits <- NULL
  top <- NULL
  
  circlizePipeline(siteFile = siteFile,
                   zoom.size = zoom.size,
                   label = label,
                   PV.cutoff = PV.cutoff,
                   bestScore.cutoff = bestScore.cutoff,
                   bestFlank.cutoff = bestFlank.cutoff,
                   showNBS = showNBS,
                   gene.bed = gene.bed,
                   ots.bed = ots.bed,
                   realigned = realigned,
                   outFile = outFile,
                   species = species,
                   nbHits = nbHits,
                   top = top
  )
  
  outFile <- file.path("~/Research/CASTSeq/circlize/test_TALEN_old/CCR5-WT_OVL2_FINAL_circlize_NOZOOM.pdf")
  circlizePipelineNOZOOM(siteFile,
                         label = label,
                         PV.cutoff = PV.cutoff,
                         bestScore.cutoff = bestScore.cutoff,
                         bestFlank.cutoff = bestFlank.cutoff,
                         showNBS = showNBS,
                         ots.bed = ots.bed, 
                         realigned = realigned,
                         outFile = outFile,
                         species = species,
                         nbHits = nbHits,
                         top = top
  )
  
  
  
  # TEST MASAKO EMAIL 30 08 22
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/GENEWIZ_90-737856818/Masako/ELANDE-KO_T3T4_OVL1/ELANDE-KO_T3T4_OVL1_FINAL.xlsx")
  zoom.size <- 25000
  label <- FALSE
  PV.cutoff <- NULL
  bestScore.cutoff <- NULL
  bestFlank.cutoff <- 25
  showNBS <- TRUE
  gene.bed <- NULL
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/Masako/ELANE_KO/ELANE-KO-Treat3/data/ots.bed"
  realigned <- TRUE
  outFile <- file.path("~/Research/CASTSeq/circlize/Masako_300822/ELANDE-KO_T3T4_OVL1/ELANDE-KO_T3T4_OVL1_circlize_25k.pdf")
  species <- "hg38"
  nbHits <- 10
  top <- NULL
  
  circlizePipeline(siteFile = siteFile,
                   zoom.size = zoom.size,
                   label = label,
                   PV.cutoff = PV.cutoff,
                   bestScore.cutoff = bestScore.cutoff,
                   bestFlank.cutoff = bestFlank.cutoff,
                   showNBS = showNBS,
                   gene.bed = gene.bed,
                   ots.bed = ots.bed,
                   realigned = realigned,
                   outFile = outFile,
                   species = species,
                   nbHits = nbHits,
                   top = top
                   )
  
  outFile <- file.path("~/Research/CASTSeq/circlize/Masako_300822/ELANDE-KO_T3T4_OVL1/ELANDE-KO_T3T4_OVL1_circlize_25k_NOZOOM.pdf")
  circlizePipelineNOZOOM(siteFile,
                        label = label,
                        PV.cutoff = PV.cutoff,
                        bestScore.cutoff = bestScore.cutoff,
                        bestFlank.cutoff = bestFlank.cutoff,
                        showNBS = showNBS,
                        ots.bed = ots.bed, 
                        realigned = realigned,
                        outFile = outFile,
                        species = species,
                        nbHits = nbHits,
                        top = top
                        )
  
  
  # DOUBLE NICKASE
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/GENEWIZ_90-737856818/Masako/ELANDE-DN_T3T4_OVL1/ELANDE-DN_T3T4_OVL1_FINAL.xlsx")
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/Masako/ELANE_DN/ELANDE-DN-treated3/data/ots.bed"
  outFile <- file.path("~/Research/CASTSeq/circlize/Masako_300822/ELANDE-DN_T3T4_OVL1/ELANDE-DN_T3T4_OVL1_circlize_25k.pdf")
  
  
  circlizePipeline(siteFile = siteFile,
                   zoom.size = zoom.size,
                   label = label,
                   PV.cutoff = PV.cutoff,
                   bestScore.cutoff = bestScore.cutoff,
                   bestFlank.cutoff = bestFlank.cutoff,
                   showNBS = showNBS,
                   gene.bed = gene.bed,
                   ots.bed = ots.bed,
                   realigned = realigned,
                   outFile = outFile,
                   species = species,
                   nbHits = nbHits,
                   top = top
  )
  
  outFile <- file.path("~/Research/CASTSeq/circlize/Masako_300822/ELANDE-DN_T3T4_OVL1/ELANDE-DN_T3T4_OVL1_circlize_25k_NOZOOM.pdf")
  circlizePipelineNOZOOM(siteFile,
                         label = label,
                         PV.cutoff = PV.cutoff,
                         bestScore.cutoff = bestScore.cutoff,
                         bestFlank.cutoff = bestFlank.cutoff,
                         showNBS = showNBS,
                         ots.bed = ots.bed, 
                         realigned = realigned,
                         outFile = outFile,
                         species = species,
                         nbHits = nbHits,
                         top = top
  )
  
  
  # TALEN
  siteFile <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/STAT3/TALEN/STAT3-CD34-T3/OVL2_SIGNIF1/STAT3-CD34-T3_FINAL.xlsx")
  ots.bed <- "~/cluster/cluster/CASTSeq/pipelineGit/samples/GENEWIZ_90-737856818/STAT3/TALEN/STAT3-CD34-T3/data/ots.bed"
  outFile <- file.path("~/Research/CASTSeq/circlize/STAT3/STAT3-CD34-T3_circlize_25k.pdf")
  
  circlizePipeline(siteFile = siteFile,
                   zoom.size = zoom.size,
                   label = label,
                   PV.cutoff = PV.cutoff,
                   bestScore.cutoff = bestScore.cutoff,
                   bestFlank.cutoff = bestFlank.cutoff,
                   showNBS = showNBS,
                   gene.bed = gene.bed,
                   ots.bed = ots.bed,
                   realigned = realigned,
                   outFile = outFile,
                   species = species,
                   nbHits = nbHits,
                   top = top
  )
  
  outFile <- file.path("~/Research/CASTSeq/circlize/STAT3/STAT3-CD34-T3_circlize_25k_NOZOOM.pdf")
  circlizePipelineNOZOOM(siteFile,
                         label = label,
                         PV.cutoff = PV.cutoff,
                         bestScore.cutoff = bestScore.cutoff,
                         bestFlank.cutoff = bestFlank.cutoff,
                         showNBS = showNBS,
                         ots.bed = ots.bed, 
                         realigned = realigned,
                         outFile = outFile,
                         species = species,
                         nbHits = nbHits,
                         top = top
  )
}

