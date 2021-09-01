#library(openxlsx)
#library(ggplot2)
#library(data.table)



getCategory <- function(x)
{
return(c("nb=0" = sum(x == 0) / length(x) * 100,
		 "nb=1" = sum(x == 1) / length(x) * 100,
		 "nb=2" = sum(x == 2) / length(x) * 100,
		 "nb=3" = sum(x == 3) / length(x) * 100,
		 "nb=4" = sum(x == 4) / length(x) * 100,
		 "nb>=5" = sum(x >= 5) / length(x) * 100))

}


piePlot <- function(inputF, outputF, hits = NULL, score = NULL, pv = NULL)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	if(!is.null(hits)) readMat <- readMat[readMat$hits > hits, ]
	if(!is.null(score)) readMat <- readMat[readMat$score > score, ]
	if(!is.null(pv)) readMat <- readMat[readMat$adj.pvalue < pv, ]

	if(nrow(readMat)>1){
		groups <- unique(readMat$group)
	
		nbMissList <- lapply(groups, function(i) readMat$nb.mism[readMat$group == i])
		nbIndelList <- lapply(groups, function(i) readMat$length.ins[readMat$group == i] + readMat$length.del[readMat$group == i])
		nbAllList <- lapply(groups, function(i) readMat$nb.mism[readMat$group == i] + readMat$length.ins[readMat$group == i] + readMat$length.del[readMat$group == i])

		catList.miss <- lapply(nbMissList, getCategory)
		catList.indel <- lapply(nbIndelList, getCategory)
		catList.all <- lapply(nbAllList, getCategory)

		ggmatList.miss <- lapply(1:length(groups), function(i)
			data.frame(Nb = catList.miss[[i]], Label = names(catList.miss[[i]]), Group = groups[i], Category = "mismatch")
			)
		ggmatList.indel <- lapply(1:length(groups), function(i)
			data.frame(Nb = catList.indel[[i]], Label = names(catList.indel[[i]]), Group = groups[i], Category = "indel")
			)	
		ggmatList.all <- lapply(1:length(groups), function(i)
			data.frame(Nb = catList.all[[i]], Label = names(catList.all[[i]]), Group = groups[i], Category = "mismatch and indel")
			)

		ggmat <- do.call(rbind, c(ggmatList.miss, ggmatList.indel))
		ggmat$Label = factor(ggmat$Label, levels = rev(c("nb=0", "nb=1", "nb=2", "nb=3", "nb=4", "nb>=5")))
		ggmat$Group = factor(ggmat$Group, levels = c("off.target", "hom.recomb", "CBS"))

		p <- ggplot(ggmat, aes(x="", y=Nb, fill=Label))
		p <- p + geom_bar(width = 1, stat = "identity")
		p <- p + coord_polar("y", start=0)
		p <- p + facet_grid(facets=Category ~ Group)
		p <- p + scale_fill_grey() + theme_minimal()
		p <- p + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
	
		pdf(outputF)
		plot(p)
		dev.off()
	}
}


getRegionMatrix <- function(rmin, rmax, by = 4)
{
	firstCol <- seq(rmin, rmax, by = by)
	if(firstCol[length(firstCol)] == rmax){
		firstCol <- firstCol[-length(firstCol)]
		secondCol <- c(firstCol[-1] - 1, rmax)
		}else(secondCol <- c(firstCol[-1] - 1, rmax))
	return(matrix(c(firstCol, secondCol), ncol = 2, byrow = FALSE))
}

getLong <- function(short)
{
	v2 <- short[,2]
	v2 <- v2*2
	v1 <- v2 + 1
	v1 <- c(1, v1[-length(v1)])
	return(data.frame(v1,v2))
}

getNBmm <- function(ref, v)
{
	if(ref == "A") return(sum(v == "C") + sum(v == "T") + sum(v == "G"))
	else if(ref == "C") return(sum(v == "A") + sum(v == "T") + sum(v == "G"))
	else if(ref == "G") return(sum(v == "A") + sum(v == "C") + sum(v == "T"))
	else if(ref == "T") return(sum(v == "A") + sum(v == "C") + sum(v == "G"))
	else if(ref == "R") return(sum(v == "C") + sum(v == "T"))
	else if(ref == "Y") return(sum(v == "A") + sum(v == "G"))
	else if(ref == "S") return(sum(v == "A") + sum(v == "T"))
	else if(ref == "W") return(sum(v == "G") + sum(v == "C"))
	else if(ref == "K") return(sum(v == "C") + sum(v == "A"))
	else if(ref == "M") return(sum(v == "G") + sum(v == "T"))
	else if(ref == "B") return(sum(v == "A"))
	else if(ref == "D") return(sum(v == "C"))
	else if(ref == "H") return(sum(v == "G"))
	else if(ref == "V") return(sum(v == "T"))
	else if(ref == "N") return(0)
}

pcBarplot <- function(inputF, PAMlength = 3, binLength = 4)
{
	# REGIONS SHOULD BE MOFIFIED ACCORDING TO GUIDE SEQUENCE !!!
	heatMat <- read.xlsx(inputF, sheet = 1)
	heatMat <- heatMat[, -1] 
	
	refSeq <- heatMat[1, ]
	heatMat <- heatMat[-1, ]
	
	# CUT PAM
	pamMat <- heatMat[, seq(ncol(heatMat)+1 - 2*PAMlength, ncol(heatMat))]
	heatMat <- heatMat[, -seq(ncol(heatMat)+1 - 2*PAMlength, ncol(heatMat))]
	
	# DEFINE REGIONS
	regions.short <- getRegionMatrix(1, ncol(heatMat) / 2, binLength)
	#regions <- getRegionMatrix(1, ncol(heatMat), binLength * 2)
	regions <- getLong(regions.short)

	nbMiss <- c()
	nbIndel <- c()
	for(i in 1:nrow(regions))
		{
		idx.start <- regions[i,1]
		idx.end <- regions[i,2]
		heatMat.sub <- heatMat[, idx.start:idx.end]
	
		nbMiss <- c(nbMiss, sum(heatMat.sub == "A") + sum(heatMat.sub == "T") + sum(heatMat.sub == "G") + sum(heatMat.sub == "C"))
		nbIndel <- c(nbIndel, sum(heatMat.sub == 1) + sum(heatMat.sub == -1))
		}
	# Get percentage
	nbMiss <- nbMiss / (nrow(heatMat)*ceiling((regions[,2] - regions[,1]) / 2)) * 100
	nbIndel <- nbIndel / (nrow(heatMat)*ceiling((regions[,2] - regions[,1]) / 2)) * 100

	# PAM
	pamSeq <- refSeq[seq(1, ncol(refSeq), by = 2)]
	pamSeq <- as.character(pamSeq[, (length(pamSeq)-length(PAMlength)-1):length(pamSeq)])
	pamMat.mm <- pamMat[, seq(1, ncol(pamMat), by = 2)]

	nbMiss.pam <- 0
	for(i in 1:length(pamSeq))
		{
		refAA <- pamSeq[i]
		nbMiss.pam <- nbMiss.pam + getNBmm(refAA, pamMat.mm[, i])
		}
		
	nbMiss.pam <- nbMiss.pam / (nrow(pamMat.mm)*PAMlength) * 100

	nbIndel.pam <- sum(pamMat == 1) + sum(pamMat == -1)
	nbIndel.pam <- nbIndel.pam / (nrow(pamMat) * PAMlength) * 100

	nbMiss <- c(nbMiss, nbMiss.pam)
	nbIndel <- c(nbIndel, nbIndel.pam)
	
	# LABEL AND PLOT
	regions.label <- apply(regions.short, 1, paste, collapse = ":")
	regions.label <- c(regions.label, "PAM")

	ggmat <- data.frame(mismatch = nbMiss, indel = nbIndel, regions = regions.label)
	ggmat <- melt(ggmat)

	ggmat$regions <- factor(ggmat$regions, levels = regions.label)

	p <- ggplot(ggmat, aes(x=regions, y=value, fill=variable))
	p <- p + geom_bar(width = 1, stat = "identity")
	p <- p + scale_fill_grey() + theme_minimal()
	p <- p + theme(text = element_text(size=17))
	p <- p + scale_y_continuous(name="percentage")
	p <- p + theme(axis.title.x=element_blank())

	pdf(gsub(".xlsx", "_pcBarplot.pdf", inputF))
	plot(p)
	dev.off()
	
}










############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


if(FALSE)
{


##########################################################################################
# PIE CHART V2


donuts <- function(x, group = 1, labels = NA, col = NULL, radius = c(.7, 1))
{
	#' x      numeric vector for each slice
	#' group  vector identifying the group for each slice
	#' labels vector of labels for individual slices
	#' col    colors for each group
	#' radius radius for inner and outer pie (usually in [0,1])
 	group <- rep_len(group, length(x))
  	ug  <- unique(group)
  	tbl <- table(group)[order(ug)]

  	col <- if (is.null(col))
    	seq_along(ug) else rep_len(col, length(ug))
  	col.main <- Map(rep, col[seq_along(tbl)], tbl)
  	col.sub  <- lapply(col.main, function(x) {
    	al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    Vectorize(adjustcolor)(x, alpha.f = al)
  	})

  plot.new()

  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],
      col = unlist(col.sub), labels = labels,
      clockwise = TRUE)

  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],
      col = unlist(col.main), labels = NA,
      clockwise = TRUE)
}

par(mfrow = c(1,2), mar = c(0,4,0,4))
with(browsers,
     donuts(share, browser, sprintf('%s: %s%%', version, share),
            col = c('cyan2','red','orange','green','dodgerblue2'))
)

with(mtcars,
     donuts(mpg, interaction(gear, cyl), rownames(mtcars))
)

setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/pipeline/results/G3_all_w250/"))
readMat <- read.xlsx("G3_all_w250_aln_stat_FLANK_GROUP.xlsx", sheet=1)

groups <- unique(readMat$group)

nbMissList <- lapply(groups, function(i) readMat$nb.mism[readMat$group == i])
nbIndelList <- lapply(groups, function(i) readMat$length.ins[readMat$group == i] + readMat$length.del[readMat$group == i])
	
for(i in 1:length(groups))
	{
	nbMissList[[i]][nbMissList[[i]] > 5] <- 5# set limit to 5
	nbIndelList[[i]][nbIndelList[[i]] > 5] <- 5# set limit to 5
	}	
	
intList <- lapply(1:length(groups), function(i)
	as.character(interaction(nbMissList[[i]], nbIndelList[[i]]))
	)	

intList.tb <- 3

ggmatList <- list()
for(i in 1:length(groups))
	{
	gnames <- strsplit(names(intList.tb[[i]]), split = "\\.")
	gnames <- unlist(lapply(gnames, function(i) i[1]))
	
	ggmatList[[i]] <- data.frame(name = names(intList.tb[[i]]), mismatch = gnames, nb = as.numeric(intList.tb[[i]]))
	
	}
	
with(a,
     donuts(nb, mismatch, name)
)


##########################################################################################
# BARPLOT

inputF <- file.path("~/cluster/master/offTargets/Giando/pipeline/FANCF_decoy/results/guide_aln/",
	"FANCF-decoy_S1_L001_w250_aln_heatmap_clust1_pv0.05.xlsx")
inputF <- file.path("~/cluster/master/offTargets/Giando/pipeline/OT_rep1/results/guide_aln/OT2-1_S4_L001_w250_aln_heatmap_clust1_pv0.05.xlsx")	
PAMlength <- 3
binLength <- 4

}