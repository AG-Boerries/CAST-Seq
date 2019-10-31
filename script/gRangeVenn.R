library(ggplot2)
library(data.table)
library(openxlsx)
library(ChIPpeakAnno)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

saveVennDeploy <- function(venn, outFile)
{
	venn_info <- data.frame(Group = names(venn), Sheet = paste("venn_", seq(1, length(venn)), sep = ""))
	venn_info <- list(Info = venn_info)

	names(venn) <- paste0("venn_", 1:length(venn))

	write.xlsx(c(venn_info, venn), outFile)	
}


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

# LOAD SITE FILES
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline"))
siteFiles <- list.files(pattern = "_SCORE.xlsx", recursive = TRUE)
names(siteFiles) <- basename(siteFiles)
names(siteFiles) <- gsub("_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx", "", names(siteFiles))

siteFiles <- siteFiles[-grep("^UT-", names(siteFiles))]# remove untreated samples

siteFiles.sub <- siteFiles[c("G3-WT-d4_S2_L001_w250", "FANCF-decoy_S1_L001_w250", "VEGFA-decoy_S5_L001_w250", "399WtD4_S5_L001_w250")]#for testing

siteList <- lapply(siteFiles.sub, read.xlsx, sheet = 1)

# SELECT CBS SITES
siteList <- lapply(siteList, function(i) i[i$group == "CBS", c("chromosome", "start", "end", "read", "collapseCluster", "read.ctl",
		"collapseCluster.ctl", "width.raw", "width", "OddRatio", "pvalue",
		"adj.pvalue", "artificial", "group")])

siteList.gr <- lapply(siteList, makeGRangesFromDataFrame, ignore.strand = TRUE, keep.extra.columns = TRUE)


# VENN ANALYSIS
ol <- findOverlapsOfPeaks(siteList.gr, maxgap=2000)
count <- makeVennDiagram(ol)$vennCounts
count <- apply(count, 2, as.numeric)
count.sum <- rowSums(count[, -ncol(count)])
count <- cbind(count, nb.samples = count.sum)
count <- count[order(count.sum), ]


setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/Venn"))
write.xlsx(count, "venn_count.xlsx")

ol.sites <- lapply(ol$peaklist, as.data.frame)
ol.sites <- lapply(ol.sites, function(i) i[,-ncol(i)])
names(ol.sites) <- gsub("///", "_AND_", names(ol.sites))
saveVennDeploy(ol.sites, "venn.xlsx")


##########################################################################################
# G3 WT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline"))
siteFiles.sub <- siteFiles[c("G3-WT-d1_S1_L001_w250", "G3-WT-d4_S2_L001_w250", "G3-WT-d14_S3_L001_w250")]#for testing
siteList <- lapply(siteFiles.sub, read.xlsx, sheet = 1)

# SELECT CBS SITES
siteList <- lapply(siteList, function(i) i[i$group == "CBS", c("chromosome", "start", "end", "read", "collapseCluster", "read.ctl",
		"collapseCluster.ctl", "width.raw", "width", "OddRatio", "pvalue",
		"adj.pvalue", "artificial", "group")])

siteList.gr <- lapply(siteList, makeGRangesFromDataFrame, ignore.strand = TRUE, keep.extra.columns = TRUE)


# VENN ANALYSIS
ol <- findOverlapsOfPeaks(siteList.gr, maxgap=2000)
count <- makeVennDiagram(ol)$vennCounts
count <- apply(count, 2, as.numeric)
count.sum <- rowSums(count[, -ncol(count)])
count <- cbind(count, nb.samples = count.sum)
count <- count[order(count.sum), ]


setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/Venn"))
write.xlsx(count, "G3_WT_venn_count.xlsx")

ol.sites <- lapply(ol$peaklist, as.data.frame)
ol.sites <- lapply(ol.sites, function(i) i[,-ncol(i)])
names(ol.sites) <- gsub("///", "_AND_", names(ol.sites))
saveVennDeploy(ol.sites, "G3_WT_venn.xlsx")



##########################################################################################
# OT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline"))
siteFiles.sub <- siteFiles[c("OT2-1_S4_L001_w250", "OT2-2_S5_L001_w250", "OT2-3_S6_L001_w250")]#for testing
siteList <- lapply(siteFiles.sub, read.xlsx, sheet = 1)

# SELECT CBS SITES
siteList <- lapply(siteList, function(i) i[, c("chromosome", "start", "end", "read", "collapseCluster", "read.ctl",
		"collapseCluster.ctl", "width.raw", "width", "OddRatio", "pvalue",
		"adj.pvalue", "artificial", "group")])

siteList.gr <- lapply(siteList, makeGRangesFromDataFrame, ignore.strand = TRUE, keep.extra.columns = TRUE)


# VENN ANALYSIS
ol <- findOverlapsOfPeaks(siteList.gr, maxgap=2000)
count <- makeVennDiagram(ol)$vennCounts
count <- apply(count, 2, as.numeric)
count.sum <- rowSums(count[, -ncol(count)])
count <- cbind(count, nb.samples = count.sum)
count <- count[order(count.sum), ]


setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/Venn"))
write.xlsx(count, "OT_venn_count.xlsx")

ol.sites <- lapply(ol$peaklist, as.data.frame)
ol.sites <- lapply(ol.sites, function(i) i[,-ncol(i)])
names(ol.sites) <- gsub("///", "_AND_", names(ol.sites))
saveVennDeploy(ol.sites, "OT_venn.xlsx")




