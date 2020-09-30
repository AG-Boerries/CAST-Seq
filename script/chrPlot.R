#library(karyoploteR)
#library(openxlsx)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

chrPlot <- function(inputF, outputF, hits = NULL, score = NULL, pv = NULL)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	if(!is.null(hits)) readMat <- readMat[readMat$hits > hits, ]
	if(!is.null(score)) readMat <- readMat[readMat$score > score, ]
	if(!is.null(pv)) readMat <- readMat[readMat$adj.pvalue < pv, ]

	if(nrow(readMat)>1){
		bed <- readMat[, c("chromosome", "start", "end", "group")]
		bed <- bed[order(bed$group), ]
		groups <- bed$group

		# DF TO GR
		bed.gr <- makeGRangesFromDataFrame(bed, seqnames.field = "chromosome",
			start.field = "start", end.field = "end", ignore.strand = TRUE,
			keep.extra.columns = TRUE)
		
		# DEFINE COLORS
		chrColor <- rep(NA, length(groups))
		chrColor[groups == "OMT"] <- "red"	
		chrColor[groups == "HMT"] <- "blue"	
		chrColor[groups == "NBS"] <- "grey50"		
	
		#plotDefaultPlotParams(plot.type=1)
		plot.params <- getDefaultPlotParams(plot.type=1)
		plot.params$data1height <- 2
		plot.params$ideogramheight <- 20
	
		pdf(outputF)
		kp <- plotKaryotype(plot.params = plot.params, labels.plotter = NULL)
		kpAddChromosomeNames(kp, cex = 1.5)
		kpPlotRegions(kp, data=extendRegions(bed.gr, extend.end = 30e5), col= chrColor,
			avoid.overlapping = FALSE, r0 =-35, r1=1
			)
		dev.off()
	}
	
}

chrPlotAside <- function(inputF, outputP, hits = NULL, score = NULL, pv = NULL)
{
	readMat <- read.xlsx(inputF, sheet = 1)
	if(!is.null(hits)) readMat <- readMat[readMat$hits > hits, ]
	if(!is.null(score)) readMat <- readMat[readMat$score > score, ]
	if(!is.null(pv)) readMat <- readMat[readMat$adj.pvalue < pv, ]

	if(nrow(readMat)>1){
		bed <- readMat[, c("chromosome", "start", "end", "group")]
		bed <- bed[order(bed$group), ]
		groups <- bed$group
	
		groups.color <- c("OMT" = "red", "HMT" = "blue", "NBS" = "grey")
		
		#plotDefaultPlotParams(plot.type=1)
		plot.params <- getDefaultPlotParams(plot.type=1)
		plot.params$data1height <- 2
		plot.params$ideogramheight <- 30	
		
		for(g in groups)
			{
			bed.sub <- bed[bed$group == g, ]
			bed.sub.gr <- makeGRangesFromDataFrame(bed.sub, seqnames.field = "chromosome",
			start.field = "start", end.field = "end", ignore.strand = TRUE,
			keep.extra.columns = TRUE)
		
			pdf(paste0(outputP, "_", g, ".pdf"))
			kp <- plotKaryotype(plot.params = plot.params, labels.plotter = NULL)
			kpAddChromosomeNames(kp, cex = 1.5)
			kpPlotRegions(kp, data=extendRegions(bed.sub.gr, extend.end = 30e5), col= groups.color[g],
				avoid.overlapping = FALSE, r0 =-35, r1=1
				)
			dev.off()
			}
	}
}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

if(FALSE)
{
# TOY
inputF <- file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D1/results/guide_aln/",
	"G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx")
clusters <- 1
score = NULL
pv <- 0.05



chrDir <- file.path("~/cluster/master/offTargets/Giando/chrPlot/")

# G3 WT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D1/results/overlap_aln/"))
chrPlot("G3_WT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "G3_WT_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("G3_WT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "G3_WT_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)

# G3 HiFi
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_HiFi_D1/results/overlap_aln/"))
chrPlot("G3_HiFi_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "G3_HiFi_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("G3_HiFi_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "G3_HiFi_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)

# 399 WT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/WT399_D1/results/overlap_aln/"))
chrPlot("WT399_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "WT399_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("WT399_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "WT399_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)


# 399 HiFi
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/HiFi399_D1/results/overlap_aln/"))
chrPlot("HiFi399_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "HiFi399_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("HiFi399_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "HiFi399_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)


# FANCF
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/FANCF_decoy_rep1/results/overlap_aln/"))
chrPlot("FANCF_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "FANCF_decoy_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("FANCF_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "FANCF_decoy_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)

# VEGFA
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/VEGFA_decoy_rep1/results/overlap_aln/"))
chrPlot("VEGFA_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "VEGFA_decoy_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("VEGFA_decoy_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "VEGFA_decoy_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)

# OT
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/OT6/results/overlap_aln/"))
chrPlot("OT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "OT_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("OT_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "OT_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)


# HBB
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/HBBTal_D1/results/overlap_aln/"))
chrPlot("HBBTal_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
	file.path(chrDir, "HBBTal_ovl_w250_aln_chrPlot.pdf"), clusters = NULL, score = NULL, pv = NULL)
#chrPlotAside("HBBTal_ovl_w250_aln_stat_FLANK_GROUP_GENES_ONCO.xlsx",
#	file.path(chrDir, "HBBTal_ovl_w250_aln_chrPlot"), clusters = NULL, score = NULL, pv = NULL)





##########################################################################################
# OLD



# G3 WT D1
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D1/results/guide_aln/"))
chrPlot("G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-WT-d1_S1_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("G3-WT-d1_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-WT-d1_S1_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
	
# G3 WT D4
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D4/results/guide_aln/"))
chrPlot("G3-WT-d4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-WT-d4_S2_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("G3-WT-d4_S2_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-WT-d4_S2_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
# G3 WT D14
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_WT_D14/results/guide_aln/"))
chrPlot("G3-WT-d14_S3_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-WT-d14_S3_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("G3-WT-d14_S3_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-WT-d14_S3_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
		
# G3 HiFi D1
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_HiFi_D1/results/guide_aln/"))
chrPlot("G3-Hifi-d1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-Hifi-d1_S4_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("G3-Hifi-d1_S4_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-Hifi-d1_S4_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
	
# G3 HiFi D4
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_HiFi_D4/results/guide_aln/"))
chrPlot("G3-Hifi-d4_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-Hifi-d4_S5_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("G3-Hifi-d4_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-Hifi-d4_S5_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
# G3 HiFi D14
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/G3_HiFi_D14/results/guide_aln/"))
chrPlot("G3-Hifi-d14_S6_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-Hifi-d14_S6_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("G3-Hifi-d14_S6_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "G3-Hifi-d14_S6_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
	
# FANCF decoy
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/FANCF_decoy/results/guide_aln/"))
chrPlot("FANCF-decoy_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "FANCF-decoy_S1_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("FANCF-decoy_S1_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "FANCF-decoy_S1_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
	
# VEGFA decoy
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/VEGFA_decoy/results/guide_aln/"))
chrPlot("VEGFA-decoy_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "VEGFA-decoy_S5_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("VEGFA-decoy_S5_L001_w250_aln_stat_FLANK_GROUP_GENES_ONCO_SCORE.xlsx",
	file.path(chrDir, "VEGFA-decoy_S5_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
		
# HBBTal D1
setwd(file.path("~/cluster/master/offTargets/Giando/pipeline/HBBTal_D1_LEFT/results/guide_aln/"))
chrPlot("HBBTalD1_S7_L001_w250_aln_stat_FLANK_GROUP_LEFT_RIGHT.xlsx",
	file.path(chrDir, "HBBTalD1_S7_L001_w250_aln_clust1_pv0.05_chrPlot.pdf"), clusters = 1, score = NULL, pv = 0.05)
chrPlotAside("HBBTalD1_S7_L001_w250_aln_stat_FLANK_GROUP_LEFT_RIGHT.xlsx",
	file.path(chrDir, "HBBTalD1_S7_L001_w250_aln_clust1_pv0.05_chrPlot"), clusters = 1, score = NULL, pv = 0.05)
		
	
	

}