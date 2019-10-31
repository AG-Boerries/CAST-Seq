library(openxlsx)


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

getCluster <- function(coord, bedMat, th){
	bedMat.sub <- bedMat[bedMat$chromosome == coord$chr &
						 bedMat$start >= coord$start &
						 bedMat$end <= coord$end, ]
	
	nbReads <- sum(bedMat.sub$read)
	nbCluster <- 0
	
	for(delta in bedMat.sub$delta){
		if(delta <= th) nbCluster <- nbCluster + 1
	}

	return(c(read = nbReads, collapseCluster = nbCluster))
}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


##########################################################################################
# G3 HiFi

# DEFINE CCR2 / CCR5 Coordinates
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters/coord")
coordMat <- read.xlsx("G3_HiFi_CCR2_CCR5.xlsx", sheet = 1)


# LOAD DELTA BED FILES
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D1/results/guide_aln/",
	"G3-Hifi-d1_S4_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D4/results/guide_aln/",
	"G3-Hifi-d4_S5_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D4_rep2/results/guide_aln/",
	"G3-d4-HF_S5_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D14/results/guide_aln/",
	"G3-Hifi-d14_S6_L001_Alignment_delta.bed")
	
fList <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList <- lapply(fList, read.delim, stringsAsFactors = FALSE)


# UT samples
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D1/results/guide_aln/",
	"UT-G3-d1_S3_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D4/results/guide_aln/",
	"UT-G3-d4_S4_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D4_rep2/results/guide_aln/",
	"G3-d4-UT_S4_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_HiFi_D14/results/guide_aln/",
	"UT-G3-d14_S5_L001_Alignment_delta.bed")
	
fList.ut <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList.ut <- lapply(fList.ut, read.delim, stringsAsFactors = FALSE)


clusterList <- lapply(1:nrow(coordMat), function(i){
	clusterMat <- do.call(rbind, lapply(bedList, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat <- t(cbind(t(clusterMat), SUM = colSums(clusterMat)))
	
	clusterMat.ut <- do.call(rbind, lapply(bedList.ut, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat.ut <- t(cbind(t(clusterMat.ut), SUM = colSums(clusterMat.ut)))
	colnames(clusterMat.ut) <- c("read.ctl", "collapseCluster.ctl")
	
	cbind(clusterMat, clusterMat.ut)
})
names(clusterList) <- coordMat$Name

setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters")
write.xlsx(clusterList, "G3_HiFi_CCR2_CCR5_manual_clusters.xlsx", row.names = TRUE)


##########################################################################################
# G3 WT

# DEFINE CCR2 / CCR5 Coordinates
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters/coord")
coordMat <- read.xlsx("G3_WT_CCR2_CCR5.xlsx", sheet = 1)

# LOAD DELTA BED FILES
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results/guide_aln/",
	"G3-WT-d1_S1_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D4/results/guide_aln/",
	"G3-WT-d4_S2_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D4_rep2/results/guide_aln/",
	"G3-d4-WT_S6_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D14/results/guide_aln/",
	"G3-WT-d14_S3_L001_Alignment_delta.bed")
	
fList <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList <- lapply(fList, read.delim, stringsAsFactors = FALSE)


# UT samples
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D1/results/guide_aln/",
	"UT-G3-d1_S3_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D4/results/guide_aln/",
	"UT-G3-d4_S4_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D4_rep2/results/guide_aln/",
	"G3-d4-UT_S4_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/G3_WT_D14/results/guide_aln/",
	"UT-G3-d14_S5_L001_Alignment_delta.bed")
	
fList.ut <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList.ut <- lapply(fList.ut, read.delim, stringsAsFactors = FALSE)


clusterList <- lapply(1:nrow(coordMat), function(i){
	clusterMat <- do.call(rbind, lapply(bedList, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat <- t(cbind(t(clusterMat), SUM = colSums(clusterMat)))
	
	clusterMat.ut <- do.call(rbind, lapply(bedList.ut, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat.ut <- t(cbind(t(clusterMat.ut), SUM = colSums(clusterMat.ut)))
	colnames(clusterMat.ut) <- c("read.ctl", "collapseCluster.ctl")
	
	cbind(clusterMat, clusterMat.ut)
})
names(clusterList) <- coordMat$Name

setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters")
write.xlsx(clusterList, "G3_WT_CCR2_CCR5_manual_clusters.xlsx", row.names = TRUE)



##########################################################################################
# HiFi 399

# DEFINE CCR2 / CCR5 Coordinates
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters/coord")
coordMat <- read.xlsx("HiFi399_CCR2_CCR5.xlsx", sheet = 1)

# LOAD DELTA BED FILES
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results/guide_aln/",
	"399HifiD1_S1_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results/guide_aln/",
	"399HifiD4_S2_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4_rep2/results/guide_aln/",
	"399-d4-HF_S2_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D14/results/guide_aln/",
	"399HifiD14_S3_L001_Alignment_delta.bed")
	
fList <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList <- lapply(fList, read.delim, stringsAsFactors = FALSE)


# UT samples
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D1/results/guide_aln/",
	"UT-399-d1_S1_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4/results/guide_aln/",
	"UT-399-d4_S2_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D4_rep2/results/guide_aln/",
	"399-d4-UT_S1_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HiFi399_D14/results/guide_aln/",
	"UT-399-d14_S7_L001_Alignment_delta.bed")
	
fList.ut <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList.ut <- lapply(fList.ut, read.delim, stringsAsFactors = FALSE)


clusterList <- lapply(1:nrow(coordMat), function(i){
	clusterMat <- do.call(rbind, lapply(bedList, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat <- t(cbind(t(clusterMat), SUM = colSums(clusterMat)))
	
	clusterMat.ut <- do.call(rbind, lapply(bedList.ut, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat.ut <- t(cbind(t(clusterMat.ut), SUM = colSums(clusterMat.ut)))
	colnames(clusterMat.ut) <- c("read.ctl", "collapseCluster.ctl")
	
	cbind(clusterMat, clusterMat.ut)
})
names(clusterList) <- coordMat$Name

setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters")
write.xlsx(clusterList, "HiFi399_CCR2_CCR5_manual_clusters.xlsx", row.names = TRUE)





##########################################################################################
# WT 399

# DEFINE CCR2 / CCR5 Coordinates
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters/coord")
coordMat <- read.xlsx("WT399_CCR2_CCR5.xlsx", sheet = 1)

# LOAD DELTA BED FILES
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D1/results/guide_aln/",
	"399WtD1_S4_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D4/results/guide_aln/",
	"399WtD4_S5_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D4_rep2/results/guide_aln/",
	"399-d4-WT_S3_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D14/results/guide_aln/",
	"399WtD14_S6_L001_Alignment_delta.bed")
	
fList <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList <- lapply(fList, read.delim, stringsAsFactors = FALSE)


# UT samples
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D1/results/guide_aln/",
	"UT-399-d1_S1_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D4/results/guide_aln/",
	"UT-399-d4_S2_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D4_rep2/results/guide_aln/",
	"399-d4-UT_S1_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/WT399_D14/results/guide_aln/",
	"UT-399-d14_S7_L001_Alignment_delta.bed")
	
fList.ut <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList.ut <- lapply(fList.ut, read.delim, stringsAsFactors = FALSE)


clusterList <- lapply(1:nrow(coordMat), function(i){
	clusterMat <- do.call(rbind, lapply(bedList, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat <- t(cbind(t(clusterMat), SUM = colSums(clusterMat)))
	
	clusterMat.ut <- do.call(rbind, lapply(bedList.ut, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat.ut <- t(cbind(t(clusterMat.ut), SUM = colSums(clusterMat.ut)))
	colnames(clusterMat.ut) <- c("read.ctl", "collapseCluster.ctl")
	
	cbind(clusterMat, clusterMat.ut)
})
names(clusterList) <- coordMat$Name

setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters")
write.xlsx(clusterList, "WT399_CCR2_CCR5_manual_clusters.xlsx", row.names = TRUE)




##########################################################################################
# HBBTal

# DEFINE HBB / HBD Coordinates
setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters/coord")
coordMat <- read.xlsx("HBBTal_HBB_HBD.xlsx", sheet = 1)

# LOAD DELTA BED FILES
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D1/results/guide_aln/",
	"HBBTalD1_S7_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D4/results/guide_aln/",
	"HBBTalD4_S8_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D4_rep2/results/guide_aln/",
	"HBB-D4-TAL_S8_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D14/results/guide_aln/",
	"HBB-Tal-d14_S8_L001_Alignment_delta.bed")
	
fList <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList <- lapply(fList, read.delim, stringsAsFactors = FALSE)


# UT samples
f1 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D1/results/guide_aln/",
	"UT-HBB-d1_S6_L001_Alignment_delta.bed")
	
f2 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D4/results/guide_aln/",
	"UT-HBB-d4_S7_L001_Alignment_delta.bed")
		
f3 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D4_rep2/results/guide_aln/",
	"HBB-D4-UT_S7_L001_Alignment_delta.bed")
		
f4 <- file.path("~/cluster/master/offTargets/Giando/pipelineGit/samples/HBBTal_D14/results/guide_aln/",
	"UT-HBB-d14_S8_L001_Alignment_delta.bed")
	
fList.ut <- c(D1 = f1, D4 = f2, D4_rep2 = f3, D14 = f4)	
bedList.ut <- lapply(fList.ut, read.delim, stringsAsFactors = FALSE)


clusterList <- lapply(1:nrow(coordMat), function(i){
	clusterMat <- do.call(rbind, lapply(bedList, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat <- t(cbind(t(clusterMat), SUM = colSums(clusterMat)))
	
	clusterMat.ut <- do.call(rbind, lapply(bedList.ut, function(j) getCluster(coordMat[i,], j, th = 1500)))
	clusterMat.ut <- t(cbind(t(clusterMat.ut), SUM = colSums(clusterMat.ut)))
	colnames(clusterMat.ut) <- c("read.ctl", "collapseCluster.ctl")
	
	cbind(clusterMat, clusterMat.ut)
})
names(clusterList) <- coordMat$Name

setwd("/Volumes/Home/Geoffroy/offTargets/Giando/manualClusters")
write.xlsx(clusterList, "HBBTal_HBB_HBD_manual_clusters.xlsx", row.names = TRUE)


















