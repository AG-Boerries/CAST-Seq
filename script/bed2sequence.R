#library(BSgenome.Hsapiens.UCSC.hg38)# hg38!!!


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

bed2sequence <- function(bedMat, g=BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
{
	sqList <- mclapply(1:nrow(bedMat), function(i)
		getSeq(g, bedMat[i, 1], as.integer(bedMat[i, 2]), as.integer(bedMat[i, 3]), as.character = TRUE), mc.cores = NBCPU
		#getSeq(Hsapiens, i[1], as.integer(i[2]), as.integer(i[3]),as.BStringViews=TRUE)
		)
	return(sqList)
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


if(FALSE)
{

# LOAD CLUSTERS
setwd(file.path("/Volumes/Home/Geoffroy/offTargets/Giando/data"))
cl.G3 <- read.delim("G3_clusters.txt")
cl.UT <- read.delim("UT_clusters.txt")

bed <- as.data.frame(cl.G3[sample(1:nrow(cl.G3), 100), c(1:3, 8)])
colnames(bed) <- c("chr", "start", "end", "value")

sq.G3 <- bed2sequence(cl.G3)

}