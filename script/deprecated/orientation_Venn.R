library(VennDiagram)
library(openxlsx)
library(GenomicRanges)
library(venneuler)

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################







############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


############################################################################################
# CCR5 #1 

# LOAD SITES
setwd(file.path("~/Research/CASTSeq/revision/overlap/venn/input/"))
siteM <- read.xlsx("CCR5_1_summary.xlsx", sheet = 1)

# REMOVE ON TARGET
siteM <- siteM[siteM$`CAST-Seq.category` != "ON", ]

dfA <- siteM[siteM$Sample == "CCR5_1_FOR", ]
dfB <- siteM[siteM$Sample == "CCR5_1_REV", ]

#dfA$Total.Hits <- dfA$Total.Hits / sum(dfA$Total.Hits)
#dfB$Total.Hits <- dfB$Total.Hits / sum(dfB$Total.Hits)

# OVERLAP SITES

grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                start.field = "start", end.field = "end",
                                keep.extra.columns = FALSE, ignore.strand = TRUE)
grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                start.field = "start", end.field = "end",
                                keep.extra.columns = FALSE, ignore.strand = TRUE)

gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = -1)

df.ovl <- data.frame(dfA[queryHits(gr.ovl),], dfB[subjectHits(gr.ovl),])

df.ovl$Total.Hits[duplicated(queryHits(gr.ovl))] <- 0
df.ovl$Total.Hits.1[duplicated(subjectHits(gr.ovl))] <- 0
df.ovl$Total.Hits.OVL <- df.ovl$Total.Hits + df.ovl$Total.Hits.1

dfA.spe <- dfA[-queryHits(gr.ovl),]
dfB.spe <- dfB[-subjectHits(gr.ovl),]

# COMBINE 1 TO N SITES
# NOT USED


# CALCULATE HITS PER SAMPLE
nbA <- 4
nbB <- 2
nbTOT <- nbA + nbB

hpsA <- dfA.spe$Total.Hits / nbA
hpsB <- dfB.spe$Total.Hits / nbB
hpsOVL <- df.ovl$Total.Hits.OVL / nbTOT


# VENN
v <- venneuler(c("FOR"=sum(hpsA),
                 "REV"=sum(hpsB),
                 "FOR&REV"=sum(hpsOVL)))

setwd(file.path("~/Research/CASTSeq/revision/overlap/venn/"))
pdf("CCR5_1_FOR_vs_REV_Venn.pdf", width = 5, height = 5)
plot(v)
dev.off()

toxlsx <- data.frame(hitPerSample.FOR = sum(hpsA),
                     hitPerSample.REV = sum(hpsB),
                     hitPerSample.BOTH = sum(hpsOVL))
write.xlsx(toxlsx, "CCR5_1_FOR_vs_REV_Venn.xlsx")

# VENN V2
v <- venneuler(c("FOR"=sum(dfA.spe$Total.Hits),
                 "REV"=sum(dfB.spe$Total.Hits),
                 "FOR&REV"=sum(df.ovl$Total.Hits.OVL)))
plot(v)



############################################################################################
# CCR5 #2 


# LOAD SITES
setwd(file.path("~/Research/CASTSeq/revision/overlap/venn/input/"))
siteM <- read.xlsx("CCR5_2_summary.xlsx", sheet = 1)

# REMOVE ON TARGET
siteM <- siteM[siteM$`CAST-Seq.category` != "ON", ]

dfA <- siteM[siteM$Sample == "CCR5_2_FOR", ]
dfB <- siteM[siteM$Sample == "CCR5_2_REV", ]

#dfA$Total.Hits <- dfA$Total.Hits / sum(dfA$Total.Hits)
#dfB$Total.Hits <- dfB$Total.Hits / sum(dfB$Total.Hits)

# OVERLAP SITES

grA <- makeGRangesFromDataFrame(dfA, seqnames.field = "chromosome",
                                start.field = "start", end.field = "end",
                                keep.extra.columns = FALSE, ignore.strand = TRUE)
grB <- makeGRangesFromDataFrame(dfB, seqnames.field = "chromosome",
                                start.field = "start", end.field = "end",
                                keep.extra.columns = FALSE, ignore.strand = TRUE)

gr.ovl <- findOverlaps(query = grA, subject = grB, type = "any", maxgap = -1)

df.ovl <- data.frame(dfA[queryHits(gr.ovl),], dfB[subjectHits(gr.ovl),])

df.ovl$Total.Hits[duplicated(queryHits(gr.ovl))] <- 0
df.ovl$Total.Hits.1[duplicated(subjectHits(gr.ovl))] <- 0
df.ovl$Total.Hits.OVL <- df.ovl$Total.Hits + df.ovl$Total.Hits.1

dfA.spe <- dfA[-queryHits(gr.ovl),]
dfB.spe <- dfB[-subjectHits(gr.ovl),]

# COMBINE 1 TO N SITES
# NOT USED


# CALCULATE HITS PER SAMPLE
nbA <- 4
nbB <- 2
nbTOT <- nbA + nbB

hpsA <- dfA.spe$Total.Hits / nbA
hpsB <- dfB.spe$Total.Hits / nbB
hpsOVL <- df.ovl$Total.Hits.OVL / nbTOT


# VENN
v <- venneuler(c("FOR"=sum(hpsA),
                 "REV"=sum(hpsB),
                 "FOR&REV"=sum(hpsOVL)))
setwd(file.path("~/Research/CASTSeq/revision/overlap/venn/"))
pdf("CCR5_2_FOR_vs_REV_Venn.pdf", width = 5, height = 5)
plot(v)
dev.off()

toxlsx <- data.frame(hitPerSample.FOR = sum(hpsA),
                     hitPerSample.REV = sum(hpsB),
                     hitPerSample.BOTH = sum(hpsOVL))
write.xlsx(toxlsx, "CCR5_2_FOR_vs_REV_Venn.xlsx")


# VENN V2
v <- venneuler(c("FOR"=sum(dfA.spe$Total.Hits),
                 "REV"=sum(dfB.spe$Total.Hits),
                 "FOR&REV"=sum(df.ovl$Total.Hits.OVL)))
plot(v)

