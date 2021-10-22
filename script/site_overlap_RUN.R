

source(file.path("~/Research/CASTSeq/pipelineGit/script/site_overlap.R"))




################################
# NICKASE (JULIA EMAIL 23.02.21)

# COL7A1-1-g1 with COL7A1-2-g1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-1-g1/results/guide_aln/COL7A1-1-g1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-2-g1/results/guide_aln/COL7A1-2-g1_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/overlap/Nickase/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "COL7A1-1-g1_with_COL7A1-2-g1.xlsx")

makeUpset("COL7A1-1-g1_with_COL7A1-2-g1.xlsx")

# COL7A1-1-g2 with COL7A1-2-g2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-1-g2/results/guide_aln/COL7A1-1-g2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-2-g2/results/guide_aln/COL7A1-2-g2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/overlap/Nickase/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "COL7A1-1-g2_with_COL7A1-2-g2.xlsx")

makeUpset("COL7A1-1-g2_with_COL7A1-2-g2.xlsx")


# COL7A1-2-D10Ag1g2 (first vs. second dataset)
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/COL7A1_G1_G2/results_middle/guide_aln/Treated-G1-G2_S4_L001_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-2-D10Ag1g2/results/guide_aln/COL7A1-2-D10Ag1g2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("COL7A1-2-D10Ag1g2_RUN1", "COL7A1-2-D10Ag1g2_RUN2")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Nickase/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "COL7A1-2-D10Ag1g2_RUN1_with_RUN2.xlsx")

makeUpset("COL7A1-2-D10Ag1g2_RUN1_with_RUN2.xlsx")





# COL17A1-1-g1 with COL17A1-2-g1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-1-g1/results/guide_aln/COL17A1-1-g1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-2-g1/results/guide_aln/COL17A1-2-g1_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/overlap/Nickase/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "COL17A1-1-g1_with_COL17A1-2-g1.xlsx")

makeUpset("COL17A1-1-g1_with_COL17A1-2-g1.xlsx")


# COL17A1-1-g3 with COL17A1-2-g3
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-1-g3/results/guide_aln/COL17A1-1-g3_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-2-g3/results/guide_aln/COL17A1-2-g3_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/overlap/Nickase/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "COL17A1-1-g3_with_COL17A1-2-g3.xlsx")

makeUpset("COL17A1-1-g3_with_COL17A1-2-g3.xlsx")


# COL17A1-1-g1g3 with COL17A1-2-D10Ag1g3
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-1-g1g3/results/guide_aln/COL17A1-1-g1g3_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL17A1-2-D10Ag1g3/results/guide_aln/COL17A1-2-D10Ag1g3_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/overlap/Nickase/")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "COL17A1-1-g1g3_with_COL17A1-2-D10Ag1g3.xlsx")

makeUpset("COL17A1-1-g1g3_with_COL17A1-2-D10Ag1g3.xlsx")




###################################################
# RAG1 TAL 2 and TAL6 (03.03.21)

# RAG1 TAL2 run1 vs. run2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run1-TAL2/results/guide_aln/RAG1-run1-TAL2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run2-TAL2/results/guide_aln/RAG1-run2-TAL2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/RAG1_TAL2_TAL6/RAG1_TAL2_run1_run2/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "RAG1_TAL2_run1_run2.xlsx")

makeUpset("RAG1_TAL2_run1_run2.xlsx")

# RAG1 TAL6 run1 vs. run2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run1-TAL6/results/guide_aln/RAG1-run1-TAL6_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run2-TAL6/results/guide_aln/RAG1-run2-TAL6_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/RAG1_TAL2_TAL6/RAG1_TAL6_run1_run2/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "RAG1_TAL6_run1_run2.xlsx")

makeUpset("RAG1_TAL6_run1_run2.xlsx")



# RAG1 TAL2 run3 vs. run4
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run3-TAL2/results/guide_aln/RAG1-run3-TAL2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run4-TAL2/results/guide_aln/RAG1-run4-TAL2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/RAG1_TAL2_TAL6/RAG1_TAL2_run3_run4/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "RAG1_TAL2_run3_run4.xlsx")

makeUpset("RAG1_TAL2_run3_run4.xlsx")

# RAG1 TAL6 run3 vs. run4
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run3-TAL6/results/guide_aln/RAG1-run3-TAL6_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/RAG1_TAL2_TAL6/RAG1-run4-TAL6/results/guide_aln/RAG1-run4-TAL6_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/RAG1_TAL2_TAL6/RAG1_TAL6_run3_run4/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "RAG1_TAL6_run3_run4.xlsx")

makeUpset("RAG1_TAL6_run3_run4.xlsx")


########################################
# EMENDO 101


# Donor 1 Bio1 tech 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/results/guide_aln/EMD101-sample1-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-2/results/guide_aln/EMD101-sample1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/EMD101-sample1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample1.xlsx")
makeUpset("EMD101-sample1.xlsx")

# Donor 1 Bio2 tech 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample2-1/results/guide_aln/EMD101-sample2-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample2-2/results/guide_aln/EMD101-sample2-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/EMD101-sample2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample2.xlsx")
makeUpset("EMD101-sample2.xlsx")

# Donor 2 Bio1 tech 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample5-1/results/guide_aln/EMD101-sample5-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample5-2/results/guide_aln/EMD101-sample5-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/EMD101-sample5")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample5.xlsx")
makeUpset("EMD101-sample5.xlsx")

# Donor 2 Bio2 tech 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample6-1/results/guide_aln/EMD101-sample6-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample6-2/results/guide_aln/EMD101-sample6-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/EMD101-sample6")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample6.xlsx")
makeUpset("EMD101-sample6.xlsx")

# Donor 3 Bio1 tech 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample9-1/results/guide_aln/EMD101-sample9-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample9-2/results/guide_aln/EMD101-sample9-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/EMD101-sample9")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample9.xlsx")
makeUpset("EMD101-sample9.xlsx")

# Donor 3 Bio2 tech 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample10-1/results/guide_aln/EMD101-sample10-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample10-2/results/guide_aln/EMD101-sample10-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Emendo/samples_0221/EMD101-sample10")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample10.xlsx")
makeUpset("EMD101-sample10.xlsx")


# EMENDO PEAR vs. BBMERGE

# EMD101-sample1-1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/results/guide_aln/EMD101-sample1-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-1/results.save/guide_aln/EMD101-sample1-1_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample1-1_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample1-1_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample1-2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-2/results/guide_aln/EMD101-sample1-2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample1-2/results_save/guide_aln/EMD101-sample1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample1-2_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample1-2_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample2-1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample2-1/results/guide_aln/EMD101-sample2-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample2-1/results_save/guide_aln/EMD101-sample2-1_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample2-1_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample2-1_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample2-2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample2-2/results/guide_aln/EMD101-sample2-2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample2-2/results_save/guide_aln/EMD101-sample2-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample2-2_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample2-2_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample5-1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample5-1/results/guide_aln/EMD101-sample5-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample5-1/results_save/guide_aln/EMD101-sample5-1_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample5-1_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample5-1_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample5-2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample5-2/results/guide_aln/EMD101-sample5-2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample5-2/results_save/guide_aln/EMD101-sample5-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample5-2_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample5-2_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample6-1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample6-1/results/guide_aln/EMD101-sample6-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample6-1/results_save/guide_aln/EMD101-sample6-1_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample6-1_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample6-1_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample6-2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample6-2/results/guide_aln/EMD101-sample6-2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample6-2/results_save/guide_aln/EMD101-sample6-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample6-2_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample6-2_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample9-1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample9-1/results/guide_aln/EMD101-sample9-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample9-1/results_save/guide_aln/EMD101-sample9-1_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample9-1_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample9-1_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample9-2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample9-2/results/guide_aln/EMD101-sample9-2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample9-2/results_save/guide_aln/EMD101-sample9-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample9-2_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample9-2_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample10-1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample10-1/results/guide_aln/EMD101-sample10-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample10-1/results_save/guide_aln/EMD101-sample10-1_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample10-1_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample10-1_PEAR_vs_BBMERGE.xlsx")

# EMD101-sample10-2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample10-2/results/guide_aln/EMD101-sample10-2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Emendo/samples_0221/EMD101-sample10-2/results_save/guide_aln/EMD101-sample10-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("BBMERGE", "PEAR")

ovlDir <- file.path("~/Research/CASTSeq/overlap/Emendo")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "EMD101-sample10-2_PEAR_vs_BBMERGE.xlsx")
makeUpset("EMD101-sample10-2_PEAR_vs_BBMERGE.xlsx")

#################################################
# HBB Treated and i53 (Kai's email 03.03.21)


# 11 vs. 21
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_1_Treated/results/guide_aln/HBB-11_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_2_Treated/results/guide_aln/HBB-21_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/HBB-Treated_11_21")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "HBB-Treated_11_21.xlsx")
makeUpset("HBB-Treated_11_21.xlsx")


# 12 vs. 22
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_1_i53/results/guide_aln/HBB-12_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_2_i53/results/guide_aln/HBB-22_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[6])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/HBB-i53_12_22")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "HBB-i53_12_22.xlsx")
makeUpset("HBB-i53_12_22.xlsx")



#################################################
# MLL SAMPLES (Masako email 16 04 21)

# replicate 1 vs. 2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/MLL_160421/gMA003/results/guide_aln/gMA003treated_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/MLL_160421/gMA4-1/results/guide_aln/gMA4-1treated_w250_FINAL.xlsx")
)
names(siteFiles) <- c("gMA003", "gMA4-1")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/MLL_160421_OVL/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 0, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "gMA003_gMA4-1_OVL.xlsx")
makeUpset("gMA003_gMA4-1_OVL.xlsx")


#################################################
# HBB Treated and i53 (Masako's email 19.04.21)

# Treated
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_160421/1HBB-Treated/results/guide_aln/1HBB-Treated-Univ-T-RAG1BVEGFA_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_160421/2HBB-Treated/results/guide_aln/2HBB-Treated_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/HBB_160421/1HBB_2HBB-Treated/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "1HBB_2HBB-Treated.xlsx")
makeUpset("1HBB_2HBB-Treated.xlsx")


# i53
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_160421/1HBB-i53/results/guide_aln/1HBB-i53_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/HBB_160421/2HBB-i53/results/guide_aln/2HBB-i53_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/HBB_160421/1HBB_2HBB-i53/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "1HBB_2HBB-i53.xlsx")
makeUpset("1HBB_2HBB-i53.xlsx")


#################################################
# AZ_160421 (Julia's email 07.08.21)

# 8 vs. 11 (8 is the new 8 that has been re-analysed on June 22nd 2021)
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/AZ_160421/AZ-Liver8_3/results/guide_aln/AZ-Liver8_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/AZ_160421/AZ-Liver11_4/results/guide_aln/AZ-Liver11_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/AZ_160421/AZ-Liver8_11_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "AZ-Liver8_11_OVL1.xlsx")
makeUpset("AZ-Liver8_11_OVL1.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/AZ_160421/AZ-Liver8_11_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "AZ-Liver8_11_OVL2.xlsx")
makeUpset("AZ-Liver8_11_OVL2.xlsx")


# 16 vs. 18
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/AZ_160421/AZ-Liver16_19/results/guide_aln/AZ-Liver16_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/AZ_160421/AZ-Liver18_20/results/guide_aln/AZ-Liver18_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/AZ_160421/AZ-Liver16_18")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "AZ-Liver16_18.xlsx")
makeUpset("AZ-Liver16_18.xlsx")


ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/AZ_160421/AZ-Liver16_18_signif1/")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "AZ-Liver16_18_signif1.xlsx")
makeUpset("AZ-Liver16_18_signif1.xlsx")


#################################################
# TRAC (Manuel's email 25.05.21)

# KKR1 vs. KKR2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/TRAC/TRAC-KKR1/results/guide_aln/TRAC-KKR1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/TRAC/TRAC-KKR2/results/guide_aln/TRAC-KKR2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/TRAC/TRAC-KKR1_KKR2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "TRAC-KKR1_KKR2.xlsx")
makeUpset("TRAC-KKR1_KKR2.xlsx")


# KVR1 vs. KVR2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/TRAC/TRAC-KVR1/results/guide_aln/TRAC-KVR1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/TRAC/TRAC-KVR2/results/guide_aln/TRAC-KVR2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/TRAC/TRAC-KVR1_KVR2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "TRAC-KVR1_KVR2.xlsx")
makeUpset("TRAC-KVR1_KVR2.xlsx")


# WT1 vs. WT2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/TRAC/TRAC-WT1/results/guide_aln/TRAC-WT1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/TRAC/TRAC-WT2/results/guide_aln/TRAC-WT2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/TRAC/TRAC-WT1_WT2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "TRAC-WT1_WT2.xlsx")
makeUpset("TRAC-WT1_WT2.xlsx")


#################################################
# STAT3 (Manuel's email 25.05.21)

# T3
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/STAT3/STAT3-T3edited-1/results/guide_aln/STAT3-T3edited-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/STAT3/STAT3-T3edited-2/results/guide_aln/STAT3-T3edited-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/STAT3/STAT3-T3_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "STAT3-T3_OVL1.xlsx")
makeUpset("STAT3-T3_OVL1.xlsx")


# T6
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/STAT3/STAT3-T6edited-1/results/guide_aln/STAT3-T6edited-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/STAT3/STAT3-T6edited-2/results/guide_aln/STAT3-T6edited-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/STAT3/STAT3-T6_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "STAT3-T6_OVL1.xlsx")
makeUpset("STAT3-T6_OVL1.xlsx")


#################################################
# KRT9_gel

# T1 T2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9_gel/KRT9-T1gel_g1/results/guide_aln/KRT9-T1gel_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9_gel/KRT9-T2gel_g1/results/guide_aln/KRT9-T2gel_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/KRT9_gel/KRT9-T1T2gel_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "KRT9-T1T2gel.xlsx")
makeUpset("KRT9-T1T2gel.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/KRT9_gel/KRT9-T1T2gel_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "KRT9-T1T2gel.xlsx")
makeUpset("KRT9-T1T2gel.xlsx")


# T3 T4
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9_gel/KRT9-T3gel_g2/results/guide_aln/KRT9-T3gel_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_180521/KRT9_gel/KRT9-T4gel_g2/results/guide_aln/KRT9-T4gel_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[8])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/KRT9_gel/KRT9-T3T4gel_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "KRT9-T3T4gel.xlsx")
makeUpset("KRT9-T3T4gel.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_180521/KRT9_gel/KRT9-T3T4gel_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "KRT9-T3T4gel.xlsx")
makeUpset("KRT9-T3T4gel.xlsx")


#################################################
# LAMA3 22.06.21

# g1-1 vs. g1-2
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-Cas9g1-1/results/guide_aln/Lama3-Cas9g1-1_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-Cas9g1-2/results/guide_aln/Lama3-Cas9g1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-Cas9g1_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Lama3-Cas9g1_d1500_OVL2.xlsx", overwrite = TRUE)
makeUpset("Lama3-Cas9g1_d1500_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-Cas9g1_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Lama3-Cas9g1_d1500_OVL1.xlsx", overwrite = TRUE)
makeUpset("Lama3-Cas9g1_d1500_OVL1.xlsx")


# g2-1 vs. g2-2
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-Cas9g2-1/results/guide_aln/Lama3-Cas9g2-1_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-Cas9g2-2/results/guide_aln/Lama3-Cas9g2-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-Cas9g2_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Lama3-Cas9g2_d1500_OVL2.xlsx", overwrite = TRUE)
makeUpset("Lama3-Cas9g2_d1500_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-Cas9g2_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Lama3-Cas9g2_d1500_OVL1.xlsx", overwrite = TRUE)
makeUpset("Lama3-Cas9g2_d1500_OVL1.xlsx")

# D10Ag1g2-1-1 vs. D10Ag1g2-1-2
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-D10Ag1g2-1-1/results/guide_aln/Lama3-D10Ag1g2-1-1_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-D10Ag1g2-1-2/results/guide_aln/Lama3-D10Ag1g2-1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-D10Ag1g2-1_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Lama3-D10Ag1g2-1_d1500_OVL2.xlsx", overwrite = TRUE)
makeUpset("Lama3-D10Ag1g2-1_d1500_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-D10Ag1g2-1_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Lama3-D10Ag1g2-1_d1500_OVL1.xlsx", overwrite = TRUE)
makeUpset("Lama3-D10Ag1g2-1_d1500_OVL1.xlsx")

# D10Ag1g2-2-1 vs. D10Ag1g2-2-2
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-D10Ag1g2-2-1/results/guide_aln/Lama3-D10Ag1g2-2-1_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-D10Ag1g2-2-2/results/guide_aln/Lama3-D10Ag1g2-2-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-D10Ag1g2-2_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Lama3-D10Ag1g2-2_d1500_OVL2.xlsx", overwrite = TRUE)
makeUpset("Lama3-D10Ag1g2-2_d1500_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-D10Ag1g2-2_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Lama3-D10Ag1g2-2_d1500_OVL1.xlsx", overwrite = TRUE)
makeUpset("Lama3-D10Ag1g2-2_d1500_OVL1.xlsx")


# D10Ag1g2-1-1 vs. D10Ag1g2-1-2 NEW INPUT FILES
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-D10Ag1g2-1-1_NEW/results/guide_aln/Lama3-D10Ag1g2-1-1_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/LAMA3/Lama3-D10Ag1g2-1-2_NEW/results/guide_aln/Lama3-D10Ag1g2-1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-D10Ag1g2-1_NEW_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Lama3-D10Ag1g2-1_d1500_OVL2.xlsx", overwrite = TRUE)
makeUpset("Lama3-D10Ag1g2-1_d1500_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/LAMA3/Lama3-D10Ag1g2-1_NEW_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Lama3-D10Ag1g2-1_d1500_OVL1.xlsx", overwrite = TRUE)
makeUpset("Lama3-D10Ag1g2-1_d1500_OVL1.xlsx")

###########################################################################
# MASAKO (email June 28th 2021)


# HAX1 HD19 vs. HD24 vs. HD33
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/HAX1/HAX1-HD19/results/guide_aln/HAX1-HD19-treat_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/HAX1/HAX1-HD24/results/guide_aln/HAX1-HD24-treat_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/HAX1/HAX1-HD33/results/guide_aln/HAX1-HD33-treat_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/HAX1/HAX1-HD19_22_34_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "HAX1-HD19_22_34_OVL2.xlsx", overwrite = TRUE)
makeUpset("HAX1-HD19_22_34_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/HAX1/HAX1-HD19_22_34_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "HAX1-HD19_22_34_OVL1.xlsx", overwrite = TRUE)
makeUpset("HAX1-HD19_22_34_OVL1.xlsx")

# HAX1 HD19 vs. HD24 vs. HD33 (PSEUDOGENE)
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/HAX1/HAX1-HD19_pseudo/results/guide_aln/HAX1-HD19-treat_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/HAX1/HAX1-HD24_pseudo/results/guide_aln/HAX1-HD24-treat_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/HAX1/HAX1-HD33_pseudo/results/guide_aln/HAX1-HD33-treat_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[9])

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/HAX1/HAX1-HD19_22_34_pseudo_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "HAX1-HD19_22_34_pseudo_OVL2.xlsx", overwrite = TRUE)
makeUpset("HAX1-HD19_22_34_pseudo_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/HAX1/HAX1-HD19_22_34_pseudo_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "HAX1-HD19_22_34_pseudo_OVL1.xlsx", overwrite = TRUE)
makeUpset("HAX1-HD19_22_34_pseudo_OVL1.xlsx")


# Universal CASTseq CCR5 with CCR5 gRNA
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/UnivCASTseq_160421/CCR5_2/results/guide_aln/Univ-T-CCR5VEGFA_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/UniversalCASTseq/Universal-CCR5-VEGFA_with_CCR5_gRNA/results/guide_aln/Universal-CCR5-VEGFA_w250_FINAL.xlsx")
)
names(siteFiles) <- c("R1", "R2")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/UniversalCASTseq/Universal-CCR5-VEGFA_with_CCR5_gRNA_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Universal-CCR5-VEGFA_with_CCR5_gRNA_OVL2.xlsx", overwrite = TRUE)
makeUpset("Universal-CCR5-VEGFA_with_CCR5_gRNA_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/UniversalCASTseq/Universal-CCR5-VEGFA_with_CCR5_gRNA_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Universal-CCR5-VEGFA_with_CCR5_gRNA_OVL1.xlsx", overwrite = TRUE)
makeUpset("Universal-CCR5-VEGFA_with_CCR5_gRNA_OVL1.xlsx")


# Universal CASTseq CCR5 with VEGFA gRNA
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/UnivCASTseq_160421/CCR5_2_VEGFA/results/guide_aln/Univ-T-CCR5VEGFA_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_170621/UniversalCASTseq/Universal-CCR5-VEGFA_with_VEGFA_gRNA/results/guide_aln/Universal-CCR5-VEGFA_w250_FINAL.xlsx")
)
names(siteFiles) <- c("R1", "R2")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/UniversalCASTseq/Universal-CCR5-VEGFA_with_VEGFA_gRNA_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "Universal-CCR5-VEGFA_with_VEGFA_gRNA_OVL2.xlsx", overwrite = TRUE)
makeUpset("Universal-CCR5-VEGFA_with_VEGFA_gRNA_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_170621/UniversalCASTseq/Universal-CCR5-VEGFA_with_VEGFA_gRNA_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "Universal-CCR5-VEGFA_with_VEGFA_gRNA_OVL1.xlsx", overwrite = TRUE)
makeUpset("Universal-CCR5-VEGFA_with_VEGFA_gRNA_OVL1.xlsx")


###############################################################################
# MASAKO EMAIL 28.07.21

# Universal CASTseq RAG1B VEGFA with RAG1B gRNA 
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_160421/GTGAA/RAG1B_gRNA/results/guide_aln/GTGAAA_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_260721/RAG1B-VEGFA/RAG1B_gRNA/results/guide_aln/RAG1B-VEGFA-treated_w250_FINAL.xlsx")
)
names(siteFiles) <- c("R1", "R2")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_260721/RAG1B-VEGFA/RAG1B_gRNA_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "RAG1B-VEGFA-treated_RAG1B_gRNA_OVL2.xlsx", overwrite = TRUE)
makeUpset("RAG1B-VEGFA-treated_RAG1B_gRNA_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_260721/RAG1B-VEGFA/RAG1B_gRNA_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "RAG1B-VEGFA-treated_RAG1B_gRNA_OVL1.xlsx", overwrite = TRUE)
makeUpset("RAG1B-VEGFA-treated_RAG1B_gRNA_OVL1.xlsx")


# Universal CASTseq RAG1B VEGFA with VEGFA gRNA 
siteFiles <- c(file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_160421/GTGAA/VEGFA_gRNA/results/guide_aln/GTGAAA_w250_FINAL.xlsx"),
               file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples/data_260721/RAG1B-VEGFA/VEGFA_gRNA/results/guide_aln/RAG1B-VEGFA-treated_w250_FINAL.xlsx")
)
names(siteFiles) <- c("R1", "R2")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_260721/RAG1B-VEGFA/VEGFA_gRNA_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "RAG1B-VEGFA-treated_VEGFA_gRNA_OVL2.xlsx", overwrite = TRUE)
makeUpset("RAG1B-VEGFA-treated_VEGFA_gRNA_OVL2.xlsx")

ovlDir <- file.path("~/cluster/cluster/CASTSeq/pipelineGit/samples_overlap/data_260721/RAG1B-VEGFA/VEGFA_gRNA_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "RAG1B-VEGFA-treated_VEGFA_gRNA_OVL1.xlsx", overwrite = TRUE)
makeUpset("RAG1B-VEGFA-treated_VEGFA_gRNA_OVL1.xlsx")


###############################################################################
# VIVIANE EMAIL 29.07.21

# UNC-u1d1 D1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u1d1-1_D1/results/guide_aln/UNC-u1d1-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u1d1-2_D1/results/guide_aln/UNC-u1d1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("UNC-u1d1-1_D1", "UNC-u1d1-2_D1")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u1d1_D1_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "UNC-u1d1_D1_OVL2.xlsx", overwrite = TRUE)
makeUpset("UNC-u1d1_D1_OVL2.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u1d1_D1_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "UNC-u1d1_D1_OVL1.xlsx", overwrite = TRUE)
makeUpset("UNC-u1d1_D1_OVL1.xlsx")

# UNC-u1d1 U1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u1d1-1_U1/results/guide_aln/UNC-u1d1-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u1d1-2_U1/results/guide_aln/UNC-u1d1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("UNC-u1d1-1_U1", "UNC-u1d1-2_U1")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u1d1_U1_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "UNC-u1d1_U1_OVL2.xlsx", overwrite = TRUE)
makeUpset("UNC-u1d1_U1_OVL2.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u1d1_U1_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "UNC-u1d1_U1_OVL1.xlsx", overwrite = TRUE)
makeUpset("UNC-u1d1_U1_OVL1.xlsx")


# UNC-u3d1 D1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u3d1-1_D1/results/guide_aln/UNC-u3d1-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u3d1-2_D1/results/guide_aln/UNC-u3d1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("UNC-u3d1-1_D1", "UNC-u3d1-2_D1")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u3d1_D1_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "UNC-u3d1_D1_OVL2.xlsx", overwrite = TRUE)
makeUpset("UNC-u3d1_D1_OVL2.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u3d1_D1_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "UNC-u3d1_D1_OVL1.xlsx", overwrite = TRUE)
makeUpset("UNC-u3d1_D1_OVL1.xlsx")


# UNC-u3d1 U3
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u3d1-1_U3/results/guide_aln/UNC-u3d1-1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/data_170621/UNC/UNC-u3d1-2_U3/results/guide_aln/UNC-u3d1-2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("UNC-u3d1-1_U3", "UNC-u3d1-2_U3")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u3d1_U3_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "UNC-u3d1_U3_OVL2.xlsx", overwrite = TRUE)
makeUpset("UNC-u3d1_U3_OVL2.xlsx")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/data_170621/UNC/UNC-u3d1_U3_OVL1")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 1500, nb.signif = 1, NBS = TRUE)
write.xlsx(compList, "UNC-u3d1_U3_OVL1.xlsx", overwrite = TRUE)
makeUpset("UNC-u3d1_U3_OVL1.xlsx")



################################
# NICKASE OVL2 (JULIA EMAIL 12.10.21)

# COL7A1-1-g1 with COL7A1-2-g1
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-1-g1/results/guide_aln/COL7A1-1-g1_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-2-g1/results/guide_aln/COL7A1-2-g1_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Nickase_OVL2/COL7A1-g1_OVL2")
dir.create(ovlDir, recursive = TRUE)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "COL7A1-g1_OVL2.xlsx")

makeUpset("COL7A1-g1_OVL2.xlsx")

# COL7A1-1-g2 with COL7A1-2-g2
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-1-g2/results/guide_aln/COL7A1-1-g2_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-2-g2/results/guide_aln/COL7A1-2-g2_w250_FINAL.xlsx")
)
names(siteFiles) <- sapply(strsplit(siteFiles, split = "/"), function(i) i[7])

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Nickase_OVL2/COL7A1-g2_OVL2")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "COL7A1-g2_OVL2.xlsx")

makeUpset("COL7A1-g2_OVL2.xlsx")


# COL7A1-2-D10Ag1g2 (first vs. second dataset)
siteFiles <- c(file.path("~/Research/CASTSeq/pipelineGit/samples/COL7A1_G1_G2/results_middle/guide_aln/Treated-G1-G2_S4_L001_w250_FINAL.xlsx"),
               file.path("~/Research/CASTSeq/pipelineGit/samples/Nickase_Jan21/COL7A1-2-D10Ag1g2/results/guide_aln/COL7A1-2-D10Ag1g2_w250_FINAL.xlsx")
)
names(siteFiles) <- c("COL7A1-2-D10Ag1g2_RUN1", "COL7A1-2-D10Ag1g2_RUN2")

ovlDir <- file.path("~/Research/CASTSeq/pipelineGit/samples_overlap/Nickase_OVL2/COL7A1-2-D10Ag1g2_OVL2")
dir.create(ovlDir)
setwd(ovlDir)
compList <- dfComparisonList(siteFiles, names(siteFiles), width = 2000, nb.signif = 2, NBS = TRUE)
write.xlsx(compList, "COL7A1-2-D10Ag1g2_OVL2.xlsx", overwrite = TRUE)

makeUpset("COL7A1-2-D10Ag1g2_OVL2.xlsx")





