
runPipelineDOUBLENICKASE <- function()
{
  
  ##########################################################################################
  ############                           RUN PIPELINE                           ############
  ##########################################################################################
  
  ################     CHECK INPUT    ################ 
  print("################     CHECK INPUT    ################")
  
  print(paste0("sampleD: ", sampleD))

  
  ##########################################################################################
  ############                        FASTQ ALIGNMENT                         ############
  ##########################################################################################
  
  ################     FASTQ ALIGNMENT    ################ 
  print("################     FASTQ ALIGNMENT    ################")
  
  # samples.treated
  samples.treated.unique <- samples.treated[!duplicated(samples.treated)]
  labels.treated.unique <- labels.treated[!duplicated(labels.treated)]
  fastqD.treated.unique <- fastqD.treated[!duplicated(labels.treated)]
  
  lapply(1:length(samples.treated.unique), function(i){
    
    resultD.current <- file.path(sampleD, labels.treated.unique[i], "fastq_aln")
    if(!file.exists(file.path(resultD.current, paste0(samples.treated.unique[i], "_Alignment.bed")))){
      # create output directory
      dir.create(resultD.current, showWarnings = FALSE, recursive = TRUE)
      
      fastqAln(functionstring=paste0("sh ", file.path(scriptD, "fastq_aln.sh")),
               homeFolder = homeD,
               annotFolder = annotD,
               sampleFolder = sampleD,
               fastqFolder = fastqD.treated.unique[i],
               testSample = samples.treated.unique[i],
               r1 = fastqExt1,
               r2 = fastqExt2,
               outFolder = resultD.current,
               cpu = NBCPU
      )	
      
    }else{
      print(paste0("skip fastq alignment because ",
                   file.path(resultD.current, paste0(samples.treated[i], "_Alignment.bed")),
                   " bed file already exists"))
    }
  
    return(NA)
  })
  
  # samples.untreated
  samples.untreated.unique <- samples.untreated[!duplicated(samples.untreated)]
  labels.untreated.unique <- labels.untreated[!duplicated(labels.untreated)]
  fastqD.untreated.unique <- fastqD.untreated[!duplicated(labels.untreated)]
  
  lapply(1:length(samples.untreated.unique), function(i){
    
    resultD.current <- file.path(sampleD, labels.untreated.unique[i], "fastq_aln")
    if(!file.exists(file.path(resultD.current, paste0(samples.untreated.unique[i], "_Alignment.bed")))){
      # create output directory
      dir.create(resultD.current, showWarnings = FALSE, recursive = TRUE)
      
      fastqAln(functionstring=paste0("sh ", file.path(scriptD, "fastq_aln.sh")),
               homeFolder = homeD,
               annotFolder = annotD,
               sampleFolder = sampleD,
               fastqFolder = fastqD.untreated.unique[i],
               testSample = samples.untreated.unique[i],
               r1 = fastqExt1,
               r2 = fastqExt2,
               outFolder = resultD.current,
               cpu = NBCPU
      )	
      
    }else{
      print(paste0("skip fastq alignment because ",
                   file.path(resultD.current, paste0(samples.untreated.unique[i], "_Alignment.bed")),
                   " bed file already exists"))
    }
    
    return(NA)
  })

  
  ################     QC: COUNT READS    ################ 
  print("################     QC: COUNT READS    ################")
  qcDir <- file.path(sampleD, "QC")
  dir.create(qcDir, showWarnings = FALSE)
  
  qcFile <- file.path(qcDir, paste0(projectName, "alignment_QC.xlsx"))
  if(file.exists(file.path(qcDir, "alignment_QC.xlsx"))){
    print("skip fastq alignment because alignment_QC.xlsx already exists")
  }else{
    countReads(fastqFolder = fastqD,
               sampleFolder = sampleD,
               testSample = c(samples.treated.unique, samples.untreated.unique),
               testLabel = c(labels.treated.unique, labels.untreated.unique),
               ext1 = fastqExt1,
               ext2 = fastqExt2,
               outputFile = file.path(qcDir, "alignment_QC.xlsx")
               )
    
    countReadsPlot(file.path(qcDir, "alignment_QC.xlsx"),
                   file.path(qcDir, "alignment_QC.pdf"))
  }
  

  ##########################################################################################
  ############                        HITS CALCULATION                          ############
  ##########################################################################################

  lapply(1:length(samples.treated), function(i){
    
    resultD <- file.path(sampleD, labels.treated[i], "hits")
    dir.create(resultD, showWarnings = FALSE)
    
    sampleName <- samples.treated[i]
    controlName <- samples.untreated[i]
    

    ################     DELTA AND HITS    ################ 
    print("################     DELTA AND HITS    ################ ")
    print(paste0(sampleName, " vs. ", controlName))
    
    # test sample
    sampleBed <- file.path(sampleD, labels.treated[i], "fastq_aln", paste0(sampleName, "_Alignment.bed"))
    
    # control sample		   
    controlBed <- file.path(sampleD, labels.untreated[i], "fastq_aln", paste0(controlName, "_Alignment.bed"))	
    
    # run delta analysis
    getDelta(sampleBed, resultD, otsBed, otsDistance, distance = distance.cutoff)		   
    getDeltaShuffle(file.path(resultD, paste0(sampleName, "_delta.bed")), resultD, nb = 10, distance = distance.cutoff, genome.size = myGenome.size)				   
    getDelta(controlBed, resultD, otsBed, otsDistance, distance = distance.cutoff)		   
    getDeltaShuffle(file.path(resultD, paste0(controlName, "_delta.bed")), resultD, nb = 10, distance = distance.cutoff, genome.size = myGenome.size)		
    
    # run hits analysis
    getHits(file.path(resultD, paste0(sampleName, "_delta.bed")), distance.cutoff)
    getHits(file.path(resultD, paste0(controlName, "_delta.bed")), distance.cutoff)
    
    ################     CALCULATE LIBRARY SIZE    ################
    print("################     CALCULATE LIBRARY SIZE    ################")
    
    # USE RAW FASTQ FILE !!!!!!!!!
    nbReads.sample <- as.numeric(nbReadFastqgz(file.path(fastqD.treated[i], paste0(sampleName, fastqExt2))))
    if(nbReads.sample == 0) print(paste0("NO reads in ", sampleName))
    
    nbReads.control <- as.numeric(nbReadFastqgz(file.path(fastqD.untreated[i], paste0(controlName, fastqExt2))))
    if(nbReads.control == 0) print(paste0("NO reads in ", controlName))
    
    ################     TEST VS. CONTROL ENRICHMENT    ################
    print("################     TEST VS. CONTROL ENRICHMENT    ################")
    
    doEnrichment(file.path(resultD, paste0(sampleName, "_hits.bed")),
                 file.path(resultD, paste0(controlName, "_hits.bed")),
                 nbReads.sample, nbReads.control, w, myGenome.size)
    
    # ADD AVERAGE MAPQ
    addMAPQ(file.path(resultD, paste0(sampleName, "_w", w, ".xlsx")),
            file.path(sampleD, labels.treated[i], "fastq_aln", paste0(sampleName, "_AlignmentSort.bam")))
    
    
    ############           ON-TARGET READOUT             ############
    print("############    ON-TARGET READOUT    ############")
    
    ONreadout(bamFile = file.path(sampleD, labels.treated[i], "fastq_aln", paste0(sampleName, "_AlignmentSort.bam")),
              otsFile = otsBed,
              gRNA.orientation = bpOrient,
              window.size = 5000,
              sampleName = sampleName,
              outDir = resultD)
    
  })
  

  ##########################################################################################
  ############                         OVERLAP ANALYSIS                         ############
  ##########################################################################################
  print("################    OVERLAP ANALYSIS    ################")
  
  ovlD <- file.path(sampleD, paste0("OVL", nb.ovl, "_SIGNIF", nb.signif))
  dir.create(ovlD, showWarnings = FALSE)
  
  siteFiles <- sapply(1:length(labels.treated), function(i)
    file.path(sampleD, labels.treated[i], "hits", paste0(samples.treated[i], "_w", w, ".xlsx")))
  
  #compList <- dfComparisonList(siteFiles, labels.treated, width = distance.cutoff, nb.signif = 1, NBS = TRUE)
  #write.xlsx(compList, file.path(ovlD, paste0(projectName, "_RAW.xlsx")), overwrite = TRUE)
  #makeUpset(file.path(ovlD, paste0(projectName, "_RAW.xlsx")))
  
  getUnion(siteFiles, file.path(ovlD, paste0(projectName, "_RAW.xlsx")),
           nList = labels.treated, width = distance.cutoff, NBS = TRUE, adjusted = TRUE)
  # TO DO
  # makeUpset(file.path(ovlD, paste0(projectName, "_RAW.xlsx")))
  
  selectSite(file.path(ovlD, paste0(projectName, "_RAW.xlsx")),
             file.path(ovlD, paste0(projectName, "_SELECTED.xlsx")),
             nbOvl = nb.ovl, nbSignif = nb.signif)
  
  
  ##########################################################################################
  ############                    BARCODE HOPING FILTER                         ############
  ##########################################################################################
  
  ################     BARCODE HOPPING FILTER    ################ 
  print("################     BARCODE HOPPING FILTER    ################")
  
  barcodeHoppingFilter(inputFile = file.path(ovlD, paste0(projectName, "_SELECTED.xlsx")),
                       outputFile = file.path(ovlD, paste0(projectName, "_BH.xlsx")),
                       hardCoef = 2.5, softCoef = 2)
  

  ##########################################################################################
  ############                        COVERAGE PER SITE                         ############
  ##########################################################################################
  print("############    COVERAGE PER SITE    ############")
  
  covD <- file.path(sampleD, "coverage")
  dir.create(covD, showWarnings = FALSE)
  
  bamFiles <- sapply(1:length(labels.treated), function(i)
    file.path(sampleD, labels.treated[i], "fastq_aln", paste0(samples.treated[i], "_AlignmentSort.bam")))

  coverage_bin(inputFile = file.path(ovlD, paste0(projectName, "_BH.xlsx")),
               outputFile = file.path(ovlD, paste0(projectName, "_BIN.xlsx")),
               bamFiles = bamFiles,
               bamLabels = labels.treated,
               bin.size = 100,
               top = "sum",
               d = distance.cov,
               plot = FALSE,
               outDir = covD)
  
  ##########################################################################################
  ############                          RANDOM SITES                            ############
  ##########################################################################################
  print("############    RANDOM SITES    ############")
  
  # GENERATE RANDOM SEQUENCE BED
  randomD <- file.path(sampleD, "random")
  dir.create(randomD, showWarnings = FALSE)
  randomName <- paste0("random_w", distance.cov)
  
  # Check if bed file already exists, if yes skip the analysis
  if(!file.exists(file.path(randomD, paste0(randomName, "_ALN.xlsx")))){
    getRandomBed(l = distance.cov*2, n = nb.rd,
                 outFile = file.path(randomD, paste0(randomName, ".bed")),
                 opt.string = paste("-g", myGenome.size, sep = " ")
                 )
    
    # DO GUIDE ALIGNMENT ON RANDOM SEQUENCES
    getGuideAlignment2(inputF = file.path(randomD, paste0(randomName, ".bed")),
                       outputF = file.path(randomD, paste0(randomName, "_ALN.xlsx")),
                       guideLeft = refSeq1, guideRight = refSeq2,
                       gnm = GNM,
                       binCoord = FALSE
                       )
    
    # CALCULATE DISTANCE BETWEEN THE 2 gRNAs
    getgRNADistance(inputF = file.path(randomD, paste0(randomName, "_ALN.xlsx")),
                    outputF = file.path(randomD, paste0(randomName, "_ALN.xlsx")),
                    guide1 = refSeq1, guide2 = refSeq2)
    
    # ADD CUMULATIVE SCORE
    getCumulScore(inputF = file.path(randomD, paste0(randomName, "_ALN.xlsx")),
                  outputF = file.path(randomD, paste0(randomName, "_ALN.xlsx")))

  }else{print(paste0("skip random gRNA alignment because ", file.path(randomD, paste0(randomName, "_ALN.xlsx")), " already exists"))}
  
  # FLANKING REGIONS
  # Check if output file already exists, if yes skip the analysis
  if(!file.exists(file.path(randomD, paste0(randomName, "_FLANK.xlsx")))){
    if(!(is.null(flank1.sq))){
      addFlankingFromSq(file.path(randomD, paste0(randomName, "_ALN.xlsx")),
                        file.path(randomD, paste0(randomName, "_FLANK.xlsx")),
                        hom.sq = flank1.sq, hom.sq2 = flank2.sq)
    }else{
      addFlanking(file.path(randomD, paste0(randomName, "_ALN.xlsx")),
                  file.path(randomD, paste0(randomName, "_FLANK.xlsx")),
                  otsBed, flankingSize)
    }
  }else{print(paste0("skip random flanking length calculation because ", file.path(randomD, paste0(randomName, "_FLANK.xlsx")), " already exists"))}
  
  # ANNOTATE GENES
  if(!file.exists(file.path(randomD, paste0(randomName, "_GENES.xlsx")))){
    annotateGeneTALEN(file.path(randomD, paste0(randomName, "_FLANK.xlsx")),
                      file.path(randomD, paste0(randomName, "_GENES.xlsx")))
  }else{print(paste0("skip random flanking length calculation because ", file.path(randomD, paste0(randomName, "_GENES.xlsx")), " already exists"))}
  
  

  ##########################################################################################
  ############                DESIGNER NUCLEASE TREATED SAMPLE                  ############
  ##########################################################################################
  
  ################     GUIDE SEQ ALIGNMENT    ################ 
  print("################     GUIDE SEQ ALIGNMENT    ################")
  
  # DO GUIDE ALIGNMENT ON REAL SEQUENCES
  getGuideAlignment2(inputF = file.path(ovlD, paste0(projectName, "_BIN.xlsx")),
                     outputF = file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
                     guideLeft = refSeq1, guideRight = refSeq2,
                     gnm = GNM,
                     binCoord = TRUE
  )			  
  #file.remove(list.files(ovlD, pattern = "_TMP.txt", full.names = TRUE))					  

  # CALCULATE DISTANCE BETWEEN THE 2 gRNAs
  getgRNADistance(inputF = file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
                  outputF = file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
                  guide1 = refSeq1, guide2 = refSeq2)
  
  # ADD CUMULATIVE SCORE
  getCumulScore(inputF = file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
                outputF = file.path(ovlD, paste0(projectName, "_ALN.xlsx")))
  
  # PLOT GUIDE ALIGNMENT
  guidePlot1(file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
            file.path(ovlD, paste0(projectName, "_ALN_heatmap.gRNA1.pdf")),
            score = NULL, pv = NULL, ref = refSeq1)# ALL
  
  guidePlot2(file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
            file.path(ovlD, paste0(projectName, "_ALN_heatmap.gRNA2.pdf")),
            score = NULL, pv = NULL, ref = refSeq2)# ALL
  
  # LOGO PLOT
  logoPlot1(file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
           file.path(ovlD, paste0(projectName, "_ALN_logo.gRNA1.pdf")),
           score = NULL, pv = NULL, ref = refSeq1)# ALL
  
  logoPlot2(file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
            file.path(ovlD, paste0(projectName, "_ALN_logo.gRNA2.pdf")),
            score = NULL, pv = NULL, ref = refSeq2)# ALL
  
  ################     FLANKING REGIONS    ################ 
  print("################     FLANKING REGIONS    ################ ")
  # REAL SEQUENCES
  if(is.null(flank1.sq)){
    addFlanking(file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
                file.path(ovlD, paste0(projectName, "_FLANK.xlsx")),
                otsBed, flankingSize)
  }else{
    addFlankingFromSq(file.path(ovlD, paste0(projectName, "_ALN.xlsx")),
                      file.path(ovlD, paste0(projectName, "_FLANK.xlsx")),
                      flank1.sq, flank2.sq)
  }
  
  ################     DEFINE GROUPS    ################ 
  print("################     DEFINE GROUPS    ################ ")
  assignGroupsDoubleNickase(file.path(ovlD, paste0(projectName, "_FLANK.xlsx")),
                            file.path(randomD, paste0(randomName, "_FLANK.xlsx")),
                            file.path(ovlD, paste0(projectName, "_GROUP.xlsx")),
                            otsBed,
                            pv.cutoff)
  
  groupSummary(file.path(ovlD, paste0(projectName, "_GROUP.xlsx")),
               file.path(ovlD, paste0(projectName, "_group_summary.xlsx")),
               hits = NULL,
               score = NULL, pv = NULL
  )	
  
  
  ################     GENE ANNOTATION    ################ 
  print("################     GENE ANNOTATION    ################ ")
  # Annotation
  annotateGene(file.path(ovlD, paste0(projectName, "_GROUP.xlsx")),
               file.path(ovlD, paste0(projectName, "_GENES.xlsx")))
  
  # Barplot per group
  geneBarplot(file.path(ovlD, paste0(projectName, "_GENES.xlsx")))
  
  # Forest plot per group
  geneForestPlot(file.path(ovlD, paste0(projectName, "_GENES.xlsx")),
                 file.path(randomD, paste0(randomName, "_GENES.xlsx")))
  
  ################     ALL GENES WITHIN 100KB ANNOTATION    ################ 
  print("################     ALL GENES WITHIN 100KB ANNOTATION    ################")
  addGenes(file.path(ovlD, paste0(projectName, "_GENES.xlsx")),
           file.path(ovlD, paste0(projectName, "_GENES2.xlsx")),
           oncoFile = oncoEntrez,
           geneMat = geneMat,
           genes.width = 0, site.width = 100000)
  
  ################     RETURN FINAL XLSX FILE    ################ 
  print("################     RETURN FINAL XLSX FILE    ################")
  finalize(file.path(ovlD, paste0(projectName, "_GENES2.xlsx")),
           file.path(ovlD, paste0(projectName, "_FINAL.xlsx")))	

  
  ##########################################################################################
  ############                           HITS BARPLOT                           ############
  ##########################################################################################
  print("############    HITS BARPLOT    ############")
  
  hitsBarplot(file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
              pv = NULL, top = 50, showNBS = TRUE, log = TRUE,
              outName = file.path(ovlD, paste0(projectName, "_hits_barplot.pdf")))
  hitsBarplot(file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
              pv = NULL, top = 50, showNBS = FALSE, log = TRUE,
              outName = file.path(ovlD, paste0(projectName, "_hits_barplot_woNBS.pdf")))
  
  ##########################################################################################
  ############               COVERAGE PLOT (WITH CLEAVAGE SITE)                 ############
  ##########################################################################################
  print("############    COVERAGE PLOT (WITH CLEAVAGE SITE    ############")
  
  coverage_cleavage_double(inputFile = file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
                          covDir = covD,
                          bin.size = 100)


  ##########################################################################################
  ############                        CHR PLOT (CIRCLIZE)                       ############
  ##########################################################################################
  print("############    CHR PLOT (CIRCLIZE)    ############")
  
  # ZOOM
  circlizePipeline(siteFile = file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
                   zoom.size = 25000, label = FALSE, 
                   PV.cutoff = NULL,
                   bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                   showNBS = TRUE,
                   gene.bed = NULL, ots.bed = otsBed, 
                   realigned = TRUE,
                   outFile = file.path(ovlD, paste0(projectName, "_circlize_25k.pdf")),
                   species = circos.sp,
                   top = NULL)
  
  circlizePipeline(siteFile = file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
                   zoom.size = 25000, label = FALSE, 
                   PV.cutoff = NULL,
                   bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                   showNBS = FALSE,
                   gene.bed = NULL, ots.bed = otsBed, 
                   realigned = TRUE,
                   outFile = file.path(ovlD, paste0(projectName, "_circlize_25k_woNBS.pdf")),
                   species = circos.sp,
                   top = NULL)
  
  # NO ZOOM
  circlizePipelineNOZOOM(siteFile = file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
                         label = FALSE, 
                         PV.cutoff = NULL,
                         bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                         showNBS = TRUE,
                         ots.bed = otsBed, 
                         realigned = TRUE,
                         outFile = file.path(ovlD, paste0(projectName, "_circlize_25k_NOZOOM.pdf")),
                         species = circos.sp,
                         top = NULL)
  
  circlizePipelineNOZOOM(siteFile = file.path(ovlD, paste0(projectName, "_FINAL.xlsx")),
                         label = FALSE, 
                         PV.cutoff = NULL,
                         bestScore.cutoff = NULL, bestFlank.cutoff = 25,
                         showNBS = FALSE,
                         ots.bed = otsBed, 
                         realigned = TRUE,
                         outFile = file.path(ovlD, paste0(projectName, "_circlize_25k_woNBS_NOZOOM.pdf")),
                         species = circos.sp,
                         top = NULL)
  
  ##########################################################################################
  ############                        REMOVE TMP FILES                          ############
  ##########################################################################################
  
  if(rmTMP){
    print("############    REMOVE TMP FILES    ############")
    mypattern <- c("_assembled.fastq.gz$", "_Filt2.fastq.gz$", "_Filt3_neg.fastq.gz$",
                   "_Filt3_pos.fastq.gz$", "_final.fastq.gz$", "_pos_NOmatch.fastq.gz$",
                   "_pos.fastq.gz$", "_trim1.fastq.gz$", "_trim2.fastq.gz$", "_trim3.fastq.gz$",
                   "_unassembled.R1.fastq.gz$", "_unassembled.R2.fastq.gz$",
                   "_merged.fastq.gz$",
                   "_BIN.xlsx$", "_BH.xlsx$", "_ALN.xlsx$", "_ALN_PV.xlsx$",
                   "_FLANK.xlsx$", "_GROUP.xlsx$", "_GENES.xlsx$", "_GENES2.xlsx$")
    torm <- unlist(lapply(mypattern, function(mp) list.files(sampleD, pattern = mp, full.names = TRUE, recursive = TRUE)))
    torm <- torm[file.exists(torm)]
    torm <- torm[!grepl("random", torm)]
    file.remove(torm)
  }
  
  ##########################################################################################
  ############                        COMPRESS RESULTS                          ############
  ##########################################################################################
  
  if(tozip){
    print("############    COMPRESS RESULTS    ############")
    setwd(file.path(sampleD))
    zip.files <- c("coverage", "QC", "data",
                   file.path(paste0("OVL", nb.ovl, "_SIGNIF", nb.signif)), 
                   file.path(labels.treated, "hits")
    )
    utils::zip(zipfile = file.path(paste0(projectName, "_CASTseq.zip")), files = zip.files)
  }
  
  
  
}







