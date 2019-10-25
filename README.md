# CAST-Seq

TO DO

## Getting Started

**General:** CAST-Seq pipeline

**Composing of Results:**  TO DO


### Prerequisites
Requiered software and databases

1. Software

	* R (3.4.2)
	* Bbmap (38.22) from https://jgi.doe.gov/data-and-tools/bbtools/
	* Bowtie2 (2.3.4.2) from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* samtools (1.9) from http://samtools.sourceforge.net
	* bedtools (2.27.1) from https://bedtools.readthedocs.io/en/latest/
	* fastQC from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* seqtk from https://github.com/lh3/seqtk

	
3. Databases from ANNOVAR
	* bowtie2Index
	* genome fasta


4. R package
	* openxlsx
	* GenomicRanges
	* Biostrings
	* data.table
	* ggplot2
	* ggseqlogo
	* textreadr
	* parallel
	* ChIPseeker
	* org.Hs.eg.db
	* clusterProfiler
	* rtracklayer
	* biomaRt
	* tools
	* BSgenome.Hsapiens.UCSC.hg38
	* TxDb.Hsapiens.UCSC.hg38.knownGene
	* karyoploteR

5. Additional Files (provided in annotation folder)
	* hg38_TSS_TES.txt (Bed file containing TSS and TES as start and end locations respectively)
	* hg38.chrom.sizes (chromosome size file)
	* TruSeq4-PE (adapter sequences)
	* CancerGenesList_ENTREZ.txt (OncoKB cancer related genes)
	* mismatch_table.xlsx

If available, the used versions are noted.


### Additional Files
#### Needed during epigenetic profil analysis
histones

## Directory structure

```
CAST-Seq
│
├── annotations
│   ├── bowtie2Index
│   │   ├── genome.1.bt2
│   │   ├── genome.2.bt2
│   │   ├── genome.3.bt2
│   │   ├── genome.4.bt2
│   │   ├── genome.fa.fai
│   │   ├── genome.rev.1.bt2
│   │   ├── genome.rev.1.bt2
│   │   └── genome.fa
│   ├── histones
│   │   ├── H3K4me3.bed
│   │   └── ...
│   ├── CancerGenesList_ENTREZ.txt
│   ├── hg38_TSS_TES.txt
│   ├── mismatch_table.xlsx
│   ├── hg38.chrom.sizes
│   └──TruSeq4-PE.fa
│
├── samples
│   ├── XXX
│   │   ├── data
│   │   │	├── fastq
│   │   │	│	├── XXX_treated_R1_001.fastq.gz
│   │   │	│	├── XXX_treated_R2_001.fastq.gz
│   │   │	│	├── XXX_UNtreated_R1_001.fastq.gz
│   │   │	│	└── XXX_UNtreated_R2_001.fastq.gz
│   │   │	├── gRNA.fa
│   │   │	├── headTOhead.fa
│   │   │	├── linker_RC.fa
│   │   │	├── linker.fa
│   │   │	├── mispriming.fa
│   │   │	├── neg.fa
│   │   │	├── ots.bed
│   │   │	└── pos.fa
│   │  	└── results
│   │   	├── fastq_aln
│   │   	├── guide_aln
│   │   	└── random
│   ├── YYY
│   │   ├── data
│   │   │	├── fastq
│   │   │	└── ...
│   │   └── ...
│   └── ...
│
└── script   
    ├── run
    │   ├── XXX_run.R
    │   ├── YYY_run.R
    │  	└── ...
    ├── lcs.py
    ├── annotateGenes.R
    ├── bed2sequence.R
    ├── bedTools_fct.R
    └── ...

```


## Running CAST-Seq
After all tools and databases are installed and work properly, you should prepare a XXX_run.sh file as described bellow.
Then CAST-Seq pipeline can be used using this single command:

```
Rscript XXX_run.R
```

### Adjusting **XXX_run.sh**
the XXX_run.R need to be adjusted for each sample as described below.

```
##################################################################################################################
#### Parameters which have to be adjusted accoridng the the environment or the users needs


# INPUT FILE AND DIRECTORY

# SET SAMPLE DIRECTORY NAME
sampleDname <- "G3_HiFi_D1"

# SET INPUT TEST FILE
sampleName <- "G3-Hifi-d1_S4_L001"

# SET INPUT CONTROL FILE
controlName <- "UT-G3-d1_S3_L001"

# SET REFERENCE FOLDER
homeD <- "/home/gandri/offTargets/Giando/pipelineGit/"


# OTHER PARAMETERS
scriptD <- file.path(homeD, "script")
annotD <- file.path(homeD, "annotations")

sampleD <- file.path(homeD, "samples", sampleDname)

dataD <- file.path(sampleD, "data")
resultD <- file.path(sampleD, "results", "guide_aln")
dir.create(resultD, showWarnings = FALSE)

# SET GUIDE SEQ
refSeq <- toupper(as.character(readDNAStringSet(file.path(dataD, "gRNA.fa"))))

# SET ON-TARGET SITE
otsBed <- file.path(dataD, "ots.bed")
otsDistance <- 50
surrounding_size <- 20000

# SET FLANKING REGIONS (+/- flanking size)
flankingSize <- 2500

# SET NUMBER OF RANDOM SEQUENCES
nb.rd <- 10000

# SET WIDTH
w = 250

# SET DISTANCE CUTOFF (FOR CLUSTERS)
distance.cutoff <- 1500

# SET PVALUE CUTOFF
pv.cutoff <- 0.05

# SET GENOME VERSION
myGenome.size <- file.path(annotD, "hg38.chrom.sizes")

# HG38 TSS TES
geneMat <- read.delim(file.path(annotD, "hg38_TSS_TES.txt"), header = FALSE)

# SET ONCO ENTREZ LIST
oncoEntrez <- file.path(annotD, "CancerGenesList_ENTREZ.txt")
onco.width <- 3000

# SET HISTONE BROAD PEAKS FILES
histoneFiles <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = TRUE)
names(histoneFiles) <- list.files(file.path(annotD, "histones"), pattern = ".broadPeak_hg38_homer.bed$", full.names = FALSE)
names(histoneFiles) <- gsub(".broadPeak_hg38_homer.bed", "", names(histoneFiles))
```


## Example
For testing MIRACUM-Pipe the public available dataset from Texas (http://txcrb.org/data.html and https://www.nature.com/articles/sdata201610) can be used. We provide the PDF reports for the samples TCRBOA6 and TCRBOA7. Both samples were run with the default parameters.

## Authors

* Geoffroy Andrieux


## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

## Acknowledgments
We thank

* The Molecular Tumor Board Freburg Team
* The whole MIRACUM consortia
* The German Ministry of Education and Research (BMBF) for funding
* The developers from Control-FREEC for the code on CNV significance


## Cite

If you use this work please cite Turchiano. et al., (submitted)

## References


