# CAST-Seq

CAST(chromosomal aberrations analysis by single targeted LM-PCR)-Seq is a novel method capable of detecting and quantifying chromosomal aberrations derived from on- and off-target activity of CRISPR-Cas nucleases or TALEN. See [Turchiano et al.](https://www.sciencedirect.com/science/article/abs/pii/S1934590921000527?via%3Dihub) for detail information about CAST-Seq background and potential clinical application.

## Cite

If you use this work please cite [Turchiano et al., Cell Stem Cell, 2021](https://www.sciencedirect.com/science/article/abs/pii/S1934590921000527?via%3Dihub)

## Getting Started

**General:** The herein code is the official CAST-Seq bioinformatic pipeline to process fastq files generated by CAST-Seq. 

**Composing of Results:**
The results directory is divided into 3 sub-directories: fastq_aln, guide_aln and random.

1. *fastq_aln* contains all pre-processing and alignment files from fastq.gz to bam files.

2. *guide_aln* contains the post-processing files from bed files to final xlsx report.

3. *random* contains the information related to the random regions that are used for normalisation.

### Prerequisites
Requiered software and databases

1. Software

	* R (3.4.2)
	* BBmap and BBmerge (38.22) from https://jgi.doe.gov/data-and-tools/bbtools/
	* Bowtie2 (2.3.4.2) from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	* samtools (1.9) from http://samtools.sourceforge.net
	* bedtools (2.27.1) from https://bedtools.readthedocs.io/en/latest/
	* fastQC from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* seqtk from https://github.com/lh3/seqtk

	
3. Databases from ANNOVAR
	* bowtie2Index
	* genome fasta


4. R packages
	* openxlsx
	* GenomicRanges
	* Biostrings
	* data.table
	* ggplot2
	* ggseqlogo
	* textreadr
	* parallel
	* ChIPseeker
	* clusterProfiler
	* rtracklayer
	* biomaRt
	* tools
	* karyoploteR
	
	* org.Hs.eg.db
	* BSgenome.Hsapiens.UCSC.hg38
	* TxDb.Hsapiens.UCSC.hg38.knownGene

5. Additional Files (provided in annotation folder)
	* hg38_TSS_TES.txt (Bed file containing TSS and TES as start and end locations respectively)
	* chrom.sizes (chromosome size file)
	* TruSeq4-PE (adapter sequences)
	* CancerGenesList_ENTREZ.txt (OncoKB cancer related genes)

If available, the used versions are noted.


### Additional Files
#### Needed during epigenetic profil analysis
Histones marks bed files (must be stored into annotations/histones/).
Example of such file is provided for H3K4me3 in Primary hematopoietic stem cells (E035 from Roadmap Epigenomics https://egg2.wustl.edu/roadmap/web_portal/processed_data.html).

## Directory structure

```
CAST-Seq
│
├── annotations
│   └── human
│       ├── bowtie2Index
│       │   ├── genome.1.bt2
│       │   ├── genome.2.bt2
│       │   ├── genome.3.bt2
│       │   ├── genome.4.bt2
│       │   ├── genome.fa.fai
│       │   ├── genome.rev.1.bt2
│       │   ├── genome.rev.1.bt2
│       │   └── genome.fa
│       ├── histones
│       │   ├── H3K4me3.bed
│       │   └── ...
│       ├── CancerGenesList_ENTREZ.txt
│       ├── hg38_TSS_TES.txt
│       ├── chrom.sizes
│       └──TruSeq4-PE.fa
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
    │   ├── XXX.R
    │   ├── YYY.R
    │  	└── ...
    ├── lcs.py
    ├── annotateGenes.R
    ├── bed2sequence.R
    ├── bedTools_fct.R
    └── ...

```


## Running CAST-Seq
After all tools and databases are installed and work properly, the whole CAST-Seq pipeline can be executed using this single command (see example in
script/run/):

```
Rscript ./CAST-Seq.R --pipeline "crispr"\
	             --sampleDname "G3_TOY"\
		     --sampleName "G3_treated"\
		     --controlName "G3_UNtreated"\
		     --homeD "../../"
```

**--pipeline** name of the pipeline you want to use. Choose between "crispr" and "talen"<br/>
**--sampleDname** name of sample directory<br/>
**--sampleName** XXX name of test (treated) file. XXX_R1_001.fastq.gz AND XXX_R2_001.fastq.gz should exist<br/>
**--controlName** XXX name of control (untreated) file. XXX_R1_001.fastq.gz AND XXX_R2_001.fastq.gz should exist<br/>
**--homeD** name of home directory<br/>

### Parameters
Additional parameters can be changed in the command above. Here is a description of these parameters:<br/>

**--fastqD** name of directory containing the fastq files<br/>
**--grna** name of gRNA fasta (default "gRNA.fa")<br/>
**--onTarget** name of ON-target bed file (default "ots.bed")<br/>
**--otsDistance** distance (bp) from the ON-target. Reads +/- this distance will be removed (default 50)<br/>
**--surrounding_size** distance (bp) from the ON-target. Use for the scoring system (default 20000)<br/>
**--flank1** name of first flanking sequence (default "flank1.fa")<br/>
**--flank12** name of second flanking sequence (default "flank12.fa")<br/>
**--flankingSize** distance to consider for HMT (default 2500)<br/>
**--random** number of random sequences to generate (default 10000)<br/>
**--width** distance to extend the putative sites (default 250)<br/>
**--distCutoff** distance to merge hits together (default 1500)<br/>
**--pvCutoff** pvalue threshold (default 0.05)<br/>
**--scoreCutoff** gRNA alignment score threshold (default *NULL*)<br/>
**--hitsCutoff** minimum number of hits per site (default 1)<br/>
**--saveReads** should reads fastq sequences be saved (default "no")<br/>
**--species** name of sample species (default "hg") *so far only hg can be used*<br/>
**--cpu** number of CPUs (default 2) *at least 4 is advised*<br/>
**--pythonPath** python path (default "/usr/bin/python")<br/>

#### TALEN specific parameters
These parameters are only used when **--pipeline** "talen" is set.<br/>
**--grnaR** name of gRNA (RIGHT) fasta file<br/>
**--grnaL** name of gRNA (LEFT) fasta file<br/>

## Authors

* Geoffroy Andrieux
* Giandomenico Turchiano


## License
This software is under AGPL3 license.

## Acknowledgments
We thank all members of our laboratories for constructive discussions and suggestions.


## References

* Langmead, B, Salzberg, SL (2012). Fast gapped-read alignment with Bowtie 2. Nat. Methods, 9, 4:357-9.
* Quinlan, AR, Hall, IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26, 6:841-2.
* Lawrence, M, Huber, W, Pagès, H, Aboyoun, P, Carlson, M, Gentleman, R, Morgan, MT, Carey, VJ (2013). Software for computing and annotating genomic ranges. PLoS Comput. Biol., 9, 8:e1003118.
* Yu, G, Wang, LG, He, QY (2015). ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. Bioinformatics, 31, 14:2382-3.
* Durinck, S, Spellman, PT, Birney, E, Huber, W (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Nat Protoc, 4, 8:1184-91.




