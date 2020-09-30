


bedTools.2in<-function(functionstring="bedIntersect",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}



bam2bed <- function(functionstring="bedtools bamtobed", bam, bed, opt.string="")
{
	options(scipen =99)
	
	command=paste(functionstring,"-i",bam,opt.string,">",bed,sep=" ")
  	cat(command,"\n")
  	try(system(command))
}


sortBed <- function(functionstring="sortBed", bed, bed.sort)
{
	options(scipen =99)
	
	command=paste(functionstring,"-i",bed,">",bed.sort,sep=" ")
  	cat(command,"\n")
  	try(system(command))
}

sortBedOLD <- function(functionstring="sortBed", bed, bed.sort)
{
	options(scipen =99)
	
	command=paste(functionstring,"-i",bed,sep=" ")
  	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

bedCoverage <- function(functionstring="bedtools genomecov", bed, bed.cov, opt.string="")
{
	options(scipen =99)
	
	command=paste(functionstring,"-i",bed,opt.string,">",bed.cov,sep=" ")
  	cat(command,"\n")
  	try(system(command))
}



getRandomBed <- function(functionstring="bedtools random", l, n, outFile, opt.string="")
{
	options(scipen =99)
	
	command=paste(functionstring, "-l", l, "-n", n, opt.string, ">", outFile, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}

getClosest <- function(functionstring="bedtools closest", bed1, bed2, outFile, opt.string="")
{
	options(scipen =99)

	command=paste(functionstring, "-a", bed1, "-b", bed2, opt.string, ">", outFile, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}


shuffleBed <- function(functionstring="bedtools shuffle", inFile, myGenome, outFile, opt.string="")
{
	options(scipen =99)
	
	command=paste(functionstring, "-i", inFile, "-g", myGenome, opt.string, ">", outFile, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}


nbReadFastq <- function(inFile)
{
	options(scipen =99)
	command=paste0("echo $(cat ", inFile, "|wc -l)/4|bc")
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

nbReadFastqgz <- function(inFile)
{
	options(scipen =99)
	command=paste0("zcat < ", inFile, " | echo $((`wc -l`/4))")
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

nbPlusReadBam <- function(inFile)
{
	options(scipen =99)
	command=paste0("samtools view -F 0x10 ", inFile, " | cut -f1 | sort | uniq | wc -l")
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}

nbMinusReadBam <- function(inFile)
{
	options(scipen =99)
	command=paste0("samtools view -f 0x10 ", inFile, " | cut -f1 | sort | uniq | wc -l")
	cat(command,"\n")
  	try(system(command, intern = TRUE))
}




if(FALSE)
{


mybam.folder <- "/Volumes/Home/Geoffroy/offTargets/Simone/bam"
mybam <- "HEK4_250bp.bam"

mybed.folder <- gsub("bam", "bed", mybam.folder)
mybed <- gsub("bam", "bed", mybam)

# Bam 2 Bed
bam2bed(bam = file.path(mybam.folder, mybam),
		bed = file.path(mybed.folder, mybed)
		)
		
# Sort Bed
mybed.sort <- gsub(".bed", ".sort.bed", mybed)
sortBed(bed = file.path(mybed.folder, mybed),
	bed.sort = file.path(mybed.folder, mybed.sort)
	)

# Bed coverage
mycov <- paste0(mybed.sort, ".cov")
option.genome = "-bg -g ~/Programs/bedtools2/genomes/human.hg19.genome"
bedCoverage(bed = file.path(mybed.folder, mybed.sort),
	bed.cov = file.path(mybed.folder, mycov),
	opt.string = option.genome
	)





# Toy example
bedCoverage(bed = "/Volumes/Home/Geoffroy/offTargets/Simone/bed/test.bed",
	bed.cov = "/Volumes/Home/Geoffroy/offTargets/Simone/bed/test.cov.bed",
	opt.string = "-g ~/Programs/bedtools2/genomes/human.hg19.genome -bg"
	)



}