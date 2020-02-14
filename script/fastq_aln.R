

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

fastqAln <- function(functionstring="fastq_aln.sh", homeFolder, annotFolder, sampleFolder, testSample, utSample)
{
	options(scipen =99)
	
	command=paste(functionstring, homeFolder, annotFolder, sampleFolder, testSample, utSample, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}


fastqAln_2OT <- function(functionstring="fastq_aln_2OT.sh", homeFolder, annotFolder, sampleFolder, testSample, utSample)
{
	options(scipen =99)
	
	command=paste(functionstring, homeFolder, annotFolder, sampleFolder, testSample, utSample, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}
