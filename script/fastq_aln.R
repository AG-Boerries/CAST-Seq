

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

fastqAln <- function(functionstring="fastq_aln.sh", homeFolder, annotFolder, sampleFolder, fastqFolder, testSample, utSample, cpu)
{
	options(scipen =99)
	
	command=paste(functionstring, homeFolder, annotFolder, sampleFolder, fastqFolder, testSample, utSample, cpu, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}


fastqAln_2OT <- function(functionstring="fastq_aln_2OT.sh", homeFolder, annotFolder, sampleFolder, fastqFolder, testSample, utSample, cpu)
{
	options(scipen =99)
	
	command=paste(functionstring, homeFolder, annotFolder, sampleFolder, testSample, utSample, cpu, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}
