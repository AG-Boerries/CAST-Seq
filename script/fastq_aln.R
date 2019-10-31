

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

fastqAln <- function(functionstring="fastq_aln.sh", homeFolder, sampleFolder, testSample, utSample)
{
	options(scipen =99)
	
	command=paste(functionstring, homeFolder, sampleFolder, testSample, utSample, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}


fastqAln_2OT <- function(functionstring="fastq_aln_2OT.sh", homeFolder, sampleFolder, testSample, utSample)
{
	options(scipen =99)
	
	command=paste(functionstring, homeFolder, sampleFolder, testSample, utSample, sep=" ")
  	cat(command,"\n")
  	try(system(command))
}
