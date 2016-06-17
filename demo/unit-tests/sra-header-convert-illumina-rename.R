library(arkas)
library(arkasData)
library(nat.utils)
od<-getwd()
fastqPath <-system.file("extdata","fastq/SRA",package="arkasData")
fastqFile<-dir(fastqPath)
system("sudo mkdir illuminaTest && chmod -R 777 ./illuminaTest")
sraOutputDir<-paste0(getwd(),"/illuminaTest")
sraFastqHeaderToIlluminaStandard(headerFormat="SRA",fastqPath=fastqPath,fastqFile=fastqFile,sraOutputDir=sraOutputDir,fastqReadNumber=1)



fastqFileUploadToBaseSpace(illuminaDirPath="~/Documents/",illuminafastqFile=dir(sraOutputDir)
,basespaceProject="SIRT1_Chip_Seq",illuminaDirs="illuminaTest",dryRun=TRUE,paired=FALSE)
