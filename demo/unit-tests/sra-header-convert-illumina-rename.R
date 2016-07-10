---
title: "Arkas: Repetitive Elements Quantification In Much Less Time"
author:  "Timothy J. Triche, Jr, Anthony R. Colombo, Harold Pimentel"
output: 
    html_document:
      toc: true 
      number_sections: true
date: "`r format(Sys.time(), '%d %B, %Y')`"
---
          
#Introduction
Importing from SRA or ENA has fixed headers and changing these headers to be compatible with illumina fastq uploading header conventions is a trivial but useful task.

```

library(arkas)
library(arkasData)
library(nat.utils)
od<-"~/Documents/"
fastqPath <-system.file("extdata","fastq/SRA",package="arkasData")
fastqFile<-dir(fastqPath)
system(paste0("sudo mkdir ",od,"/illuminaTest"))
system(paste0("sudo chmod -R 777 ",od,"/illuminaTest"))
sraOutputDir<-paste0(getwd(),"/illuminaTest")
sraFastqHeaderToIlluminaStandard(headerFormat="SRA",fastqPath=fastqPath,fastqFile=fastqFile,sraOutputDir=sraOutputDir,fastqReadNumber=1)



fastqFileUploadToBaseSpace(illuminaDirPath="~/Documents/",illuminafastqFile=dir(sraOutputDir)
,basespaceProject="SIRT1_Chip_Seq",illuminaDirs="illuminaTest",dryRun=TRUE,paired=FALSE)

system(paste0("sudo rm -r ",od,"/illuminaTest"))

```
