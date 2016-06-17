dataIsInstalled <- suppressPackageStartupMessages(require(arkasData))
if(!dataIsInstalled){ 
  message("You need to install the arkasData package.")
  message("Try biocLite(\"ramsinghlab/arkasData\") then come back.")
  stop("Go ahead, we'll wait!")
} 
library(matrixStats)

message("Looking for Kallisto in ~/bin...")
kallisto <- paste0(path.expand("~/bin"), "/kallisto")
if(!file.exists(kallisto)) {
  message("You do not seem to have kallisto installed.  We can't proceed.")
} else {
  message("Found it, proceeding...")
}

pathBase <- system.file("extdata", "", package="arkasData")
fastaPath <- paste0(pathBase, "/fasta")
fastqPath <- paste0(pathBase, "/fastq")
samples <- c(MrN="MrN", MrT="MrT") ## normally set by appSession
fastaFiles <- c( "ERCC.fa.gz", ## spike-in controls  
                 "Homo_sapiens.RepBase.20_07.merged.fa.gz")

## build an index if it isn't already there (in artemisData, it is)
indexName <- indexKallisto(fastaFiles=fastaFiles, fastaPath=fastaPath,makeUnique=TRUE)$indexName

## run pseudoalignments 
library(parallel)
results <- mclapply(samples, 
                    runKallisto,
                    indexName=indexName,
                    fastqPath=fastqPath,
                    fastaPath=fastaPath,
                    bootstraps=100,
                    outputPath=".")
merged <- mergeKallisto(samples, outputPath=".")

## plot repeat element txn using (counts/bootstrap MADs) as "effect size"
topKbyMAD <- function(kexp, k=25) {
  tpm(kexp)[rev(order(rowMeans(counts(kexp) / mad(kexp))))[1:k],]
}
pdf("ExamplePlot.pdf")
heatmap(log1p(topKbyMAD(merged)), scale="none", 
        col=colorRampPalette(c("white","red","darkred"))(255),
        main="Repeat transcription, teratoma vs. normal")
dev.off()
