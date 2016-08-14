#' Analysis of raw transcript abundance estimates.
#' 
#' @param kexp        a KallistoExperiment or SummarizedExperiment-like object 
#' @param design      a design matrix w/contrast or coefficient to test in col2
#' @param p.cutoff    where to set the p-value cutoff for plots, etc. (0.05)
#' @param fold.cutoff where to set the log2-FC cutoff for plots, etc. (1 == 2x)
#' @param coef        which column of the design matrix to test on (2nd)
#' @param read.cutoff a cutoff to filter reads below across samples
#' @param tx_biotype  optionally restrict to one or more tx_biotype classes 
#' @param gene_biotype optionally restrict to one or more gene_biotype classes 
#' @param biotype_class optionally restrict to one or more biotype_class ...es
#' @param adjustBy either none , BH, BY, holm for limma adjust type
#' @import edgeR 
#' @import limma
#'
#' @export
transcriptWiseAnalysis <- function(kexp, design, p.cutoff=0.05, fold.cutoff=1, 
                                   coef=2,read.cutoff=1,tx_biotype=NULL, gene_biotype=NULL,
                                   biotype_class=NULL,adjustBy="holm",species=c("Homo.sapiens","Mus.musculus")){ 

  ## this is really only meant for a KallistoExperiment
  if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }

 if (is.null(design)) {
    if (!is.null(exptData(kexp)$design)) {
      design <- exptData(kexp)$design
    } else {
      stop("A design matrix must be supplied, or present in metadata.")
    }
  }

 species <- match.arg(species, c("Homo.sapiens", "Mus.musculus"))
  commonName <- switch(species,
                       Mus.musculus="mouse",
                       Homo.sapiens="human")
  adjustBy<-match.arg(adjustBy, c("none","BH","BY","holm"))
  choices<- c("holm", "BY", "BH","none")
 ranked<-data.frame(type=choices,rank=c(1,2,3,4))
  message("Fitting bundles...")
  initialRank<-ranked[which(ranked$type==adjustBy),2]

  if (all(sapply(c(tx_biotype, gene_biotype, biotype_class), is.null))) {
    res <- fitTranscripts(kexp, design, read.cutoff)
     while( initialRank <=4 ) {
     message(paste0("fitting using FDR: ",adjustBy))
     res$top <- with(res, topTable(fit, coef=2, p=p.cutoff,adjust.method=adjustBy, n=nrow(kexp)))

      if(nrow(res$top)==0){
     message(paste0("no DE found for using FDR: ",adjustBy))
     initialRank<-initialRank + 1
     adjustBy<-as.character(ranked$type[initialRank])
       }
      else{
      message(paste0("found ", nrow(res$top), " DE genes using FDR procedure ", as.character(ranked$type[initialRank]) ))
     initialRank<-5 #break the loop at first instance
       }
    }

    res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
    topTranscripts <- rownames(res$top)
    res$topTranscripts <- topTranscripts

  } else {
    keep <- seq_len(nrow(kexp))
    if (!is.null(biotype_class)) {
      keep <- intersect(keep, which(mcols(kexp)$biotype_class == biotype_class))
    }
    if (!is.null(gene_biotype)) {
      keep <- intersect(keep, which(mcols(kexp)$gene_biotype == gene_biotype))
    }
    if (!is.null(tx_biotype)) {
      keep <- intersect(keep, which(mcols(kexp)$tx_biotype == tx_biotype))
    }
    res <- fitTranscripts(kexp[keep, ], design, read.cutoff)
    res$top <- with(res, topTable(fit, coef=2, p=p.cutoff, n=nrow(kexp[keep,])))
    res$top <- res$top[ abs(res$top$logFC) >= fold.cutoff, ] ## per SEQC
    topTranscripts <- rownames(res$top)
    res$topTranscripts <- topTranscripts
    }
 
  res$biotype_class <- biotype_class
  res$gene_biotype <- gene_biotype
  res$tx_biotype <- tx_biotype

  res$limmaWithMeta<-cbind(res$top,rowRanges(kexp)[rownames(res$top)]$gene_name,rowRanges(kexp)[rownames(res$top)]$gene_id )
  colnames(res$limmaWithMeta)[ncol(res$top)+1]<-"gene.name"
  colnames(res$limmaWithMeta)[ncol(res$top)+2]<-"gene.id"

 return(res)


}
