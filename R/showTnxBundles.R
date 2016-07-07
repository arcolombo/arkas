#' Groups bundles of transcripts, discard any that represent pointless tests,
#' and shows the transcript bundles for each gene. 
#' @param kexp          A KallistoExperiment (or something very much like it)
#' @param design the design matrix to be used to fit a linear model for Diff Expr
#' @param species either human or mouse
#' @import GenomicRanges
#' @return              a matrix of summarized counts per sample bundle 
#'
#' @seealso collapseTranscripts
#'
#' @export 
showTnxBundles<-function(kexp,design,species=c("Homo.sapiens","Mus.musculus")){
#check input
 if (!is(kexp, "KallistoExperiment")) {
    message("This function is optimized for KallistoExperiment objects.")
    message("It may work for other classes, but we make no guarantees.")
  }
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))

#call twa  use holm filter
  twa<-transcriptWiseAnalysis(kexp,design, species=species,adjustBy="holm")
  bundled<-split.data.frame(twa$limmaWithMeta,as.character(twa$limmaWithMeta$gene.name))

return(bundled)
} #{{{ main
