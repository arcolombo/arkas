#' encapsulate limma/voom analysis and TMM normalization at transcript level 
#' FIXME: just farm it out to sleuth, e.g. as seen in ?artemisData::withSleuth
#' 
#' @param kexp        A KallistoExperiment
#' @param design      A model matrix 
#' @param read.cutoff Exclude transcripts where the maximum count is < this 
#' @param ...          any additional parameters
#' @return        A list with elements (design, voomed, fit)
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma ebayes
#' @importFrom limma lmFit
#' @export
fitTranscripts <- function(kexp, design, read.cutoff=1, ...) { 

  res <- list()
  filteredCounts <- collapseTranscripts(kexp, read.cutoff=read.cutoff, ...)
  dge <- DGEList(counts=filteredCounts)
  dge <- calcNormFactors(dge)
  res$design <- design 
  res$voomed <- voom(dge, res$design)
  res$fit <- eBayes(lmFit(res$voomed, res$design))
  return(res)

}
