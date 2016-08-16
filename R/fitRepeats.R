#' repeat analysis is used for analysis of repeat regions from an annotated kexp
#' @param kexp  kallisto experiment object
#' @param design design matrix for modeling linear fit
#' @param ...    any additional params
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom limma voom
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @return a limma list of linear model statistics
#' @export 
fitRepeats<-function(kexp,design,...){

    rps<-list()
    idx<-which(rowRanges(kexp)$biotype_class=="repeat")
    repFeatures<-rowRanges(kexp)[idx]
    repeatNames<-as.character(seqnames(repFeatures))
    tt<-rownames(counts(kexp)) %in% repeatNames
    repKexp<-counts(kexp)[tt,]
    rge<-DGEList(counts=repKexp)
    rge<-calcNormFactors(rge)
    rps$design<-design
    rps$voomed<-voom(rge,rps$design)
    rps$fit<-eBayes(lmFit(rps$voomed,rps$design))
    return(rps)

}#{{{Main
