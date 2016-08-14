#' RUVSeq methods allow for proper normalization by using a GLM to find the linear space of unwanted variance; it uses factor analysis to find the span of unwanted variance using ERCC expression negative controls, or in silico negative controls.  the advantage is that it does not assume a constant global normalization correction for constant assumed technical noise.  RUVg instead determines the linear space where the negative controls span. 
#' 
#' @param kexp kallisto Experiment object or something along its line
#' @param k integer, this is the k value for number of unwanted variance
#' @param spikeIns  boolean, whether ERCC spike-ins are to be used (FALSE) 
#' @param p.cutoff  numeric  a p value cutoff it should be fixed at 1 or 2 becuaes you need to find insignificant negative controls
#' @param inSilico for when spikeIns is flagged as FALSE, inSilcio must be a vector names of in silico genes which are constant across samples apriori. housekeeping genes will do fine.  the insilico vector can be derived here if it is unknown by taking the bottom quartile, bottom 10 percent ranked by P.Value, of significant genes after running a raw DE analysis.
#' @param read.cutoff , integer here we employ a read cutoff that filters out any rows where the rowSums falls under this category.  
#' @param byLevel a string character which must match the names of the meta-columns of the rowRanges(kexp), this collapses the count data by this feature term, and performs filtering
#' @param controlNumber integer, this is the amount of negative controls to use to normalize; the higher the number the more stable your normalization but higher chance of including a false negative as a negative control; the lower the number the lower the chance of false negative, but not as stable normalization.
#' @importFrom RUVSeq RUVg
#' @export 
#' @return return a list object of the selected negative controls, and RUVg normalization design matrix weights

ruvNormalization <- function(kexp, k=1, spikeIns=FALSE, p.cutoff=1, 
                             inSilico=NULL, read.cutoff=1, 
                             byLevel=c("gene_id", "tx_id"),
                             controlNumber=100) {

  if(!is(kexp, "KallistoExperiment")) {
    warning("This method only works with KallistoExperiment-like objects.")
  }

  byLevel <- match.arg(byLevel)
  collapsedKexp <- collapseBundles(kexp, 
                                   bundleID=byLevel, 
                                   read.cutoff=read.cutoff) 
  exprs <- round(collapsedKexp) # must be integers? wtf

  if(spikeIns == "TRUE") {

    #{{{
    spikes <- rownames(exprs)[grep("^ERCC",rownames(exprs))] #grabbing ERCCs
    if(length(spikes) == 0) {
      stop("ERCC spike-ins were not found... please try again..")
    }
    message("detected ERCC spike ins ...")
    # }}}

    ruvOutput <- RUVg(exprs,spikes,k=k)

  } else { 

    #inSilico must be a character vector
    if(is.null(inSilico)) { 
      message("Did not detect a vector of in silico negative controls...")
      message("Checking for a design matrix in metadata(kexp)...")
      if( is.null(metadata(kexp)$design)) { # {{{ stop
        stop("please include a vector of row names, or add a design matrix.")
      } # }}}
      message("Found design matrix, determining in silico negative controls...")
    
      if(byLevel=="gene_id") {
        # {{{ collapse by gene_id
        message("performing gene-wise-analysis with higher cutoff to search for negative controls...")
        GWA <- geneWiseAnalysis(kexp, design=metadata(kexp)$design, how="cpm",
                                p.cutoff=p.cutoff,
                                fold.cutoff=1, 
                                read.cutoff=read.cutoff,
                                fitOnly=TRUE,adjustBy="BH")
        
        idx <- rev(order(GWA$top$adj.P.Val))
        master<-GWA$top[idx,]
        rIdx<-which(GWA$top$adj.P.Val[idx]>=0.60)
             
        print(paste0("found ",length(rIdx), " with adj.P.Val greater than 0.60 used as negative in silico controls "))
        derived.inSilico <- rownames(master[rIdx,])
        derived.inSilico<-derived.inSilico[1:controlNumber]
        ruvOutput <- RUVg(exprs,derived.inSilico,k=k)
      }
    
      if(byLevel=="tx_id"){
        # {{{ collapse by tx_id
        message("performing transcript-wise-analysis...")
        #need to collapseTranscripts             
        TWA<-transcriptWiseAnalysis(kexp, design=metadata(kexp)$design,
                                    p.cutoff=p.cutoff, fold.cutoff=1,
                                    read.cutoff=read.cutoff,adjustBy="BH")  
        
        idx <- rev(order(TWA$top$adj.P.Value))
         master<-TWA$top[idx,]
        rIdx<-which(TWA$top$adj.P.Val[idx]>0.60)
        print(paste0("found ",length(rIdx), " with adj.P.Val greater than 0.60 used as negative in silico controls "))
        derived.inSilico <- rownames(master[rIdx,])
        derived.inSilico<-derived.inSilico[1:controlNumber]
        trnxExprs <- collapseTranscripts(kexp,read.cutoff=read.cutoff)
        trnxExprs <- round(trnxExprs)
        # }}}
        ruvOutput <- RUVg(trnxExprs,derived.inSilico,k=k)
      }
    } #design matrix found
  
  }
  if (!is.null(inSilico)) { 
    idx <- (rownames(assays(kexp)$est_counts) %in% inSilico)
    silico<-rownames(assays(kexp)$est_counts)[idx] #grabbing ERCCs
    if (length(silico)>=1) {
      ruvOutput<-RUVg(exprs,silico,k=k)
    }
  }

  return(list(negative.insilico.controls=derived.inSilico,RUV=ruvOutput))
}
