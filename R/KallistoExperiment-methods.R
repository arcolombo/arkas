#' @rdname tpm-methods
#' @param object a feature defined object
#' @aliases counts,tpm
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

#' covariates generic
#' @name covariates
#' @rdname covariates-methods
#' @exportMethod covariates
setGeneric("covariates", function(object) standardGeneric("covariates"))


#' covariates setter/replacement generic
#' @name covariates<-
#' @rdname covariates-methods
#' @exportMethod covariates<-
setGeneric("covariates<-", 
           function(object, value) standardGeneric("covariates<-"))

#' @rdname covariates-methods
#' @aliases covariates
setMethod("covariates", "KallistoExperiment",
          function (object) return(colData(object)))

#' @rdname covariates-methods
setMethod("pData", "KallistoExperiment",
          function (object) return(colData(object)))

#' @rdname covariates-methods
#' @param object a feature defined object
#' @param value a value to replace features
#' @aliases covariates
setReplaceMethod("covariates", "KallistoExperiment",
                 function (object, value) {
                   object <- BiocGenerics:::replaceSlots(object, colData=value)
                   msg <- SummarizedExperiment:::.valid.SummarizedExperiment0.assays_ncol(object)
                   if (!is.null(msg)) stop(msg)
                   else return(object)
                 })


#' @rdname covariates-methods
setReplaceMethod("pData", c("KallistoExperiment", "DataFrame"),
                 function (object, value) {
                   object <-BiocGenerics:::replaceSlots(object, colData=value)
                   msg <- SummarizedExperiment:::.valid.SummarizedExperiment0.assays_ncol(object)
                   if (!is.null(msg)) stop(msg)
                   else return(object)
                 })


#' features generic
#' @name features
#' @rdname features-methods
#' @exportMethod features
setGeneric("features", function(object) standardGeneric("features"))


#' method for features<-
#' @name features<-
#' @rdname features-methods
#' @exportMethod features<-
setGeneric("features<-", function(object, value) standardGeneric("features<-"))

#' @rdname features-methods
#' @aliases features
setMethod("features", "KallistoExperiment", function (object) rowRanges(object))

#' @rdname features-methods
#' @param object feature defined object
#' @param value a value defined object
#' @aliases features
setReplaceMethod("features", c("KallistoExperiment", "ANY"),
                function(object, value) {
                  object <- BiocGenerics:::replaceSlots(object,
                            rowRanges=value)
                  msg <- SummarizedExperiment:::.valid.SummarizedExperiment0.assays_nrow(object)
                  if (!is.null(msg)) stop(msg)
                  else return(object)
})

#' eff length generic
#' @name eff_length
#' @rdname eff_length-methods
#' @exportMethod eff_length
setGeneric("eff_length", function(object) standardGeneric("eff_length"))

#' @rdname eff_length-methods
#' @param object a feature defined object
#' @aliases eff_length
setMethod("eff_length", "KallistoExperiment",
          function (object) return(assays(object)$eff_length))

#' tpm generic
#' @name tpm
#' @rdname tpm-methods
#' @exportMethod tpm
setGeneric("tpm", function(object) standardGeneric("tpm"))

#' @rdname tpm-methods
#' @aliases tpm,counts
setMethod("tpm", "KallistoExperiment",
          function (object) {

          if(length(assays(object)$tpm)==0) {
            rate <- log(counts(object)) - log(eff_length(object))
            exp(rate - log(sum(exp(rate))) + log(1e6))
            }#if
        else{
        return(assays(object)$tpm)
        } #else

               })   
       
#' kallisto version accessor
#' @name kallistoVersion
#' @rdname kallistoVersion-methods
#' @exportMethod kallistoVersion
setGeneric("kallistoVersion", 
           function(object) standardGeneric("kallistoVersion"))

#' @rdname kallistoVersion-methods
#' @param object a kallisto experiment
#' @aliases kallistoVersion
setMethod("kallistoVersion", "KallistoExperiment",
          function (object) return(object@kallistoVersion))

#' transcriptome accessor 
#' @name transcriptomes
#' @rdname transcriptomes-methods
#' @exportMethod transcriptomes
setGeneric("transcriptomes", 
           function(object) standardGeneric("transcriptomes"))

#' @rdname transcriptomes-methods
#' @param object feature defined object
#' @aliases transcriptomes
setMethod("transcriptomes", "KallistoExperiment",
          function (object) return(object@transcriptomes))

setMethod("transcriptsBy", "KallistoExperiment",
          function(x, by="gene", ...) {
            if (by == "gene") { 
              split(x, mcols(x)$gene_name)
            } else { 
              return(x[mcols(x)$gene_name == by, ])
            }
          })

setMethod("mad", "KallistoExperiment", function(x) assays(x)$est_counts_mad)

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 



#' Convert a KallistoExperiment to a SummarizedExperiment without losing data
#' @name as
#' 
setAs("KallistoExperiment", "SummarizedExperiment", 
      function(from) {
        metadata(from)$transcriptomes <- transcriptomes(from)
        metadata(from)$kallistoVersion <- kallistoVersion(from)
        SummarizedExperiment(assays(from), rowRanges=rowRanges(from), 
                             colData=colData(from), metadata=metadata(from))
      })



#' Convert suitably annotated SummarizedExperiment back to a KallistoExperiment
#' @name as
#' 
setAs("SummarizedExperiment", "KallistoExperiment", 
      function(from) {
        txomes <- metadata(from)$transcriptomes
        kversion <- metadata(from)$kallistoVersion
        new("KallistoExperiment", from, 
            kallistoVersion=kversion, transcriptomes=txomes)
      })
