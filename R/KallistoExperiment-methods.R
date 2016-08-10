#' @describeIn KallistoExperiment 
#' Retrieve the estimated count matrix from a KallistoExperiment. 
#' @export
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

#' new generics for arkas covariates definition
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod covariates
setGeneric("covariates", function(object) standardGeneric("covariates"))


#' generics arkas defined covariates setter
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod covariates<-
setGeneric("covariates<-", 
           function(object, value) standardGeneric("covariates<-"))
#' @rdname KallistoExperiment-class
#' @aliases covariates,KallistoExperiment-method
#' @param object for which features are defined
#' @export
setMethod("covariates", "KallistoExperiment",
          function (object) return(colData(object)))
#' @describeIn KallistoExperiment 
#' Retrieve the pData matrix from a KallistoExperiment. 
#' @param object: something from which to retrieve covariates to Retrieve the estimated count matrix from a KallistoExperiment. 
#'
#' @export
setMethod("pData", "KallistoExperiment",
          function (object) return(colData(object)))

#' @describeIn KallistoExperiment 
#'
#' @param   object: something to which features should be assigned
#' @param   value: the features to assign (usually a GRanges or GRangesList)
#' @export
setReplaceMethod("covariates", "KallistoExperiment",
                 function (object, value) {
                   object <- BiocGenerics:::replaceSlots(object, colData=value)
                   msg <- SummarizedExperiment:::.valid.SummarizedExperiment0.assays_ncol(object)
                   if (!is.null(msg)) stop(msg)
                   else return(object)
                 })


#' @describeIn KallistoExperiment 
#' @param   object: something to which features should be assigned
#' @param   value: the features to assign (usually a GRanges or GRangesList)
#'
#' @export
setReplaceMethod("pData", c("KallistoExperiment", "DataFrame"),
                 function (object, value) {
                   object <-BiocGenerics:::replaceSlots(object, colData=value)
                   msg <- SummarizedExperiment:::.valid.SummarizedExperiment0.assays_ncol(object)
                   if (!is.null(msg)) stop(msg)
                   else return(object)
                 })

#' set in GenomicFeatures, which we have to import anyways setGeneric("features", function(object) standardGeneric("features"))
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod features<-
setGeneric("features<-", function(object, value) standardGeneric("features<-"))

#' @rdname KallistoExperiment-class
#' @aliases features,KallistoExperiment-method
#' Retrieve the per-row annotations for a KallistoExperiment. 
#' @export
setMethod("features", "KallistoExperiment", function (object) rowRanges(object))

#' @describeIn KallistoExperiment
#'
#' Assign per-row annotations to a KallistoExperiment.
#' @param   object: something to which features should be assigned
#' @param   value: the features to assign (usually a GRanges or GRangesList)
#' @export
setReplaceMethod("features", c("KallistoExperiment", "ANY"),
                function(object, value) {
                  object <- BiocGenerics:::replaceSlots(object,
                            rowRanges=value)
                  msg <- SummarizedExperiment:::.valid.SummarizedExperiment0.assays_nrow(object)
                  if (!is.null(msg)) stop(msg)
                  else return(object)
})

#' eff_length generic 
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod eff_length
setGeneric("eff_length", function(object) standardGeneric("eff_length"))

#' @rdname KallistoExperiment-class
#' @aliases eff_length,KallistoExperiment-method
#' Retrieve the matrix of effective transcript lengths from a KallistoExperiment
#' @export
setMethod("eff_length", "KallistoExperiment",
          function (object) return(assays(object)$eff_length))

#' tpm generic 
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod tpm
setGeneric("tpm", function(object) standardGeneric("tpm"))

#' @rdname KallistoExperiment-class
#' @aliases tpm,KallistoExperiment-method
#' Obtain tpm from the precomputed matrix ( computed as shown in 
#' https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ ), specifically,
#' \code{ rate <- log(counts(object)) - log(eff_length(object)); }
#' 
#' \code{ tpm <- exp(rate - log(sum(exp(rate))) + log(1e6)) } 
#' 
#' @export
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
       
#' kallistoVersion generic 
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod kallistoVersion
setGeneric("kallistoVersion", 
           function(object) standardGeneric("kallistoVersion"))

#' @rdname KallistoExperiment-class
#' @aliases kallistoVersion, KallistoExperiment-method
#' Retrieve the version of Kallisto used for alignment from a KallistoExperiment
#' @export
setMethod("kallistoVersion", "KallistoExperiment",
          function (object) return(object@kallistoVersion))

#' transcriptomes generic 
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @exportMethod transcriptomes
setGeneric("transcriptomes", 
           function(object) standardGeneric("transcriptomes"))

#' @rdname KallistoExperiment-class
#' @aliases transcriptomes,KallistoExperiment-method
#' @export
setMethod("transcriptomes", "KallistoExperiment",
          function (object) return(object@transcriptomes))

#' @describeIn KallistoExperiment
#' Fetch transcripts for a gene, or all transcripts bundled by gene.
#' @export
setMethod("transcriptsBy", "KallistoExperiment",
          function(x, by="gene", ...) {
            if (by == "gene") { 
              split(x, mcols(x)$gene_name)
            } else { 
              return(x[mcols(x)$gene_name == by, ])
            }
          })

#' @describeIn KallistoExperiment 
#'
#' Fetch the matrix of MADs for estimated counts, if bootstraps were run. 
#' @export
setMethod("mad", "KallistoExperiment", function(x) assays(x)$est_counts_mad)

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 
