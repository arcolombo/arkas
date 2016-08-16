#' @export
#' @rdname KallistoExperiment-class
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

#' @export
#' @rdname KallistoExperiment-class
setMethod("pData", "KallistoExperiment",
          function (object) return(colData(object)))

#' @export
#' @param value this is the replacement value for object
#' @aliases pData
#' @rdname KallistoExperiment-class
setReplaceMethod("pData", c("KallistoExperiment", "DataFrame"),
                 function (object, value) {
                   colData(object) <- value
                   return(object)
                 })
#' finds the effective length of a transcript entry
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @export
setGeneric("eff_length", function(object) standardGeneric("eff_length"))

#' @rdname KallistoExperiment-class
#' @aliases eff_length KallistoExperiment-method
#' @export
setMethod("eff_length", "KallistoExperiment",
          function (object) return(assays(object)$eff_length))
#' finds the transcript per million
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @export
setGeneric("tpm", function(object) standardGeneric("tpm"))

#' @rdname KallistoExperiment-class
#' @aliases tpm KallistoExperiment-method
#' @export
setMethod("tpm", "KallistoExperiment",
          function (object) return(assays(object)$tpm))   
#' finds the kallistoVersion used in quantification
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#'
#' @export
setGeneric("kallistoVersion", 
           function(object) standardGeneric("kallistoVersion"))
#' @rdname KallistoExperiment-class
#' @aliases kallistoVersion KallistoExperiment-method
#' @export
setMethod("kallistoVersion", "KallistoExperiment",
      function (object) return(object@kallistoVersion))


#' finds the transcriptomes pseudo-aligned against
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#' @param object KallistoExperiment
#' @export
setGeneric("transcriptomes", 
           function(object) standardGeneric("transcriptomes"))
#' @rdname KallistoExperiment-class
#' @aliases transcriptomes KallistoExperiment-method
#' @export
setMethod("transcriptomes", "KallistoExperiment",
          function (object) return(object@transcriptomes))

#
#setMethod("transcriptsBy", "KallistoExperiment",
 #         function(x, by="gene", ...) {
  #          if (by == "gene") { 
   #           split(x, mcols(x)$gene_name)
    #        } else { 
     #         return(x[mcols(x)$gene_name == by, ])
      #      }
       #   })

#' median abs deviation
#' @name KallistoExperiment-class
#' @rdname KallistoExperiment-class
#'
#' @export
setGeneric("mabDev", function(object) standardGeneric("mabDev"))


#' @rdname KallistoExperiment-class
#' @aliases mabDev KallistoExperiment-method
#' @export
setMethod("mabDev", "KallistoExperiment", function(object) assays(object)$est_counts_mad)

# FIXME: add method to retrieve normalization factors if ERCC spike-ins used 

#' @export
setAs("KallistoExperiment", "SummarizedExperiment", 
      function(from) {
        metanames <- names(from@metadata)
        metaorder <- c("transcriptomes", "kallistoVersion", metanames)
        from@metadata$transcriptomes <- from@transcriptomes
        from@metadata$kallistoVersion <- from@kallistoVersion
        from@metadata <- from@metadata[metaorder]
        if (!identical(colnames(from@assays$data[[1]]),rownames(from@colData))){
          for (i in names(from@assays$data)) {
            colnames(from@assays$data[[i]]) <- rownames(from@colData)
          }
        }
        SummarizedExperiment(assays=from@assays$data, 
                             rowRanges=from@rowRanges,
                             colData=from@colData, 
                             metadata=from@metadata)
      })


#' @export
setAs("SummarizedExperiment", "KallistoExperiment", 
      function(from) {
        txomes <- metadata(from)$transcriptomes
        kversion <- metadata(from)$kallistoVersion
        new("KallistoExperiment", from, 
            kallistoVersion=kversion, transcriptomes=txomes)
      })
