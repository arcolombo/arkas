#' @export
setMethod("counts", "KallistoExperiment",
          function (object) return(assays(object)$est_counts))

#' @export
setMethod("pData", "KallistoExperiment",
          function (object) return(colData(object)))

#' @export
setReplaceMethod("pData", c("KallistoExperiment", "DataFrame"),
                 function (object, value) {
                   colData(object) <- value
                   return(object)
                 })

#' @export
setGeneric("eff_length", function(object) standardGeneric("eff_length"))

#' @export
setMethod("eff_length", "KallistoExperiment",
          function (object) return(assays(object)$eff_length))

#' @export
setGeneric("tpm", function(object) standardGeneric("tpm"))

#' @export
setMethod("tpm", "KallistoExperiment",
          function (object) return(assays(object)$tpm))   
       
#' @export
setGeneric("kallistoVersion", 
           function(object) standardGeneric("kallistoVersion"))

#' @export
setMethod("kallistoVersion", "KallistoExperiment",
          function (object) return(object@kallistoVersion))

#' @export
setGeneric("transcriptomes", 
           function(object) standardGeneric("transcriptomes"))

#' @export
setMethod("transcriptomes", "KallistoExperiment",
          function (object) return(object@transcriptomes))

#' @export
setMethod("transcriptsBy", "KallistoExperiment",
          function(x, by="gene", ...) {
            if (by == "gene") { 
              split(x, mcols(x)$gene_name)
            } else { 
              return(x[mcols(x)$gene_name == by, ])
            }
          })

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
