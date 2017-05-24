#' @title converts a kallistoExperiment object to SummarizedExperiment
#' @param from must be a kexp
KEtoSE<-function(from){
  stopifnot(class(from)=="KallistoExperiment")
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
      }


