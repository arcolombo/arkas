#' @title convertes a summarizedExperiment to KallistoExperiment
SEtoKE<-function(from){
       
        txomes <- metadata(from)$transcriptomes
        kversion <- metadata(from)$kallistoVersion
        new("KallistoExperiment", from,
            kallistoVersion=kversion, transcriptomes=txomes)
      }
