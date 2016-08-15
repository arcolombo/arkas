#' A SummarizedExperiment subclass that stores multiple Kallisto runs 
#' 
#'
#' @slot transcriptomes   Transcriptomes against which reads were pseudoaligned
#' @slot kallistoVersion  The version of Kallisto used to pseudoalign the reads
#' @name KallistoExperiment-class
#' @section Slots: 
#'  \describe{
#'    \item{\code{transcriptomes}:}{ character class \code{"character"}, containing data from transcriptomes}
#'    \item{\code{kallistoVersion}:}{Object of class \code{"character"}, containing data that needs to go in kallistoVersion.}
#'  }  
#' 
#'  
#' @rdname KallistoExperiment-class  
#' @export 
setClass("KallistoExperiment",
         representation(transcriptomes="character", 
                        kallistoVersion="character"),
              contains="RangedSummarizedExperiment")

.checkAssayNames <- function (object, names) { # {{{
  if (!all(names %in% names(assays(object, withDimnames = FALSE)))) {
    return(sprintf("object of class '%s' needs assays with names '%s'", 
                   class(object), paste0(names, collapse = ", ")))
  } else {
    NULL
  }
} # }}}

setValidity("KallistoExperiment", function(object) { # {{{
  msg <- validMsg(NULL, NULL)
  msg <- .checkAssayNames(object, c("est_counts", "eff_length"))
  if (is.null(msg)) TRUE else msg
}) # }}}
