#' A SummarizedExperiment subclass that stores multiple Kallisto runs 
#' 
#'
#' @slot transcriptomes   Transcriptomes against which reads were pseudoaligned
#' @slot kallistoVersion  The version of Kallisto used to pseudoalign the reads
#' @param object a feature defined object
#' @param names an object with names defined.
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
