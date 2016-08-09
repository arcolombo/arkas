#' given a vector or matrix of est_counts and the same for eff_len, compute TPM
#' 
#' @param est_counts    estimated counts, usually from an HDF5 file
#' @param eff_len       effective transcipt lengths, usually from an HDF5 file
#'
#' @return              either a vector or a matrix of TPMs (matches the args)
#' 
#' @export
calcTpm <- function(est_counts, eff_len) { 
  if (is(est_counts, "matrix") || is(eff_len, "matrix")) {
    stopifnot(ncol(est_counts) == ncol(eff_len))
    stopifnot(nrow(est_counts) == nrow(eff_len))
    stopifnot(identical(colnames(est_counts), colnames(eff_len)))
    stopifnot(identical(rownames(est_counts), rownames(eff_len)))
    tpm <- matrix(0, nrow=nrow(est_counts), ncol=ncol(est_counts),
                  dimnames=list(rownames(est_counts), colnames(est_counts)))
    for (j in colnames(est_counts)) {
      tpm[, j] <- est_counts[, j] / eff_len[, j]
      totalMass <- sum(tpm[, j])
      tpm[, j] <- (tpm[, j] / totalMass) * 1e6
    }
    return(tpm)
  } else { # two vectors 
    stopifnot(length(est_counts) == length(eff_len))
    if (any(eff_len < 1)) message("You have transcripts w/effective length < 1")
    tpm <- (est_counts / eff_len)
    totalMass <- sum(tpm)
    (tpm / totalMass) * 1e6
  }
}
