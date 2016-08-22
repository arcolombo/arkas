#' get the running version of Kallisto in the default path (if there is one)
#' 
#' @return  a string (the version of Kallisto that was found)
#' @importFrom TxDbLite strpop
#' @export
getKallistoVersion <- function() {
  strpop(system2("kallisto", args="version", stdout=TRUE))
}
