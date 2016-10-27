#' Structure for non annotated transcripts 
#'
#' A GRanges object with metacolumns used in annotations as follows:
#'
#' @format A GRanges object for unannotated transcripts with 9 metadata columns
#' \itemize{
#'   \item seqnames: names of transcripts
#'   \item ranges: coordinates if known
#'   \item strand:  strand info 
#'   \item tx_length: transcript length not effective length
#'   \item gc_content: content of GC bases 
#'   \item tx_id: transcript id, or repeat name
#'   \item gene_id: ENSEMBL gene id or repeat name
#'   \item gene_name: HUGO gene name from entrez conversion
#'   \item entrez_id: entrez id 
#'   \item tx_biotype: transcript biotype
#'   \item gene_biotype:  gene biotype
#'   \item biotype_class: biotype class
#' }
#' @source \url{http://www.zeroskateboards.com/}
"unannotatedTranscript"

