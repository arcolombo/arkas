#' groundwork for collapsing rows of a RangedSummarizedExperiment assay by TSS
#' 
#' @param kexp    the RSE, in this case a KallistoExperiment
#' @param groups  groups of chromosomes to keep (default: auto+sex+circular)
#' @param slop    how wide a range to reduce TSSes over for grouping (200bp) 
#' 
#' @return  the TSS for transcripts on chr1:22/X/Y/M to collapse by TSS,
#'          or NA for transcripts originating anywhere else (ERCC, repeats, &c)
#'
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb seqlevelsInGroup seqnames
#' @importFrom IRanges promoters reduce findOverlaps
#' @export
collapseByTss <- function(kexp, groups=c("auto","sex","circular"), slop=200) {

  chroms <- seqlevels(kexp)
  style <- seqlevelsStyle(kexp)
  if (any(grepl("Mmusculus", transcriptomes(kexp)))) species <- "Mus musculus"
  if (any(grepl("Hsapiens", transcriptomes(kexp)))) species <- "Homo sapiens"
  keepchroms <- do.call(c, 
                        lapply(groups, function(group) 
                               seqlevelsInGroup(chroms, group, species, style)))
  bundleable <- as.character(seqnames(kexp)) %in% keepchroms

  bundles <- rowRanges(promoters(kexp), up=slop, down=round(slop / 10))
  if (slop > 0) {
    # aggregate by reduced TSS 
    sloppyTSS <- reduce(bundles)
    tssMap <- findOverlaps(bundles, sloppyTSS)
    bundles$sloppyTSS <- as.character(sloppyTSS[subjectHits(tssMap)])
  } else { 
    bundles$sloppyTSS <- as.character(bundles)
  } 
  bundleIDs <- bundles$sloppyTSS
  is.na(bundleIDs) <- !bundleable
  return(bundleIDs)

}

