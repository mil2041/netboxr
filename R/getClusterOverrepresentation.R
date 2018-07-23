#' Run an enrichment analysis for a set of clustered data
#'
#' @param clusters a set of clusters as a list of vectors
#' @param gmtFile a path to a GMT file where entry IDs match up the IDs in the vectors
#' @param totalEntriesNum a number to be used as a total entry count for phyper,
#'   if NULL this will automatically be calculated as the number of unique entries
#'   across all clusters
#' @param verbose a boolean whether to display debugging information
#'
#' @return a data.frame with the following entries:
#'
#' \itemize{
#'   \item{entriesInGmtNum} {TODO}
#'   \item{gmtEntryNum} {TODO}
#' }
#'
#' @examples
#' #TODO
#'
#' @concept netboxr
#' @export
#' @importFrom paxtoolsr readGmt
getClusterOverrepresentation <- function(clusters, gmtFile, totalEntriesNum=NULL, verbose=FALSE) {
  gmtList <- readGmt(gmtFile)

  for(listEntry in 1:length(gmtList)) {
    gmtList[[listEntry]] <- gmtList[[listEntry]][which(gmtList[[listEntry]] != "" &
                                                       gmtList[[listEntry]] != "-")]
  }

  pvals <- NULL

  entriesInGmtNumVec <- NULL
  gmtEntryNumVec <- NULL
  totalClusterEntriesNumVec <- NULL
  clusterEntryNumVec <- NULL
  clusterIdxVec <- NULL
  gmtEntryIdxVec <- NULL
  gmtEntryNameVec <- NULL
  entriesInGmtVec <- NULL
  pvalVec <- NULL

  # totalEntriesNum is needed for phyper so it must be calculated prior to the
  # iteration.
  if (is.null(totalEntriesNum)){
    totalEntriesNum <- length(unique(c(clusters, recursive=TRUE))) # could be changed to all netbox input genes
  }

  for(clusterIdx in 1:length(clusters)) {
    for(gmtEntryIdx in 1:length(gmtList)) {
      entriesInGmt <- intersect(gmtList[[gmtEntryIdx]],
                                    clusters[[clusterIdx]])
      entriesInGmtNum <- length(entriesInGmt)

      gmtEntryName <- names(gmtList)[gmtEntryIdx]

      gmtEntryNum <- length(gmtList[[gmtEntryIdx]])

      clusterEntryNum <- length(clusters[[clusterIdx]])

      pval <- phyper(entriesInGmtNum-1,
                     gmtEntryNum,
                     totalEntriesNum-gmtEntryNum,
                     clusterEntryNum,
                     lower.tail=FALSE)

      pvals <- c(pvals, pval)

      # Store results in individual vectors
      entriesInGmtNumVec <- c(entriesInGmtNumVec, entriesInGmtNum)
      gmtEntryNumVec <- c(gmtEntryNumVec, gmtEntryNum)
      totalClusterEntriesNumVec <- c(totalClusterEntriesNumVec, totalEntriesNum)
      clusterEntryNumVec <- c(clusterEntryNumVec, clusterEntryNum)
      clusterIdxVec <- c(clusterIdxVec, clusterIdx)
      gmtEntryIdxVec <- c(gmtEntryIdxVec, gmtEntryIdx)
      gmtEntryNameVec <- c(gmtEntryNameVec, gmtEntryName)
      entriesInGmtVec <- c(entriesInGmtVec, paste(entriesInGmt, collapse=" "))
      pvalVec <- c(pvalVec, pval)

      if(verbose) {
        cat("GN: ", entriesInGmtNum,
            " PN: ", gmtEntryNum,
            " TG: ", totalEntriesNum,
            " CG: ", clusterEntryNum,
            " P: ", pval,
            " CI: ", clusterIdx,
            " PI: ", gmtEntryIdx,
            " P: ", gmtEntryName,
            " G: ", entriesInGmt,
            "\n")
      }
    }
  }

  pvalAdjust <- p.adjust(pvals, "fdr")

  results <- data.frame(entriesInGmtNum=entriesInGmtNumVec,
                        gmtEntryNum=gmtEntryNumVec,
                        totalClusterEntriesNum=totalClusterEntriesNumVec,
                        clusterEntryNum=clusterEntryNumVec,
                        clusterIdx=clusterIdxVec,
                        gmtEntryIdx=gmtEntryIdxVec,
                        gmtEntryName=gmtEntryNameVec,
                        entriesInGmt=entriesInGmtVec,
                        pval=pvalVec,
                        pvalAdjust=pvalAdjust)

  #save(results, file=output_file)

  return(results)
}
